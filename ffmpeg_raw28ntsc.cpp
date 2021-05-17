
// NTS: This is not like modern "posterize" filters where the pixels are quantizied to N levels then scaled out to 0..255
//      That requires a multiply/divide per pixel. Think old-school hardware where such operations were too expensive.
//      The "posterize" we emulate here is more the type where you run the video through an ADC, truncate the least significant
//      bits, then run back through a DAC on the other side (well within the realm of 1980s/1990s hardware)

#define __STDC_CONSTANT_MACROS
#define __STDC_LIMIT_MACROS

#include <sys/types.h>
#include <signal.h>
#include <stdint.h>
#include <assert.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <fcntl.h>
#include <math.h>

extern "C" {
#include <libavutil/opt.h>
#include <libavutil/avutil.h>
#include <libavutil/pixfmt.h>
#include <libavutil/pixdesc.h>
#include <libavutil/samplefmt.h>
#include <libavutil/pixelutils.h>

#include <libavcodec/avcodec.h>
#include <libavcodec/version.h>

#include <libavformat/avformat.h>
#include <libavformat/avio.h>
#include <libavformat/version.h>

#include <libswscale/swscale.h>
#include <libswscale/version.h>

#include <libswresample/swresample.h>
#include <libswresample/version.h>
}

using namespace std;

#include <map>
#include <list>
#include <string>
#include <vector>
#include <stdexcept>
#include <algorithm>

/* return a floating point value specifying what to scale the sample
 * value by to reduce it from full volume to dB decibels */
double dBFS(double dB)
{
	/* 10 ^ (dB / 20),
	   based on reversing the formula for converting samples to decibels:
	   dB = 20.0 * log10(sample);
	   where "sample" is -1.0 <= x <= 1.0 */
	return pow(10.0,dB / 20.0);
}

/* attenuate a sample value by this many dBFS */
/* so if you want to reduce it by 20dBFS you pass -20 as dB */
double attenuate_dBFS(double sample,double dB)
{
	return sample * dBFS(dB);
}

/* opposite: convert sample to decibels */
double dBFS_measure(double sample) {
	return 20.0 * log10(sample);
}

// lowpass filter
// you can make it a highpass filter by applying a lowpass then subtracting from source.
class LowpassFilter {
public:
	LowpassFilter() : timeInterval(0), cutoff(0), alpha(0), prev(0), tau(0) {
	}
	void setFilter(const double rate/*sample rate of audio*/,const double hz/*cutoff*/) {
#ifndef M_PI
#error your math.h does not include M_PI constant
#endif
		timeInterval = 1.0 / rate;
		tau = 1 / (hz * 2 * M_PI);
		cutoff = hz;
		alpha = timeInterval / (tau + timeInterval);
	}
	void resetFilter(const double val=0) {
		prev = val;
	}
	double lowpass(const double sample) {
		const double stage1 = sample * alpha;
		const double stage2 = prev - (prev * alpha); /* NTS: Instead of prev * (1.0 - alpha) */
		return (prev = (stage1 + stage2)); /* prev = stage1+stage2 then return prev */
	}
	double highpass(const double sample) {
		const double stage1 = sample * alpha;
		const double stage2 = prev - (prev * alpha); /* NTS: Instead of prev * (1.0 - alpha) */
		return sample - (prev = (stage1 + stage2)); /* prev = stage1+stage2 then return (sample - prev) */
	}
public:
	double			timeInterval;
	double			cutoff;
	double			alpha; /* timeInterval / (tau + timeInterval) */
	double			prev;
	double			tau;
};

class HiLoPair {
public:
	LowpassFilter		hi,lo;	// highpass, lowpass
public:
	void setFilter(const double rate/*sample rate of audio*/,const double low_hz,const double high_hz) {
		lo.setFilter(rate,low_hz);
		hi.setFilter(rate,high_hz);
	}
	double filter(const double sample) {
		return hi.highpass(lo.lowpass(sample)); /* first lowpass, then highpass */
	}
};

class HiLoPass : public vector<HiLoPair> { // all passes, one sample of one channel
public:
	HiLoPass() : vector() { }
public:
	void setFilter(const double rate/*sample rate of audio*/,const double low_hz,const double high_hz) {
		for (size_t i=0;i < size();i++) (*this)[i].setFilter(rate,low_hz,high_hz);
	}
	double filter(double sample) {
		for (size_t i=0;i < size();i++) sample = (*this)[i].lo.lowpass(sample);
		for (size_t i=0;i < size();i++) sample = (*this)[i].hi.highpass(sample);
		return sample;
	}
	void init(const unsigned int passes) {
		clear();
		resize(passes);
		assert(size() >= passes);
	}
};

class HiLoSample : public vector<HiLoPass> { // all passes, all channels of one sample period
public:
	HiLoSample() : vector() { }
public:
	void init(const unsigned int channels,const unsigned int passes) {
		clear();
		resize(channels);
		assert(size() >= channels);
		for (size_t i=0;i < size();i++) (*this)[i].init(passes);
	}
	void setFilter(const double rate/*sample rate of audio*/,const double low_hz,const double high_hz) {
		for (size_t i=0;i < size();i++) (*this)[i].setFilter(rate,low_hz,high_hz);
	}
};

class HiLoComboPass {
public:
	HiLoComboPass() : passes(0), channels(0), rate(0), low_cutoff(0), high_cutoff(0) {
	}
	~HiLoComboPass() {
		clear();
	}
	void setChannels(const size_t _channels) {
		if (channels != _channels) {
			clear();
			channels = _channels;
		}
	}
	void setCutoff(const double _low_cutoff,const double _high_cutoff) {
		if (low_cutoff != _low_cutoff || high_cutoff != _high_cutoff) {
			clear();
			low_cutoff = _low_cutoff;
			high_cutoff = _high_cutoff;
		}
	}
	void setRate(const double _rate) {
		if (rate != _rate) {
			clear();
			rate = _rate;
		}
	}
	void setPasses(const size_t _passes) {
		if (passes != _passes) {
			clear();
			passes = _passes;
		}
	}
	void clear() {
		audiostate.clear();
	}
	void init() {
		clear();
		if (channels == 0 || passes == 0 || rate == 0 || low_cutoff == 0 || high_cutoff == 0) return;
		audiostate.init(channels,passes);
		audiostate.setFilter(rate,low_cutoff,high_cutoff);
	}
public:
	double		rate;
	size_t		passes;
	size_t		channels;
	double		low_cutoff;
	double		high_cutoff;
	HiLoSample	audiostate;
};

std::list<string>           src_composite;

unsigned long long          src_byte_counter = 0; /* at beginning of input buffer */
int                         src_fd = -1;

bool            mark_sync = false;
bool            disable_sync = false;
bool            disable_wp_equ = false;
bool            show_subcarrier = false;
bool            use_422_colorspace = true;
bool            disable_subcarrier = false;
bool            disable_equalization = false;
AVRational	output_field_rate = { 60000, 1001 };	// NTSC 60Hz default
int		output_width = 720;
int		output_height = 480;
bool        input_ntsc = false;
bool		output_ntsc = true;	// NTSC color subcarrier emulation
bool		output_pal = false;	// PAL color subcarrier emulation
int		output_audio_channels = 2;	// VHS stereo (set to 1 for mono)
int		output_audio_rate = 44100;	// VHS Hi-Fi goes up to 20KHz
string	sratep = "ntsc28";

static double           subcarrier_freq = 0;
static double           sample_rate = 0;
static double           one_frame_time = 0;
static double           one_scanline_time = 0;
static unsigned int     one_scanline_raw_length = 0;
static double           one_scanline_width = 0;
static double           one_scanline_width_err = 0;

void NTSCAnyMHz(const char *str) {
	sample_rate =				atof(str);
}

void NTSC28MHz() {
	sample_rate =				((315000000.00 * 8.0) / 88.00);        /* 315/88 * MHz = about 28.636363 MHz */
}

void Do40MHz() {
	sample_rate =				40000000.00;                           /* 40MHz */
}

void compute_NTSC() {
	subcarrier_freq =			315000000.00 / 88.00;                   /* 315/88MHz or about 3.5795454...MHz */
	one_frame_time =			sample_rate / (30000.00 / 1001.00);  /* 30000/1001 = about 29.97 */
	one_scanline_time =			one_frame_time / 525.00;            /* one scanline */
	one_scanline_raw_length =	(unsigned int)(one_scanline_time + 0.5);
	one_scanline_width =		one_scanline_raw_length;
	one_scanline_width_err =	0;
}

signed int                                  int_scanline[4096];
signed int                                  int_chroma[4096];
signed int                                  int_luma2[4096];
signed int                                  int_luma[4096];
unsigned char                               read_tmp[4096];

struct oneprocsamp {
    uint8_t                                 raw;
    uint8_t                                 hsync_dc_raw;
};

std::vector<oneprocsamp>                    input_samples;
std::vector<oneprocsamp>::iterator          input_samples_read,input_samples_end;

unsigned long long total_count_src() {
    return src_byte_counter + (input_samples_read - input_samples.begin());
}

size_t count_src() {
    assert(input_samples_read <= input_samples_end);
    return (input_samples_end - input_samples_read);
}

void empty_src() {
    src_byte_counter = total_count_src();
    input_samples_read = input_samples_end = input_samples.begin();
}

void flush_src() {
    assert(input_samples_read <= input_samples_end);
    if (input_samples_read != input_samples.begin()) {
        src_byte_counter = total_count_src();
        size_t move = input_samples_read - input_samples.begin();
        assert(move != 0);
        size_t todo = input_samples_end - input_samples_read;
        if (todo > 0) memmove(&(*input_samples.begin()),&(*input_samples_read),todo*sizeof(oneprocsamp));
        input_samples_read -= move;
        assert(input_samples_read == input_samples.begin());
        input_samples_end -= move;
    }
}

void do_filter_new_input(oneprocsamp *samples,int count);

void refill_src() {
    if (src_fd >= 0) {
        assert(input_samples_read <= input_samples_end);
        while (input_samples_end < input_samples.end()) {
            size_t todo = input_samples.end() - input_samples_end;
            if (todo > sizeof(read_tmp)) todo = sizeof(read_tmp);
            assert(todo <= input_samples.size());
            assert(todo <= sizeof(read_tmp));
            int rd = read(src_fd,read_tmp,todo);
            if (rd > 0) {
                assert(rd <= todo);
                for (int x=0;x < rd;x++)
                    input_samples_end[x].raw = read_tmp[x];

                do_filter_new_input(&(*input_samples_end),rd);
                assert((size_t)rd <= todo);
                input_samples_end += (size_t)rd;
                assert(input_samples_end <= input_samples.end());
            }
            if (rd < todo) break;
        }
    }
}

void lazy_flush_src() {
    if (input_samples_read > (input_samples.begin()+(input_samples.size()/2u))) flush_src();
    refill_src();
}

void rewind_src() {
    if (src_fd >= 0) lseek(src_fd,0,SEEK_SET);
}

bool open_src() {
    while (src_fd < 0) {
        if (src_composite.empty()) return false;

        std::string path = src_composite.front();
        src_composite.pop_front();

        if (path == "-") { // STDIN
            src_fd = dup(0/*STDIN*/);
            close(0/*STDIN*/);
            if (src_fd < 0) return false;
        }
        else {
            src_fd = open(path.c_str(),O_RDONLY);
        }
    }

    input_samples.resize(one_scanline_raw_length*2048);
    input_samples_read = input_samples_end = input_samples.begin();
    return true;
}

void close_src() {
    if (src_fd >= 0) {
        close(src_fd);
        src_fd = -1;
    }
}

#define RGBTRIPLET(r,g,b)       (((uint32_t)(r) << (uint32_t)16) + ((uint32_t)(g) << (uint32_t)8) + ((uint32_t)(b) << (uint32_t)0))

AVFormatContext*	        output_avfmt = NULL;
AVStream*		            output_avstream_video = NULL;	// do not free
AVCodecContext*		        output_avstream_video_codec_context = NULL; // do not free
AVFrame*		            output_avstream_video_frame = NULL;         // ARGB
AVFrame*		            output_avstream_video_encode_frame = NULL;  // 4:2:2 or 4:2:0
struct SwsContext*          output_avstream_video_resampler = NULL;

std::string                 output_file;

volatile int DIE = 0;

void sigma(int x) {
	if (++DIE >= 20) abort();
}

void preset_PAL() {
	output_field_rate.num = 25;
	output_field_rate.den = 1;
	output_height = 576;
	output_width = 720;
	output_pal = true;
	output_ntsc = false;
}

void preset_NTSC() {
	output_field_rate.num = 60000;
	output_field_rate.den = 1001;
	output_height = 262;
	output_width = (one_scanline_raw_length + 1) & (~1);
	output_pal = false;
	output_ntsc = true;
}

void preset_720p60() {
	output_field_rate.num = 30000;
	output_field_rate.den = 1001;
	output_height = 720;
	output_width = 1280;
	output_pal = false;
	output_ntsc = true;
}

void preset_1080p60() {
	output_field_rate.num = 30000;
	output_field_rate.den = 1001;
	output_height = 1080;
	output_width = 1920;
	output_pal = false;
	output_ntsc = true;
}

static void help(const char *arg0) {
	fprintf(stderr,"%s [options]\n",arg0);
	fprintf(stderr," -i <input file>               you can specify more than one input file, in order of layering\n");
	fprintf(stderr," -o <output file>\n");
	fprintf(stderr," -s <rate>                     ntsc28, 40mhz\n");
}

static int parse_argv(int argc,char **argv) {
	const char *a;
	int i;

	for (i=1;i < argc;) {
		a = argv[i++];

		if (*a == '-') {
			do { a++; } while (*a == '-');

			if (!strcmp(a,"h") || !strcmp(a,"help")) {
				help(argv[0]);
				return 1;
            }
            else if (!strcmp(a,"marksig")) {
                mark_sync = true;
            }
            else if (!strcmp(a,"noequ")) {
                disable_equalization = true;
            }
            else if (!strcmp(a,"nowequ")) {
                disable_wp_equ = true;
            }
            else if (!strcmp(a,"nosig")) {
                disable_sync = true;
            }
            else if (!strcmp(a,"nosc")) {
                disable_subcarrier = true;
            }
            else if (!strcmp(a,"showsc")) {
                show_subcarrier = true;
            }
            else if (!strcmp(a,"s")) {
                a = argv[i++];
                if (a == NULL) return 1;
                sratep = a;
            }
            else if (!strcmp(a,"width")) {
                a = argv[i++];
                if (a == NULL) return 1;
                output_width = (int)strtoul(a,NULL,0);
                if (output_width < 32) return 1;
            }
            else if (!strcmp(a,"i")) {
                a = argv[i++];
                if (a == NULL) return 1;
                src_composite.push_back(a);
            }
            else if (!strcmp(a,"o")) {
                a = argv[i++];
                if (a == NULL) return 1;
                output_file = a;
            }
            else if (!strcmp(a,"422")) {
                use_422_colorspace = true;
            }
            else if (!strcmp(a,"420")) {
                use_422_colorspace = false;
            }
            else if (!strcmp(a,"inntsc")) {
                input_ntsc = true;
            }
			else {
				fprintf(stderr,"Unknown switch '%s'\n",a);
				return 1;
			}
		}
		else {
			fprintf(stderr,"Unhandled arg '%s'\n",a);
			return 1;
		}
	}

    if (output_file.empty()) {
        fprintf(stderr,"No output file specified\n");
        return 1;
    }
    if (src_composite.empty()) {
        fprintf(stderr,"No input file specified\n");
        return 1;
    }

	return 0;
}

void output_frame(AVFrame *frame,unsigned long long field_number) {
	int gotit = 0;
	AVPacket pkt;

	av_init_packet(&pkt);
	if (av_new_packet(&pkt,50000000/8) < 0) {
		fprintf(stderr,"Failed to alloc vid packet\n");
		return;
	}

	frame->key_frame = (field_number % (15ULL * 2ULL)) == 0 ? 1 : 0;

    {
		frame->interlaced_frame = 0;
		frame->pts = field_number;
		pkt.pts = field_number;
		pkt.dts = field_number;
	}

	fprintf(stderr,"\x0D" "Output field %llu ",field_number); fflush(stderr);
    if (avcodec_encode_video2(output_avstream_video_codec_context,&pkt,frame,&gotit) == 0) {
        if (gotit) {
            pkt.stream_index = output_avstream_video->index;
            av_packet_rescale_ts(&pkt,output_avstream_video_codec_context->time_base,output_avstream_video->time_base);

            if (av_interleaved_write_frame(output_avfmt,&pkt) < 0)
                fprintf(stderr,"AV write frame failed video\n");
        }
    }

	av_packet_unref(&pkt);
}

// From cxadc you can't assume specific sample values when detecting sync pulses,
// so the first step is to "DC normalize" the samples using the hsync pulse detection.
#define                     hsync_dc_detect_passes     (3)
LowpassFilter               hsync_dc_detect[3];
double                      hsync_dc_level = 128.0;
vector<uint8_t>             hsync_dc_detect_delay;
vector<uint8_t>::iterator   hsync_dc_detect_delay_i;

uint8_t                     sync_threshhold = (uint8_t)(192 * 0.25 * 0.5);
uint8_t                     blank_level = (uint8_t)0;
uint8_t                     white_level = (uint8_t)192;

oneprocsamp hsync_dc_proc(oneprocsamp v) {
    double lv = v.raw;

    for (size_t i=0;i < hsync_dc_detect_passes;i++)
        lv = hsync_dc_detect[i].lowpass(lv);

    if (hsync_dc_level > lv) {
        hsync_dc_level = lv; // lowpass filter already smooths it out
    }
    else {
        const double a = 1.0 / (one_frame_time * 0.6);
        hsync_dc_level = (hsync_dc_level * (1.0 - a)) + (lv * a);
    }

    {
        int x = v.raw - hsync_dc_level;
        if (x < 0) x = 0;
        if (x > 255) x = 255;
        v.raw = (uint8_t)x;

        assert(hsync_dc_detect_delay_i >= hsync_dc_detect_delay.begin());
        assert(hsync_dc_detect_delay_i < hsync_dc_detect_delay.end());
        uint8_t ov = *hsync_dc_detect_delay_i;
        *hsync_dc_detect_delay_i = v.raw;
        v.raw = ov;
        hsync_dc_detect_delay_i++;
        if (hsync_dc_detect_delay_i >= hsync_dc_detect_delay.end())
            hsync_dc_detect_delay_i  = hsync_dc_detect_delay.begin();
    }

    {
        int x = lv - hsync_dc_level;
        if (x < 0) x = 0;
        if (x > 255) x = 255;
        v.hsync_dc_raw = (uint8_t)x;
    }

    if (mark_sync && v.hsync_dc_raw < sync_threshhold)
        v.raw = 255;

    return v;
}

void do_filter_new_input(oneprocsamp *samples,int count) {
    for (int i=0;i < count;i++)
        samples[i] = hsync_dc_proc(samples[i]);
}

// This code assumes ARGB and the frame match resolution/
void composite_layer(AVFrame *dstframe,unsigned int field,unsigned long long fieldno) {
    double sx,sy,tx,ty;
    unsigned int dx,dy;
    unsigned int ystep;
    double dot_radius;
    unsigned int x,y;
    double sigscalxy;
    double frame_t;
    uint32_t rgba;
    double signal;
    double t;

    if (dstframe == NULL) return;
    if (dstframe->data[0] == NULL) return;
    if (dstframe->linesize[0] < (dstframe->width*4)) return; // ARGB

    double fieldt = ((double)fieldno * output_field_rate.den) / output_field_rate.num;
    double filet = (double)total_count_src() / sample_rate;

    /* look for vsync. We're using the hsync_dc_raw which is the lowpass filtered version of the video signal. */
    one_scanline_width_err = 0;
    lazy_flush_src();
    refill_src();

    if (!disable_sync) {
        vector<oneprocsamp>::iterator i = input_samples_read,last_pulse = input_samples_read;
        int vsb_count = 0;

        while (i < input_samples_end) {
            while (i < input_samples_end && (*i).hsync_dc_raw >= sync_threshhold) i++;
            vector<oneprocsamp>::iterator si = i;
            while (i < input_samples_end && (*i).hsync_dc_raw < sync_threshhold) i++;
            vector<oneprocsamp>::iterator ei = i;

            size_t synclen = (size_t)(ei-si);

            if (synclen >= (int)(one_scanline_raw_length * 0.3)) { /* vertical sync pulse (0.5H - 0.07H) */
                i = si + (int)(one_scanline_raw_length * 0.3);
                if (i < ei) i = ei;
                vsb_count++;
            }
            else if (synclen >= (int)(one_scanline_raw_length * 0.06)) { /* hsync pulse */
                if (vsb_count >= (3*3)) {
                    input_samples_read = si + (synclen / 2); /* use the center of the pulse for stability */
                    i = si + (int)(one_scanline_raw_length * 0.3);
                    if (i < ei) i = ei;
                    break;
                }
            }
            else if (synclen >= (int)(one_scanline_raw_length * 0.02)) { /* equalization pulse */
                i = si + (int)(one_scanline_raw_length * 0.3);
                if (i < ei) i = ei;
                vsb_count++;

                /* use the values to calibrate the black level */
                {
                    vector<oneprocsamp>::iterator j = si;
                    int mina = 0,mind = 0;
                    int maxa = 0,maxd = 0;

                    while (j < i) {
                        if ((*j).hsync_dc_raw >= sync_threshhold) {
                            maxa += (*j).hsync_dc_raw;
                            maxd++;
                        }
                        else {
                            mina += (*j).hsync_dc_raw;
                            mind++;
                        }

                        j++;
                    }

                    if (mind > 0) mina /= mind;
                    if (maxd > 0) maxa /= maxd;

                    int nwhite = (uint8_t)min(max((int)(maxa + ((maxa - mina) / 0.25)),maxa+1),240);
                    white_level = ((white_level * 7) + nwhite + 4) / 8;

                    int nblack = maxa;
                    blank_level = ((blank_level * 7) + nblack + 4) / 8;
                }
            }

            last_pulse = si;
        }
    }

    {
        vector<oneprocsamp>::iterator input_scan = input_samples_read;
        vector<oneprocsamp>::iterator input_start = input_samples_read;

        /* render normally */
        for (y=0;y < dstframe->height && (input_scan+(one_scanline_raw_length*2)) < input_samples_end;y++) {
            /* use interpolation because our concept of "one scanline" isn't exactly an integer.
             * without interpolation, adjustment can be a bit jagged. */
            {
                int a = (int)floor(one_scanline_width_err * 256);
                if (a < 0) a = 0;
                if (a > 256) a = 256;
                for (x=0;x < (one_scanline_raw_length+16);x++) {
                    int v = ((input_scan[x].raw * (256 - a)) + (input_scan[x+1].raw * a)) >> 8;
                    if (!disable_equalization) {
                        v -= blank_level;
                        if (!disable_wp_equ) {
                            v = (v * 255) / (white_level - blank_level);
                        }
                    }
                    if (v < 0) v = 0;
                    if (v > 255) v = 255;
                    int_scanline[x] = (uint8_t)v;
                }
            }

            if (disable_subcarrier) {
                for (x=0;x < (one_scanline_raw_length+16);x++) {
                    int_luma[x] = int_scanline[x];
                    int_chroma[x] = 0;
                }
            }
            else {
                /* 28.6MHz is exactly 8x the chroma subcarrier.
                 * So instead of complex filtering, we can just average the scanline with itself delayed 4 (half of 8) samples
                 * and make use of destructive interference to filter out the chroma subcarrier. This would not work if using,
                 * say, 40MHz. Then lowpass a bit more to remove remaining subcarrier.
                 *
                 * Once luma is determined, subtract from original to get chroma subcarrier.
                 * Note that some of the edge detail in luma will also end up in the chroma subcarrier.
                 * We'll lowpass the decoded I and Q later to help filter that out, but it is the reason
                 * fine details have color artifacts with composite video. */
                for (x=0;x < (one_scanline_raw_length+16-4);x++)
                    int_luma2[x] = (int_scanline[x] + int_scanline[x+4] + 1) / 2;
                for (x=0;x < (one_scanline_raw_length+16-4);x++)
                    int_luma[x] = (int_luma2[x] + int_luma2[x+1] + int_luma2[x+2] + int_luma2[x+3] + 2) / 4;
                for (x=0;x < (one_scanline_raw_length+16-4);x++)
                    int_chroma[x] = int_scanline[x] - int_luma[x];
            }

            uint32_t *dst = (uint32_t*)(dstframe->data[0] + (dstframe->linesize[0] * y));
            for (x=0;x < dstframe->width;x++) {
                int r,g,b;
                int Y = int_luma[x];

                if (show_subcarrier)
                    Y = int_chroma[x] + 128;

                r = g = b = Y;
                if (r < 0) r = 0;
                if (g < 0) g = 0;
                if (b < 0) b = 0;
                if (r > 255) r = 255;
                if (g > 255) g = 255;
                if (b > 255) b = 255;

                dst[x] = RGBTRIPLET(r,g,b);
            }

            {
                unsigned int adj = floor(one_scanline_width);
                one_scanline_width_err += one_scanline_width - adj;
                if (one_scanline_width_err >= 1.0) {
                    one_scanline_width_err -= 1.0;
                    adj++;
                }
                input_scan += adj;
                if (input_scan > input_samples_end)
                    input_scan = input_samples_end;
            }

            if (!disable_sync) {
                vector<oneprocsamp>::iterator i = input_scan;
                int vsb_count = 0;

                if (i > input_samples_read) {
                    size_t avail = (size_t)(i-input_samples_read);
                    if (avail >= (one_scanline_raw_length * 0.1))
                        avail  = (one_scanline_raw_length * 0.1);

                    i -= avail;
                }

                while (i < input_samples_end) {
                    while (i < input_samples_end && (*i).hsync_dc_raw >= sync_threshhold) i++;
                    vector<oneprocsamp>::iterator si = i;
                    while (i < input_samples_end && (*i).hsync_dc_raw < sync_threshhold) i++;
                    vector<oneprocsamp>::iterator ei = i;

                    size_t synclen = (size_t)(ei-si);

                    if (synclen >= (int)(one_scanline_raw_length * 0.3)) { /* vertical sync pulse (0.5H - 0.07H) */
                        i = si + (int)(one_scanline_raw_length * 0.3);
                        if (i < ei) i = ei;
                        vsb_count++;
                    }
                    else if (synclen >= (int)(one_scanline_raw_length * 0.06)) { /* hsync pulse */
                        input_scan = si + (synclen / 2); /* use the center of the pulse for stability */
                        break;
                    }
                    else if (synclen >= (int)(one_scanline_raw_length * 0.02)) { /* equalization pulse */
                        i = si + (int)(one_scanline_raw_length * 0.3);
                        if (i < ei) i = ei;
                        vsb_count++;
                    }

                    if (vsb_count >= (3*3)) {
                        y = INT_MAX;
                        break;
                    }
                }
            }
        }

        if (disable_sync)
            input_samples_read = input_scan;

        {
            vector<oneprocsamp>::iterator input_should = input_start + (one_scanline_raw_length * 240);

            if (input_should > input_samples_end)
                input_should = input_samples_end;

            if (input_samples_read < input_should)
                input_samples_read = input_should;
        }
    }
}

int main(int argc,char **argv) {
    if (parse_argv(argc,argv))
		return 1;

	av_register_all();
	avformat_network_init();
	avcodec_register_all();

    /* open output file */
	assert(output_avfmt == NULL);
	if (avformat_alloc_output_context2(&output_avfmt,NULL,NULL,output_file.c_str()) < 0) {
		fprintf(stderr,"Failed to open output file\n");
		return 1;
	}

    if (sratep == "ntsc28")
        NTSC28MHz();
    else if (sratep == "40mhz")
        Do40MHz();
    else if (!sratep.empty() && isdigit(sratep[0]))
        NTSCAnyMHz(sratep.c_str());
    else {
        fprintf(stderr,"Unknown -s preset '%s'\n",sratep.c_str());
        NTSC28MHz();
    }
    compute_NTSC();
    preset_NTSC();

    fprintf(stderr,"Subcarrier:             %.3f\n",subcarrier_freq);
    fprintf(stderr,"Sample rate:            %.3f\n",sample_rate);
    fprintf(stderr,"One frame duration:     %.3f (%.3fHz)\n",one_frame_time,sample_rate / one_frame_time);
    fprintf(stderr,"One field duration:     %.3f (%.3fHz)\n",one_frame_time / 2.0,sample_rate / (one_frame_time / 2.0));
    fprintf(stderr,"One scanline duration:  %.3f (%.3fHz)\n",one_scanline_time,sample_rate / one_scanline_time);
    fprintf(stderr,"Raw render to:          %u\n",one_scanline_raw_length);

    /* raw samples need to be delayed to match hsync detect, lowpass filtering has a delay effect */
    hsync_dc_detect_delay.resize((size_t)((one_scanline_time * 0.075 * 0.75) * 0.5));
    hsync_dc_detect_delay_i = hsync_dc_detect_delay.begin();
    /* set up lowpass filter */
    for (size_t i=0;i < hsync_dc_detect_passes;i++) {
        hsync_dc_detect[i].setFilter(sample_rate,sample_rate / (one_scanline_time * 0.075 * 0.75));
        for (size_t j=0;j < one_frame_time;j++) hsync_dc_detect[i].lowpass(128);
    }

    if (!open_src()) {
        fprintf(stderr,"Failed to open src\n");
        return 1;
    }

	{
		output_avstream_video = avformat_new_stream(output_avfmt, NULL);
		if (output_avstream_video == NULL) {
			fprintf(stderr,"Unable to create output video stream\n");
			return 1;
		}

		output_avstream_video_codec_context = output_avstream_video->codec;
		if (output_avstream_video_codec_context == NULL) {
			fprintf(stderr,"Output stream video no codec context?\n");
			return 1;
		}

		// FIXME: How do I get FFMPEG to write raw YUV 4:2:2?
		avcodec_get_context_defaults3(output_avstream_video_codec_context,avcodec_find_encoder(AV_CODEC_ID_H264));
		output_avstream_video_codec_context->width = output_width;
		output_avstream_video_codec_context->height = output_height;
		output_avstream_video_codec_context->sample_aspect_ratio = (AVRational){output_height*4, output_width*3};
		output_avstream_video_codec_context->pix_fmt = use_422_colorspace ? AV_PIX_FMT_YUV422P : AV_PIX_FMT_YUV420P;
		output_avstream_video_codec_context->gop_size = 15;
		output_avstream_video_codec_context->max_b_frames = 0;
		output_avstream_video_codec_context->time_base = (AVRational){output_field_rate.den, output_field_rate.num};
        output_avstream_video_codec_context->bit_rate = 15000000;

		output_avstream_video->time_base = output_avstream_video_codec_context->time_base;
		if (output_avfmt->oformat->flags & AVFMT_GLOBALHEADER)
			output_avstream_video_codec_context->flags |= AV_CODEC_FLAG_GLOBAL_HEADER;

		if (avcodec_open2(output_avstream_video_codec_context,avcodec_find_encoder(AV_CODEC_ID_H264),NULL) < 0) {
			fprintf(stderr,"Output stream cannot open codec\n");
			return 1;
		}
	}

	if (!(output_avfmt->oformat->flags & AVFMT_NOFILE)) {
		if (avio_open(&output_avfmt->pb, output_file.c_str(), AVIO_FLAG_WRITE) < 0) {
			fprintf(stderr,"Output file cannot open file\n");
			return 1;
		}
	}

	if (avformat_write_header(output_avfmt,NULL) < 0) {
		fprintf(stderr,"Failed to write header\n");
		return 1;
	}

	/* soft break on CTRL+C */
	signal(SIGINT,sigma);
	signal(SIGHUP,sigma);
	signal(SIGQUIT,sigma);
	signal(SIGTERM,sigma);

    /* prepare video encoding */
    output_avstream_video_frame = av_frame_alloc();
    if (output_avstream_video_frame == NULL) {
        fprintf(stderr,"Failed to alloc video frame\n");
        return 1;
    }
    output_avstream_video_frame->format = AV_PIX_FMT_BGRA;
    output_avstream_video_frame->height = output_height;
    output_avstream_video_frame->width = output_width;
    if (av_frame_get_buffer(output_avstream_video_frame,64) < 0) {
        fprintf(stderr,"Failed to alloc render frame\n");
        return 1;
    }

    {
        output_avstream_video_encode_frame = av_frame_alloc();
        if (output_avstream_video_encode_frame == NULL) {
            fprintf(stderr,"Failed to alloc video frame3\n");
            return 1;
        }
        av_frame_set_colorspace(output_avstream_video_encode_frame,AVCOL_SPC_SMPTE170M);
        av_frame_set_color_range(output_avstream_video_encode_frame,AVCOL_RANGE_MPEG);
        output_avstream_video_encode_frame->format = output_avstream_video_codec_context->pix_fmt;
        output_avstream_video_encode_frame->height = output_height;
        output_avstream_video_encode_frame->width = output_width;
        if (av_frame_get_buffer(output_avstream_video_encode_frame,64) < 0) {
            fprintf(stderr,"Failed to alloc render frame2\n");
            return 1;
        }
    }

    if (output_avstream_video_resampler == NULL) {
        output_avstream_video_resampler = sws_getContext(
                // source
                output_avstream_video_frame->width,
                output_avstream_video_frame->height,
                (AVPixelFormat)output_avstream_video_frame->format,
                // dest
                output_avstream_video_encode_frame->width,
                output_avstream_video_encode_frame->height,
                (AVPixelFormat)output_avstream_video_encode_frame->format,
                // opt
                SWS_BILINEAR, NULL, NULL, NULL);
        if (output_avstream_video_resampler == NULL) {
            fprintf(stderr,"Failed to alloc ARGB -> codec converter\n");
            return 1;
        }
    }

    /* run all inputs and render to output, until done */
    {
        bool eof,copyaud;
        signed long long current=0;

        do {
            if (DIE) break;

            refill_src();
            if (count_src() < (one_scanline_raw_length*256)) {
                close_src();
                if (!open_src()) break;
            }

            memset(output_avstream_video_frame->data[0],0,output_avstream_video_frame->linesize[0]*output_avstream_video_frame->height);

            // composite the layer, keying against the color. all code assumes ARGB
            composite_layer(output_avstream_video_frame,(current & 1) ^ 1,current);

            // convert ARGB to whatever the codec demands, and encode
            output_avstream_video_encode_frame->pts = output_avstream_video_frame->pts;
            output_avstream_video_encode_frame->pkt_pts = output_avstream_video_frame->pkt_pts;
            output_avstream_video_encode_frame->pkt_dts = output_avstream_video_frame->pkt_dts;
            output_avstream_video_encode_frame->top_field_first = output_avstream_video_frame->top_field_first;
            output_avstream_video_encode_frame->interlaced_frame = output_avstream_video_frame->interlaced_frame;

            if (sws_scale(output_avstream_video_resampler,
                        // source
                        output_avstream_video_frame->data,
                        output_avstream_video_frame->linesize,
                        0,output_avstream_video_frame->height,
                        // dest
                        output_avstream_video_encode_frame->data,
                        output_avstream_video_encode_frame->linesize) <= 0)
                fprintf(stderr,"WARNING: sws_scale failed\n");

            output_frame(output_avstream_video_encode_frame,current);
            current++;
        } while (!eof);
    }

    /* flush encoder delay */
    do {
        AVPacket pkt;
        int gotit=0;

        av_init_packet(&pkt);
        if (av_new_packet(&pkt,50000000/8) < 0) break;

        if (avcodec_encode_video2(output_avstream_video_codec_context,&pkt,NULL,&gotit) == 0) {
            if (gotit) {
                pkt.stream_index = output_avstream_video->index;
                av_packet_rescale_ts(&pkt,output_avstream_video_codec_context->time_base,output_avstream_video->time_base);

                if (av_interleaved_write_frame(output_avfmt,&pkt) < 0)
                    fprintf(stderr,"AV write frame failed video\n");
            }
        }

        av_packet_unref(&pkt);
        if (!gotit) break;
    } while (1);

    /* close output */
    if (output_avstream_video_resampler != NULL) {
        sws_freeContext(output_avstream_video_resampler);
        output_avstream_video_resampler = NULL;
    }
    if (output_avstream_video_encode_frame != NULL)
		av_frame_free(&output_avstream_video_encode_frame);
	if (output_avstream_video_frame != NULL)
		av_frame_free(&output_avstream_video_frame);
	av_write_trailer(output_avfmt);
	if (output_avfmt != NULL && !(output_avfmt->oformat->flags & AVFMT_NOFILE))
		avio_closep(&output_avfmt->pb);
	avformat_free_context(output_avfmt);
    close_src();

	return 0;
}

