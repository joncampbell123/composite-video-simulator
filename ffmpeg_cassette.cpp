// recreate hiss, warble, L-R temporal distortions as if recorded to audio cassette

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
#include <string>
#include <vector>

volatile int DIE = 0;

void sigma(int x) {
	if (++DIE >= 20) abort();
}

int		audio_stream_index = 0;

string		input_file;
string		output_file;

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

HiLoComboPass		audio_hilopass;

// preemphsis emuluation
LowpassFilter		audio_linear_preemphasis_pre[2];
LowpassFilter		audio_linear_preemphasis_post[2];

AVFormatContext*	input_avfmt = NULL;
AVStream*		input_avstream_audio = NULL;	// do not free
AVCodecContext*		input_avstream_audio_codec_context = NULL; // do not free
AVFrame*		input_avstream_audio_frame = NULL;

int             input_avstream_audio_resampler_rate = -1;
int             input_avstream_audio_resampler_channels = -1;

struct SwrContext*	input_avstream_audio_resampler = NULL;

AVFormatContext*	output_avfmt = NULL;
AVStream*		output_avstream_audio = NULL;	// do not free
AVCodecContext*		output_avstream_audio_codec_context = NULL; // do not free

double          transcode_start = -1;
double          transcode_end = -1;
double          transcode_dur = -1;

int		output_audio_channels = 2;	// VHS stereo (set to 1 for mono)
int		output_audio_rate = 44100;	// VHS Hi-Fi goes up to 20KHz
double		output_audio_hiss_db = -45; // FIXME: guess
double		output_audio_highpass = 20; // highpass to filter out below 20Hz
double		output_audio_lowpass = 20000; // lowpass to filter out above 20KHz
// NTS:
//   VHS Hi-Fi: 20Hz - 20KHz                  (70dBFS S/N)
//   VHS SP:    100Hz - 10KHz                 (42dBFS S/N)
//   VHS LP:    100Hz - 7KHz (right??)        (42dBFS S/N)
//   VHS EP:    100Hz - 4KHz                  (42dBFS S/N)
bool		emulating_preemphasis = false;		// emulate preemphasis
bool		emulating_deemphasis = false;		// emulate deemphasis

int		output_audio_hiss_level = 0; // out of 10000

enum {
	VHS_SP=0,
	VHS_LP,
	VHS_EP
};

int		output_vhs_tape_speed = VHS_SP;

static inline int clampu8(const int x) {
	if (x > 255)
		return 255;
	else if (x < 0)
		return 0;

	return x;
}

static inline int clips16(const int x) {
	if (x < -32768)
		return -32768;
	else if (x > 32767)
		return 32767;

	return x;
}

static unsigned long long audio_proc_count = 0;

class ConvolutionMap {
public:
    ConvolutionMap() : length(0), map(NULL), multiply(NULL) {
    }
    ~ConvolutionMap() {
        freemap();
    }
    bool allocmap(const size_t len) {
        if (map != NULL && len == length)
            return true;

        freemap();
        length = len;
        if (length == 0)
            return true;

        map = new double[length];
        memset(map,0,sizeof(double) * length);
        multiply = new double[length];
        memset(multiply,0,sizeof(double) * length);
        return true;
    }
public:
    void freemap(void) {
        if (map) delete[] map;
        map = NULL;

        if (multiply) delete[] multiply;
        multiply = NULL;
    }
    double calc(double s) {
        double r = 0;
        size_t i;

        for (i=0;(i+1) < length;i++) map[i] = map[i+1];
        map[i] = s;

        for (i=0;i < length;i++)
            r += map[i] * multiply[i];

        return r;
    }
public:
    size_t                  length;
    double*                 map;
    double*                 multiply;
};

ConvolutionMap          audio_conv[2];
double                  lr_delay = 2;               // part of head tilt, as a consequence of storing stereo left + right on separate halves of the tape
double                  head_tilt = 0.2;            // everyone's a little out of alignment
double                  head_tilt_waver = 0.5;      // and variation in tape speed changes it over time
double                  head_tilt_final = 0;

bool                    mono_downmix = false;

void composite_audio_process(int16_t *audio,unsigned int samples) { // number of channels = output_audio_channels, sample rate = output_audio_rate. audio is interleaved.
	assert(audio_hilopass.audiostate.size() >= output_audio_channels);

    if (audio_conv[0].map == NULL) {
        audio_conv[0].allocmap((int)floor(fabs(head_tilt * 2) + fabs(head_tilt * 3) + 7.5));
        audio_conv[1].allocmap((int)floor(fabs(head_tilt * 2) + fabs(head_tilt * 3) + 7.5));
    }

	for (unsigned int s=0;s < samples;s++,audio += output_audio_channels) {
        {
            double t = (double)audio_proc_count / output_audio_rate;
            head_tilt_final = (head_tilt_waver * sin(t * M_PI * 2 * 1.5)) + head_tilt;
            lr_delay = head_tilt_final * 1.5;

            {
                double mid = lr_delay + ((double)audio_conv[0].length / 2);

                for (size_t i=0;i < audio_conv[0].length;i++) {
                    double d = ((double)i - mid) / (fabs(head_tilt_final) + 1.0);
                    d = 1.0 - fabs(d); // FIXME: sinc would be more appropriate?
                    if (d < 0) d = 0;
                    d /= fabs(head_tilt_final) + 1.0;
                    audio_conv[0].multiply[i] = d;
                }
            }

            {
                double mid = -lr_delay + ((double)audio_conv[1].length / 2);

                for (size_t i=0;i < audio_conv[1].length;i++) {
                    double d = ((double)i - mid) / (fabs(head_tilt_final) + 1.0);
                    d = 1.0 - fabs(d); // FIXME: sinc would be more appropriate?
                    if (d < 0) d = 0;
                    d /= fabs(head_tilt_final) + 1.0;
                    audio_conv[1].multiply[i] = d;
                }
            }
        }

		for (unsigned int c=0;c < output_audio_channels;c++) {
			double s;

			s = (double)audio[c] / 32768;

			/* lowpass filter */
			s = audio_hilopass.audiostate[c].filter(s);

			/* preemphasis */
			if (emulating_preemphasis) {
				for (unsigned int i=0;i < output_audio_channels;i++) {
					s = s + audio_linear_preemphasis_pre[i].highpass(s);
				}
			}

			/* analog limiting (when the signal is too loud) */
			if (s > 1.0)
				s = 1.0;
			else if (s < -1.0)
				s = -1.0;

			/* hiss */
			if (output_audio_hiss_level != 0)
				s += ((double)(((int)((unsigned int)rand() % ((output_audio_hiss_level * 2) + 1))) - output_audio_hiss_level)) / 20000;

            /* convolution */
            s = audio_conv[c].calc(s);

			/* deemphasis */
			if (emulating_deemphasis) {
				for (unsigned int i=0;i < output_audio_channels;i++) {
					s = audio_linear_preemphasis_post[i].lowpass(s);
				}
			}

			audio[c] = clips16(s * 32768);
		}

        if (mono_downmix)
            audio[0] = audio[1] = (audio[0] + audio[1]) / 2;

		audio_proc_count++;
	}
}

static void help(const char *arg0) {
	fprintf(stderr,"%s [options]\n",arg0);
	fprintf(stderr," -i <input file>\n");
	fprintf(stderr," -o <output file>\n");
	fprintf(stderr," -preemphasis <0|1>        Enable preemphasis emulation\n");
	fprintf(stderr," -deemphasis <0|1>         Enable deepmhasis emulation\n");
	fprintf(stderr," -audio-hiss <-120..0>     Audio hiss in decibels (0=100%)\n");
	fprintf(stderr," -a <n>                    Pick the n'th audio stream\n");
	fprintf(stderr," -an                       Don't render any audio stream\n");
    fprintf(stderr," -ss <t>                   Start transcoding from t seconds\n");
    fprintf(stderr," -se <t>                   Stop transcoding at t seconds\n");
    fprintf(stderr," -t <t>                    Transcode only t seconds\n");
    fprintf(stderr," -low <n>                  Lowpass frequency\n");
    fprintf(stderr," -high <n>                 Highpass frequency\n");
    fprintf(stderr," -headalign <n>            Head misalignment (0 = perfectly aligned)\n");
    fprintf(stderr," -headalignwaver <n>       Head misalignment wavering (0 no waver)\n");
    fprintf(stderr," -mono                     Mono playback\n");
    fprintf(stderr," -preset <x>               Preset to use (0, 1, 2...)\n");
	fprintf(stderr,"\n");
	fprintf(stderr," Output file will be up/down converted to 720x480 (NTSC 29.97fps) or 720x576 (PAL 25fps).\n");
	fprintf(stderr," Output will be rendered as interlaced video.\n");
}

static int parse_argv(int argc,char **argv) {
	const char *a;
	int i;

    // default to modern "good" cassette brands.
    output_audio_highpass = 20; // highpass to filter out below 20Hz
    output_audio_lowpass = 20000; // lowpass to filter out above 20KHz
    output_audio_channels = 2;

	for (i=1;i < argc;) {
		a = argv[i++];

		if (*a == '-') {
			do { a++; } while (*a == '-');

			if (!strcmp(a,"h") || !strcmp(a,"help")) {
				help(argv[0]);
				return 1;
			}
            else if (!strcmp(a,"mono")) {
                mono_downmix = 1;
            }
            else if (!strcmp(a,"headalign")) {
                a = argv[i++];
                if (a == NULL) return 1;
                head_tilt = atoi(a);
            }
            else if (!strcmp(a,"headalignwaver")) {
                a = argv[i++];
                if (a == NULL) return 1;
                head_tilt_waver = atoi(a);
            }
            else if (!strcmp(a,"low")) {
                a = argv[i++];
                if (a == NULL) return 1;
                output_audio_lowpass = (int)strtoul(a,NULL,0);
            }
            else if (!strcmp(a,"high")) {
                a = argv[i++];
                if (a == NULL) return 1;
                output_audio_highpass = (int)strtoul(a,NULL,0);
            }
            else if (!strcmp(a,"ss")) {
                transcode_start = atof(argv[i++]);
            }
            else if (!strcmp(a,"se")) {
                transcode_end = atof(argv[i++]);
            }
            else if (!strcmp(a,"t")) {
                transcode_dur = atof(argv[i++]);
            }
			else if (!strcmp(a,"a")) {
				audio_stream_index = atoi(argv[i++]);
			}
			else if (!strcmp(a,"an")) {
				audio_stream_index = -1;
			}
			else if (!strcmp(a,"audio-hiss")) {
				output_audio_hiss_db = atof(argv[i++]);
			}
			else if (!strcmp(a,"preemphasis")) {
				int x = atoi(argv[i++]);
				emulating_preemphasis = (x > 0);
			}
			else if (!strcmp(a,"deemphasis")) {
				int x = atoi(argv[i++]);
				emulating_deemphasis = (x > 0);
			}
			else if (!strcmp(a,"i")) {
				input_file = argv[i++];
			}
			else if (!strcmp(a,"o")) {
				output_file = argv[i++];
			}
            else if (!strcmp(a,"preset")) {
                a = argv[i++];
                if (a == NULL) return 1;

                int preset = atoi(a);

                switch (preset) {
                    case 0: /* old tape, lower grade recording, with some misalignment, somewhat muffled */
                        output_audio_lowpass = 16000;
                        output_audio_highpass = 100;
                        head_tilt_waver = 0.55;
                        head_tilt = 3.5;
                        break;
                    case 1: /* older! (best used with mono) */
                        output_audio_lowpass = 14000;
                        output_audio_highpass = 100;
                        head_tilt_waver = 0.6;
                        head_tilt = 6;
                        break;
                    case 2: /* voice grade crap (best used with mono) */
                        output_audio_lowpass = 10000;
                        output_audio_highpass = 100;
                        head_tilt_waver = 0.5;
                        head_tilt = 3;
                        break;
                    case 3: /* recorded on badly aligned deck */
                        output_audio_lowpass = 16000;
                        output_audio_highpass = 20;
                        head_tilt_waver = 0.75;
                        head_tilt = 10;
                        break;
                    case 4: /* good recording, good deck */
                        output_audio_lowpass = 16000;
                        output_audio_highpass = 20;
                        head_tilt_waver = 0.25;
                        head_tilt = 1.1;
                        break;
                    default:
                        fprintf(stderr,"Unknown preset\n");
                        return 1;
                };
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

    if (transcode_start >= 0 && transcode_end >= 0)
        transcode_dur = transcode_end - transcode_start;
    if (transcode_start < 0)
        transcode_start = 0;
    if (transcode_end < 0 && transcode_dur >= 0)
        transcode_end = transcode_start + transcode_dur;
    if (transcode_start >= 0 && transcode_end >= 0 && transcode_start >= transcode_end) {
        fprintf(stderr,"nothing to transcode\n");
        return 1;
    }

    fprintf(stderr,"Transcoding from %.2f to %.2f\n",
        transcode_start,transcode_end);

	output_audio_hiss_level = dBFS(output_audio_hiss_db) * 5000;

	if (input_file.empty() || output_file.empty()) {
		fprintf(stderr,"You must specify an input and output file (-i and -o).\n");
		return 1;
	}

	return 0;
}

struct AVDelayedFrameInfo {
    AVDelayedFrameInfo() : duration(0) {
    }
    unsigned int        duration;
};

std::map<unsigned long long,AVDelayedFrameInfo> AVDelayed;

uint8_t **audio_dst_data = NULL;
int audio_dst_data_alloc_samples = 0;
int audio_dst_data_linesize = 0;
int audio_dst_data_samples = 0;

bool do_audio_decode_and_render(AVPacket &pkt,unsigned long long &audio_sample) {
    int got_frame = 0;

    if (avcodec_decode_audio4(input_avstream_audio_codec_context,input_avstream_audio_frame,&got_frame,&pkt) >= 0) {
        if (got_frame != 0 && input_avstream_audio_frame->nb_samples != 0) {
            unsigned long long tgt_sample = input_avstream_audio_frame->pts;
            if (tgt_sample == AV_NOPTS_VALUE) tgt_sample = pkt.pts;

            if (tgt_sample == AV_NOPTS_VALUE)
                tgt_sample = audio_sample; // don't want me to guess? give me PTS timestamps then!
            else {
                if ((signed long long)tgt_sample < 0LL) tgt_sample = 0LL;

                // deal with imperfections, prevent them from making an unstable frame rate
                signed long long d = (signed long long)tgt_sample - (signed long long)audio_sample;

                if (llabs(d) < (output_audio_rate/30) && tgt_sample < audio_sample)
                    tgt_sample = audio_sample;
            }

            if (input_avstream_audio_resampler != NULL) {
                if (input_avstream_audio_resampler_rate != input_avstream_audio_codec_context->sample_rate ||
                        input_avstream_audio_resampler_channels != input_avstream_audio_codec_context->channels) {
                    fprintf(stderr,"Audio format changed\n");
                    swr_free(&input_avstream_audio_resampler);
                }
            }

            if (input_avstream_audio_resampler == NULL) {
                input_avstream_audio_resampler = swr_alloc();
                av_opt_set_int(input_avstream_audio_resampler, "in_channel_count", input_avstream_audio_codec_context->channels, 0); // FIXME: FFMPEG should document this!!
                av_opt_set_int(input_avstream_audio_resampler, "out_channel_count", output_avstream_audio_codec_context->channels, 0); // FIXME: FFMPEG should document this!!
                av_opt_set_int(input_avstream_audio_resampler, "in_channel_layout", input_avstream_audio_codec_context->channel_layout, 0);
                av_opt_set_int(input_avstream_audio_resampler, "out_channel_layout", output_avstream_audio_codec_context->channel_layout, 0);
                av_opt_set_int(input_avstream_audio_resampler, "in_sample_rate", input_avstream_audio_codec_context->sample_rate, 0);
                av_opt_set_int(input_avstream_audio_resampler, "out_sample_rate", output_avstream_audio_codec_context->sample_rate, 0);
                av_opt_set_sample_fmt(input_avstream_audio_resampler, "in_sample_fmt", input_avstream_audio_codec_context->sample_fmt, 0);
                av_opt_set_sample_fmt(input_avstream_audio_resampler, "out_sample_fmt", output_avstream_audio_codec_context->sample_fmt, 0);
                if (swr_init(input_avstream_audio_resampler) < 0) {
                    fprintf(stderr,"Failed to init audio resampler\n");
                    swr_free(&input_avstream_audio_resampler);
                    return got_frame;
                }
                input_avstream_audio_resampler_rate = input_avstream_audio_codec_context->sample_rate;
                input_avstream_audio_resampler_channels = input_avstream_audio_codec_context->channels;

                if (audio_dst_data != NULL) {
                    av_freep(&audio_dst_data[0]); // NTS: Why??
                    av_freep(&audio_dst_data);
                }

                audio_dst_data_alloc_samples = 0;
                fprintf(stderr,"Audio resampler init %uHz -> %uHz\n",
                        input_avstream_audio_codec_context->sample_rate,
                        output_avstream_audio_codec_context->sample_rate);
            }

            audio_dst_data_samples = av_rescale_rnd(
                    swr_get_delay(input_avstream_audio_resampler, input_avstream_audio_frame->sample_rate) + input_avstream_audio_frame->nb_samples,
                    output_avstream_audio_codec_context->sample_rate, input_avstream_audio_frame->sample_rate, AV_ROUND_UP);

            if (audio_dst_data == NULL || audio_dst_data_samples > audio_dst_data_alloc_samples) {
                if (audio_dst_data != NULL) {
                    av_freep(&audio_dst_data[0]); // NTS: Why??
                    av_freep(&audio_dst_data);
                }

                audio_dst_data_alloc_samples = 0;
                fprintf(stderr,"Allocating audio buffer %u samples\n",(unsigned int)audio_dst_data_samples);
                if (av_samples_alloc_array_and_samples(&audio_dst_data,&audio_dst_data_linesize,
                            output_avstream_audio_codec_context->channels,audio_dst_data_samples,
                            output_avstream_audio_codec_context->sample_fmt, 0) >= 0) {
                    audio_dst_data_alloc_samples = audio_dst_data_samples;
                }
                else {
                    fprintf(stderr,"Failure to allocate audio buffer\n");
                    audio_dst_data_alloc_samples = 0;
                }
            }

            /* pad-fill */
            while (audio_sample < tgt_sample) {
                unsigned long long out_samples = tgt_sample - audio_sample;

                if (out_samples > output_audio_rate)
                    out_samples = output_audio_rate;

                AVPacket dstpkt;
                av_init_packet(&dstpkt);
                if (av_new_packet(&dstpkt,out_samples * 2 * output_audio_channels) >= 0) { // NTS: Will reset fields too!
                    assert(dstpkt.data != NULL);
                    assert(dstpkt.size >= (out_samples * 2 * output_audio_channels));
                    memset(dstpkt.data,0,out_samples * 2 * output_audio_channels);
                }
                dstpkt.pts = audio_sample;
                dstpkt.dts = audio_sample;
                dstpkt.stream_index = output_avstream_audio->index;
                av_packet_rescale_ts(&dstpkt,output_avstream_audio_codec_context->time_base,output_avstream_audio->time_base);
                if (av_interleaved_write_frame(output_avfmt,&dstpkt) < 0)
                    fprintf(stderr,"Failed to write frame\n");
                av_packet_unref(&dstpkt);

                fprintf(stderr,"Pad fill %llu samples\n",out_samples);
                audio_sample += out_samples;
            }

            if (audio_dst_data != NULL && tgt_sample >= audio_sample) {
                int out_samples;

                if ((out_samples=swr_convert(input_avstream_audio_resampler,audio_dst_data,audio_dst_data_samples,
                                (const uint8_t**)input_avstream_audio_frame->data,input_avstream_audio_frame->nb_samples)) > 0) {
                    // PROCESS THE AUDIO. At this point by design the code can assume S16LE (16-bit PCM interleaved)
                    composite_audio_process((int16_t*)audio_dst_data[0],out_samples);
                    // write it out. TODO: At some point, support conversion to whatever the codec needs and then convert to it.
                    // that way we can render directly to MP4 our VHS emulation.
                    AVPacket dstpkt;
                    av_init_packet(&dstpkt);
                    if (av_new_packet(&dstpkt,out_samples * 2 * output_audio_channels) >= 0) { // NTS: Will reset fields too!
                        assert(dstpkt.data != NULL);
                        assert(dstpkt.size >= (out_samples * 2 * output_audio_channels));
                        memcpy(dstpkt.data,audio_dst_data[0],out_samples * 2 * output_audio_channels);
                    }
                    dstpkt.pts = audio_sample;
                    dstpkt.dts = audio_sample;
                    dstpkt.stream_index = output_avstream_audio->index;
                    av_packet_rescale_ts(&dstpkt,output_avstream_audio_codec_context->time_base,output_avstream_audio->time_base);
                    if (av_interleaved_write_frame(output_avfmt,&dstpkt) < 0)
                        fprintf(stderr,"Failed to write frame\n");
                    av_packet_unref(&dstpkt);

                    audio_sample += out_samples;
                }
                else if (out_samples < 0) {
                    fprintf(stderr,"Failed to resample audio\n");
                }
            }
        }
    }
    else {
        fprintf(stderr,"No audio decoded\n");
    }

    return (got_frame != 0);
}

int main(int argc,char **argv) {
	if (parse_argv(argc,argv))
		return 1;

	av_register_all();
	avformat_network_init();
	avcodec_register_all();

	assert(input_avfmt == NULL);
	if (avformat_open_input(&input_avfmt,input_file.c_str(),NULL,NULL) < 0) {
		fprintf(stderr,"Failed to open input file\n");
		return 1;
	}

	if (avformat_find_stream_info(input_avfmt,NULL) < 0)
		fprintf(stderr,"WARNING: Did not find stream info on input\n");

	/* scan streams for one video, one audio */
	{
		size_t i;
		AVStream *is;
		int ac=0,vc=0;
		AVCodecContext *isctx;

		fprintf(stderr,"Input format: %u streams found\n",input_avfmt->nb_streams);
		for (i=0;i < (size_t)input_avfmt->nb_streams;i++) {
			is = input_avfmt->streams[i];
			if (is == NULL) continue;

			isctx = is->codec;
			if (isctx == NULL) continue;

			if (isctx->codec_type == AVMEDIA_TYPE_AUDIO) {
				if (input_avstream_audio == NULL && ac == audio_stream_index) {
					if (avcodec_open2(isctx,avcodec_find_decoder(isctx->codec_id),NULL) >= 0) {
						input_avstream_audio = is;
						input_avstream_audio_codec_context = isctx;
						fprintf(stderr,"Found audio stream idx=%zu %u-channel %uHz\n",
							i,
							input_avstream_audio_codec_context->channels,
							input_avstream_audio_codec_context->sample_rate);
					}
					else {
						fprintf(stderr,"Found audio stream but not able to decode\n");
					}
				}

				ac++;
			}
		}

		if (input_avstream_audio == NULL) {
			fprintf(stderr,"Neither video nor audio found\n");
			return 1;
		}
	}

	assert(output_avfmt == NULL);
	if (avformat_alloc_output_context2(&output_avfmt,NULL,NULL,output_file.c_str()) < 0) {
		fprintf(stderr,"Failed to open output file\n");
		return 1;
	}

	if (input_avstream_audio != NULL) {
		output_avstream_audio = avformat_new_stream(output_avfmt, NULL);
		if (output_avstream_audio == NULL) {
			fprintf(stderr,"Unable to create output audio stream\n");
			return 1;
		}

		output_avstream_audio_codec_context = output_avstream_audio->codec;
		if (output_avstream_audio_codec_context == NULL) {
			fprintf(stderr,"Output stream audio no codec context?\n");
			return 1;
		}

		if (output_audio_channels == 2)
			output_avstream_audio_codec_context->channel_layout = AV_CH_LAYOUT_STEREO;
		else
			output_avstream_audio_codec_context->channel_layout = AV_CH_LAYOUT_MONO;

		output_avstream_audio_codec_context->sample_rate = output_audio_rate;
		output_avstream_audio_codec_context->channels = output_audio_channels;
		output_avstream_audio_codec_context->sample_fmt = AV_SAMPLE_FMT_S16;
		output_avstream_audio_codec_context->time_base = (AVRational){1, output_audio_rate};
		output_avstream_audio->time_base = output_avstream_audio_codec_context->time_base;

		if (output_avfmt->oformat->flags & AVFMT_GLOBALHEADER)
			output_avstream_audio_codec_context->flags |= AV_CODEC_FLAG_GLOBAL_HEADER;

		if (avcodec_open2(output_avstream_audio_codec_context,avcodec_find_encoder(AV_CODEC_ID_PCM_S16LE),NULL) < 0) {
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

	/* prepare audio filtering */
	audio_hilopass.setChannels(output_audio_channels);
	audio_hilopass.setRate(output_audio_rate);
	audio_hilopass.setCutoff(output_audio_lowpass,output_audio_highpass); // hey, our filters aren't perfect
	audio_hilopass.setPasses(6);
	audio_hilopass.init();

	if (emulating_preemphasis) {
        for (unsigned int i=0;i < output_audio_channels;i++) {
            audio_linear_preemphasis_pre[i].setFilter(output_audio_rate,4000/*FIXME: Guess! Also let user set this.*/);
        }
    }
	if (emulating_deemphasis) {
        for (unsigned int i=0;i < output_audio_channels;i++) {
            audio_linear_preemphasis_post[i].setFilter(output_audio_rate,4000/*FIXME: Guess! Also let user set this.*/);
        }
    }

	/* prepare audio decoding */
	if (input_avstream_audio != NULL) {
		input_avstream_audio_frame = av_frame_alloc();
		if (input_avstream_audio_frame == NULL) {
			fprintf(stderr,"Failed to alloc audio frame\n");
			return 1;
		}
	}

	// PARSE
	{
        unsigned long long av_frame_counter = 0;
        unsigned long long audio_sample = 0;
		unsigned long long video_field = 0;
        double adj_time = 0;
        int got_frame = 0;
        double t,pt = -1;
		AVPacket pkt;

		av_init_packet(&pkt);
		while (av_read_frame(input_avfmt,&pkt) >= 0) {
			if (DIE != 0) break;

            // ugh... this can happen if the source is an AVI file
            if (pkt.pts == AV_NOPTS_VALUE) pkt.pts = pkt.dts;

            /* track time and keep things monotonic for our code */
            if (pkt.stream_index < input_avfmt->nb_streams) {
                if (pkt.pts != AV_NOPTS_VALUE) {
                    t = pkt.pts * av_q2d(input_avfmt->streams[pkt.stream_index]->time_base);
                    if (transcode_end >= 0 && t >= transcode_end)
                        break;

                    if (t < transcode_start) {
                        av_packet_unref(&pkt);
                        av_init_packet(&pkt);
                        continue;
                    }

                    if (pt < 0)
                        adj_time = -t;
                    else if ((t+1.5) < pt) { // time code jumps backwards (1.5 is safe for DVD timecode resets)
                        adj_time += pt - t;
                        fprintf(stderr,"Time code jump backwards %.6f->%.6f. adj_time=%.6f\n",pt,t,adj_time);
                    }
                    else if (t > (pt+5)) { // time code jumps forwards
                        adj_time += pt - t;
                        fprintf(stderr,"Time code jump forwards %.6f->%.6f. adj_time=%.6f\n",pt,t,adj_time);
                    }

                    pt = t;
                }

                if (pt < 0) {
                    av_packet_unref(&pkt);
                    av_init_packet(&pkt);
                    continue;
                }

                if (pkt.pts != AV_NOPTS_VALUE) {
                    pkt.pts += (adj_time * input_avfmt->streams[pkt.stream_index]->time_base.den) /
                        input_avfmt->streams[pkt.stream_index]->time_base.num;
                }

                if (pkt.dts != AV_NOPTS_VALUE) {
                    pkt.dts += (adj_time * input_avfmt->streams[pkt.stream_index]->time_base.den) /
                        input_avfmt->streams[pkt.stream_index]->time_base.num;
                }
            }

			if (input_avstream_audio != NULL && pkt.stream_index == input_avstream_audio->index) {
				av_packet_rescale_ts(&pkt,input_avstream_audio->time_base,output_avstream_audio->time_base);
                do_audio_decode_and_render(/*&*/pkt,/*&*/audio_sample);
			}

			av_packet_unref(&pkt);
			av_init_packet(&pkt);
		}

        av_packet_unref(&pkt);
        av_init_packet(&pkt);

        /* the encoder usually has a delay.
         * we need the encoder to flush those frames out. */
        {
            do {
                if (DIE != 0) break;

                pkt.size = 0;
                pkt.data = NULL;
                if (input_avstream_audio != NULL)
                    got_frame |= do_audio_decode_and_render(/*&*/pkt,/*&*/audio_sample) ? 1 : 0;
            } while (got_frame);
        }

		if (audio_dst_data != NULL) {
			av_freep(&audio_dst_data[0]); // NTS: Why??
			av_freep(&audio_dst_data);
		}
	}

	if (input_avstream_audio_frame != NULL)
		av_frame_free(&input_avstream_audio_frame);
	audio_hilopass.clear();
	av_write_trailer(output_avfmt);
	if (input_avstream_audio_resampler != NULL)
		swr_free(&input_avstream_audio_resampler);
    if (input_avstream_audio_codec_context != NULL) {
        avcodec_close(input_avstream_audio_codec_context);
        input_avstream_audio_codec_context = NULL;
        input_avstream_audio = NULL;
    }
	if (output_avfmt != NULL && !(output_avfmt->oformat->flags & AVFMT_NOFILE))
		avio_closep(&output_avfmt->pb);
	avformat_free_context(output_avfmt);
	avformat_close_input(&input_avfmt);
	return 0;
}

