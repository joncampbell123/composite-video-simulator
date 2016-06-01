// what to do next:
//
// VHS Hi-Fi "buzz"
//
// Video decoding (scaling to 720x480 or 720x576, frame/field rate conversion to interlaced video, encode)
// When done, all video will be rendered as if 720x480 29.97fps interlaced or 720x576 25fps interlaced.
//
// Fake NTSC color subcarrier generation (from chroma planes into luma).
//
// Rendering of color back out from luma plane fake subcarrier (NTSC artifacts and all)
//
// If emulating PAL, fake subcarrier needs to alternate color phase every other field scanline.
//
// Luma smearing/sharpen filtering to emulate analog degredation and "sharpening" in analog equipment to compensate.
//
// Additional luma smearing/sharpening as a side effect of preemphasis and modulation onto VHS tape.
//
// Video noise.
//
// VHS audio linear track video+audio crosstalk (ever notice you can almost hear the video signal in the audio?) and hiss.
//
// Luma/chroma instability and noise.
//
// Fake macrovision "darkening" of the top of the picture, and
//
// Fake macrovision "bend" at the top of the picture.
//
// Analog video "banding", slight but noticeable bright/dark bands that vary according to video content.
//
// Minor VHS noise (occasional white specks)
//
// Tracking noise (the staticky band at the bottom of the screen).
//
// After processing audio, allow encode through non-PCM codec and write to output.
//
// After processing video, allow encode through non-uncompressed codec and write to output.
//
// FFMPEG-like options to specify which audio track to use (so DVD rips don't accidentally use the non-English tracks).

#define __STDC_CONSTANT_MACROS

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

#include <string>
#include <vector>

volatile int DIE = 0;

void sigma(int x) {
	if (++DIE >= 20) abort();
}

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
	void resetFilter() {
		prev = 0;
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
AVStream*		input_avstream_video = NULL;	// do not free
AVCodecContext*		input_avstream_video_codec_context = NULL; // do not free
AVFrame*		input_avstream_audio_frame = NULL;
AVFrame*		input_avstream_video_frame = NULL;

struct SwrContext*	input_avstream_audio_resampler = NULL;
struct SwsContext*	input_avstream_video_resampler = NULL;

AVFormatContext*	output_avfmt = NULL;
AVStream*		output_avstream_audio = NULL;	// do not free
AVCodecContext*		output_avstream_audio_codec_context = NULL; // do not free
AVStream*		output_avstream_video = NULL;	// do not free
AVCodecContext*		output_avstream_video_codec_context = NULL; // do not free
AVFrame*		output_avstream_video_input_frame = NULL;
AVFrame*		output_avstream_video_frame = NULL;

int		video_chroma_noise = 0;
int		video_noise = 2;
int		subcarrier_amplitude = 80;
AVRational	output_field_rate = { 60000, 1001 };	// NTSC 60Hz default
int		output_width = 720;
int		output_height = 480;
bool		output_ntsc = true;	// NTSC color subcarrier emulation
bool		output_pal = false;	// PAL color subcarrier emulation
int		output_audio_channels = 2;	// VHS stereo (set to 1 for mono)
int		output_audio_rate = 44100;	// VHS Hi-Fi goes up to 20KHz
double		output_audio_hiss_db = -72;
double		output_audio_linear_buzz = -42;	// how loud the "buzz" is audible in dBFS (S/N). Ever notice on old VHS tapes (prior to Hi-Fi) you can almost hear the video signal sync pulses in the audio?
double		output_audio_highpass = 20; // highpass to filter out below 20Hz
double		output_audio_lowpass = 20000; // lowpass to filter out above 20KHz
// NTS:
//   VHS Hi-Fi: 20Hz - 20KHz                  (70dBFS S/N)
//   VHS SP:    100Hz - 10KHz                 (42dBFS S/N)
//   VHS LP:    100Hz - 7KHz (right??)        (42dBFS S/N)
//   VHS EP:    100Hz - 4KHz                  (42dBFS S/N)
bool		output_vhs_hifi = true;
bool		output_vhs_linear_stereo = false; // not common
bool		output_vhs_linear_audio = false; // if true (non Hi-Fi) then we emulate hiss and noise of linear VHS tracks including the video sync pulses audible in the audio.
bool		emulating_vhs = false;
bool		emulating_preemphasis = true;		// emulate preemphasis
bool		emulating_deemphasis = true;		// emulate deemphasis
bool		nocolor_subcarrier = false;		// if set, emulate subcarrier but do not decode back to color (debug)
bool		nocolor_subcarrier_after_yc_sep = false;// if set, separate luma-chroma but do not decode back to color (debug)

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

void composite_audio_process(int16_t *audio,unsigned int samples) { // number of channels = output_audio_channels, sample rate = output_audio_rate. audio is interleaved.
	assert(audio_hilopass.audiostate.size() >= output_audio_channels);

	for (unsigned int s=0;s < samples;s++,audio += output_audio_channels) {
		for (unsigned int c=0;c < output_audio_channels;c++) {
			double s;

			s = (double)audio[c] / 32768;
			s = audio_hilopass.audiostate[c].filter(s);

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

			if (emulating_deemphasis) {
				for (unsigned int i=0;i < output_audio_channels;i++) {
					s = audio_linear_preemphasis_post[i].lowpass(s);
				}
			}

			audio[c] = clips16(s * 32768);
		}
	}
}

void composite_video_process(AVFrame *dst,unsigned int field,unsigned long long fieldno) {
	unsigned int x,y;

	{ /* lowpass the chroma more. composite video does not allocate as much bandwidth to color as luma. */
		for (unsigned int p=1;p <= 2;p++) {
			for (y=field;y < dst->height;y += 2) {
				unsigned char *s = dst->data[p] + (y * dst->linesize[p]);
				unsigned char pppppppppc = 128,ppppppppc = 128,pppppppc = 128,ppppppc = 128,pppppc = 128,ppppc = 128,pppc,ppc,pc,c;

				ppppc = s[0];
				pppc = s[1];
				ppc = s[2];
				pc = s[3];
				for (x=0;x < (dst->width/2)/*4:2:2*/;x++) {
					c = s[x+4];
					s[x] = (pppppppppc + ppppppppc + pppppppc + pppppppc + ppppppc + ppppppc +
						pppppc + pppppc + ppppc + ppppc + pppc + pppc + ppc + ppc + pc + c + 8) >> 4;
					pppppppppc = ppppppppc;
					ppppppppc = pppppppc;
					pppppppc = ppppppc;
					ppppppc = pppppc;
					pppppc = ppppc;
					ppppc = pppc;
					pppc = ppc;
					ppc = pc;
					pc = c;
				}
			}
		}
	}

	/* render the chroma into the luma as a fake NTSC color subcarrier */
	{
		for (y=field;y < dst->height;y += 2) {
			unsigned char *Y = dst->data[0] + (y * dst->linesize[0]);
			unsigned char *U = dst->data[1] + (y * dst->linesize[1]);
			unsigned char *V = dst->data[2] + (y * dst->linesize[2]);
			unsigned int xc = dst->width;

			if (((fieldno^y)>>1)&1) {
				Y += 2;
				xc -= 2;
			}

			/* remember: this code assumes 4:2:2 */
			/* NTS: the subcarrier is two sine waves superimposed on top of each other, 90 degrees apart */
			for (x=0;x < xc;x += 4,Y += 4,U += 2,V += 2) {
				Y[0] = clampu8(Y[0] + ((((int)U[0] - 128) * subcarrier_amplitude) / 50));
				Y[1] = clampu8(Y[1] + ((((int)V[0] - 128) * subcarrier_amplitude) / 50));
				Y[2] = clampu8(Y[2] - ((((int)U[1] - 128) * subcarrier_amplitude) / 50));
				Y[3] = clampu8(Y[3] - ((((int)V[1] - 128) * subcarrier_amplitude) / 50));

				if (nocolor_subcarrier)
					U[0] = V[0] = U[1] = V[1] = 128;
			}
		}
	}

	/* add video noise */
	if (video_noise != 0) {
		int noise = 0,noise_mod = (video_noise * 255) / 100;

		for (y=field;y < dst->height;y += 2) {
			unsigned char *Y = dst->data[0] + (y * dst->linesize[0]);

			for (x=0;x < dst->width;x++) {
				Y[x] = clampu8(Y[x] + noise);
				noise += ((int)((unsigned int)rand() % ((video_noise*2)+1))) - video_noise;
				noise /= 2;
			}
		}
	}

	/* filter subcarrier back out, use result to emulate NTSC luma-chroma artifacts */
	if (!nocolor_subcarrier) {
		unsigned char chroma[dst->width]; // WARNING: This is more GCC-specific C++ than normal

		for (y=field;y < dst->height;y += 2) {
			unsigned char *Y = dst->data[0] + (y * dst->linesize[0]);
			unsigned char *U = dst->data[1] + (y * dst->linesize[1]);
			unsigned char *V = dst->data[2] + (y * dst->linesize[2]);
			unsigned char delay[4] = {16,16,16,16};
			unsigned int sum = 16 * (4 - 2);
			unsigned char c;

			// precharge by 2 pixels to center box blur
			delay[2] = Y[0]; sum += delay[2];
			delay[3] = Y[1]; sum += delay[3];
			for (x=0;x < dst->width;x++) {
				c = Y[x+2];
				sum -= delay[0];
				for (unsigned int j=0;j < (4-1);j++) delay[j] = delay[j+1];
				delay[3] = c;
				sum += delay[3];
				Y[x] = sum / 4;
				chroma[x] = clampu8(c + 128 - Y[x]);

				if (nocolor_subcarrier_after_yc_sep) {
					// debug option to SHOW what we got after filtering
					Y[x] = chroma[x];
					U[x/2] = V[x/2] = 128;
				}
			}

			if (!nocolor_subcarrier_after_yc_sep) {
				unsigned int xi = 0;
				if (((fieldno^y)>>1)&1) xi = 2;

				for (x=xi;x < dst->width;x += 4) { // flip the part of the sine wave that would correspond to negative U and V values
					chroma[x+2] = 255 - chroma[x+2];
					chroma[x+3] = 255 - chroma[x+3];
				}

				for (x=0;x < dst->width;x++) {
					chroma[x] = clampu8(((((int)chroma[x] - 128) * 50) / subcarrier_amplitude) + 128);
				}

				/* decode the color right back out from the subcarrier we generated */
				for (x=0;x < (dst->width/2);x++) {
					U[x] = 255 - chroma[(x*2)+0];
					V[x] = 255 - chroma[(x*2)+1];
				}
			}
		}
	}

	/* add video noise */
	if (video_chroma_noise != 0) {
		int noiseU = 0,noiseV = 0,noise_mod = (video_chroma_noise * 255) / 100;

		for (y=field;y < dst->height;y += 2) {
			unsigned char *U = dst->data[1] + (y * dst->linesize[1]);
			unsigned char *V = dst->data[2] + (y * dst->linesize[2]);

			for (x=0;x < (dst->width/2);x++) {
				U[x] = clampu8(U[x] + noiseU);
				V[x] = clampu8(V[x] + noiseV);
				noiseU += ((int)((unsigned int)rand() % ((video_chroma_noise*2)+1))) - video_chroma_noise;
				noiseU /= 2;
				noiseV += ((int)((unsigned int)rand() % ((video_chroma_noise*2)+1))) - video_chroma_noise;
				noiseV /= 2;
			}
		}
	}

	// NTS: At this point, the video best resembles what you'd get from a typical DVD player's composite video output.
	//      Slightly blurry, some color artifacts, and edges will have that "buzz" effect, but still a good picture.
}

void render_field(AVFrame *dst,AVFrame *src,unsigned int field,unsigned long long field_number) {
	unsigned int y,sy,sy2,syf;

	// NTS: dst is 4:2:2 of output_width x output_height
	//      src is 4:2:2 of output_width x source frame height
	//
	//      the reason we do that is so that swscale handles horizontal scaling, then we handle
	//      vertical scaling in the way we need to in order to render interlaced video properly.
	//
	//      this code renders only one field or the other at a time.
	for (y=field;y < dst->height;y += 2) {
		sy = (y * 0x100 * src->height) / dst->height;
		syf = sy & 0xFF;
		sy >>= 8;

		if (src->interlaced_frame) {
			unsigned int which_field = src->top_field_first ? 0/*top*/ : 1/*bottom*/;
			unsigned int src_pts = src->pkt_pts;
			unsigned long long pts_delta = field_number - src_pts;

			if (pts_delta >= ((unsigned long long)input_avstream_video_codec_context->ticks_per_frame / 2ULL))
				which_field ^= 1;

			if (which_field == 0) { // make it even. do not interpolate if first even line of the pair.
				sy++; // but shift up the frame 1 line
				if (!(sy & 1U)) syf = 0;
				else sy--;
			}
			else {
				if (!(sy & 1U)) { // make it odd. do not interpolate if first odd line of the pair.
					syf = 0;
					sy++;
				}
			}

			if (sy >= (src->height - 2)) {
				sy = src->height - 2;
				syf = 0;
			}
			sy2 = sy + 2;
		}
		else {
			if (sy >= (src->height - 1)) {
				sy = src->height - 1;
				syf = 0;
			}
			sy2 = sy + 1;
		}

		if (syf == 0) {
			for (unsigned int p=0;p < 3;p++) {
				unsigned char *s = src->data[p] + (src->linesize[p] * sy);
				unsigned char *d = dst->data[p] + (dst->linesize[p] * y);
				memcpy(d,s,src->linesize[p]);
			}
		}
		else {
			for (unsigned int p=0;p < 3;p++) {
				unsigned char *s1 = src->data[p] + (src->linesize[p] * sy);
				unsigned char *s2 = src->data[p] + (src->linesize[p] * sy);
				unsigned char *d = dst->data[p] + (dst->linesize[p] * y);

				for (unsigned int x=0;x < src->linesize[p];x++)
					d[x] = s1[x] + ((unsigned char)((((int)s2[x] - (int)s1[x]) * (int)syf) >> (int)8));
			}
		}
	}
}

void output_frame(AVFrame *frame,unsigned long long field_number) {
	int gotit = 0;
	AVPacket pkt;

	av_init_packet(&pkt);
	if (av_new_packet(&pkt,50000000/8) < 0) {
		fprintf(stderr,"Failed to alloc vid packet\n");
		return;
	}

	frame->interlaced_frame = 1;
	frame->top_field_first = 0; // NTSC is bottom field first
	frame->pts = field_number / 2ULL;
	pkt.pts = field_number / 2ULL;
	pkt.dts = field_number / 2ULL;
	frame->key_frame = (field_number % (15ULL * 2ULL)) == 0 ? 1 : 0;

	fprintf(stderr,"Output field %llu\n",field_number);
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

void preset_PAL() {
	output_field_rate.num = 50;
	output_field_rate.den = 1;
	output_height = 576;
	output_width = 720;
	output_pal = true;
	output_ntsc = false;
}

void preset_NTSC() {
	output_field_rate.num = 60000;
	output_field_rate.den = 1001;
	output_height = 480;
	output_width = 720;
	output_pal = false;
	output_ntsc = true;
}
	
static void help(const char *arg0) {
	fprintf(stderr,"%s [options]\n",arg0);
	fprintf(stderr," -i <input file>\n");
	fprintf(stderr," -o <output file>\n");
	fprintf(stderr," -tvstd <pal|ntsc>\n");
	fprintf(stderr," -vhs                      Emulation of VHS artifacts\n");
	fprintf(stderr," -vhs-hifi <0|1>           (default on)\n");
	fprintf(stderr," -vhs-speed <ep|lp|sp>     (default sp)\n");
	fprintf(stderr," -preemphasis <0|1>        Enable preemphasis emulation\n");
	fprintf(stderr," -deemphasis <0|1>         Enable deepmhasis emulation\n");
	fprintf(stderr," -nocolor-subcarrier       Emulate color subcarrier but do not decode back (debug)\n");
	fprintf(stderr," -nocolor-subcarrier-after-yc-sep Emulate Y/C subcarrier separation but do not decode back (debug)\n");
	fprintf(stderr," -subcarrier-amp <0...100> Subcarrier amplitude (0 to 100 percent of luma)\n");
	fprintf(stderr," -noise <0..100>           Noise amplitude\n");
	fprintf(stderr," -chroma-noise <0..100>    Chroma noise amplitude\n");
	fprintf(stderr," -audio-hiss <-120..0>     Audio hiss in decibels (0=100%)\n");
	fprintf(stderr,"\n");
	fprintf(stderr," Output file will be up/down converted to 720x480 (NTSC 29.97fps) or 720x576 (PAL 25fps).\n");
	fprintf(stderr," Output will be rendered as interlaced video.\n");
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
			else if (!strcmp(a,"audio-hiss")) {
				output_audio_hiss_db = atof(argv[i++]);
			}
			else if (!strcmp(a,"chroma-noise")) {
				int x = atoi(argv[i++]);
				video_chroma_noise = x;
			}
			else if (!strcmp(a,"noise")) {
				int x = atoi(argv[i++]);
				video_noise = x;
			}
			else if (!strcmp(a,"subcarrier-amp")) {
				int x = atoi(argv[i++]);
				subcarrier_amplitude = x;
			}
			else if (!strcmp(a,"nocolor-subcarrier")) {
				nocolor_subcarrier = true;
			}
			else if (!strcmp(a,"nocolor-subcarrier-after-yc-sep")) {
				nocolor_subcarrier_after_yc_sep = true;
			}
			else if (!strcmp(a,"vhs")) {
				emulating_vhs = true;
				emulating_preemphasis = false; // no preemphasis by default
				emulating_deemphasis = false; // no preemphasis by default
				output_audio_hiss_db = -70;
				video_chroma_noise = 6;
				video_noise = 4; // VHS is a bit noisy
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
			else if (!strcmp(a,"vhs-speed")) {
				a = argv[i++];

				emulating_vhs = true;			// implies -vhs
				if (!strcmp(a,"ep")) {
					output_vhs_tape_speed = VHS_EP;
					video_chroma_noise = 8;
					video_noise = 6;
				}
				else if (!strcmp(a,"lp")) {
					output_vhs_tape_speed = VHS_LP;
					video_chroma_noise = 7;
					video_noise = 5;
				}
				else if (!strcmp(a,"sp")) {
					output_vhs_tape_speed = VHS_SP;
					video_chroma_noise = 6;
					video_noise = 4;
				}
				else {
					fprintf(stderr,"Unknown vhs tape speed '%s'\n",a);
					return 1;
				}
			}
			else if (!strcmp(a,"vhs-hifi")) {
				int x = atoi(argv[i++]);
				output_vhs_hifi = (x > 0);
				output_vhs_linear_audio = !output_vhs_hifi;
				emulating_vhs = true;			// implies -vhs
				if (output_vhs_hifi) {
					emulating_preemphasis = true;
					emulating_deemphasis = true;
					output_audio_hiss_db = -70;
				}
				else {
					output_audio_hiss_db = -42;
				}
			}
			else if (!strcmp(a,"tvstd")) {
				a = argv[i++];

				if (!strcmp(a,"pal")) {
					preset_PAL();
				}
				else if (!strcmp(a,"ntsc")) {
					preset_NTSC();
				}
				else {
					fprintf(stderr,"Unknown tv std '%s'\n",a);
					return 1;
				}
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

	if (emulating_vhs) {
		if (output_vhs_hifi) {
			output_audio_highpass = 20; // highpass to filter out below 20Hz
			output_audio_lowpass = 20000; // lowpass to filter out above 20KHz
			output_audio_channels = 2;
		}
		else if (output_vhs_linear_audio) {
			switch (output_vhs_tape_speed) {
				case VHS_SP:
					output_audio_highpass = 100; // highpass to filter out below 100Hz
					output_audio_lowpass = 10000; // lowpass to filter out above 10KHz
					break;
				case VHS_LP:
					output_audio_highpass = 100; // highpass to filter out below 100Hz
					output_audio_lowpass = 7000; // lowpass to filter out above 7KHz
					break;
				case VHS_EP:
					output_audio_highpass = 100; // highpass to filter out below 100Hz
					output_audio_lowpass = 4000; // lowpass to filter out above 4KHz
					break;
			}

			if (!output_vhs_linear_stereo)
				output_audio_channels = 1;
			else
				output_audio_channels = 2;
		}
	}
	else {
		// not emulating VHS
		output_audio_highpass = 20; // highpass to filter out below 20Hz
		output_audio_lowpass = 20000; // lowpass to filter out above 20KHz
		output_audio_channels = 2;
	}

	output_audio_hiss_level = dBFS(output_audio_hiss_db) * 10000;

	if (input_file.empty() || output_file.empty()) {
		fprintf(stderr,"You must specify an input and output file (-i and -o).\n");
		return 1;
	}

	return 0;
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
		AVCodecContext *isctx;

		fprintf(stderr,"Input format: %u streams found\n",input_avfmt->nb_streams);
		for (i=0;i < (size_t)input_avfmt->nb_streams;i++) {
			is = input_avfmt->streams[i];
			if (is == NULL) continue;

			isctx = is->codec;
			if (isctx == NULL) continue;

			if (isctx->codec_type == AVMEDIA_TYPE_AUDIO) {
				if (input_avstream_audio == NULL) {
					if (avcodec_open2(isctx,avcodec_find_decoder(isctx->codec_id),NULL) >= 0) {
						input_avstream_audio = is;
						input_avstream_audio_codec_context = isctx;
						fprintf(stderr,"Found audio stream idx=%zu\n",i);
					}
					else {
						fprintf(stderr,"Found audio stream but not able to decode\n");
					}
				}
			}
			else if (isctx->codec_type == AVMEDIA_TYPE_VIDEO) {
				if (avcodec_open2(isctx,avcodec_find_decoder(isctx->codec_id),NULL) >= 0) {
					input_avstream_video = is;
					input_avstream_video_codec_context = isctx;
					fprintf(stderr,"Found video stream idx=%zu\n",i);
				}
				else {
					fprintf(stderr,"Found video stream but not able to decode\n");
				}
			}
		}

		if (input_avstream_video == NULL && input_avstream_audio == NULL) {
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

		input_avstream_audio_resampler = swr_alloc();
		av_opt_set_int(input_avstream_audio_resampler, "in_channel_layout", input_avstream_audio_codec_context->channel_layout, 0);
		av_opt_set_int(input_avstream_audio_resampler, "out_channel_layout", output_avstream_audio_codec_context->channel_layout, 0);
		av_opt_set_int(input_avstream_audio_resampler, "in_sample_rate", input_avstream_audio_codec_context->sample_rate, 0);
		av_opt_set_int(input_avstream_audio_resampler, "out_sample_rate", output_avstream_audio_codec_context->sample_rate, 0);
		av_opt_set_sample_fmt(input_avstream_audio_resampler, "in_sample_fmt", input_avstream_audio_codec_context->sample_fmt, 0);
		av_opt_set_sample_fmt(input_avstream_audio_resampler, "out_sample_fmt", output_avstream_audio_codec_context->sample_fmt, 0);
		if (swr_init(input_avstream_audio_resampler) < 0) {
			fprintf(stderr,"Failed to init audio resampler\n");
			return 1;
		}

		fprintf(stderr,"Audio resampler init %uHz -> %uHz\n",
			input_avstream_audio_codec_context->sample_rate,
			output_avstream_audio_codec_context->sample_rate);
	}

	if (input_avstream_video != NULL) {
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
		output_avstream_video_codec_context->pix_fmt = AV_PIX_FMT_YUV422P;
		output_avstream_video_codec_context->time_base = (AVRational){output_field_rate.den, (output_field_rate.num/2)};
		output_avstream_video_codec_context->gop_size = 15;
		output_avstream_video_codec_context->max_b_frames = 0;
		output_avstream_video->time_base = output_avstream_video_codec_context->time_base;

		if (output_avfmt->oformat->flags & AVFMT_GLOBALHEADER)
			output_avstream_video_codec_context->flags |= AV_CODEC_FLAG_GLOBAL_HEADER;

		// our output is interlaced
		output_avstream_video_codec_context->flags |= CODEC_FLAG_INTERLACED_DCT;

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

	/* prepare audio filtering */
	audio_hilopass.setChannels(output_audio_channels);
	audio_hilopass.setRate(output_audio_rate);
	audio_hilopass.setCutoff(output_audio_lowpass,output_audio_highpass); // hey, our filters aren't perfect
	audio_hilopass.setPasses(6);
	audio_hilopass.init();

	// TODO: VHS Hi-Fi is also documented to use 2:1 companding when recording, which we do not yet emulate

	if (emulating_preemphasis) {
		if (output_vhs_hifi) {
			for (unsigned int i=0;i < output_audio_channels;i++) {
				audio_linear_preemphasis_pre[i].setFilter(output_audio_rate,16000/*FIXME: Guess! Also let user set this.*/);
			}
		}
		else {
			for (unsigned int i=0;i < output_audio_channels;i++) {
				audio_linear_preemphasis_pre[i].setFilter(output_audio_rate,8000/*FIXME: Guess! Also let user set this.*/);
			}
		}
	}
	if (emulating_deemphasis) {
		if (output_vhs_hifi) {
			for (unsigned int i=0;i < output_audio_channels;i++) {
				audio_linear_preemphasis_post[i].setFilter(output_audio_rate,16000/*FIXME: Guess! Also let user set this.*/);
			}
		}
		else {
			for (unsigned int i=0;i < output_audio_channels;i++) {
				audio_linear_preemphasis_post[i].setFilter(output_audio_rate,8000/*FIXME: Guess! Also let user set this.*/);
			}
		}
	}

	/* prepare audio decoding */
	input_avstream_audio_frame = av_frame_alloc();
	if (input_avstream_audio_frame == NULL) {
		fprintf(stderr,"Failed to alloc audio frame\n");
		return 1;
	}

	/* prepare video decoding */
	input_avstream_video_frame = av_frame_alloc();
	if (input_avstream_video_frame == NULL) {
		fprintf(stderr,"Failed to alloc video frame\n");
		return 1;
	}

	/* prepare video encoding */
	output_avstream_video_frame = av_frame_alloc();
	if (output_avstream_video_frame == NULL) {
		fprintf(stderr,"Failed to alloc video frame\n");
		return 1;
	}
	av_frame_set_colorspace(output_avstream_video_frame,AVCOL_SPC_SMPTE170M);
	av_frame_set_color_range(output_avstream_video_frame,AVCOL_RANGE_MPEG);
	output_avstream_video_frame->format = output_avstream_video_codec_context->pix_fmt;
	output_avstream_video_frame->height = output_height;
	output_avstream_video_frame->width = output_width;
	if (av_frame_get_buffer(output_avstream_video_frame,64) < 0) {
		fprintf(stderr,"Failed to alloc render frame\n");
		return 1;
	}

	// PARSE
	{
		uint8_t **dst_data = NULL;
		AVPixelFormat input_avstream_video_resampler_format = AV_PIX_FMT_NONE;
		int input_avstream_video_resampler_height = -1;
		int input_avstream_video_resampler_width = -1;
		unsigned long long video_field = 0;
		int dst_data_alloc_samples = 0;
		int dst_data_linesize = 0;
		int dst_data_samples = 0;
		int got_frame = 0;
		AVPacket pkt;

		av_init_packet(&pkt);
		while (av_read_frame(input_avfmt,&pkt) >= 0) {
			if (DIE != 0) break;

			if (input_avstream_audio != NULL && pkt.stream_index == input_avstream_audio->index) {
				av_packet_rescale_ts(&pkt,input_avstream_audio->time_base,output_avstream_audio->time_base);
				if (avcodec_decode_audio4(input_avstream_audio_codec_context,input_avstream_audio_frame,&got_frame,&pkt) >= 0) {
					if (got_frame != 0 && input_avstream_audio_frame->nb_samples != 0) {
						dst_data_samples = av_rescale_rnd(
							swr_get_delay(input_avstream_audio_resampler, input_avstream_audio_frame->sample_rate) + input_avstream_audio_frame->nb_samples,
							output_avstream_audio_codec_context->sample_rate, input_avstream_audio_frame->sample_rate, AV_ROUND_UP);

						if (dst_data == NULL || dst_data_samples > dst_data_alloc_samples) {
							if (dst_data != NULL) {
								av_freep(&dst_data[0]); // NTS: Why??
								av_freep(&dst_data);
							}

							dst_data_alloc_samples = 0;
							fprintf(stderr,"Allocating audio buffer %u samples\n",(unsigned int)dst_data_samples);
							if (av_samples_alloc_array_and_samples(&dst_data,&dst_data_linesize,
								output_avstream_audio_codec_context->channels,dst_data_samples,
								output_avstream_audio_codec_context->sample_fmt, 0) >= 0) {
								dst_data_alloc_samples = dst_data_samples;
							}
							else {
								fprintf(stderr,"Failure to allocate audio buffer\n");
								dst_data_alloc_samples = 0;
							}
						}

						if (dst_data != NULL) {
							int out_samples;

							if ((out_samples=swr_convert(input_avstream_audio_resampler,dst_data,dst_data_samples,
								(const uint8_t**)input_avstream_audio_frame->data,input_avstream_audio_frame->nb_samples)) > 0) {
								// PROCESS THE AUDIO. At this point by design the code can assume S16LE (16-bit PCM interleaved)
								composite_audio_process((int16_t*)dst_data[0],out_samples);
								// write it out. TODO: At some point, support conversion to whatever the codec needs and then convert to it.
								// that way we can render directly to MP4 our VHS emulation.
								AVPacket dstpkt;
								av_init_packet(&dstpkt);
								if (av_new_packet(&dstpkt,out_samples * 2 * output_audio_channels) >= 0) { // NTS: Will reset fields too!
									assert(dstpkt.data != NULL);
									assert(dstpkt.size >= (out_samples * 2 * output_audio_channels));
									memcpy(dstpkt.data,dst_data[0],out_samples * 2 * output_audio_channels);
								}
								av_packet_copy_props(&dstpkt,&pkt);
								dstpkt.stream_index = output_avstream_audio->index;
								if (av_interleaved_write_frame(output_avfmt,&dstpkt) < 0)
									fprintf(stderr,"Failed to write frame\n");
								av_packet_unref(&dstpkt);
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
			}
			else if (input_avstream_video != NULL && pkt.stream_index == input_avstream_video->index) {
				AVRational m = (AVRational){output_field_rate.den, output_field_rate.num};
				av_packet_rescale_ts(&pkt,input_avstream_video->time_base,m); // convert to FIELD number

				if (avcodec_decode_video2(input_avstream_video_codec_context,input_avstream_video_frame,&got_frame,&pkt) >= 0) {
					if (got_frame != 0 && input_avstream_video_frame->width > 0 && input_avstream_video_frame->height > 0) {
						unsigned long long tgt_field = input_avstream_video_frame->pkt_pts;
						if (tgt_field == AV_NOPTS_VALUE) tgt_field = input_avstream_video_frame->pkt_dts;
						if (tgt_field == AV_NOPTS_VALUE) tgt_field = pkt.dts;

						if (output_avstream_video_input_frame != NULL) {
							if (output_avstream_video_input_frame->height != input_avstream_video_frame->height) {
								if (output_avstream_video_input_frame != NULL)
									av_frame_free(&output_avstream_video_input_frame);
							}
						}

						if (output_avstream_video_input_frame != NULL) {
							while (video_field < tgt_field) {
								render_field(output_avstream_video_frame,output_avstream_video_input_frame,(int)(video_field & 1ULL) ^ 1/*bottom field first*/,video_field);
								composite_video_process(output_avstream_video_frame,(int)(video_field & 1ULL) ^ 1/*bottom field first*/,video_field);
								if ((video_field & 1ULL)) output_frame(output_avstream_video_frame,video_field - 1ULL);
								video_field++;
							}
						}

						if (output_avstream_video_input_frame == NULL) {
							fprintf(stderr,"New input frame\n");
							output_avstream_video_input_frame = av_frame_alloc();
							if (output_avstream_video_input_frame == NULL) {
								fprintf(stderr,"Failed to alloc video frame\n");
								return 1;
							}
							av_frame_set_colorspace(output_avstream_video_input_frame,AVCOL_SPC_SMPTE170M);
							av_frame_set_color_range(output_avstream_video_input_frame,AVCOL_RANGE_MPEG);
							output_avstream_video_input_frame->format = output_avstream_video_codec_context->pix_fmt;
							output_avstream_video_input_frame->height = input_avstream_video_frame->height;
							output_avstream_video_input_frame->width = output_width;
							if (av_frame_get_buffer(output_avstream_video_input_frame,64) < 0) {
								fprintf(stderr,"Failed to alloc render frame\n");
								return 1;
							}
						}

						if (input_avstream_video_resampler != NULL) { // pixel format change or width/height change = free resampler and reinit
							if (input_avstream_video_resampler_format != input_avstream_video_frame->format ||
								input_avstream_video_resampler_width != input_avstream_video_frame->width ||
								input_avstream_video_resampler_height != input_avstream_video_frame->height) {
								sws_freeContext(input_avstream_video_resampler);
								input_avstream_video_resampler = NULL;
							}
						}

						if (input_avstream_video_resampler == NULL) {
							input_avstream_video_resampler = sws_getContext(
								// source
								input_avstream_video_frame->width,
								input_avstream_video_frame->height,
								(AVPixelFormat)input_avstream_video_frame->format,
								// dest
								output_avstream_video_codec_context->width,
								input_avstream_video_frame->height,
								output_avstream_video_codec_context->pix_fmt,
								// opt
								SWS_BILINEAR, NULL, NULL, NULL);

							if (input_avstream_video_resampler != NULL) {
								fprintf(stderr,"sws_getContext new context\n");
								input_avstream_video_resampler_format = (AVPixelFormat)input_avstream_video_frame->format;
								input_avstream_video_resampler_width = input_avstream_video_frame->width;
								input_avstream_video_resampler_height = input_avstream_video_frame->height;
							}
							else {
								fprintf(stderr,"sws_getContext fail\n");
							}
						}

						if (input_avstream_video_resampler != NULL) {
							output_avstream_video_input_frame->pts = input_avstream_video_frame->pts;
							output_avstream_video_input_frame->pkt_pts = input_avstream_video_frame->pkt_pts;
							output_avstream_video_input_frame->pkt_dts = input_avstream_video_frame->pkt_dts;
							output_avstream_video_input_frame->top_field_first = input_avstream_video_frame->top_field_first;
							output_avstream_video_input_frame->interlaced_frame = input_avstream_video_frame->interlaced_frame;

							if (sws_scale(input_avstream_video_resampler,
								// source
								input_avstream_video_frame->data,
								input_avstream_video_frame->linesize,
								0,input_avstream_video_frame->height,
								// dest
								output_avstream_video_input_frame->data,
								output_avstream_video_input_frame->linesize) <= 0)
								fprintf(stderr,"WARNING: sws_scale failed\n");

							while (video_field < tgt_field) {
								render_field(output_avstream_video_frame,output_avstream_video_input_frame,(int)(video_field & 1ULL) ^ 1/*bottom field first*/,video_field);
								composite_video_process(output_avstream_video_frame,(int)(video_field & 1ULL) ^ 1/*bottom field first*/,video_field);
								if ((video_field & 1ULL)) output_frame(output_avstream_video_frame,video_field - 1ULL);
								video_field++;
							}
						}
					}
				}
				else {
					fprintf(stderr,"No video decoded\n");
				}
			}

			av_packet_unref(&pkt);
			av_init_packet(&pkt);
		}

		if (dst_data != NULL) {
			av_freep(&dst_data[0]); // NTS: Why??
			av_freep(&dst_data);
		}
	}

	if (output_avstream_video_input_frame != NULL)
		av_frame_free(&output_avstream_video_input_frame);
	if (output_avstream_video_frame != NULL)
		av_frame_free(&output_avstream_video_frame);
	if (input_avstream_video_frame != NULL)
		av_frame_free(&input_avstream_video_frame);
	if (input_avstream_audio_frame != NULL)
		av_frame_free(&input_avstream_audio_frame);
	audio_hilopass.clear();
	av_write_trailer(output_avfmt);
	if (input_avstream_video_resampler != NULL) {
		sws_freeContext(input_avstream_video_resampler);
		input_avstream_video_resampler = NULL;
	}
	if (input_avstream_audio_resampler != NULL)
		swr_free(&input_avstream_audio_resampler);
	if (output_avfmt != NULL && !(output_avfmt->oformat->flags & AVFMT_NOFILE))
		avio_closep(&output_avfmt->pb);
	avformat_free_context(output_avfmt);
	avformat_close_input(&input_avfmt);
	return 0;
}

