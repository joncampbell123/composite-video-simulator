// what to do next:
//
// If emulating PAL, fake subcarrier needs to alternate color phase every other field scanline.
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

int		audio_stream_index = 0;
int		video_stream_index = 0;

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
AVFrame*		output_avstream_video_input_frame = NULL;   // 4:2:2
AVFrame*		output_avstream_video_frame = NULL;         // 4:2:2

AVFrame*		output_avstream_video_bob_frame = NULL;     // 4:2:0 or 4:2:2

bool            use_422_colorspace = false; // I would default this to true but Adobe Premiere Pro apparently can't handle 4:2:2 H.264 >:(

double			composite_preemphasis = 0;	// analog artifacts related to anything that affects the raw composite signal i.e. CATV modulation
double			composite_preemphasis_cut = 1000000;

double			vhs_out_sharpen = 1.5;
double			vhs_out_sharpen_chroma = 0.85;

bool			vhs_head_switching = false;
double			vhs_head_switching_phase = 1.0 - ((4.5+0.01/*slight error, like most VHS tapes*/) / 262.5); // 4 scanlines NTSC up from vsync
double			vhs_head_switching_phase_noise = (((1.0 / 300)/*slight error, like most VHS tapes*/) / 262.5); // 1/300th of a scanline

int		video_yc_recombine = 0;			// additional Y/C combine/sep phases (testing)
int		video_color_fields = 4;			// NTSC color framing
int		video_chroma_noise = 0;
int		video_chroma_phase_noise = 0;
int		video_chroma_loss = 0;
int		video_noise = 2;
int		subcarrier_amplitude = 50;
int		subcarrier_amplitude_back = 50;
bool		output_video_as_interlaced = false;	// render as 480i (half field rate). else, render at field rate with bob filter
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
double		vhs_linear_high_boost = 0.25;
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
bool		vhs_chroma_vert_blend = true;		// if set, and VHS, blend vertically the chroma scanlines (as the VHS format does)
bool		vhs_svideo_out = false;			// if not set, and VHS, video is recombined as if composite out on VCR
bool        enable_composite_emulation = true; // if not set, video goes straight back out to the encoder.

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

/* render the chroma into the luma as a fake NTSC color subcarrier */
void composite_video_yuv_to_ntsc(AVFrame *dst,unsigned int field,unsigned long long fieldno,const int subcarrier_amplitude) {
	unsigned int x,y;

	{ /* lowpass the chroma more. composite video does not allocate as much bandwidth to color as luma. */
		for (unsigned int p=1;p <= 2;p++) {
			for (y=field;y < dst->height;y += 2) {
				unsigned char *P = dst->data[p] + (y * dst->linesize[p]);
				LowpassFilter lp[3];
				LowpassFilter hp;
				double cutoff;
				int delay;
				double s;

				if (output_ntsc) {
					// NTSC YIQ bandwidth: I=1.3MHz Q=0.6MHz
					cutoff = (p == 1) ? 1300000 : 600000;
					delay = (p == 1) ? 2 : 4;
				}
				else {
					// PAL: R-Y and B-Y are 1.3MHz
					cutoff = 1300000;
					delay = 2;
				}

				hp.setFilter((315000000.00 * 4) / (88 * 2),cutoff/2); // 315/88 Mhz rate * 4 (divide by 2 for 4:2:2)  vs 600KHz cutoff
				hp.resetFilter(128);
				for (unsigned int f=0;f < 3;f++) {
					lp[f].setFilter((315000000.00 * 4) / (88 * 2),cutoff); // 315/88 Mhz rate * 4 (divide by 2 for 4:2:2)  vs 600KHz cutoff
					lp[f].resetFilter(128);
				}

				for (x=0;x < (dst->width/2)/*4:2:2*/;x++) {
					s = P[x];
					s += hp.highpass(s);
					for (unsigned int f=0;f < 3;f++) s = lp[f].lowpass(s);
					if (x >= delay) P[x-delay] = clampu8(s);
				}
			}
		}
	}

	for (y=field;y < dst->height;y += 2) {
		static const int8_t Umult[4] = { 1, 0,-1, 0 };
		static const int8_t Vmult[4] = { 0, 1, 0,-1 };
		unsigned char *Y = dst->data[0] + (y * dst->linesize[0]);
		unsigned char *U = dst->data[1] + (y * dst->linesize[1]);
		unsigned char *V = dst->data[2] + (y * dst->linesize[2]);
		unsigned int xc = dst->width;
		unsigned int xi;

		if (output_ntsc) { // NTSC 2 color frames long
			xi = (fieldno + (y >> 1)) & 3;
		}
		else/*PAL*/ {
			// FIXME: Is this right?
			xi = (fieldno + y) & 3;
		}

		/* remember: this code assumes 4:2:2 */
		/* NTS: the subcarrier is two sine waves superimposed on top of each other, 90 degrees apart */
		for (x=0;x < xc;x += 2,Y += 2,U++,V++) {
			for (unsigned int sx=0;sx < 2;sx++) {
				unsigned int sxi = xi+x+sx;
				int chroma;

				chroma  = ((int)U[0] - 128) * subcarrier_amplitude * Umult[sxi&3];
				chroma += ((int)V[0] - 128) * subcarrier_amplitude * Vmult[sxi&3];
				Y[sx] = clampu8(Y[sx] + (chroma / 50));
			}

			if (nocolor_subcarrier)
				U[0] = V[0] = 128;
		}
	}
}

/* filter subcarrier back out, use result to emulate NTSC luma-chroma artifacts */
void composite_ntsc_to_yuv(AVFrame *dst,unsigned int field,unsigned long long fieldno,const int subcarrier_amplitude_back) {
	unsigned char chroma[dst->width]; // WARNING: This is more GCC-specific C++ than normal
	unsigned int x,y;

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

			if (output_ntsc) { // NTSC 2 color frames long
				xi = (fieldno + (y >> 1)) & 3;
			}
			else/*PAL*/ {
				// FIXME: Is this right?
				xi = (fieldno + y) & 3;
			}

			for (x=((4-xi)&3);x < dst->width;x += 4) { // flip the part of the sine wave that would correspond to negative U and V values
				chroma[x+2] = 255 - chroma[x+2];
				chroma[x+3] = 255 - chroma[x+3];
			}

			for (x=0;x < dst->width;x++) {
				chroma[x] = clampu8(((((int)chroma[x] - 128) * 50) / subcarrier_amplitude_back) + 128);
			}

			/* decode the color right back out from the subcarrier we generated */
			if (xi & 1) {
				for (x=0;x < (dst->width/2);x++) {
					U[x] = 255 - chroma[(x*2)+1];
					V[x] = 255 - chroma[(x*2)+0];
				}
			}
			else {
				for (x=0;x < (dst->width/2);x++) {
					U[x] = 255 - chroma[(x*2)+0];
					V[x] = 255 - chroma[(x*2)+1];
				}
			}
		}
	}

	{ /* lowpass the chroma more. composite video does not allocate as much bandwidth to color as luma. */
		for (unsigned int p=1;p <= 2;p++) {
			for (y=field;y < dst->height;y += 2) {
				unsigned char *P = dst->data[p] + (y * dst->linesize[p]);
				LowpassFilter lp[3];
				double cutoff;
				int delay;
				double s;

				if (output_ntsc) {
					// NTSC YIQ bandwidth: I=1.3MHz Q=0.6MHz
					cutoff = (p == 1) ? 1300000 : 600000;
					delay = (p == 1) ? 2 : 4;
				}
				else {
					// PAL: R-Y and B-Y are 1.3MHz
					cutoff = 1300000;
					delay = 2;
				}

				for (unsigned int f=0;f < 3;f++) {
					lp[f].setFilter((315000000.00 * 4) / (88 * 2),cutoff); // 315/88 Mhz rate * 4 (divide by 2 for 4:2:2)  vs 600KHz cutoff
					lp[f].resetFilter(128);
				}

				for (x=0;x < (dst->width/2)/*4:2:2*/;x++) {
					s = P[x];
					for (unsigned int f=0;f < 3;f++) s = lp[f].lowpass(s);
					if (x >= delay) P[x-delay] = clampu8(s);
				}
			}
		}
	}
}

static unsigned long long audio_proc_count = 0;
static LowpassFilter audio_post_vhs_boost[2];

void composite_audio_process(int16_t *audio,unsigned int samples) { // number of channels = output_audio_channels, sample rate = output_audio_rate. audio is interleaved.
	assert(audio_hilopass.audiostate.size() >= output_audio_channels);
	double linear_buzz = dBFS(output_audio_linear_buzz);
	double hsync_hz = output_ntsc ? /*NTSC*/15734 : /*PAL*/15625;
	int vsync_lines = output_ntsc ? /*NTSC*/525 : /*PAL*/625;
	int vpulse_end = output_ntsc ? /*NTSC*/10 : /*PAL*/12;
	double hpulse_end = output_ntsc ? /*NTSC*/(hsync_hz * (4.7/*us*/ / 1000000)) : /*PAL*/(hsync_hz * (4.0/*us*/ / 1000000));

	for (unsigned int s=0;s < samples;s++,audio += output_audio_channels) {
		for (unsigned int c=0;c < output_audio_channels;c++) {
			double s;

			s = (double)audio[c] / 32768;

			/* that faint "buzzing" noise on linear tracks because of audio/video crosstalk */
			if (!output_vhs_hifi && linear_buzz > 0.000000001) {
				const unsigned int oversample = 16;
				for (unsigned int oi=0;oi < oversample;oi++) {
					double t = ((((double)audio_proc_count * oversample) + oi) * hsync_hz) / output_audio_rate / oversample;
					double hpos = fmod(t,1.0);
					int vline = (int)fmod(floor(t + 0.0001/*fudge*/ - hpos),(double)vsync_lines / 2);
					bool pulse = false;

					if (hpos < hpulse_end)
						pulse = true; // HSYNC
					if (vline < vpulse_end)
						pulse = true; // VSYNC

					if (pulse)
						s -= linear_buzz / oversample / 2;
				}
			}

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

			/* some VCRs (at least mine) will boost higher frequencies if playing linear tracks */
			if (!output_vhs_hifi && vhs_linear_high_boost > 0)
				s += audio_post_vhs_boost[c].highpass(s) * vhs_linear_high_boost;

			/* deemphasis */
			if (emulating_deemphasis) {
				for (unsigned int i=0;i < output_audio_channels;i++) {
					s = audio_linear_preemphasis_post[i].lowpass(s);
				}
			}

			audio[c] = clips16(s * 32768);
		}

		audio_proc_count++;
	}
}

void composite_video_process(AVFrame *dst,unsigned int field,unsigned long long fieldno) {
	unsigned int x,y;

	composite_video_yuv_to_ntsc(dst,field,fieldno,subcarrier_amplitude);

	/* video composite preemphasis */
	if (composite_preemphasis != 0 && composite_preemphasis_cut > 0) {
		for (y=field;y < dst->height;y += 2) {
			unsigned char *Y = dst->data[0] + (y * dst->linesize[0]);
			LowpassFilter pre;
			double s;

			pre.setFilter((315000000.00 * 4) / 88,composite_preemphasis_cut); // 315/88 Mhz rate * 4  vs 1.0MHz cutoff
			pre.resetFilter(16);
			for (x=0;x < dst->width;x++) {
				s = Y[x];
				s += pre.highpass(s) * composite_preemphasis;
				Y[x] = clampu8(s);
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

	// VHS head switching noise
	if (vhs_head_switching) {
		unsigned int twidth = dst->width + (dst->width / 10);
		unsigned int tx,x,p,x2,shy=0;
		double noise = 0;
		int shif,ishif,y;
		double t;

		if (vhs_head_switching_phase_noise != 0) {
			unsigned int x = (unsigned int)rand() * (unsigned int)rand() * (unsigned int)rand() * (unsigned int)rand();
			x %= 2000000000U;
			noise = ((double)x / 1000000000U) - 1.0;
			noise *= vhs_head_switching_phase_noise;
		}

		if (output_ntsc)
			t = twidth * 262.5;
		else
			t = twidth * 312.5;

		p = (unsigned int)(fmod(vhs_head_switching_phase + noise,1.0) * t);
		x = p % (unsigned int)twidth;
		y = ((p / (unsigned int)twidth) * 2) + field;

		if (output_ntsc)
			y -= (262 - 240) * 2;
		else
			y -= (312 - 288) * 2;

		tx = x;
		if (x >= (twidth/2))
			ishif = x - twidth;
		else
			ishif = x;

		shif = 0;
		while (y < dst->height) {
			if (y >= 0) {
				unsigned char *Y = dst->data[0] + (y * dst->linesize[0]);

				if (shif != 0) {
					char tmp[twidth];

					/* WARNING: This is not 100% accurate. On real VHS you'd see the line shifted over and the next line's contents after hsync. */

					/* luma. the chroma subcarrier is there, so this is all we have to do. */
					x2 = (tx + twidth + (unsigned int)shif) % (unsigned int)twidth;
					memset(tmp,16,sizeof(tmp));
					memcpy(tmp,Y,dst->width);
					for (x=tx;x < dst->width;x++) {
						Y[x] = tmp[x2];
						if ((++x2) == twidth) x2 = 0;
					}
				}
			}

			if (shy == 0)
				shif = ishif;
			else
				shif = (shif * 7) / 8;

			tx = 0;
			y += 2;
			shy++;
		}
	}

	if (!nocolor_subcarrier)
		composite_ntsc_to_yuv(dst,field,fieldno,subcarrier_amplitude_back);

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
	if (video_chroma_phase_noise != 0) {
		int noise = 0,noise_mod = (video_chroma_noise * 255) / 100;
		double pi,u,v,u_,v_;

		for (y=field;y < dst->height;y += 2) {
			unsigned char *U = dst->data[1] + (y * dst->linesize[1]);
			unsigned char *V = dst->data[2] + (y * dst->linesize[2]);

			noise += ((int)((unsigned int)rand() % ((video_chroma_phase_noise*2)+1))) - video_chroma_phase_noise;
			noise /= 2;
			pi = ((double)noise * M_PI) / 100;

			for (x=0;x < (dst->width/2);x++) {
				u = (int)U[x] - 128; // think of 'u' as x-coord
				v = (int)V[x] - 128; // and 'v' as y-coord

				// then this 2D rotation then makes more sense
				u_ = (u * cos(pi)) - (u * sin(pi));
				v_ = (v * cos(pi)) + (v * sin(pi));

				// put it back
				U[x] = clampu8(u_ + 128);
				V[x] = clampu8(v_ + 128);
			}
		}
	}

	// NTS: At this point, the video best resembles what you'd get from a typical DVD player's composite video output.
	//      Slightly blurry, some color artifacts, and edges will have that "buzz" effect, but still a good picture.

	if (emulating_vhs) {
		double luma_cut,chroma_cut;
		int chroma_delay;

		switch (output_vhs_tape_speed) {
			case VHS_SP:
				luma_cut = 2400000; // 3.0MHz x 80%
				chroma_cut = 320000; // 400KHz x 80%
				chroma_delay = 4;
				break;
			case VHS_LP:
				luma_cut = 1900000; // ..
				chroma_cut = 300000; // 375KHz x 80%
				chroma_delay = 5;
				break;
			case VHS_EP:
				luma_cut = 1400000; // ..
				chroma_cut = 280000; // 350KHz x 80%
				chroma_delay = 6;
				break;
			default:
				abort();
		};

		// luma lowpass
		for (y=field;y < dst->height;y += 2) {
			unsigned char *Y = dst->data[0] + (y * dst->linesize[0]);
			LowpassFilter lp[3];
			LowpassFilter pre;
			double s;

			for (unsigned int f=0;f < 3;f++) {
				lp[f].setFilter((315000000.00 * 4) / 88,luma_cut); // 315/88 Mhz rate * 4  vs 3.0MHz cutoff
				lp[f].resetFilter(16);
			}
			pre.setFilter((315000000.00 * 4) / 88,luma_cut); // 315/88 Mhz rate * 4  vs 1.0MHz cutoff
			pre.resetFilter(16);
			for (x=0;x < dst->width;x++) {
				s = Y[x];
				for (unsigned int f=0;f < 3;f++) s = lp[f].lowpass(s);
				s += pre.highpass(s) * 1.6;
				Y[x] = clampu8(s);
			}
		}

		// chroma lowpass
		for (y=field;y < dst->height;y += 2) {
			unsigned char *U = dst->data[1] + (y * dst->linesize[1]);
			unsigned char *V = dst->data[2] + (y * dst->linesize[2]);
			LowpassFilter lpU[3],lpV[3];
			double s;

			for (unsigned int f=0;f < 3;f++) {
				lpU[f].setFilter((315000000.00 * 4) / (88 * 2/*4:2:2*/),chroma_cut); // 315/88 Mhz rate * 4 (divide by 2 for 4:2:2) vs 400KHz cutoff
				lpU[f].resetFilter(128);
				lpV[f].setFilter((315000000.00 * 4) / (88 * 2/*4:2:2*/),chroma_cut); // 315/88 Mhz rate * 4 (divide by 2 for 4:2:2) vs 400KHz cutoff
				lpV[f].resetFilter(128);
			}
			for (x=0;x < (dst->width/2);x++) {
				s = U[x];
				for (unsigned int f=0;f < 3;f++) s = lpU[f].lowpass(s);
				if (x >= chroma_delay) U[x-chroma_delay] = clampu8(s);

				s = V[x];
				for (unsigned int f=0;f < 3;f++) s = lpV[f].lowpass(s);
				if (x >= chroma_delay) V[x-chroma_delay] = clampu8(s);
			}
		}

		// VHS decks also vertically smear the chroma subcarrier using a delay line
		// to add the previous line's color subcarrier to the current line's color subcarrier.
		// note that phase changes in NTSC are compensated for by the VHS deck to make the
		// phase line up per scanline (else summing the previous line's carrier would
		// cancel it out).
		if (vhs_chroma_vert_blend && output_ntsc) {
			unsigned char delayU[dst->width/2];
			unsigned char delayV[dst->width/2];

			memset(delayU,128,dst->width/2);
			memset(delayV,128,dst->width/2);
			for (y=(field+2);y < dst->height;y += 2) {
				unsigned char *U = dst->data[1] + (y * dst->linesize[1]);
				unsigned char *V = dst->data[2] + (y * dst->linesize[2]);
				unsigned char cU,cV;

				for (x=0;x < (dst->width/2);x++) {
					cU = U[x];
					cV = V[x];
					U[x] = (delayU[x]+cU+1)>>1;
					V[x] = (delayV[x]+cV+1)>>1;
					delayU[x] = cU;
					delayV[x] = cV;
				}
			}
		}

		// VHS decks tend to sharpen the picture on playback
		if (true/*TODO make option*/) {
			// luma
			for (y=field;y < dst->height;y += 2) {
				unsigned char *Y = dst->data[0] + (y * dst->linesize[0]);
				LowpassFilter lp[3];
				double s,ts;

				for (unsigned int f=0;f < 3;f++) {
					lp[f].setFilter((315000000.00 * 4) / 88,luma_cut*2); // 315/88 Mhz rate * 4  vs 3.0MHz cutoff
					lp[f].resetFilter(16);
				}
				for (x=0;x < dst->width;x++) {
					s = ts = Y[x];
					for (unsigned int f=0;f < 3;f++) ts = lp[f].lowpass(ts);
					Y[x] = clampu8(s + ((s - ts) * vhs_out_sharpen));
				}
			}

			// chroma
			for (y=field;y < dst->height;y += 2) {
				unsigned char *U = dst->data[1] + (y * dst->linesize[1]);
				unsigned char *V = dst->data[2] + (y * dst->linesize[2]);
				LowpassFilter lpU[3],lpV[3];
				double s,ts;

				for (unsigned int f=0;f < 3;f++) {
					lpU[f].setFilter((315000000.00 * 4) / (88 * 2/*4:2:2*/),chroma_cut*2); // 315/88 Mhz rate * 4 (divide by 2 for 4:2:2) vs 400KHz cutoff
					lpU[f].resetFilter(128);
					lpV[f].setFilter((315000000.00 * 4) / (88 * 2/*4:2:2*/),chroma_cut*2); // 315/88 Mhz rate * 4 (divide by 2 for 4:2:2) vs 400KHz cutoff
					lpV[f].resetFilter(128);
				}
				for (x=0;x < (dst->width/2);x++) {
					s = ts = U[x];
					for (unsigned int f=0;f < 3;f++) ts = lpU[f].lowpass(ts);
					U[x] = clampu8(s + ((s - ts) * vhs_out_sharpen_chroma));

					s = ts = V[x];
					for (unsigned int f=0;f < 3;f++) ts = lpV[f].lowpass(ts);
					V[x] = clampu8(s + ((s - ts) * vhs_out_sharpen_chroma));
				}
			}
		}

		if (!vhs_svideo_out) {
			composite_video_yuv_to_ntsc(dst,field,fieldno,subcarrier_amplitude);
			composite_ntsc_to_yuv(dst,field,fieldno,subcarrier_amplitude);
		}
	}

	if (video_chroma_loss != 0) {
		for (y=field;y < dst->height;y += 2) {
			unsigned char *U = dst->data[1] + (y * dst->linesize[1]);
			unsigned char *V = dst->data[2] + (y * dst->linesize[2]);

			if ((((unsigned int)rand())%100000) < video_chroma_loss) {
				memset(U,128,dst->width/2);
				memset(V,128,dst->width/2);
			}
		}
	}

	for (int i=0;i < video_yc_recombine;i++) {
		composite_video_yuv_to_ntsc(dst,field,fieldno,subcarrier_amplitude);
		composite_ntsc_to_yuv(dst,field,fieldno,subcarrier_amplitude);
	}
}

void render_field(AVFrame *dst,AVFrame *src,unsigned int field,unsigned long long field_number) {
	unsigned int y,sy,sy2,syf,csy,csy2,csyf;
    unsigned int chroma_height;

    if (output_avstream_video_input_frame->format == AV_PIX_FMT_YUV420P)
        chroma_height = src->height >> 1;
    else
        chroma_height = src->height;

	// NTS: dst is 4:2:2 of output_width x output_height
	//      src is 4:2:2 of output_width x source frame height
    //
    //      if the video source is 4:2:0 then src is 4:2:0
	//
	//      the reason we do that is so that swscale handles horizontal scaling, then we handle
	//      vertical scaling in the way we need to in order to render interlaced video properly.
	//
	//      this code renders only one field or the other at a time.
	for (y=field;y < dst->height;y += 2) {
		sy = (y * 0x100 * src->height) / dst->height;
		syf = sy & 0xFF;
		sy >>= 8;

        csy = sy;
        csyf = syf;
        if (output_avstream_video_input_frame->format == AV_PIX_FMT_YUV420P) {
            if (!(csy&1)) csyf = 0;
            csy >>= 1;
        }

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

			if (which_field == 0) { // make it even. do not interpolate if first even line of the pair.
				csy++; // but shift up the frame 1 line
				if (!(csy & 1U)) csyf = 0;
				else csy--;
			}
			else {
				if (!(csy & 1U)) { // make it odd. do not interpolate if first odd line of the pair.
					csyf = 0;
					csy++;
				}
			}

			if (sy >= (src->height - 2)) {
				sy = src->height - 2;
				syf = 0;
			}
			sy2 = sy + 2;

			if (csy >= (chroma_height - 2)) {
				csy = chroma_height - 2;
				csyf = 0;
			}
			csy2 = csy + 1;
		}
		else {
			if (sy >= (src->height - 1)) {
				sy = src->height - 1;
				syf = 0;
			}
			sy2 = sy + 1;

			if (csy >= (chroma_height - 1)) {
				csy = chroma_height - 1;
				csyf = 0;
			}
			csy2 = csy + 1;
		}

        if (output_avstream_video_input_frame->format == AV_PIX_FMT_YUV420P) {
            for (unsigned int p=0;p < 3;p++) {
                unsigned char *s1 = src->data[p] + (src->linesize[p] * sy);
                unsigned char *s2 = src->data[p] + (src->linesize[p] * sy2);
                unsigned char *d = dst->data[p] + (dst->linesize[p] * y);

                if (syf == 0) {
                    memcpy(d,s1,src->linesize[p]);
                }
                else {
                    for (unsigned int x=0;x < src->linesize[p];x++)
                        d[x] = s1[x] + ((unsigned char)((((int)s2[x] - (int)s1[x]) * (int)syf) >> (int)8));
                }

                if (p == 0) {
                    sy = csy;
                    sy2 = csy2;
                    syf = csyf;
                }
            }
        }
        else {
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
                    unsigned char *s2 = src->data[p] + (src->linesize[p] * sy2);
                    unsigned char *d = dst->data[p] + (dst->linesize[p] * y);

                    for (unsigned int x=0;x < src->linesize[p];x++)
                        d[x] = s1[x] + ((unsigned char)((((int)s2[x] - (int)s1[x]) * (int)syf) >> (int)8));
                }
            }
        }
	}
}

void output_frame(AVFrame *frame,unsigned long long field_number,unsigned int field) {
	int gotit = 0;
	AVPacket pkt;

	av_init_packet(&pkt);
	if (av_new_packet(&pkt,50000000/8) < 0) {
		fprintf(stderr,"Failed to alloc vid packet\n");
		return;
	}

	frame->key_frame = (field_number % (15ULL * 2ULL)) == 0 ? 1 : 0;

	if (output_video_as_interlaced) {
		frame->interlaced_frame = 1;
		frame->top_field_first = (field == 0)?1:0;
		frame->pts = field_number / 2ULL;
		pkt.pts = field_number / 2ULL;
		pkt.dts = field_number / 2ULL;
	}
	else {
		frame->interlaced_frame = 0;
		frame->pts = field_number;
		pkt.pts = field_number;
		pkt.dts = field_number;
	}

	fprintf(stderr,"Output field %llu\n",field_number);
	if (output_video_as_interlaced && use_422_colorspace) { // 4:2:2 interlaced = use as-is
        if (avcodec_encode_video2(output_avstream_video_codec_context,&pkt,frame,&gotit) == 0) {
            if (gotit) {
                pkt.stream_index = output_avstream_video->index;
                av_packet_rescale_ts(&pkt,output_avstream_video_codec_context->time_base,output_avstream_video->time_base);

                if (av_interleaved_write_frame(output_avfmt,&pkt) < 0)
                    fprintf(stderr,"AV write frame failed video\n");
            }
        }
    }
	else {
		output_avstream_video_bob_frame->interlaced_frame = frame->interlaced_frame;
		output_avstream_video_bob_frame->top_field_first = frame->top_field_first;
		output_avstream_video_bob_frame->pts = frame->pts;

		assert(frame->height <= output_avstream_video_bob_frame->height);
		assert(frame->width <= output_avstream_video_bob_frame->width);

        if (use_422_colorspace) {
            for (unsigned int y=0;y < frame->height;y++) {
                unsigned int sy;

                if (field)
                    sy = (y|1); // 1, 1, 3, 3, 5, 5, ....
                else
                    sy = (y+1)&(~1); // 0, 2, 2, 4, 4, 6, 6, ...

                if (sy >= frame->height)
                    sy -= 2;

                memcpy(output_avstream_video_bob_frame->data[0] + (output_avstream_video_bob_frame->linesize[0] * y),
                        frame->data[0] + (frame->linesize[0] * sy),frame->width);//luma

                for (unsigned int p=1;p <= 2;p++) {
                    memcpy(output_avstream_video_bob_frame->data[p] + (output_avstream_video_bob_frame->linesize[p] * y),
                            frame->data[p] + (frame->linesize[p] * sy),frame->width/2);//chroma
                }
            }
        }
        else { // 4:2:0
            for (unsigned int y=0;y < frame->height;y++) {
                unsigned int sy;

                if (output_video_as_interlaced)
                    sy = y;
                else if (field)
                    sy = (y|1); // 1, 1, 3, 3, 5, 5, ....
                else
                    sy = (y+1)&(~1); // 0, 2, 2, 4, 4, 6, 6, ...

                if (sy >= frame->height)
                    sy -= 2;

                memcpy(output_avstream_video_bob_frame->data[0] + (output_avstream_video_bob_frame->linesize[0] * y),
                        frame->data[0] + (frame->linesize[0] * sy),frame->width);//luma

                if (output_video_as_interlaced) {
                    if ((y&2) == 0) {
                        unsigned int cy = (y & 1) + ((y & (~3)) >> 1);

                        for (unsigned int p=1;p <= 2;p++) {
                            memcpy(output_avstream_video_bob_frame->data[p] + (output_avstream_video_bob_frame->linesize[p] * cy),
                                    frame->data[p] + (frame->linesize[p] * sy),frame->width/2);//chroma
                        }
                    }
                }
                else {
                    if ((y&1) == 0) {
                        unsigned int cy = y >> 1;

                        for (unsigned int p=1;p <= 2;p++) {
                            memcpy(output_avstream_video_bob_frame->data[p] + (output_avstream_video_bob_frame->linesize[p] * cy),
                                    frame->data[p] + (frame->linesize[p] * sy),frame->width/2);//chroma
                        }
                    }
                }
            }
        }

		if (avcodec_encode_video2(output_avstream_video_codec_context,&pkt,output_avstream_video_bob_frame,&gotit) == 0) {
			if (gotit) {
				pkt.stream_index = output_avstream_video->index;
				av_packet_rescale_ts(&pkt,output_avstream_video_codec_context->time_base,output_avstream_video->time_base);

				if (av_interleaved_write_frame(output_avfmt,&pkt) < 0)
					fprintf(stderr,"AV write frame failed video\n");
			}
		}
	}

	av_packet_unref(&pkt);
}

void preset_PAL() {
	output_field_rate.num = 50;
	output_field_rate.den = 1;
	video_color_fields = 8;
	output_height = 576;
	output_width = 720;
	output_pal = true;
	output_ntsc = false;
}

void preset_NTSC() {
	output_field_rate.num = 60000;
	output_field_rate.den = 1001;
	video_color_fields = 4;
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
	fprintf(stderr," -vhs-linear-video-crosstalk <x> Emulate video crosstalk in audio. Loudness in dBFS (0=100%)\n");
	fprintf(stderr," -chroma-phase-noise <x>   Chroma phase noise (0...100)\n");
	fprintf(stderr," -vhs-chroma-vblend <0|1>  Vertically blend chroma scanlines (as VHS format does)\n");
	fprintf(stderr," -vhs-svideo <0|1>         Render VHS as if S-Video (luma and chroma separate out of VHS)\n");
	fprintf(stderr," -yc-recomb <n>            Recombine Y/C n-times\n");
	fprintf(stderr," -a <n>                    Pick the n'th audio stream\n");
	fprintf(stderr," -an                       Don't render any audio stream\n");
	fprintf(stderr," -v <n>                    Pick the n'th video stream\n");
	fprintf(stderr," -vn                       Don't render any video stream\n");
	fprintf(stderr," -comp-pre <s>             Composite preemphasis scale\n");
	fprintf(stderr," -comp-cut <f>             Composite preemphasis freq\n");
	fprintf(stderr," -comp-catv                Composite preemphasis preset, as if CATV #1\n");
	fprintf(stderr," -comp-catv2               Composite preemphasis preset, as if CATV #2\n");
	fprintf(stderr," -comp-catv3               Composite preemphasis preset, as if CATV #3\n");
	fprintf(stderr," -vi                       Render video at frame rate, interlaced\n");
	fprintf(stderr," -vp                       Render video at field rate, progressive (with bob filter)\n");
	fprintf(stderr," -chroma-dropout <x>       Chroma scanline dropouts (0...10000)\n");
	fprintf(stderr," -vhs-linear-high-boost <x> Boost high frequencies in VHS audio (linear tracks)\n");
	fprintf(stderr," -vhs-head-switching <0|1> Enable/disable VHS head switching emulation\n");
	fprintf(stderr," -vhs-head-switching-point <x> Head switching point (0....1)\n");
	fprintf(stderr," -vhs-head-switching-noise-level <x> Head switching noise (variation)\n");
    fprintf(stderr," -422                      Render in 4:2:2 colorspace\n");
    fprintf(stderr," -420                      Render in 4:2:0 colorspace (default)\n"); // dammit Premiere >:(
    fprintf(stderr," -nocomp                   Don't apply emulation, just transcode\n");
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
            else if (!strcmp(a,"nocomp")) {
                enable_composite_emulation = false;
            }
            else if (!strcmp(a,"422")) {
                use_422_colorspace = true;
            }
            else if (!strcmp(a,"420")) {
                use_422_colorspace = false;
            }
			else if (!strcmp(a,"a")) {
				audio_stream_index = atoi(argv[i++]);
			}
			else if (!strcmp(a,"v")) {
				video_stream_index = atoi(argv[i++]);
			}
			else if (!strcmp(a,"an")) {
				audio_stream_index = -1;
			}
			else if (!strcmp(a,"vn")) {
				video_stream_index = -1;
			}
			else if (!strcmp(a,"vi")) {
				output_video_as_interlaced = true;
			}
			else if (!strcmp(a,"vp")) {
				output_video_as_interlaced = false;
			}
			else if (!strcmp(a,"vhs-head-switching-point")) {
				vhs_head_switching_phase = atof(argv[i++]);
			}
			else if (!strcmp(a,"vhs-head-switching-noise-level")) {
				vhs_head_switching_phase_noise = atof(argv[i++]);
			}
			else if (!strcmp(a,"vhs-head-switching")) {
				int x = atoi(argv[i++]);
				vhs_head_switching = (x > 0)?true:false;
			}
			else if (!strcmp(a,"vhs-linear-high-boost")) {
				vhs_linear_high_boost = atof(argv[i++]);
			}
			else if (!strcmp(a,"comp-pre")) {
				composite_preemphasis = atof(argv[i++]);
			}
			else if (!strcmp(a,"comp-cut")) {
				composite_preemphasis_cut = atof(argv[i++]);
			}
			else if (!strcmp(a,"comp-catv")) {
				composite_preemphasis = 1.5;
				composite_preemphasis_cut = 315000000 / 88 / 2;
				video_chroma_phase_noise = 2;
			}
			else if (!strcmp(a,"comp-catv2")) {
				composite_preemphasis = 2.5;
				composite_preemphasis_cut = 315000000 / 88 / 2;
				video_chroma_phase_noise = 4;
			}
			else if (!strcmp(a,"comp-catv3")) {
				composite_preemphasis = 4;
				composite_preemphasis_cut = 315000000 / 88 / 2;
				video_chroma_phase_noise = 6;
			}
			else if (!strcmp(a,"vhs-linear-video-crosstalk")) {
				output_audio_linear_buzz = atof(argv[i++]);
			}
			else if (!strcmp(a,"chroma-phase-noise")) {
				int x = atoi(argv[i++]);
				video_chroma_phase_noise = x;
			}
			else if (!strcmp(a,"yc-recomb")) {
				video_yc_recombine = atof(argv[i++]);
			}
			else if (!strcmp(a,"audio-hiss")) {
				output_audio_hiss_db = atof(argv[i++]);
			}
			else if (!strcmp(a,"vhs-svideo")) {
				int x = atoi(argv[i++]);
				vhs_svideo_out = (x > 0);
			}
			else if (!strcmp(a,"vhs-chroma-vblend")) {
				int x = atoi(argv[i++]);
				vhs_chroma_vert_blend = (x > 0);
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
				subcarrier_amplitude_back = x;
			}
			else if (!strcmp(a,"nocolor-subcarrier")) {
				nocolor_subcarrier = true;
			}
			else if (!strcmp(a,"nocolor-subcarrier-after-yc-sep")) {
				nocolor_subcarrier_after_yc_sep = true;
			}
			else if (!strcmp(a,"chroma-dropout")) {
				int x = atoi(argv[i++]);
				video_chroma_loss = x;
			}
			else if (!strcmp(a,"vhs")) {
				emulating_vhs = true;
				vhs_head_switching = true;
				emulating_preemphasis = false; // no preemphasis by default
				emulating_deemphasis = false; // no preemphasis by default
				output_audio_hiss_db = -70;
				video_chroma_phase_noise = 4;
				video_chroma_noise = 16;
				video_chroma_loss = 4;
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
					video_chroma_phase_noise = 6;
					video_chroma_noise = 22;
					video_chroma_loss = 8;
					video_noise = 6;
				}
				else if (!strcmp(a,"lp")) {
					output_vhs_tape_speed = VHS_LP;
					video_chroma_phase_noise = 5;
					video_chroma_noise = 19;
					video_chroma_loss = 6;
					video_noise = 5;
				}
				else if (!strcmp(a,"sp")) {
					output_vhs_tape_speed = VHS_SP;
					video_chroma_phase_noise = 4;
					video_chroma_noise = 16;
					video_chroma_loss = 4;
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

	if (composite_preemphasis != 0)
		subcarrier_amplitude_back += (50 * composite_preemphasis) / 4;

	output_audio_hiss_level = dBFS(output_audio_hiss_db) * 5000;

	fprintf(stderr,"VHS head switching point: %.6f\n",vhs_head_switching_phase);
	fprintf(stderr,"VHS head switching noise: %.6f\n",vhs_head_switching_phase_noise);
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
			else if (isctx->codec_type == AVMEDIA_TYPE_VIDEO) {
				if (input_avstream_video == NULL && vc == video_stream_index) {
					if (avcodec_open2(isctx,avcodec_find_decoder(isctx->codec_id),NULL) >= 0) {
						input_avstream_video = is;
						input_avstream_video_codec_context = isctx;
						fprintf(stderr,"Found video stream idx=%zu\n",i);
					}
					else {
						fprintf(stderr,"Found video stream but not able to decode\n");
					}
				}

				vc++;
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
		output_avstream_video_codec_context->pix_fmt = use_422_colorspace ? AV_PIX_FMT_YUV422P : AV_PIX_FMT_YUV420P;
		output_avstream_video_codec_context->gop_size = 15;
		output_avstream_video_codec_context->max_b_frames = 0;

		if (output_video_as_interlaced)
			output_avstream_video_codec_context->time_base = (AVRational){output_field_rate.den, (output_field_rate.num/2)};
		else
			output_avstream_video_codec_context->time_base = (AVRational){output_field_rate.den, output_field_rate.num};

		output_avstream_video->time_base = output_avstream_video_codec_context->time_base;
		if (output_avfmt->oformat->flags & AVFMT_GLOBALHEADER)
			output_avstream_video_codec_context->flags |= AV_CODEC_FLAG_GLOBAL_HEADER;

		// our output is interlaced
		if (output_video_as_interlaced)
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

	/* high boost on playback */
	for (unsigned int i=0;i < 2;i++)
		audio_post_vhs_boost[i].setFilter(output_audio_rate,10000);

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
	if (input_avstream_audio != NULL) {
		input_avstream_audio_frame = av_frame_alloc();
		if (input_avstream_audio_frame == NULL) {
			fprintf(stderr,"Failed to alloc audio frame\n");
			return 1;
		}
	}

	/* prepare video decoding */
	if (input_avstream_video != NULL) {
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
		output_avstream_video_frame->format = AV_PIX_FMT_YUV422P;
		output_avstream_video_frame->height = output_height;
		output_avstream_video_frame->width = output_width;
		if (av_frame_get_buffer(output_avstream_video_frame,64) < 0) {
			fprintf(stderr,"Failed to alloc render frame\n");
			return 1;
		}

		if (!output_video_as_interlaced || !use_422_colorspace) {
			output_avstream_video_bob_frame = av_frame_alloc();
			if (output_avstream_video_bob_frame == NULL) {
				fprintf(stderr,"Failed to alloc video frame2\n");
				return 1;
			}
			av_frame_set_colorspace(output_avstream_video_bob_frame,AVCOL_SPC_SMPTE170M);
			av_frame_set_color_range(output_avstream_video_bob_frame,AVCOL_RANGE_MPEG);
			output_avstream_video_bob_frame->format = output_avstream_video_codec_context->pix_fmt;
			output_avstream_video_bob_frame->height = output_height;
			output_avstream_video_bob_frame->width = output_width;
			if (av_frame_get_buffer(output_avstream_video_bob_frame,64) < 0) {
				fprintf(stderr,"Failed to alloc render frame2\n");
				return 1;
			}
		}
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

                                if (enable_composite_emulation)
                                    composite_video_process(output_avstream_video_frame,(int)(video_field & 1ULL) ^ 1/*bottom field first*/,video_field);

                                if (output_video_as_interlaced) {
                                    if ((video_field & 1ULL)) output_frame(output_avstream_video_frame,video_field - 1ULL,(int)((video_field - 1ULL) & 1ULL) ^ 1/*bottom field first*/);
                                }
                                else {
                                    output_frame(output_avstream_video_frame,video_field,(int)(video_field & 1ULL) ^ 1/*bottom field first*/);
                                }

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

                            // HACK: libswscale does NOT do proper 4:2:0 to 4:2:2 interlaced conversion.
                            //       so if the source is 4:2:0 then we want upconversion to 4:2:0 and
                            //       we'll do it properly ourselves. using libswscale means we end up
                            //       with the chroma fields backwards.
                            switch (input_avstream_video_frame->format) {
                                case AV_PIX_FMT_YUV420P:
                                case AV_PIX_FMT_YUVJ420P:
                                    output_avstream_video_input_frame->format = AV_PIX_FMT_YUV420P;
                                    break;
                                default:
                                    output_avstream_video_input_frame->format = AV_PIX_FMT_YUV422P;
                                    break;
                            }

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
								output_avstream_video_input_frame->width,
								output_avstream_video_input_frame->height,
								(AVPixelFormat)output_avstream_video_input_frame->format,
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

								if (output_video_as_interlaced) {
									if ((video_field & 1ULL)) output_frame(output_avstream_video_frame,video_field - 1ULL,(int)((video_field - 1ULL) & 1ULL) ^ 1/*bottom field first*/);
								}
								else {
									output_frame(output_avstream_video_frame,video_field,(int)(video_field & 1ULL) ^ 1/*bottom field first*/);
								}

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
	if (output_avstream_video_bob_frame != NULL)
		av_frame_free(&output_avstream_video_bob_frame);
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

