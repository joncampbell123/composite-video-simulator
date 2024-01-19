
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
#include <stdexcept>

double          gamma_correction = -1;

int             underscan = 0;

bool            use_422_colorspace = false;
AVRational	output_field_rate = { 60000, 1001 };	// NTSC 60Hz default
int		output_width = -1;
int		output_height = -1;
int		output_ar_n = 1,output_ar_d = 1;

#define RGBTRIPLET(r,g,b)           (((uint32_t)(r) << (uint32_t)16) + ((uint32_t)(g) << (uint32_t)8) + ((uint32_t)(b) << (uint32_t)0))

AVFormatContext*                    output_avfmt = NULL;
AVStream*                           output_avstream_video = NULL;	// do not free
AVCodecContext*	                    output_avstream_video_codec_context = NULL; // do not free
AVFrame*                            output_avstream_video_frame = NULL;         // ARGB
AVFrame*                            output_avstream_video_encode_frame = NULL;  // 4:2:2 or 4:2:0
struct SwsContext*                  output_avstream_video_resampler = NULL;

class InputFile {
	public:
		InputFile() {
			input_avfmt = NULL;
			input_avstream_video = NULL;
			input_avstream_video_frame = NULL;
			input_avstream_video_frame_rgb = NULL;
			input_avstream_video_resampler = NULL;
			input_avstream_video_codec_context = NULL;
			next_pts = next_dts = -1LL;
			avpkt_valid = false;
			eof_stream = false;
			eof = false;
		}
		~InputFile() {
			close_input();
		}
	public:
		double video_frame_to_output_f(void) {
			if (input_avstream_video_frame != NULL) {
				if (input_avstream_video_frame->pts != AV_NOPTS_VALUE) {
					double n = input_avstream_video_frame->pts;

					n *= (signed long long)input_avstream_video->time_base.num * (signed long long)output_field_rate.num;
					n /= (signed long long)input_avstream_video->time_base.den * (signed long long)output_field_rate.den;

					return n;
				}
			}

			return AV_NOPTS_VALUE;
		}
		double video_frame_rgb_to_output_f(void) {
			if (input_avstream_video_frame_rgb != NULL) {
				if (input_avstream_video_frame_rgb->pts != AV_NOPTS_VALUE) {
					double n = input_avstream_video_frame_rgb->pts;

					n *= (signed long long)input_avstream_video->time_base.num * (signed long long)output_field_rate.num;
					n /= (signed long long)input_avstream_video->time_base.den * (signed long long)output_field_rate.den;

					return n;
				}
			}

			return AV_NOPTS_VALUE;
		}
		void reset_on_dup(void) {
			path.clear();
		}
		bool open_input(void) {
			if (input_avfmt == NULL) {
				if (avformat_open_input(&input_avfmt,path.c_str(),NULL,NULL) < 0) {
					fprintf(stderr,"Failed to open input file\n");
					close_input();
					return false;
				}

				if (avformat_find_stream_info(input_avfmt,NULL) < 0)
					fprintf(stderr,"WARNING: Did not find stream info on input\n");

				/* scan streams for one video, one audio */
				{
					size_t i;
					AVStream *is;
					int ac=0,vc=0;
					AVCodecParameters *ispar;

					fprintf(stderr,"Input format: %u streams found\n",input_avfmt->nb_streams);
					for (i=0;i < (size_t)input_avfmt->nb_streams;i++) {
						is = input_avfmt->streams[i];
						if (is == NULL) continue;

						ispar = is->codecpar;
						if (ispar == NULL) continue;

						if (ispar->codec_type == AVMEDIA_TYPE_VIDEO) {
							if (input_avstream_video == NULL && vc == 0) {
								if ((input_avstream_video_codec_context=avcodec_alloc_context3(avcodec_find_decoder(ispar->codec_id))) != NULL) {
									if (avcodec_parameters_to_context(input_avstream_video_codec_context,ispar) < 0)
										fprintf(stderr,"WARNING: parameters to context failed\n");

									if (avcodec_open2(input_avstream_video_codec_context,avcodec_find_decoder(ispar->codec_id),NULL) >= 0) {
										input_avstream_video = is;
										fprintf(stderr,"Found video stream idx=%zu\n",i);
									}
									else {
										fprintf(stderr,"Found video stream but not able to decode\n");
										avcodec_free_context(&input_avstream_video_codec_context);
									}
								}
							}

							vc++;
						}
					}

					if (input_avstream_video == NULL) {
						fprintf(stderr,"Video not found\n");
						close_input();
						return 1;
					}
				}
			}

			/* prepare video decoding */
			if (input_avstream_video != NULL) {
				input_avstream_video_frame = av_frame_alloc();
				if (input_avstream_video_frame == NULL) {
					fprintf(stderr,"Failed to alloc video frame\n");
					close_input();
					return 1;
				}
			}

			input_avstream_video_resampler_format = AV_PIX_FMT_NONE;
			input_avstream_video_resampler_height = -1;
			input_avstream_video_resampler_width = -1;
			eof_stream = false;
			got_video = false;
			adj_time = 0;
			t = pt = -1;
			eof = false;
			avpkt_init();
			next_pts = next_dts = -1LL;
			return (input_avfmt != NULL);
		}

		uint32_t *copy_rgba(const AVFrame * const src) {
			assert(src != NULL);
			assert(src->data[0] != NULL);
			assert(src->linesize[0] != 0);
			assert(src->height != 0);

			assert(src->linesize[0] >= (src->width * 4));

			uint32_t *r = (uint32_t*)(new uint8_t[src->linesize[0] * src->height]);
			memcpy(r,src->data[0],src->linesize[0] * src->height);

			return r;
		}

		bool next_packet(void) {
			if (eof) return false;
			if (input_avfmt == NULL) return false;

			do {
				if (eof_stream) break;
				avpkt_release();
				avpkt_init();
				if (av_read_frame(input_avfmt,avpkt) < 0) {
					eof_stream = true;
					return false;
				}
				if (avpkt->stream_index >= input_avfmt->nb_streams)
					continue;

				// ugh... this can happen if the source is an AVI file
				if (avpkt->pts == AV_NOPTS_VALUE) avpkt->pts = avpkt->dts;

				/* track time and keep things monotonic for our code */
				if (avpkt->pts != AV_NOPTS_VALUE) {
					t = avpkt->pts * av_q2d(input_avfmt->streams[avpkt->stream_index]->time_base);

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

				if (pt < 0)
					continue;

				if (avpkt->pts != AV_NOPTS_VALUE) {
					avpkt->pts += (adj_time * input_avfmt->streams[avpkt->stream_index]->time_base.den) /
						input_avfmt->streams[avpkt->stream_index]->time_base.num;
				}

				if (avpkt->dts != AV_NOPTS_VALUE) {
					avpkt->dts += (adj_time * input_avfmt->streams[avpkt->stream_index]->time_base.den) /
						input_avfmt->streams[avpkt->stream_index]->time_base.num;
				}

				got_video = false;
				if (input_avstream_video != NULL && avpkt->stream_index == input_avstream_video->index) {
					if (got_video) fprintf(stderr,"Video content lost\n");
					//				AVRational m = (AVRational){output_field_rate.den, output_field_rate.num};
					//				av_packet_rescale_ts(avpkt,input_avstream_video->time_base,m); // convert to FIELD number
					handle_frame(/*&*/(*avpkt)); // will set got_video
					break;
				}

				avpkt_release();
			} while (1);

			if (eof_stream) {
				avpkt_release();
				handle_frame(); // will set got_video
				if (!got_video) eof = true;
				else fprintf(stderr,"Got latent frame\n");
			}

			return true;
		}
		void frame_copy_scale(void) {
			if (input_avstream_video_frame_rgb == NULL) {
				fprintf(stderr,"New input frame\n");
				input_avstream_video_frame_rgb = av_frame_alloc();
				if (input_avstream_video_frame_rgb == NULL) {
					fprintf(stderr,"Failed to alloc video frame\n");
					return;
				}

				input_avstream_video_frame_rgb->format = AV_PIX_FMT_BGRA;
				input_avstream_video_frame_rgb->height = output_height;
				input_avstream_video_frame_rgb->width = output_width;
				if (av_frame_get_buffer(input_avstream_video_frame_rgb,64) < 0) {
					fprintf(stderr,"Failed to alloc render frame\n");
					return;
				}
				memset(input_avstream_video_frame_rgb->data[0],0,input_avstream_video_frame_rgb->linesize[0]*input_avstream_video_frame_rgb->height);

				fprintf(stderr,"RGB is %d, %d\n",input_avstream_video_frame_rgb->width,input_avstream_video_frame_rgb->height);
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
				int final_w,final_h;

				final_w = (input_avstream_video_frame_rgb->width * (100 - underscan)) / 100;
				final_h = (input_avstream_video_frame_rgb->height * (100 - underscan)) / 100;

				if (final_w < 1) final_w = 1;
				if (final_h < 1) final_h = 1;

				input_avstream_video_resampler = sws_getContext(
						// source
						input_avstream_video_frame->width,
						input_avstream_video_frame->height,
						(AVPixelFormat)input_avstream_video_frame->format,
						// dest
						final_w,// input_avstream_video_frame_rgb->width,
						final_h,// input_avstream_video_frame_rgb->height,
						(AVPixelFormat)input_avstream_video_frame_rgb->format,
						// opt
						SWS_BILINEAR, NULL, NULL, NULL);

				if (input_avstream_video_resampler != NULL) {
					fprintf(stderr,"sws_getContext new context\n");
					input_avstream_video_resampler_format = (AVPixelFormat)input_avstream_video_frame->format;
					input_avstream_video_resampler_width = input_avstream_video_frame->width;
					input_avstream_video_resampler_height = input_avstream_video_frame->height;
					input_avstream_video_resampler_x = (input_avstream_video_frame_rgb->width - final_w) / 2;
					input_avstream_video_resampler_y = (input_avstream_video_frame_rgb->height - final_h) / 2;
					assert(input_avstream_video_resampler_x >= 0);
					assert(input_avstream_video_resampler_y >= 0);
					fprintf(stderr,"dst %d, %d\n",input_avstream_video_frame_rgb->width,input_avstream_video_frame_rgb->height);
					fprintf(stderr,"ofs %d, %d\n",input_avstream_video_resampler_x,input_avstream_video_resampler_y);
				}
				else {
					fprintf(stderr,"sws_getContext fail\n");
				}
			}

			if (input_avstream_video_resampler != NULL) {
				input_avstream_video_frame_rgb->pts = input_avstream_video_frame->pts;
				input_avstream_video_frame_rgb->pkt_dts = input_avstream_video_frame->pkt_dts;
				input_avstream_video_frame_rgb->top_field_first = input_avstream_video_frame->top_field_first;
				input_avstream_video_frame_rgb->interlaced_frame = input_avstream_video_frame->interlaced_frame;

				unsigned char *dst_planes[8] = {NULL};

				dst_planes[0]  = input_avstream_video_frame_rgb->data[0];
				dst_planes[0] += input_avstream_video_resampler_y * input_avstream_video_frame_rgb->linesize[0];
				dst_planes[0] += input_avstream_video_resampler_x * 4;

				if (sws_scale(input_avstream_video_resampler,
							// source
							input_avstream_video_frame->data,
							input_avstream_video_frame->linesize,
							0,input_avstream_video_frame->height,
							// dest
							dst_planes,//input_avstream_video_frame_rgb->data,
							input_avstream_video_frame_rgb->linesize) <= 0)
					fprintf(stderr,"WARNING: sws_scale failed\n");
			}
		}
		void handle_frame(void) {
			avcodec_send_packet(input_avstream_video_codec_context,NULL);

			if (avcodec_receive_frame(input_avstream_video_codec_context,input_avstream_video_frame) >= 0) {
				got_video = true;
			}
			else {
				got_video = false;
				fprintf(stderr,"No video decoded\n");
			}
		}
		void handle_frame(AVPacket &pkt) {
			avcodec_send_packet(input_avstream_video_codec_context,&pkt);

			if (avcodec_receive_frame(input_avstream_video_codec_context,input_avstream_video_frame) >= 0) {
				got_video = true;
			}
			else {
				got_video = false;
				fprintf(stderr,"No video decoded\n");
			}
		}
		void avpkt_init(void) {
			if (!avpkt_valid) {
				avpkt_valid = true;
				avpkt = av_packet_alloc();
			}
		}
		void avpkt_release(void) {
			if (avpkt_valid) {
				avpkt_valid = false;
				av_packet_free(&avpkt);
			}
			got_video = false;
		}
		void close_input(void) {
			eof = true;
			avpkt_release();
			if (input_avstream_video_codec_context != NULL) {
				avcodec_close(input_avstream_video_codec_context);
				avcodec_free_context(&input_avstream_video_codec_context);
				assert(input_avstream_video_codec_context == NULL);
				input_avstream_video = NULL;
			}

			if (input_avstream_video_frame != NULL)
				av_frame_free(&input_avstream_video_frame);
			if (input_avstream_video_frame_rgb != NULL)
				av_frame_free(&input_avstream_video_frame_rgb);

			if (input_avstream_video_resampler != NULL) {
				sws_freeContext(input_avstream_video_resampler);
				input_avstream_video_resampler = NULL;
			}

			avformat_close_input(&input_avfmt);
		}
	public:
		std::string             path;
		uint32_t                color;
		bool                    eof;
		bool                    eof_stream;
		bool                    got_video;
	public:
		AVFormatContext*        input_avfmt;
		AVStream*               input_avstream_video;	            // do not free
		AVCodecContext*         input_avstream_video_codec_context;
		AVFrame*                input_avstream_video_frame;
		AVFrame*                input_avstream_video_frame_rgb;
		struct SwsContext*      input_avstream_video_resampler;
		AVPixelFormat           input_avstream_video_resampler_format;
		int                     input_avstream_video_resampler_height;
		int                     input_avstream_video_resampler_width;
		int                     input_avstream_video_resampler_y;
		int                     input_avstream_video_resampler_x;
		signed long long        next_pts;
		signed long long        next_dts;
		AVPacket*               avpkt = NULL;
		bool                    avpkt_valid;
		double                  adj_time;
		double                  t,pt;
};

InputFile                   input_file;
std::string                 output_file;

volatile int DIE = 0;

void sigma(int x) {
	if (++DIE >= 20) abort();
}

void preset_NTSC() {
	output_field_rate.num = 60000;
	output_field_rate.den = 1001;
}

static void help(const char *arg0) {
	fprintf(stderr,"%s [options]\n",arg0);
	fprintf(stderr," -i <input file>               you can specify more than one input file, in order of layering\n");
	fprintf(stderr," -o <output file>\n");
	fprintf(stderr," -or <frame rate>\n");
	fprintf(stderr," -width <x>\n");
	fprintf(stderr," -height <x>\n");
	fprintf(stderr," -fa <x>                       Interpolate alternate frames\n");
	fprintf(stderr," -gamma <x>                    Interpolate with gamma correction (number, ntsc, vga)\n");
	fprintf(stderr," -underscan <x>                Underscan the image during rendering\n");
	fprintf(stderr," -422                          Render to 4:2:2 colorspace\n");
	fprintf(stderr," -420                          Render to 4:2:0 colorspace\n");
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
			else if (!strcmp(a,"width")) {
				a = argv[i++];
				if (a == NULL) return 1;
				output_width = (int)strtoul(a,NULL,0);
				if (output_width < 32) return 1;
			}
			else if (!strcmp(a,"height")) {
				a = argv[i++];
				if (a == NULL) return 1;
				output_height = (int)strtoul(a,NULL,0);
				if (output_height < 32) return 1;
			}
			else if (!strcmp(a,"gamma")) {
				a = argv[i++];
				if (a == NULL) return 1;

				if (isdigit(*a))
					gamma_correction = atof(a);
				else if (!strcmp(a,"vga") || !strcmp(a,"ntsc"))
					gamma_correction = 2.2;
			}
			else if (!strcmp(a,"i")) {
				a = argv[i++];
				if (a == NULL) return 1;
				input_file.path = a;
			}
			else if (!strcmp(a,"or")) {
				a = argv[i++];
				if (a == NULL) return 1;

				int d = 1;
				double n = strtof(a,(char**)(&a));
				if (*a == ':' || *a == '/' || *a == '\\') {
					a++;
					d = strtoul(a,(char**)(&a),10);
					if (d < 1) d = 1;
				}

				if (n < 0) n = 0;

				/* this code can cause problems below 5fps */
				if ((n/d) < 5) {
					n = 5;
					d = 1;
				}

				if (d > 1) {
					output_field_rate.num = (long)floor(n + 0.5);
					output_field_rate.den = (long)d;
				}
				else {
					output_field_rate.num = (long)floor((n * 10000) + 0.5);
					output_field_rate.den = (long)10000;
				}
			}
			else if (!strcmp(a,"o")) {
				a = argv[i++];
				if (a == NULL) return 1;
				output_file = a;
			}
			else if (!strcmp(a,"underscan")) {
				a = argv[i++];
				if (a == NULL) return 1;
				underscan = atoi(a);
				if (underscan < 0) underscan = 0;
				if (underscan > 99) underscan = 99;
			}
			else if (!strcmp(a,"422")) {
				use_422_colorspace = true;
			}
			else if (!strcmp(a,"420")) {
				use_422_colorspace = false;
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
	if (input_file.path.empty()) {
		fprintf(stderr,"No input files specified\n");
		return 1;
	}

	return 0;
}

void output_frame(AVFrame *frame,unsigned long long field_number) {
	int gotit = 0;
	AVPacket* pkt = av_packet_alloc();

	frame->key_frame = (field_number % (15ULL * 2ULL)) == 0 ? 1 : 0;

	{
		frame->interlaced_frame = 0;
		frame->pts = field_number;
		pkt->pts = field_number;
		pkt->dts = field_number;
	}

	fprintf(stderr,"\x0D" "Output field %llu ",field_number); fflush(stderr);

	avcodec_send_frame(output_avstream_video_codec_context,frame);
	while (avcodec_receive_packet(output_avstream_video_codec_context,pkt) >= 0) {
		pkt->stream_index = output_avstream_video->index;
		av_packet_rescale_ts(pkt,output_avstream_video_codec_context->time_base,output_avstream_video->time_base);

		if (av_interleaved_write_frame(output_avfmt,pkt) < 0)
			fprintf(stderr,"AV write frame failed video\n");
	}

	av_packet_free(&pkt);
}

// This code assumes ARGB and the frame match resolution/
void composite_layer(AVFrame *dstframe,AVFrame *srcframe,InputFile &inputfile) {
	uint32_t *dscan,*sscan;
	unsigned int x,y;
	unsigned int shr;

	if (dstframe == NULL || srcframe == NULL) return;
	if (dstframe->data[0] == NULL || srcframe->data[0] == 0) return;
	if (dstframe->linesize[0] < (dstframe->width*4)) return; // ARGB
	if (srcframe->linesize[0] < (srcframe->width*4)) return; // ARGB
	if (dstframe->width != srcframe->width) return;
	if (dstframe->height != srcframe->height) return;

	for (y=0;y < dstframe->height;y++) {
		sscan = (uint32_t*)(srcframe->data[0] + (srcframe->linesize[0] * y));
		dscan = (uint32_t*)(dstframe->data[0] + (dstframe->linesize[0] * y));
		for (x=0;x < dstframe->width;x++,dscan++,sscan++) {
			*dscan = *sscan;
		}
	}
}

int clamp255(int x) {
	if (x > 255)
		return 255;
	if (x < 0)
		return 0;
	return x;
}

double gamma_dec(double x) {
	return pow(x,1.0 / gamma_correction);
}

unsigned long gamma_dec16_table[256];
unsigned long gamma_enc16_table[8192 + 1];

bool gamma16_init = false;

void gamma16_do_init(void);

unsigned long gamma_dec16(unsigned long x) {
	if (!gamma16_init) gamma16_do_init();

	if (x > 255u) x = 255u;

	return gamma_dec16_table[x];
}

double gamma_enc(double x) {
	return pow(x,gamma_correction);
}

unsigned long gamma_enc16(unsigned long x) {
	if (!gamma16_init) gamma16_do_init();

	if (x > 8192u) x = 8192u;

	return gamma_enc16_table[x];
}

void gamma16_do_init(void) {
	gamma16_init = true;

	for (unsigned int i=0;i < 256;i++)
		gamma_dec16_table[i] = (unsigned long)(gamma_dec(i / 255.0) * 8192);

	for (unsigned int i=0;i <= 8192;i++)
		gamma_enc16_table[i] = (unsigned long)(gamma_enc(i / 8192.0) * 255);
}

int main(int argc,char **argv) {
	preset_NTSC();
	if (parse_argv(argc,argv))
		return 1;

	/* open all input files */
	if (!input_file.open_input()) {
		fprintf(stderr,"Failed to open %s\n",input_file.path.c_str());
		return 1;
	}

	/* pick output */
	if (output_width < 1 || output_height < 1) {
		if (input_file.input_avstream_video_codec_context != NULL) {
			output_width = input_file.input_avstream_video_codec_context->width;
			output_height = input_file.input_avstream_video_codec_context->height;
			output_ar_n = input_file.input_avstream_video_codec_context->sample_aspect_ratio.num;
			output_ar_d = input_file.input_avstream_video_codec_context->sample_aspect_ratio.den;
		}
	}
	fprintf(stderr,"Output frame: %d x %d with %d:%d PAR\n",output_width,output_height,output_ar_n,output_ar_d);

	/* no decision, no frame */
	if (output_width < 16 || output_height < 16) {
		fprintf(stderr,"None or invalid output dimensions\n");
		return 1;
	}

	/* open output file */
	assert(output_avfmt == NULL);
	if (avformat_alloc_output_context2(&output_avfmt,NULL,NULL,output_file.c_str()) < 0) {
		fprintf(stderr,"Failed to open output file\n");
		return 1;
	}

	{
		AVDictionary *opt_dict = NULL;

		output_avstream_video = avformat_new_stream(output_avfmt, NULL);
		if (output_avstream_video == NULL) {
			fprintf(stderr,"Unable to create output video stream\n");
			return 1;
		}

		output_avstream_video_codec_context = avcodec_alloc_context3(avcodec_find_encoder(AV_CODEC_ID_H264));
		if (output_avstream_video_codec_context == NULL) {
			fprintf(stderr,"Output stream video no codec context?\n");
			return 1;
		}

		output_avstream_video_codec_context->width = output_width;
		output_avstream_video_codec_context->height = output_height;
		output_avstream_video_codec_context->sample_aspect_ratio = (AVRational){output_ar_n,output_ar_d};
		output_avstream_video_codec_context->pix_fmt = use_422_colorspace ? AV_PIX_FMT_YUV422P : AV_PIX_FMT_YUV420P;
		output_avstream_video_codec_context->gop_size = 15;
		output_avstream_video_codec_context->max_b_frames = 0;
		output_avstream_video_codec_context->time_base = (AVRational){output_field_rate.den, output_field_rate.num};

		av_dict_set(&opt_dict,"crf","16",0);
		av_dict_set(&opt_dict,"crf_max","16",0);
		av_dict_set(&opt_dict,"preset","superfast",0);

		output_avstream_video->time_base = output_avstream_video_codec_context->time_base;
		if (output_avfmt->oformat->flags & AVFMT_GLOBALHEADER)
			output_avstream_video_codec_context->flags |= AV_CODEC_FLAG_GLOBAL_HEADER;

		if (avcodec_open2(output_avstream_video_codec_context,avcodec_find_encoder(AV_CODEC_ID_H264),&opt_dict) < 0) {
			fprintf(stderr,"Output stream cannot open codec\n");
			return 1;
		}

		if (avcodec_parameters_from_context(output_avstream_video->codecpar,output_avstream_video_codec_context) < 0)
			fprintf(stderr,"WARNING: parameters from context failed\n");

		av_dict_free(&opt_dict);
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
		output_avstream_video_encode_frame->colorspace = AVCOL_SPC_SMPTE170M;
		output_avstream_video_encode_frame->color_range = AVCOL_RANGE_MPEG;
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
		signed long long outbase=0;

		memset(output_avstream_video_frame->data[0],0x00,output_avstream_video_frame->linesize[0] * output_avstream_video_frame->height);

		{
			signed long long current=0;

			while (!input_file.eof && !DIE) {
				if (!input_file.got_video)
					input_file.next_packet();
				else
					break;
			}

			if (input_file.input_avstream_video_frame != NULL && input_file.got_video) {
				input_file.frame_copy_scale();
				input_file.got_video = false;
			}

			std::vector<uint32_t*> frames; /* they're all the same RGBA frame, so just store a pointer */
			std::vector<double> frame_t;

			if (input_file.input_avstream_video_frame_rgb != NULL) {
				frame_t.push_back(input_file.video_frame_rgb_to_output_f());
				frames.push_back(input_file.copy_rgba(input_file.input_avstream_video_frame_rgb));
			}

			while (!DIE) {
				size_t cutoff = 0;

				while (!input_file.eof && !DIE && input_file.video_frame_to_output_f() < (current + 30LL)) {
					input_file.next_packet();

					if (input_file.input_avstream_video_frame != NULL && input_file.got_video) {
						input_file.frame_copy_scale();
						input_file.got_video = false;

						if (input_file.input_avstream_video_frame_rgb != NULL) {
							frame_t.push_back(input_file.video_frame_rgb_to_output_f());
							frames.push_back(input_file.copy_rgba(input_file.input_avstream_video_frame_rgb));
						}
					}
				}

				if (input_file.eof &&
					(input_file.video_frame_to_output_f() < -1000/*AV_NOPTS_VALUE*/ ||
					 current > (unsigned long long)ceil(input_file.video_frame_to_output_f())))
					break;

				for (size_t i=0;(i+1ul) < frames.size();i++) {
					double bt = frame_t[i];
					double et = frame_t[i+1];

					if (i != 0) {
						if ((et + 2.0) < current) {
							cutoff = i;
						}
					}
				}

				long *lframe = new long[output_width*output_height*3];

				if (gamma_correction > 1) {
					for (unsigned int y=0;y < output_height;y++) {
						unsigned char *outframe = (unsigned char*)(output_avstream_video_frame->data[0] + (y * (output_avstream_video_frame->linesize[0])));
						unsigned char *inframe = ((unsigned char*)frames[0] + (y * (input_file.input_avstream_video_frame_rgb->linesize[0])));
						long *longframe = lframe + (y * output_width * 3);

						for (unsigned int x=0;x < output_width;x++) {
							longframe[x*3+0] = gamma_dec16(inframe[x*4+0]) << 16ul;
							longframe[x*3+1] = gamma_dec16(inframe[x*4+1]) << 16ul;
							longframe[x*3+2] = gamma_dec16(inframe[x*4+2]) << 16ul;
						}
					}
				}
				else {
					for (unsigned int y=0;y < output_height;y++) {
						unsigned char *outframe = (unsigned char*)(output_avstream_video_frame->data[0] + (y * (output_avstream_video_frame->linesize[0])));
						unsigned char *inframe = ((unsigned char*)frames[0] + (y * (input_file.input_avstream_video_frame_rgb->linesize[0])));
						long *longframe = lframe + (y * output_width * 3);

						for (unsigned int x=0;x < output_width;x++) {
							longframe[x*3+0] = inframe[x*4+0] << 16ul;
							longframe[x*3+1] = inframe[x*4+1] << 16ul;
							longframe[x*3+2] = inframe[x*4+2] << 16ul;
						}
					}
				}

				long scaleto = gamma_correction > 1 ? (0x10000l * 8192l) : 0x10000l;
				long minv = (scaleto * 6l) / 10l;
				long maxv = (scaleto * 4l) / 10l;

				{
					unsigned int minx = (output_width*15)/100;
					unsigned int maxx = (output_width*90)/100;
					unsigned int miny = (output_height*10)/100;
					unsigned int maxy = (output_height*90)/100;

					for (unsigned int y=miny;(y+16) < maxy;y += 16) {
						for (unsigned int x=minx;(x+16) < maxx;x += 16) {
							long gr = 0;

							for (unsigned int sy=0;sy < 16;sy++) {
								long *longframe = lframe + ((y+sy) * output_width * 3);
								for (unsigned int sx=0;sx < 16;sx++) {
									/* BGR */
									gr +=	((longframe[(x+sx)*3+0]/*B*/ * 11l) +
										 (longframe[(x+sx)*3+1]/*G*/ * 59l) +
										 (longframe[(x+sx)*3+2]/*R*/ * 30l) + 50l) / 100l;
								}
							}

							gr = (gr + 128) / 256;

							if (minv > gr)
								minv = gr;
							if (maxv < gr)
								maxv = gr;
						}
					}
				}

				if (minv == maxv) maxv++;

				{
					const long dist = maxv - minv;
					minv -= dist / 10;
					maxv += dist / 25;
				}

//				fprintf(stderr,"\nmin=%.15f max=%.15f\n",(double)minv/scaleto,(double)maxv/scaleto);

				for (unsigned int y=0;y < output_height;y++) {
					long *longframe = lframe + (y * output_width * 3);
					for (unsigned int x=0;x < output_width*3;x++) {
						long long v = (((long long)(longframe[x] - minv)) * (long long)scaleto) / ((long long)(maxv - minv));
						if (v < -0x7FFFFFFFl) v = -0x7FFFFFFFl;
						if (v >  0x7FFFFFFFl) v =  0x7FFFFFFFl;
						longframe[x] = (long)v;
					}
				}

				if (gamma_correction > 1) {
					for (unsigned int y=0;y < output_height;y++) {
						unsigned char *outframe = (unsigned char*)(output_avstream_video_frame->data[0] + (y * (output_avstream_video_frame->linesize[0])));
						unsigned char *inframe = ((unsigned char*)frames[0] + (y * (input_file.input_avstream_video_frame_rgb->linesize[0])));
						long *longframe = lframe + (y * output_width * 3);

						for (unsigned int x=0;x < output_width;x++) {
							outframe[x*4+0] = clamp255(gamma_enc16(std::max(0l,longframe[x*3+0] >> 16l)));
							outframe[x*4+1] = clamp255(gamma_enc16(std::max(0l,longframe[x*3+1] >> 16l)));
							outframe[x*4+2] = clamp255(gamma_enc16(std::max(0l,longframe[x*3+2] >> 16l)));
							outframe[x*4+3] = 0xFF;
						}
					}
				}
				else {
					for (unsigned int y=0;y < output_height;y++) {
						unsigned char *outframe = (unsigned char*)(output_avstream_video_frame->data[0] + (y * (output_avstream_video_frame->linesize[0])));
						unsigned char *inframe = ((unsigned char*)frames[0] + (y * (input_file.input_avstream_video_frame_rgb->linesize[0])));
						long *longframe = lframe + (y * output_width * 3);

						for (unsigned int x=0;x < output_width;x++) {
							outframe[x*4+0] = clamp255(std::max(0l,longframe[x*3+0] >> 16l));
							outframe[x*4+1] = clamp255(std::max(0l,longframe[x*3+1] >> 16l));
							outframe[x*4+2] = clamp255(std::max(0l,longframe[x*3+2] >> 16l));
							outframe[x*4+3] = 0xFF;
						}
					}
				}

				delete[] lframe;

				output_avstream_video_frame->pts = current;
				output_avstream_video_frame->pkt_dts = current;
				output_avstream_video_frame->top_field_first = 0;
				output_avstream_video_frame->interlaced_frame = 0;

				// convert ARGB to whatever the codec demands, and encode
				output_avstream_video_encode_frame->pts = output_avstream_video_frame->pts;
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

				if (cutoff > 0) {
					assert(frame_t.size() > cutoff);
					assert(frames.size() > cutoff);

					for (size_t i=0;i < cutoff;i++) {
						if (frames[i] != NULL) {
							delete[] frames[i];
							frames[i] = NULL;
						}
					}

					frame_t.erase(frame_t.begin(),frame_t.begin()+cutoff);
					frames.erase(frames.begin(),frames.begin()+cutoff);
				}
			}

			for (size_t i=0;i < frames.size();i++) {
				if (frames[i] != NULL) {
					delete[] frames[i];
					frames[i] = NULL;
				}
			}
			frames.clear();
			frame_t.clear();

			outbase += current;
		}
	}

	/* flush encoder delay */
	fprintf(stderr,"Flushing delayed frames\n");
	{
		AVPacket *pkt = av_packet_alloc();

		avcodec_send_frame(output_avstream_video_codec_context,NULL);
		while (avcodec_receive_packet(output_avstream_video_codec_context,pkt) >= 0) {
			pkt->stream_index = output_avstream_video->index;
			av_packet_rescale_ts(pkt,output_avstream_video_codec_context->time_base,output_avstream_video->time_base);

			if (av_interleaved_write_frame(output_avfmt,pkt) < 0)
				fprintf(stderr,"AV write frame failed video\n");
		}

		av_packet_free(&pkt);
	}
	fprintf(stderr,"Flushing delayed frames--done\n");

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

	avcodec_free_context(&output_avstream_video_codec_context);
	avformat_free_context(output_avfmt);

	/* close all */
	input_file.close_input();

	return 0;
}

