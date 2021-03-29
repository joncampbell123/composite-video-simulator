
// NTS: This is not like modern "posterize" filters where the pixels are quantizied to N levels then scaled out to 0..255
//      That requires a multiply/divide per pixel. Think old-school hardware where such operations were too expensive.
//      The "posterize" we emulate here is more the type where you run the video through an ADC, truncate the least significant
//      bits, then run back through a DAC on the other side (well within the realm of 1980s/1990s hardware)

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

#include <map>
#include <list>
#include <string>
#include <vector>
#include <stdexcept>

std::list<string>           src_composite;

int                         src_fd = -1;

bool            use_422_colorspace = true;
AVRational	output_field_rate = { 60000, 1001 };	// NTSC 60Hz default
int		output_width = 720;
int		output_height = 480;
bool        input_ntsc = false;
bool		output_ntsc = true;	// NTSC color subcarrier emulation
bool		output_pal = false;	// PAL color subcarrier emulation
int		output_audio_channels = 2;	// VHS stereo (set to 1 for mono)
int		output_audio_rate = 44100;	// VHS Hi-Fi goes up to 20KHz

static const double         subcarrier_freq = (315000000.00        / 88.00);        /* 315/88 MHz = about 3.57954 MHz */
static const double         sample_rate =    ((315000000.00 * 8.0) / 88.00);        /* 315/88 * MHz = about 28.636363 MHz */
static const double         one_frame_time =   sample_rate / (30000.00 / 1001.00);  /* 30000/1001 = about 29.97 */
static const double         one_scanline_time = one_frame_time / 525.00;            /* one scanline */
static const unsigned int   one_scanline_raw_length = (unsigned int)(one_scanline_time + 0.5);

double                      one_scanline_width = one_scanline_time + 0.0495;
double                      one_scanline_width_err = 0;

unsigned char               int_scanline[4096];

std::vector<uint8_t>                input_samples;
std::vector<uint8_t>::iterator      input_samples_read,input_samples_end;

size_t count_src() {
    assert(input_samples_read <= input_samples_end);
    return (input_samples_end - input_samples_read);
}

void empty_src() {
    input_samples_read = input_samples_end = input_samples.begin();
}

void flush_src() {
    assert(input_samples_read <= input_samples_end);
    if (input_samples_read != input_samples.begin()) {
        size_t move = input_samples_read - input_samples.begin();
        assert(move != 0);
        size_t todo = input_samples_end - input_samples_read;
        if (todo > 0) memmove(&(*input_samples.begin()),&(*input_samples_read),todo);
        input_samples_read -= move;
        assert(input_samples_read == input_samples.begin());
        input_samples_end -= move;
    }
}

void refill_src() {
    if (src_fd >= 0) {
        assert(input_samples_read <= input_samples_end);
        if (input_samples_end < input_samples.end()) {
            size_t todo = input_samples.end() - input_samples_end;
            assert(todo <= input_samples.size());
            int rd = read(src_fd,&(*input_samples_end),todo);
            if (rd > 0) {
                assert((size_t)rd <= todo);
                input_samples_end += (size_t)rd;
                assert(input_samples_end <= input_samples.end());
            }
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

        src_fd = open(path.c_str(),O_RDONLY);
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
	output_field_rate.num = 30000;
	output_field_rate.den = 1001;
	output_height = 525;
	output_width = 910;
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
			else if (!strcmp(a,"tvstd")) {
				a = argv[i++];

				if (!strcmp(a,"pal")) {
					preset_PAL();
				}
				else if (!strcmp(a,"ntsc")) {
					preset_NTSC();
				}
                else if (!strcmp(a,"720p60")) {
                    preset_720p60();
                }
                else if (!strcmp(a,"1080p60")) {
                    preset_1080p60();
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

    for (y=0;y < dstframe->height;y++) {
        lazy_flush_src();
        if (count_src() < (one_scanline_raw_length*2)) {
            empty_src();
            break;
        }

        {
            int a = (int)floor(one_scanline_width_err * 256);
            if (a < 0) a = 0;
            if (a > 256) a = 256;
            for (x=0;x < (one_scanline_raw_length+16);x++)
                int_scanline[x] = (unsigned char)(((input_samples_read[x] * (256 - a)) + (input_samples_read[x+1] * a)) >> 8);
        }

        uint32_t *dst = (uint32_t*)(dstframe->data[0] + (dstframe->linesize[0] * y));
        for (x=0;x < dstframe->width;x++) {
            size_t si = x * 2;
            int r,g,b;
            int Y = int_scanline[si];

            r = g = b = Y;
            if (r < 0) r = 0;
            if (g < 0) g = 0;
            if (b < 0) b = 0;
            if (r > 255) r = 255;
            if (g > 255) g = 255;
            if (b > 255) b = 255;

            dst[x] = (r + (g << 8) + (b << 16) + 0xFF000000);
        }

        {
            unsigned int adj = floor(one_scanline_width);
            one_scanline_width_err += one_scanline_width - adj;
            if (one_scanline_width_err >= 1.0) {
                one_scanline_width_err -= 1.0;
                adj++;
            }
            input_samples_read += adj;
        }
    }
}

int main(int argc,char **argv) {
    preset_NTSC();
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

    fprintf(stderr,"Subcarrier:             %.3f\n",subcarrier_freq);
    fprintf(stderr,"Sample rate:            %.3f\n",sample_rate);
    fprintf(stderr,"One frame duration:     %.3f\n",one_frame_time);
    fprintf(stderr,"One scanline duration:  %.3f\n",one_scanline_time);
    fprintf(stderr,"Raw render to:          %u\n",one_scanline_raw_length);

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
            if (count_src() < one_scanline_raw_length) {
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

