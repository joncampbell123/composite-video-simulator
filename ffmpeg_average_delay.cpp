
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
#include <string>
#include <vector>
#include <stdexcept>

bool            use_422_colorspace = false; // I would default this to true but Adobe Premiere Pro apparently can't handle 4:2:2 H.264 >:(
AVRational	output_field_rate = { 60000, 1001 };	// NTSC 60Hz default
int		output_width = 720;
int		output_height = 480;
bool		output_ntsc = true;	// NTSC color subcarrier emulation
bool		output_pal = false;	// PAL color subcarrier emulation
int		output_audio_channels = 2;	// VHS stereo (set to 1 for mono)
int		output_audio_rate = 44100;	// VHS Hi-Fi goes up to 20KHz

#define RGBTRIPLET(r,g,b)       (((uint32_t)(r) << (uint32_t)16) + ((uint32_t)(g) << (uint32_t)8) + ((uint32_t)(b) << (uint32_t)0))

AVFormatContext*	        output_avfmt = NULL;
AVStream*		            output_avstream_audio = NULL;	// do not free
AVCodecContext*		        output_avstream_audio_codec_context = NULL; // do not free
AVStream*		            output_avstream_video = NULL;	// do not free
AVCodecContext*		        output_avstream_video_codec_context = NULL; // do not free
AVFrame*		            output_avstream_video_frame = NULL;         // ARGB
AVFrame*		            output_avstream_video_encode_frame = NULL;  // 4:2:2 or 4:2:0
struct SwsContext*          output_avstream_video_resampler = NULL;

class InputFile {
public:
    InputFile() : newlevel(128) {
        input_avfmt = NULL;
        audio_dst_data = NULL;
        input_avstream_audio = NULL;
        input_avstream_audio_frame = NULL;
        input_avstream_audio_resampler = NULL;
        input_avstream_audio_codec_context = NULL;
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
                AVCodecContext *isctx;

                fprintf(stderr,"Input format: %u streams found\n",input_avfmt->nb_streams);
                for (i=0;i < (size_t)input_avfmt->nb_streams;i++) {
                    is = input_avfmt->streams[i];
                    if (is == NULL) continue;

                    isctx = is->codec;
                    if (isctx == NULL) continue;

                    if (isctx->codec_type == AVMEDIA_TYPE_AUDIO) {
                        if (input_avstream_audio == NULL && ac == 0) {
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
                        if (input_avstream_video == NULL && vc == 0) {
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
                    close_input();
                    return 1;
                }
            }
        }

        /* prepare audio decoding */
        if (input_avstream_audio != NULL) {
            input_avstream_audio_frame = av_frame_alloc();
            if (input_avstream_audio_frame == NULL) {
                fprintf(stderr,"Failed to alloc audio frame\n");
                close_input();
                return 1;
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

            /* prepare video encoding */
            input_avstream_video_frame_rgb = av_frame_alloc();
            if (input_avstream_video_frame_rgb == NULL) {
                fprintf(stderr,"Failed to alloc video frame\n");
                close_input();
                return 1;
            }
            input_avstream_video_frame_rgb->format = AV_PIX_FMT_BGRA;
            input_avstream_video_frame_rgb->height = output_height;
            input_avstream_video_frame_rgb->width = output_width;
            if (av_frame_get_buffer(input_avstream_video_frame_rgb,64) < 0) {
                fprintf(stderr,"Failed to alloc render frame\n");
                close_input();
                return 1;
            }
            memset(input_avstream_video_frame_rgb->data[0],0,input_avstream_video_frame_rgb->linesize[0]*input_avstream_video_frame_rgb->height);
        }

        input_avstream_video_resampler_format = AV_PIX_FMT_NONE;
        input_avstream_video_resampler_height = -1;
        input_avstream_video_resampler_width = -1;
        last_written_sample = 0;
        audio_dst_data_out_audio_sample = 0;
        audio_sample = 0;
        audio_dst_data = NULL;
        audio_dst_data_alloc_samples = 0;
        audio_dst_data_linesize = 0;
        audio_dst_data_samples = 0;
        audio_dst_data_out_samples = 0;
        input_avstream_audio_resampler_channels = -1;
        input_avstream_audio_resampler_rate = -1;
        eof_stream = false;
        got_audio = false;
        got_video = false;
        adj_time = 0;
        t = pt = -1;
        eof = false;
        avpkt_init();
        next_pts = next_dts = -1LL;
        return (input_avfmt != NULL);
    }
    bool next_packet(void) {
        if (eof) return false;
        if (input_avfmt == NULL) return false;

        do {
            if (eof_stream) break;
            avpkt_release();
            avpkt_init();
            if (av_read_frame(input_avfmt,&avpkt) < 0) {
                eof_stream = true;
                return false;
            }
            if (avpkt.stream_index >= input_avfmt->nb_streams)
                continue;

            /* track time and keep things monotonic for our code */
            if (avpkt.pts != AV_NOPTS_VALUE) {
                t = avpkt.pts * av_q2d(input_avfmt->streams[avpkt.stream_index]->time_base);

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

            if (avpkt.pts != AV_NOPTS_VALUE) {
                avpkt.pts += (adj_time * input_avfmt->streams[avpkt.stream_index]->time_base.den) /
                    input_avfmt->streams[avpkt.stream_index]->time_base.num;
            }

            if (avpkt.dts != AV_NOPTS_VALUE) {
                avpkt.dts += (adj_time * input_avfmt->streams[avpkt.stream_index]->time_base.den) /
                    input_avfmt->streams[avpkt.stream_index]->time_base.num;
            }

            got_audio = false;
            got_video = false;
			if (input_avstream_audio != NULL && avpkt.stream_index == input_avstream_audio->index) {
                if (got_audio) fprintf(stderr,"Audio content lost\n");
				av_packet_rescale_ts(&avpkt,input_avstream_audio->time_base,output_avstream_audio->time_base);
                handle_audio(/*&*/avpkt);
                got_audio = true;
                break;
			}
			else if (input_avstream_video != NULL && avpkt.stream_index == input_avstream_video->index) {
                if (got_video) fprintf(stderr,"Video content lost\n");
				AVRational m = (AVRational){output_field_rate.den, output_field_rate.num};
				av_packet_rescale_ts(&avpkt,input_avstream_video->time_base,m); // convert to FIELD number

                // ugh... this can happen if the source is an AVI file
                if (avpkt.pts == AV_NOPTS_VALUE) avpkt.pts = avpkt.dts;

                handle_frame(/*&*/avpkt); // will set got_video
                break;
			}

            avpkt_release();
        } while (1);

        if (eof_stream) {
            avpkt_release();
            avpkt.size = 0;
            avpkt.data = NULL;
            handle_frame(/*&*/avpkt); // will set got_video
            if (!got_video) eof = true;
            else fprintf(stderr,"Got latent frame\n");
        }

        return true;
    }
    void handle_audio(AVPacket &pkt) {
        int got_frame = 0;

        if (avcodec_decode_audio4(input_avstream_audio_codec_context,input_avstream_audio_frame,&got_frame,&pkt) >= 0) {
            if (got_frame != 0 && input_avstream_audio_frame->nb_samples != 0) {
                if (input_avstream_audio_frame->pts == AV_NOPTS_VALUE)
                    input_avstream_audio_frame->pts = pkt.pts;

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
                        return;
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

                if (audio_dst_data != NULL) {
                    if ((audio_dst_data_out_samples=swr_convert(input_avstream_audio_resampler,audio_dst_data,audio_dst_data_samples,
                                    (const uint8_t**)input_avstream_audio_frame->data,input_avstream_audio_frame->nb_samples)) > 0) {
                    }
                    else if (audio_dst_data_out_samples < 0) {
                        fprintf(stderr,"Failed to resample audio\n");
                    }

                    audio_dst_data_out_audio_sample = audio_sample;
                    audio_sample += audio_dst_data_out_samples;
                }
            }
        }
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
                    input_avstream_video_frame_rgb->width,
                    input_avstream_video_frame_rgb->height,
                    (AVPixelFormat)input_avstream_video_frame_rgb->format,
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
            input_avstream_video_frame_rgb->pts = input_avstream_video_frame->pts;
            input_avstream_video_frame_rgb->pkt_pts = input_avstream_video_frame->pkt_pts;
            input_avstream_video_frame_rgb->pkt_dts = input_avstream_video_frame->pkt_dts;
            input_avstream_video_frame_rgb->top_field_first = input_avstream_video_frame->top_field_first;
            input_avstream_video_frame_rgb->interlaced_frame = input_avstream_video_frame->interlaced_frame;

            if (sws_scale(input_avstream_video_resampler,
                        // source
                        input_avstream_video_frame->data,
                        input_avstream_video_frame->linesize,
                        0,input_avstream_video_frame->height,
                        // dest
                        input_avstream_video_frame_rgb->data,
                        input_avstream_video_frame_rgb->linesize) <= 0)
                fprintf(stderr,"WARNING: sws_scale failed\n");
        }
    }
    void handle_frame(AVPacket &pkt) {
        int got_frame = 0;

        if (avcodec_decode_video2(input_avstream_video_codec_context,input_avstream_video_frame,&got_frame,&pkt) >= 0) {
            if (got_frame != 0 && input_avstream_video_frame->width > 0 && input_avstream_video_frame->height > 0) {
                got_video = true;
            }
        }
        else {
            fprintf(stderr,"No video decoded\n");
        }
    }
    void avpkt_init(void) {
        if (!avpkt_valid) {
            avpkt_valid = true;
            av_init_packet(&avpkt);
        }
    }
    void avpkt_release(void) {
        if (avpkt_valid) {
            avpkt_valid = false;
            av_packet_unref(&avpkt);
        }
        got_audio = false;
        got_video = false;
    }
    void close_input(void) {
        eof = true;
        avpkt_release();
        if (input_avstream_audio_codec_context != NULL) {
            avcodec_close(input_avstream_audio_codec_context);
            input_avstream_audio_codec_context = NULL;
            input_avstream_audio = NULL;
        }
        if (input_avstream_video_codec_context != NULL) {
            avcodec_close(input_avstream_video_codec_context);
            input_avstream_video_codec_context = NULL;
            input_avstream_video = NULL;
        }

        if (input_avstream_audio_frame != NULL)
            av_frame_free(&input_avstream_audio_frame);
        if (input_avstream_video_frame != NULL)
            av_frame_free(&input_avstream_video_frame);
        if (input_avstream_video_frame_rgb != NULL)
            av_frame_free(&input_avstream_video_frame_rgb);

        if (input_avstream_audio_resampler != NULL)
            swr_free(&input_avstream_audio_resampler);
        if (input_avstream_video_resampler != NULL) {
            sws_freeContext(input_avstream_video_resampler);
            input_avstream_video_resampler = NULL;
        }

		if (audio_dst_data != NULL) {
			av_freep(&audio_dst_data[0]); // NTS: Why??
			av_freep(&audio_dst_data);
		}

        input_avstream_audio_resampler_channels = -1;
        input_avstream_audio_resampler_rate = -1;
        avformat_close_input(&input_avfmt);
    }
public:
    std::string             path;
    uint32_t                color;
    int                     newlevel;
    bool                    eof;
    bool                    eof_stream;
    bool                    got_audio;
    bool                    got_video;
public:
    unsigned long long      last_written_sample;
    unsigned long long      audio_sample;
    uint8_t**               audio_dst_data;
    int                     audio_dst_data_alloc_samples;
    int                     audio_dst_data_linesize;
    int                     audio_dst_data_samples;
    int                     audio_dst_data_out_samples;
    unsigned long long      audio_dst_data_out_audio_sample;
    int                     input_avstream_audio_resampler_rate;
    int                     input_avstream_audio_resampler_channels;
    AVFormatContext*        input_avfmt;
    AVStream*               input_avstream_audio;	            // do not free
    AVCodecContext*         input_avstream_audio_codec_context; // do not free
    AVFrame*                input_avstream_audio_frame;
    AVStream*               input_avstream_video;	            // do not free
    AVCodecContext*         input_avstream_video_codec_context; // do not free
    AVFrame*		        input_avstream_video_frame;
    AVFrame*		        input_avstream_video_frame_rgb;
    struct SwrContext*      input_avstream_audio_resampler;
    struct SwsContext*	    input_avstream_video_resampler;
    AVPixelFormat           input_avstream_video_resampler_format;
    int                     input_avstream_video_resampler_height;
    int                     input_avstream_video_resampler_width;
    signed long long        next_pts;
    signed long long        next_dts;
    AVPacket                avpkt;
    bool                    avpkt_valid;
    double                  adj_time;
    double                  t,pt;
};

std::vector<InputFile>      input_files;
std::string                 output_file;

InputFile &current_input_file(void) {
    if (input_files.empty()) {
        std::string what = "input files empty";
        throw std::out_of_range(/*&*/what);
    }

    return *(input_files.rbegin()); /* last one */
}

InputFile &new_input_file(void) {
    if (!input_files.empty()) {
        /* copy the last one, except for some fields */
        {
            InputFile &last = current_input_file();
            input_files.push_back(last);
        }
        {
            InputFile &last = current_input_file();
            last.reset_on_dup();
        }
    }
    else {
        /* make a new one with defaults */
        input_files.push_back(InputFile());
    }

    return current_input_file();
}

volatile int DIE = 0;

void sigma(int x) {
	if (++DIE >= 20) abort();
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
	fprintf(stderr," -i <input file>               you can specify more than one input file, in order of layering\n");
	fprintf(stderr," -o <output file>\n");
    fprintf(stderr," -n <n>                        New content averaging level (256=100% 0=0%)\n");
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
            else if (!strcmp(a,"n")) {
                a = argv[i++];
                if (a == NULL) return 1;
                current_input_file().newlevel = (int)strtoul(a,NULL,0);
            }
            else if (!strcmp(a,"i")) {
                a = argv[i++];
                if (a == NULL) return 1;
                new_input_file().path = a;
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

    if (output_file.empty()) {
        fprintf(stderr,"No output file specified\n");
        return 1;
    }
    if (input_files.empty()) {
        fprintf(stderr,"No input files specified\n");
        return 1;
    }

	return 0;
}

void process_audio(InputFile &fin) {
    if (fin.audio_dst_data == NULL || fin.audio_dst_data_out_samples == 0)
        return;

    // we don't do anything with audio in this filter
}

void write_out_audio(InputFile &fin) {
    if (fin.audio_dst_data == NULL || fin.audio_dst_data_out_samples == 0)
        return;

    /* pad-fill */
    while (fin.last_written_sample < fin.audio_dst_data_out_audio_sample) {
        unsigned long long out_samples = fin.audio_dst_data_out_audio_sample - fin.last_written_sample;

        if (out_samples > output_audio_rate)
            out_samples = output_audio_rate;

        AVPacket dstpkt;
        av_init_packet(&dstpkt);
        if (av_new_packet(&dstpkt,out_samples * 2 * output_audio_channels) >= 0) { // NTS: Will reset fields too!
            assert(dstpkt.data != NULL);
            assert(dstpkt.size >= (out_samples * 2 * output_audio_channels));
            memset(dstpkt.data,0,out_samples * 2 * output_audio_channels);
        }
        dstpkt.pts = fin.last_written_sample;
        dstpkt.dts = fin.last_written_sample;
        dstpkt.stream_index = output_avstream_audio->index;
        av_packet_rescale_ts(&dstpkt,output_avstream_audio_codec_context->time_base,output_avstream_audio->time_base);
        if (av_interleaved_write_frame(output_avfmt,&dstpkt) < 0)
            fprintf(stderr,"Failed to write frame\n");
        av_packet_unref(&dstpkt);

        fprintf(stderr,"Pad fill %llu samples\n",out_samples);
        fin.last_written_sample += out_samples;
    }

    // write it out. TODO: At some point, support conversion to whatever the codec needs and then convert to it.
    // that way we can render directly to MP4 our VHS emulation.
    AVPacket dstpkt;
    av_init_packet(&dstpkt);
    if (av_new_packet(&dstpkt,fin.audio_dst_data_out_samples * 2 * output_audio_channels) >= 0) { // NTS: Will reset fields too!
        assert(dstpkt.data != NULL);
        assert(dstpkt.size >= (fin.audio_dst_data_out_samples * 2 * output_audio_channels));
        memcpy(dstpkt.data,fin.audio_dst_data[0],fin.audio_dst_data_out_samples * 2 * output_audio_channels);
    }
    dstpkt.pts = fin.audio_dst_data_out_audio_sample;
    dstpkt.dts = fin.audio_dst_data_out_audio_sample;
    dstpkt.stream_index = output_avstream_audio->index;
    av_packet_rescale_ts(&dstpkt,output_avstream_audio_codec_context->time_base,output_avstream_audio->time_base);
    if (av_interleaved_write_frame(output_avfmt,&dstpkt) < 0)
        fprintf(stderr,"Failed to write frame\n");
    av_packet_unref(&dstpkt);

    fin.audio_sample = fin.last_written_sample = fin.audio_dst_data_out_audio_sample + fin.audio_dst_data_out_samples;
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
void composite_layer(AVFrame *dstframe,AVFrame *srcframe,InputFile &inputfile,unsigned long long field) {
    uint32_t *dscan,*sscan;
    unsigned int r,g,b;
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
            r  = (((*sscan >> 16UL) & 0xFF) * inputfile.newlevel);
            r += (((*dscan >> 16UL) & 0xFF) * (256 - inputfile.newlevel));
            r += ((((x^y)+field)&3)*255)/3;
            r >>= 8;

            g  = (((*sscan >>  8UL) & 0xFF) * inputfile.newlevel);
            g += (((*dscan >>  8UL) & 0xFF) * (256 - inputfile.newlevel));
            g += ((((x^y)+field)&3)*255)/3;
            g >>= 8;

            b  = (((*sscan >>  0UL) & 0xFF) * inputfile.newlevel);
            b += (((*dscan >>  0UL) & 0xFF) * (256 - inputfile.newlevel));
            b += ((((x^y)+field)&3)*255)/3;
            b >>= 8;

            *dscan = (r << 16) + (g << 8) + b;
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

    /* open all input files */
    for (std::vector<InputFile>::iterator i=input_files.begin();i!=input_files.end();i++) {
        if (!(*i).open_input()) {
            fprintf(stderr,"Failed to open %s\n",(*i).path.c_str());
            return 1;
        }
    }

    /* open output file */
	assert(output_avfmt == NULL);
	if (avformat_alloc_output_context2(&output_avfmt,NULL,NULL,output_file.c_str()) < 0) {
		fprintf(stderr,"Failed to open output file\n");
		return 1;
	}

	{
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
        signed long long upto=0;
        signed long long current=0;

        do {
            if (DIE) break;

            eof = true;
            copyaud = false;
            for (std::vector<InputFile>::iterator i=input_files.begin();i!=input_files.end();i++) {
                if ((*i).eof == false) {
                    eof = false;
                    if (!((*i).got_audio) && !((*i).got_video))
                        (*i).next_packet();

                    if ((*i).got_audio) {
                        /* we don't do anything with audio, but we do copy through the first input file's audio */
                        if (!copyaud) {
                            copyaud = true;
                            process_audio((*i)); // NTS: We do chroma/color keying, we don't do anything with audio
                            write_out_audio((*i));
                        }
                        (*i).got_audio = false;
                    }
                }
            }

            upto = -1LL;
            for (std::vector<InputFile>::iterator i=input_files.begin();i!=input_files.end();i++) {
                if ((*i).eof == false) {
                    if ((*i).input_avstream_video_frame != NULL) {
                        if ((*i).got_video) {
                            if ((*i).input_avstream_video_frame->pkt_pts != AV_NOPTS_VALUE) {
                                if (upto == (-1LL) || upto > (*i).input_avstream_video_frame->pkt_pts)
                                    upto = (*i).input_avstream_video_frame->pkt_pts;
                            }

                            if ((*i).input_avstream_video_frame->pkt_pts == AV_NOPTS_VALUE || current >= (*i).input_avstream_video_frame->pkt_pts) {
                                (*i).frame_copy_scale();
                                (*i).got_video = false;
                            }
                        }
                        else {
                            (*i).got_video = false;
                            upto = current;
                        }
                    }
                    else {
                        (*i).got_video = false;
                        upto = current;
                    }
                }
                else {
                    if ((*i).got_video) {
                        (*i).frame_copy_scale();
                        (*i).got_video = false;
                    }
                }
            }

            while (current < upto) {
                for (std::vector<InputFile>::iterator i=input_files.begin();i!=input_files.end();i++) {
                    if ((*i).eof == false) {
                        if ((*i).input_avstream_video_frame != NULL) {
                            if ((*i).got_video) {
                                if ((*i).input_avstream_video_frame->pkt_pts == AV_NOPTS_VALUE || current >= (*i).input_avstream_video_frame->pkt_pts) {
                                    (*i).frame_copy_scale();
                                    (*i).got_video = false;
                                }
                            }
                            else {
                                (*i).got_video = false;
                            }
                        }
                        else {
                            (*i).got_video = false;
                            upto = current;
                        }
                    }
                    else {
                        if ((*i).got_video) {
                            (*i).frame_copy_scale();
                            (*i).got_video = false;
                        }
                    }

                    // composite the layer, keying against the color. all code assumes ARGB
                    composite_layer(output_avstream_video_frame,(*i).input_avstream_video_frame_rgb,*i,current);
                }

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
            }
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

    /* close all */
    for (std::vector<InputFile>::iterator i=input_files.begin();i!=input_files.end();i++)
        (*i).close_input();

	return 0;
}

