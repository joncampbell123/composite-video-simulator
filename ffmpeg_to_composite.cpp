
#define __STDC_CONSTANT_MACROS

#include <sys/types.h>
#include <stdint.h>
#include <assert.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <fcntl.h>

extern "C" {
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

string		input_file;
string		output_file;

AVFormatContext*	input_avfmt = NULL;
AVFormatContext*	output_avfmt = NULL;

AVRational	output_field_rate = { 60000, 1001 };	// NTSC 60Hz default
int		output_width = 720;
int		output_height = 480;
bool		output_ntsc = true;	// NTSC color subcarrier emulation
bool		output_pal = false;	// PAL color subcarrier emulation

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
			else if (!strcmp(a,"i")) {
				input_file = argv[i++];
			}
			else if (!strcmp(a,"o")) {
				output_file = argv[i++];
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

	assert(output_avfmt == NULL);
	if (avformat_alloc_output_context2(&output_avfmt,NULL,NULL,output_file.c_str()) < 0) {
		fprintf(stderr,"Failed to open output file\n");
		return 1;
	}

	if (avformat_write_header(output_avfmt,NULL) < 0) {
		fprintf(stderr,"Failed to write header\n");
		return 1;
	}

	av_write_trailer(output_avfmt);
	avformat_free_context(output_avfmt);
	av_free(output_avfmt);
	avformat_close_input(&input_avfmt);
	return 0;
}

