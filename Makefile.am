SUBDIRS =
ACLOCAL_AMFLAGS = -I m4
top = $(shell pwd)

bin_PROGRAMS = \
	ffmpeg_to_composite \
	ffmpeg_colorkey \
	ffmpeg_posterize

ffmpeg_to_composite_SOURCES = ffmpeg_to_composite.cpp
ffmpeg_to_composite_CXXFLAGS = $(AVCODEC_CFLAGS) $(AVFORMAT_CFLAGS) $(AVUTIL_CFLAGS) $(SWSCALE_CFLAGS) $(SWRESAMPLE_CFLAGS)
ffmpeg_to_composite_LDADD = $(AVCODEC_LIBS) $(AVFORMAT_LIBS) $(AVUTIL_LIBS) $(SWSCALE_LIBS) $(SWRESAMPLE_LIBS)

ffmpeg_colorkey_SOURCES = ffmpeg_colorkey.cpp
ffmpeg_colorkey_CXXFLAGS = $(AVCODEC_CFLAGS) $(AVFORMAT_CFLAGS) $(AVUTIL_CFLAGS) $(SWSCALE_CFLAGS) $(SWRESAMPLE_CFLAGS)
ffmpeg_colorkey_LDADD = $(AVCODEC_LIBS) $(AVFORMAT_LIBS) $(AVUTIL_LIBS) $(SWSCALE_LIBS) $(SWRESAMPLE_LIBS)

ffmpeg_posterize_SOURCES = ffmpeg_posterize.cpp
ffmpeg_posterize_CXXFLAGS = $(AVCODEC_CFLAGS) $(AVFORMAT_CFLAGS) $(AVUTIL_CFLAGS) $(SWSCALE_CFLAGS) $(SWRESAMPLE_CFLAGS)
ffmpeg_posterize_LDADD = $(AVCODEC_LIBS) $(AVFORMAT_LIBS) $(AVUTIL_LIBS) $(SWSCALE_LIBS) $(SWRESAMPLE_LIBS)

