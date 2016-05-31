SUBDIRS =
ACLOCAL_AMFLAGS = -I m4
top = $(shell pwd)

bin_PROGRAMS = \
	ffmpeg_to_composite

ffmpeg_to_composite_SOURCES = ffmpeg_to_composite.cpp
ffmpeg_to_composite_CXXFLAGS = $(AVCODEC_CFLAGS) $(AVFORMAT_CFLAGS) $(AVUTIL_CFLAGS) $(SWSCALE_CFLAGS) $(SWRESAMPLE_CFLAGS)
ffmpeg_to_composite_LDADD = $(AVCODEC_LIBS) $(AVFORMAT_LIBS) $(AVUTIL_LIBS) $(SWSCALE_LIBS) $(SWRESAMPLE_LIBS)

