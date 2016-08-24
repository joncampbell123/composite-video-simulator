SUBDIRS =
ACLOCAL_AMFLAGS = -I m4
top = $(shell pwd)

bin_PROGRAMS = \
	ffmpeg_to_composite \
	ffmpeg_colorkey \
	ffmpeg_posterize \
	ffmpeg_average_delay \
	ffmpeg_colormap \
	ffmpeg_scanimate \
	ffmpeg_ntsc

ffmpeg_to_composite_SOURCES = ffmpeg_to_composite.cpp
ffmpeg_to_composite_CXXFLAGS = $(AVCODEC_CFLAGS) $(AVFORMAT_CFLAGS) $(AVUTIL_CFLAGS) $(SWSCALE_CFLAGS) $(SWRESAMPLE_CFLAGS)
ffmpeg_to_composite_LDADD = $(AVCODEC_LIBS) $(AVFORMAT_LIBS) $(AVUTIL_LIBS) $(SWSCALE_LIBS) $(SWRESAMPLE_LIBS)

ffmpeg_scanimate_SOURCES = ffmpeg_scanimate.cpp
ffmpeg_scanimate_CXXFLAGS = $(AVCODEC_CFLAGS) $(AVFORMAT_CFLAGS) $(AVUTIL_CFLAGS) $(SWSCALE_CFLAGS) $(SWRESAMPLE_CFLAGS)
ffmpeg_scanimate_LDADD = $(AVCODEC_LIBS) $(AVFORMAT_LIBS) $(AVUTIL_LIBS) $(SWSCALE_LIBS) $(SWRESAMPLE_LIBS)

ffmpeg_colorkey_SOURCES = ffmpeg_colorkey.cpp
ffmpeg_colorkey_CXXFLAGS = $(AVCODEC_CFLAGS) $(AVFORMAT_CFLAGS) $(AVUTIL_CFLAGS) $(SWSCALE_CFLAGS) $(SWRESAMPLE_CFLAGS)
ffmpeg_colorkey_LDADD = $(AVCODEC_LIBS) $(AVFORMAT_LIBS) $(AVUTIL_LIBS) $(SWSCALE_LIBS) $(SWRESAMPLE_LIBS)

ffmpeg_posterize_SOURCES = ffmpeg_posterize.cpp
ffmpeg_posterize_CXXFLAGS = $(AVCODEC_CFLAGS) $(AVFORMAT_CFLAGS) $(AVUTIL_CFLAGS) $(SWSCALE_CFLAGS) $(SWRESAMPLE_CFLAGS)
ffmpeg_posterize_LDADD = $(AVCODEC_LIBS) $(AVFORMAT_LIBS) $(AVUTIL_LIBS) $(SWSCALE_LIBS) $(SWRESAMPLE_LIBS)

ffmpeg_average_delay_SOURCES = ffmpeg_average_delay.cpp
ffmpeg_average_delay_CXXFLAGS = $(AVCODEC_CFLAGS) $(AVFORMAT_CFLAGS) $(AVUTIL_CFLAGS) $(SWSCALE_CFLAGS) $(SWRESAMPLE_CFLAGS)
ffmpeg_average_delay_LDADD = $(AVCODEC_LIBS) $(AVFORMAT_LIBS) $(AVUTIL_LIBS) $(SWSCALE_LIBS) $(SWRESAMPLE_LIBS)

ffmpeg_colormap_SOURCES = ffmpeg_colormap.cpp
ffmpeg_colormap_CXXFLAGS = $(AVCODEC_CFLAGS) $(AVFORMAT_CFLAGS) $(AVUTIL_CFLAGS) $(SWSCALE_CFLAGS) $(SWRESAMPLE_CFLAGS)
ffmpeg_colormap_LDADD = $(AVCODEC_LIBS) $(AVFORMAT_LIBS) $(AVUTIL_LIBS) $(SWSCALE_LIBS) $(SWRESAMPLE_LIBS)

ffmpeg_ntsc_SOURCES = ffmpeg_ntsc.cpp
ffmpeg_ntsc_CXXFLAGS = $(AVCODEC_CFLAGS) $(AVFORMAT_CFLAGS) $(AVUTIL_CFLAGS) $(SWSCALE_CFLAGS) $(SWRESAMPLE_CFLAGS)
ffmpeg_ntsc_LDADD = $(AVCODEC_LIBS) $(AVFORMAT_LIBS) $(AVUTIL_LIBS) $(SWSCALE_LIBS) $(SWRESAMPLE_LIBS)

