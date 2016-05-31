#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.59]) # Sorry troy
AC_INIT([composite-video-simulator], [1.0], [jonathan@castus.tv])
AM_INIT_AUTOMAKE([composite-video-simulator],[1.0])
AC_CONFIG_SRCDIR([configure.ac])
AC_CONFIG_HEADER([config.h])
AC_CONFIG_MACRO_DIR([m4])
AC_PREFIX_DEFAULT([/usr])

# Checks for programs.
AC_PROG_CC
AC_PROG_CXX
AM_PROG_CC_C_O

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([fcntl.h stdint.h stdlib.h string.h sys/socket.h unistd.h])

# we require pkg-config
PKG_PROG_PKG_CONFIG

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_INT64_T
AC_TYPE_UINT32_T
AC_TYPE_UINT8_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_MEMCMP
AC_CHECK_FUNCS([memset])

# FFMPEG
ISP_PKG_CHECK(AVCODEC,				libavcodec >= 55.39.101)
ISP_PKG_CHECK(AVFORMAT,				libavformat >= 55.19.104)
ISP_PKG_CHECK(SWSCALE,				libswscale >= 4.0.100)
ISP_PKG_CHECK(AVUTIL,				libavutil >= 52.48.101)
ISP_PKG_CHECK(SWRESAMPLE,			libswresample >= 2.0.101)

# variables for multi-target
AM_CONDITIONAL(WIN32,false)
AM_CONDITIONAL(LINUX,true)

AC_CONFIG_FILES([Makefile])
AC_OUTPUT

# where are we?
top="`pwd`"

