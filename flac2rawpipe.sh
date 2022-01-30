#!/bin/bash
ffmpeg -i "$1" -ss "$2" -acodec pcm_u8 -f u8 - 2>/dev/null </dev/null

