#!/bin/bash
# OK do it
echo >NEWS
echo >AUTHORS
echo >ChangeLog
mkdir -p m4 || exit 1
(aclocal && autoheader && automake --add-missing && autoconf) || exit 1
mkdir -p m4 || exit 1

