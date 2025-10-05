#!/bin/bash

for dir in */
do
   cd "$dir"
   touch configure.ac aclocal.m4 configure Makefile.am Makefile.in
   make
   cd ..
done

./dune-common/bin/dunecontrol --opts=configure.opts all
