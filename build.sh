#!/bin/sh
#########################################################################
#    Copyright 2020 Anton Douginets
#    This file is part of solensim.
#
#    solensim is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    solensim is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with solensim.  If not, see <https://www.gnu.org/licenses/>.
#########################################################################

# Release .zip "build" script

# clear previous zip
rm solensim.zip

# clean up
rm -rf src/__pycache__
rm -rf src/mcode/__pycache__
rm -rf src/pycode/__pycache__
rm src/octave-workspace
rm src/mcode/octave-workspace

# pack .zip
mkdir solensim
cp -r src/* solensim/
cp LICENSE solensim/
cp readme.txt solensim/

# pack docs into zip
#cp tex/doc/doc.pdf solensim/solensim.pdf

# pack matlab code separately, for now:
mv solensim/mcode/Magnetic_Bodge.m solensim/m_script.m
rm -rf solensim/mcode

zip -r solensim.zip solensim

# clean tmp build dir
rm -rf solensim
