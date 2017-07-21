# SVM with stochastic gradient

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111, USA

#define HINGELOSS 1
#define SMOOTHHINGELOSS 2
#define SQUAREDHINGELOSS 3
#define LOGLOSS 10
#define LOGLOSSMARGIN 11

CXX=g++
OPTS=-g -O3 -Wno-deprecated -static -DLOSS=1 #-DDEBUGON
#OPTS=-g -Wno-deprecated -DLOSS=1 #-DDEBUGON
#OPTS=-g -Wno-deprecated -static -pg
CXXFLAGS= ${OPTS} -Wall
LIBS = -lz -lm

PROGRAMS = svmsgdnspdk 

all: ${PROGRAMS}

clean:
	-rm ${PROGRAMS} 2>/dev/null
	-rm *.o 2>/dev/null

svmsgdnspdk: svmsgdnspdk.o vectors.o gzstream.o GraphClass.o Histogram.o NSPDK_FeatureGenerator.o vectors.o FlagsService.o 
	${CXX} ${CXXFLAGS} -o $@ svmsgdnspdk.o vectors.o gzstream.o BaseGraphClass.o GraphClass.o Histogram.o NSPDK_FeatureGenerator.o FlagsService.o ${LIBS}

svmsgdnspdk.o: svmsgdnspdk.cpp vectors.h  gzstream.h 
	${CXX} ${CXXFLAGS} -c -o $@ svmsgdnspdk.cpp

vectors.o: vectors.cpp vectors.h  wrapper.h
	${CXX} ${CXXFLAGS} -c -o $@ vectors.cpp

gzstream.o: gzstream.cpp gzstream.h 
	${CXX} ${CXXFLAGS} -c -o $@ gzstream.cpp

FlagsService.o:FlagsService.h

NSPDK_FeatureGenerator.o:NSPDK_FeatureGenerator.h GraphClass.o 

GraphClass.o:GraphClass.h BaseGraphClass.o 

BaseGraphClass.o:BaseGraphClass.h 

Histogram.o:Histogram.h Utility.h
