#    Makefile from rmap software
#
#    Copyright (C) 2008 Cold Spring Harbor Laboratory, 
#                       University of Southern California and
#                       Andrew D. Smith
#
#    Authors: Andrew D. Smith
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

ifndef RMAP
$(error Must define RMAP variable)
endif

HEADERS = $(shell ls *.hpp)
SOURCES = $(shell ls *.cpp)

OBJS = $(subst .cpp,.o,$(SOURCES))

CXX = g++
CFLAGS = -Wall -fPIC -fmessage-length=50
OPTFLAGS = -O3
DEBUGFLAGS = -g

ifeq "$(shell uname)" "Darwin"
CFLAGS += -arch x86_64
endif

ifdef DEBUG
CFLAGS += $(DEBUGFLAGS)
endif

ifdef OPT
CFLAGS += $(OPTFLAGS)
endif

all: $(OBJS)

%.o: %.cpp %.hpp
	$(CXX) $(CFLAGS) -c -o $@ $<

clean:
	@-rm -f $(PROGS) *.o *.so *.a *~

.PHONY: clean
