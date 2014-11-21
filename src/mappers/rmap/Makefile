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

ifndef SRC_ROOT
SRC_ROOT=../..
endif

ifndef SMITHLAB_CPP
SMITHLAB_CPP=../../smithlab_cpp
endif

PROGS = rmap rmap-pe

CXX = g++
CXXFLAGS = -Wall -fPIC -fmessage-length=50
OPTFLAGS = -O3
DEBUGFLAGS = -g
COMMON_DIR = ../../common
TEST_DIR = $(SRC_ROOT)/test


# Flags for Mavericks
ifeq "$(shell uname)" "Darwin"
CXXFLAGS += -arch x86_64
ifeq "$(shell if [ `sysctl -n kern.osrelease | cut -d . -f 1` -ge 13 ];\
              then echo 'true'; fi)" "true"
CXXFLAGS += -stdlib=libstdc++
endif
endif

ifdef DEBUG
CXXFLAGS += $(DEBUGFLAGS)
endif

ifdef OPT
CXXFLAGS += $(OPTFLAGS)
endif


all:	$(PROGS)

%.o: %.cpp %.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ $< -I$(COMMON_DIR) -I$(SMITHLAB_CPP)

%: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^ -I$(COMMON_DIR) -I$(SMITHLAB_CPP)

rmap rmap-pe: $(addprefix $(COMMON_DIR)/, SeedMaker.o FastRead.o load_reads.o \
	clip_adaptor_from_reads.o) \
	$(addprefix $(SMITHLAB_CPP)/, GenomicRegion.o smithlab_os.o smithlab_utils.o \
	OptionParser.o QualityScore.o)

rmap-pe: $(addprefix $(SMITHLAB_CPP)/, MappedRead.o)

install: all
	@mkdir -p $(SRC_ROOT)/bin
	@install -m 755 $(PROGS) $(SRC_ROOT)/bin

test_%:	%
	@$(TEST_DIR)/$@ $(TEST_DIR)

test:	$(addprefix test_, $(PROGS))

clean:
	@-rm -f $(PROGS) *.o *.so *.a *~

.PHONY: clean
