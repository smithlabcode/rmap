ifndef SRC_ROOT
SRC_ROOT=../..
endif

ifndef SMITHLAB_CPP
SMITHLAB_CPP=../../smithlab_cpp
endif

PROGS = build

CXX = g++
CFLAGS = -Wall -fPIC -fmessage-length=50
OPTFLAGS = -O3
DEBUGFLAGS = -g
COMMON_DIR = ../../common
TEST_DIR = $(SRC_ROOT)/test

ifeq "$(shell uname)" "Darwin"
CFLAGS += -arch x86_64
endif

ifdef DEBUG
CFLAGS += $(DEBUGFLAGS)
endif

ifdef OPT
CFLAGS += $(OPTFLAGS)
endif

SRC = reference.cpp index.cpp
OBJ = $(patsubst %.cpp,%.o,$(SRC))

all:	$(PROGS)

%.o: %.cpp %.hpp
	$(CXX) $(CFLAGS) -c -o $@ $< -I$(SMITHLAB_CPP)

build : $(OBJ) $(addprefix $(SMITHLAB_CPP)/, GenomicRegion.o smithlab_os.o smithlab_utils.o \
		OptionParser.o QualityScore.o)


install: all
	@mkdir -p $(SRC_ROOT)/bin
	@install -m 755 $(PROGS) $(SRC_ROOT)/bin

%: %.cpp
	$(CXX) $(CFLAGS) -o $@ $^ -I$(SMITHLAB_CPP)

test_%:	%
	@$(TEST_DIR)/$@ $(TEST_DIR)

test:	$(addprefix test_, $(PROGS))

clean:
	@-rm -f $(PROGS) *.o *.so *.a *~

.PHONY: clean
