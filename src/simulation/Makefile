#    Makefile from rmap software
#
#    Copyright (C) 2009 University of Southern California and
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

ifndef SMITHLAB_CPP
SMITHLAB_CPP=../smithlab_cpp/
endif

app_subdirs=simreads simreadsbs simreadspe

all:
	@for i in $(app_subdirs); do \
		make -C $${i} RMAP=$(RMAP) SMITHLAB_CPP=$(SMITHLAB_CPP) OPT=1; \
	done;

install:
	@for i in $(app_subdirs); do \
		make -C $${i} RMAP=$(RMAP) SMITHLAB_CPP=$(SMITHLAB_CPP) OPT=1 install; \
	done;

test:
	@for i in $(app_subdirs); do \
		make -C $${i} RMAP=$(RMAP) SMITHLAB_CPP=$(SMITHLAB_CPP) test; \
	done;

clean:
	@for i in $(app_subdirs); do \
		make -C $${i} RMAP=$(RMAP) SMITHLAB_CPP=$(SMITHLAB_CPP) clean; \
	done;
