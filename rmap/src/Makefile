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

app_subdirs=mappers simulation

all:
	@make -C common RMAP=$(RMAP) OPT=1
	@for i in $(app_subdirs); do \
		make -C $${i} RMAP=$(RMAP) OPT=1; \
	done;

install:
	@make -C common RMAP=$(RMAP) OPT=1
	@for i in $(app_subdirs); do \
		make -C $${i} RMAP=$(RMAP) OPT=1 install; \
	done;

test:
	@for i in $(app_subdirs); do \
		make -C $${i} RMAP=$(RMAP) test; \
	done;

clean:
	@make -C common RMAP=$(RMAP) OPT=1 clean
	@for i in $(app_subdirs); do \
		make -C $${i} RMAP=$(RMAP) OPT=1 clean; \
	done;
