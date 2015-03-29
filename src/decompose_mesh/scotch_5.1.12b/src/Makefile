## Copyright 2004,2007-2011 ENSEIRB, INRIA & CNRS
##
## This file is part of the Scotch software package for static mapping,
## graph partitioning and sparse matrix ordering.
##
## This software is governed by the CeCILL-C license under French law
## and abiding by the rules of distribution of free software. You can
## use, modify and/or redistribute the software under the terms of the
## CeCILL-C license as circulated by CEA, CNRS and INRIA at the following
## URL: "http://www.cecill.info".
## 
## As a counterpart to the access to the source code and rights to copy,
## modify and redistribute granted by the license, users are provided
## only with a limited warranty and the software's author, the holder of
## the economic rights, and the successive licensors have only limited
## liability.
## 
## In this respect, the user's attention is drawn to the risks associated
## with loading, using, modifying and/or developing or reproducing the
## software by the user in light of its specific status of free software,
## that may mean that it is complicated to manipulate, and that also
## therefore means that it is reserved for developers and experienced
## professionals having in-depth computer knowledge. Users are therefore
## encouraged to load and test the software's suitability as regards
## their requirements in conditions enabling the security of their
## systems and/or data to be ensured and, more generally, to use and
## operate it in the same conditions as regards security.
## 
## The fact that you are presently reading this means that you have had
## knowledge of the CeCILL-C license and that you accept its terms.
##

VERSION	= 5
RELEASE = 1
PATCHLEVEL = 12
export VERSION RELEASE PATCHLEVEL

.PHONY				:	clean	default	install	realclean	required	scotch

default				:	scotch

required			:	Makefile.inc	../bin	../include	../lib

include Makefile.inc

prefix		?= /usr/local
bindir		?= $(prefix)/bin
includedir	?= $(prefix)/include
libdir		?= $(prefix)/lib
datarootdir	?= $(prefix)/share
mandir		?= $(datarootdir)/man

../bin				:
					-$(MKDIR) ../bin

../include			:
					-$(MKDIR) ../include

../lib				:
					-$(MKDIR) ../lib

$(bindir)			:
					-$(MKDIR) $(bindir)

$(datarootdir)			:
					-$(MKDIR) $(datarootdir)

$(includedir)			:
					-$(MKDIR) $(includedir)

$(libdir)			:
					-$(MKDIR) $(libdir)

$(mandir)			:	$(datarootdir)
					-$(MKDIR) $(mandir)

$(mandir)/man1			:	$(mandir)
					-$(MKDIR) $(mandir)/man1

scotch				:	required
					(cd libscotch ;      $(MAKE) scotch ; $(MAKE) install)
					(cd scotch ;         $(MAKE) scotch ; $(MAKE) install)

install				:	required	$(bindir)	$(includedir)	$(libdir)	$(mandir)/man1
					-$(CP) -f ../bin/[agm]*$(EXE) $(bindir)
					-$(CP) -f ../include/*scotch*.h $(includedir)
					-$(CP) -f ../lib/*scotch*$(LIB) $(libdir)
					-$(CP) -Rf ../man/* $(mandir)

clean				:	required
					(cd libscotch ;      $(MAKE) clean)
					(cd scotch ;         $(MAKE) clean)

realclean			:	required
					(cd libscotch ;      $(MAKE) realclean)
					(cd scotch ;         $(MAKE) realclean)
					-$(RM) ../bin/* ../include/* ../lib/*
