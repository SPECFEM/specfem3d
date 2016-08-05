## * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
##* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
##=======================================================================
##Copyright (C) 2010-2013 Leonardo A. BAUTISTA GOMEZ
##This program is free software; you can redistribute it and/or modify
##it under the terms of the GNU General Public License (GPL) as published
##of the License, or (at your option) any later version.
##
##This program is distributed in the hope that it will be useful,
##but WITHOUT ANY WARRANTY; without even the implied warranty of
##MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##GNU General Public License for more details.
##
##To read the license please visit http://www.gnu.org/copyleft/gpl.html
##=======================================================================

##=======================================================================
##   PLEASE SET THESE VARIABLES BEFORE COMPILING
##=======================================================================

FTIPATH		= /path/to/fti/install/directory
MPIPATH		= /path/to/mpi/directory

##=======================================================================
##   DIRECTORY TREE
##=======================================================================

LIB 		= lib
OBJ		= obj
SRC		= src
DOC		= doc
INC		= include

##=======================================================================
##   COMPILERS
##=======================================================================

CC 		= gcc
MPICC 		= mpicc

##=======================================================================
##   FLAGS
##=======================================================================

FTIFLAGS	= -fPIC -g -Iinclude/ -c
MPIFLAGS	= -I$(MPIPATH)

##=======================================================================
##   TARGETS
##=======================================================================


OBJS		= $(OBJ)/galois.o $(OBJ)/jerasure.o $(OBJ)/dictionary.o $(OBJ)/iniparser.o \
		$(OBJ)/recover.o $(OBJ)/tools.o $(OBJ)/checkpoint.o $(OBJ)/topo.o \
		$(OBJ)/api.o

SHARED		= libfti.so
STATIC		= libfti.a

all: 		$(SHARED) $(STATIC)

doc:
		doxygen $(DOC)/Doxyfile

$(OBJ)/%.o:	$(SRC)/%.c
		$(MPICC) $(FTIFLAGS) $< -o $@

$(SHARED):	$(OBJS)
		$(CC) $(MPIFLAGS) -shared -o $(LIB)/$(SHARED) $(OBJS) -lc

$(STATIC):	$(OBJS)
		$(AR) -cvq $(LIB)/$(STATIC) $(OBJS)
install:
		if [ ! -d "$(FTIPATH)/FTI" ]; then mkdir $(FTIPATH)/FTI; fi
		if [ ! -d "$(FTIPATH)/FTI/$(LIB)" ]; then mkdir $(FTIPATH)/FTI/$(LIB); fi
		if [ ! -d "$(FTIPATH)/FTI/$(INC)" ]; then mkdir $(FTIPATH)/FTI/$(INC); fi
		rm -f $(FTIPATH)/FTI/$(LIB)/* $(FTIPATH)/FTI/$(INC)/*
		cp $(INC)/* $(FTIPATH)/FTI/$(INC)/
		cp $(LIB)/* $(FTIPATH)/FTI/$(LIB)/

uninstall:
		rm -f $(FTIPATH)/FTI/$(LIB)/* $(FTIPATH)/FTI/$(INC)/*
		if [ -d "$(FTIPATH)/FTI/$(LIB)" ]; then rmdir $(FTIPATH)/FTI/$(LIB); fi
		if [ -d "$(FTIPATH)/FTI/$(INC)" ]; then rmdir $(FTIPATH)/FTI/$(INC); fi
		if [ -d "$(FTIPATH)/FTI" ]; then rmdir $(FTIPATH)/FTI; fi

clean:
		rm -f $(OBJ)/* $(LIB)/*

.PHONY:		$(SHARED) $(STATIC) doc install uninstall clean


