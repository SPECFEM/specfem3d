## * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
##* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
##=======================================================================
##Copyright (C) 2010-2011 Leonardo BAUTISTA GOMEZ
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
##
## File          : Makefile
## Created on    : 27 Sep 2010
## Author        : Leonardo BAUTISTA G <leobago@matsulab.is.titech.ac.jp>
##
## Last modified : 27 Sep 2010 (12:20:49 PM)
## Author        : Leonardo BAUTISTA G <leobago@matsulab.is.titech.ac.jp>
## Description   : .
##

NBPR		= 24
BIN 		= bin
LIB 		= lib
OBJ		= obj
SRC		= src
SPT		= scripts
INC		= include
CC 		= gcc
MPICC 		= mpicc
MPIRUN 		= mpirun
INCS		= -I/usr/lib/openmpi/include/ -Iinclude/
MPIFLAGS	= -g -fPIC $(INCS)
PYFLAGS 	= -lpython2.4 -lm -lutil -ldl -lpthread
OBJS		= $(OBJ)/galois.o $(OBJ)/enc.o $(OBJ)/dec.o $(OBJ)/tools.o $(OBJ)/topo.o

SHARED		= libfti.so.1

all: 		$(SHARED) prog

$(OBJ)/%.o:	$(SRC)/%.c
	$(MPICC) $(MPIFLAGS) -c $< -o $@

$(SHARED):	$(OBJS)
	rm -f libfti.so $(SHARED)
	$(CC) $(INCS) $(PYFLAGS) -shared -o $@ $(OBJS)
	ln -s $(SHARED) libfti.so

prog:	$(SHARED) $(SRC)/program.c
	$(MPICC) $(INCS) -L$(PWD) -lfti $(PYFLAGS) -o $(BIN)/prog $(SRC)/program.c

# Test creating 100MB local checkpoint files
test1:
	$(MPIRUN) --hostfile Hostfile -np 24 $(BIN)/prog 400000000 50;

clean:
	rm -f  $(OBJ)/* $(BIN)/* $(LIB)/*

.PHONY:		$(SHARED) prog test1 clean
