#!/usr/bin/perl

#
#  replace "use mpi" with include "mpif.h" in all source files;
#  also swap it with the "implicit none" line, which must be located after
#  in the case of "use mpi" but before in the case of include "mpif.h"
#
#  Author : Dimitri Komatitsch, EPS - Harvard University, January 1998 and CNRS Marseille, France, July 2013
#

#
# first clean trailing white spaces in all f90 files in the src/ sub-directories
#

      @objects = `ls src/*/*.f90 src/*/*.F90`;

      foreach $name (@objects) {
            chop $name;
# change tabs to white spaces
            system("expand -2 < $name > _____dummy08_____");
            $f90name = $name;
            print STDOUT "Cleaning trailing white spaces in $f90name ...\n";

            open(FILE_INPUT,"<_____dummy08_____");
            open(FILEF90,">$f90name");

# open the input f90 file
      while($line = <FILE_INPUT>) {

# suppress trailing white spaces and carriage return
      $line =~ s/\s*$//;

      print FILEF90 "$line\n";

      }

            close(FILE_INPUT);
            close(FILEF90);

      }

#
# then perform the replacement in all f90 and F90 files in the src/ sub-directories
#

      @objects = `ls src/*/*.f90 src/*/*.F90`;

      foreach $name (@objects) {
            chop $name;
# change tabs to white spaces
            system("expand -2 < $name > _____dummy08_____");
            $f90name = $name;
            print STDOUT "Replacing 'use mpi' (if any) in $f90name ...\n";

            open(FILE_INPUT,"<_____dummy08_____");
            open(FILEF90,">$f90name");

# to read the whole file in the variable instead of a single line
      undef $/;

# read the whole input file
      $whole_file = <FILE_INPUT>;

# make the replacement
      $whole_file =~ s/\n\s*use mpi\s*\n*\s*implicit none\s*\n/\n\n  implicit none\n\n  include 'mpif.h'\n\n/og;

      print FILEF90 "$whole_file";

            close(FILE_INPUT);
            close(FILEF90);

      }

            system("rm -f _____dummy08_____");

