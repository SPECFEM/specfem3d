#!/usr/bin/perl


################
################ run this script from the root directory of the package
################


# List of source code beautifiers found by Dimitri, June 2014:
# -----------------------------------------------------------

# This Perl script (below), which we developed ourselves and cannot hurt anyway, and which BuildBot runs on the three source codes.

# Basic, but replaces all tabs, removes useless white spaces, switches all comparison operators to their modern syntax
# (changing .GE. to >= and so on), and a few other useful things.
# It does not analyze indenting though.

# For C files we could use "A-Style" instead ( http://astyle.sourceforge.net/ ) as we mentioned (I am not sure if it
# handles CUDA in addition to C; it does NOT support Fortran for sure).

# --------------------

# For Fortran:

# https://www.freelancer.com/projects/Python-Fortran/Fortran-source-code-beautifier-Python.html
# Quote from that page: "NOTE: as example of beautiful Fortran code refer to http://www.cp2k.org/. By the way, as part of
# CP2K there is Python normalizer which performs good Fortran beautification and can be used as a baseline."
# http://www.cp2k.org
# http://www.cp2k.org/dev:codingconventions

# It contains a Python script called prettify.py, which calls another Python program called normalizeFortranFile.py.
# We could try it (I am not sure I like their idea of converting all Fortran keywords to uppercase, but I think there is
# a flag to turn that off).
# Not clear if it takes care of indenting, it probably does but we should check (if not then it is probably not very
# useful for SPECFEM because my simple Perl script below is then probably sufficient).

# --------------------

# http://dev.eclipse.org/mhonarc/lists/photran/msg01761.html
# http://dev.eclipse.org/mhonarc/lists/photran/msg01767.html
# http://www.ifremer.fr/ditigo/molagnon/fortran90/

# Was nice, but written in 2001 and unsupported since then (although the language has not changed that much since then,
# thus it may not be a problem; but Fortran2003 has a few extra keywords, and it seems this package only supports F90 but
# not F95, which we use a lot).
# I would tend to favor the Python program above because it is actively maintained and was written much more recently.
# And it is used routinely by CP2K developers, CP2K being a well-known open-source project for molecular dynamics.

# --------------------

# http://ray.met.fsu.edu/~bret/fpret.html

# Was nice, but written in 1995 and only supports Fortran77, thus cannot be used in our case.

# --------------------

# Conclusions:

# - BuildBot already automatically runs this Perl script on every pull request.

# - see if calling prettify.py from CP2K (with upper case conversion turned off?) for our Fortran files
# could be a good idea.


# --------------------
# --------------------
# --------------------


#
#  Clean spaces, tabs and other non-standard or obsolete things in Fortran files
#
#  Author : Dimitri Komatitsch, EPS - Harvard University, USA, January 1998
#

#
# read and clean all Fortran files in the current directory and subdirectories
#

#     @objects = `ls *.f90 *.F90 *.h *.h.in *.fh */*.f90 */*.F90 */*.h */*.h.in */*.fh */*/*.f90 */*/*.F90 */*/*.h */*/*.h.in */*/*.fh */*/*/*.f90 */*/*/*.F90 */*/*/*.h */*/*/*.h.in */*/*/*.fh`;
# when using this "find" command from Perl we need to use \\ instead of \ below otherwise Perl tries to interpret it
      @objects = `find . -name '.git' -prune -o -name 'm4' -prune -o -path './utils/ADJOINT_TOMOGRAPHY_TOOLS/flexwin' -prune -o -type f -regextype posix-extended -regex '.*\\.(fh|f90|F90|h|h\\.in)' -print`;

      foreach $name (@objects) {
            chop $name;
# change tabs to white spaces
            system("expand -2 < $name > _____temp08_____");
            $f90name = $name;
            print STDOUT "Cleaning $f90name ...\n";

            open(FILE_INPUT,"<_____temp08_____");
            open(FILEF90,">$f90name");

# open the input f90 file
      while($line = <FILE_INPUT>) {

# suppress trailing white spaces and carriage return
      $line =~ s/\s*$//;

# use new syntax of comparison operators, ignoring case in starting pattern (useful in case of mixed case)
      $line =~ s#\.le\.#<=#ogi;
      $line =~ s#\.ge\.#>=#ogi;
      $line =~ s#\.lt\.#<#ogi;
      $line =~ s#\.gt\.#>#ogi;
      $line =~ s#\.eq\.#==#ogi;
      $line =~ s#\.ne\.#/=#ogi;

# switch to lowercase for comparison operators
      $line =~ s#\.and\.#\.and\.#ogi;
      $line =~ s#\.or\.#\.or\.#ogi;
      $line =~ s#\.not\.#\.not\.#ogi;
      $line =~ s#\.true\.#\.true\.#ogi;
      $line =~ s#\.false\.#\.false\.#ogi;

#
# known problem: makes the changes also in constant strings, and not only
# in the code (which is dangerous, but really easier to program...)
#
# DK DK this could be dangerous if these words appear in strings or print statements
      $line =~ s#end\s*if#endif#ogi;
      $line =~ s#end\s*do#enddo#ogi;
      $line =~ s#elseif#else if#ogi;
      $line =~ s#use\s*::\s*mpi#use mpi#ogi;
      $line =~ s#enddo_LOOP_IJK#ENDDO_LOOP_IJK#ogi;
      $line =~ s#print\*#print \*#ogi;

      print FILEF90 "$line\n";

      }

            close(FILE_INPUT);
            close(FILEF90);

      }

            system("rm -f _____temp08_____");

################################################################################################
################################################################################################
################################################################################################

#
#  Clean spaces in different ASCII files
#
#  Author : Dimitri Komatitsch, EPS - Harvard University, USA, January 1998
#

#
# read and clean all these files in the current directory and subdirectories
#

#     @objects = `ls *.c *.cu *.h *.h.in *.fh */*.c */*.cu */*.h */*.h.in */*.fh */*/*.c */*/*.cu */*/*.h */*/*.h.in */*/*.fh */*/*/*.c */*/*/*.cu */*/*/*.h */*/*/*.h.in */*/*/*.fh`;
# when using this "find" command from Perl we need to use \\ instead of \ below otherwise Perl tries to interpret it
# purposely excluding Python files from this list, since (I think) Python can use tabs for indentation (?) and thus they should not be converted to two spaces (?)
      @objects = `find . -type f -name \\*Par_file\\* -print -o -name '.git' -prune -o -name 'm4' -prune -o -path './utils/ADJOINT_TOMOGRAPHY_TOOLS/flexwin' -prune -o -type f -regextype posix-extended -regex '.*\\.(bash|c|csh|cu|fh|f90|F90|h|h\\.in|pl|tex|txt|sh)' -print`;

      foreach $name (@objects) {
            chop $name;
# change tabs to white spaces
            system("expand -2 < $name > _____temp08_____");
            $cname = $name;
            print STDOUT "Cleaning $cname ...\n";

            open(FILEDEPART,"<_____temp08_____");
            open(FILEC,">$cname");

# open the input C file
      while($line = <FILEDEPART>) {

# suppress trailing white spaces and carriage return
      $line =~ s/\s*$//;

      print FILEC "$line\n";

      }

            close(FILEDEPART);
            close(FILEC);

      }

            system("rm -f _____temp08_____");

################################################################################################
################################################################################################
################################################################################################

#
#  Clean all accented letters and non-ASCII characters to remain portable
#
#  Authors : David Luet, Princeton University, USA and Dimitri Komatitsch, CNRS, France, February 2015; "find" command from Elliott Sales de Andrade
#

#     @objects = `ls *.txt *.c *.cu *.h *.h.in *.fh */*.c */*.cu */*.h */*.h.in */*.fh */*/*.c */*/*.cu */*/*.h */*/*.h.in */*/*.fh */*/*/*.c */*/*/*.cu */*/*/*.h */*/*/*.h.in */*/*/*.fh *.f90 *.F90 *.h *.h.in *.fh */*.f90 */*.F90 */*.h */*.h.in */*.fh */*/*.f90 */*/*.F90 */*/*.h */*/*.h.in */*/*.fh */*/*/*.f90 */*/*/*.F90 */*/*/*.h */*/*/*.h.in */*/*/*.fh */*.txt */*/*.txt */*/*/*.txt */*.tex */*/*.tex */*/*/*.tex *.sh */*.sh */*/*.sh */*/*/*.sh *.csh */*.csh */*/*.csh */*/*/*.csh *.bash */*.bash */*/*.bash */*/*/*.bash *.pl */*.pl */*/*.pl */*/*/*.pl *.py */*.py */*/*.py */*/*/*.py`;
# when using this "find" command from Perl we need to use \\ instead of \ below otherwise Perl tries to interpret it
      @objects = `find . -type f -name \\*Par_file\\* -print -o -name '.git' -prune -o -name 'm4' -prune -o -path './utils/ADJOINT_TOMOGRAPHY_TOOLS/flexwin' -prune -o -type f -regextype posix-extended -regex '.*\\.(bash|c|csh|cu|fh|f90|F90|h|h\\.in|pl|py|tex|txt|sh)' -print`;

      system("rm -f _____temp08_____ _____temp09_____");

      foreach $name (@objects) {
            chop $name;
            print STDOUT "Cleaning $name ...\n";

# this preserves file permissions of executable files
# the second "cp" is there just to be safe: if "iconv" is not installed on a system, then the last "cp", if executed, will copy the file identical to itself
            system("cp -p $name _____temp08_____ ; cp -p $name _____temp09_____ ; iconv -f utf-8 -t ASCII//TRANSLIT _____temp08_____ -o _____temp09_____ ; cp _____temp09_____ $name");

      }

      system("rm -f _____temp08_____ _____temp09_____");

