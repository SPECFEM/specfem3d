#!/usr/bin/perl

#
# script to beautify a source code (in Fortran + C + CUDA), remove extra white spaces,
# and convert accented letters to non-accented to maximize portability.
#
# Dimitri Komatitsch, CNRS Marseille, France, 2015.
#

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

use File::Basename;

# when using this "find" command from Perl we need to use \\ instead of \ below otherwise Perl tries to interpret it
@objects = `find . -name '.git' -prune -o -name 'm4' -prune -o -path './utils/ADJOINT_TOMOGRAPHY_TOOLS/flexwin' -prune -o -type f -regextype posix-extended -regex '.*\\.(fh|f90|F90|h\\.in|fh\\.in)' -print`;

foreach $name (@objects) {
  chop $name;
# change tabs to white spaces
  system("expand -2 < $name > _____temp08_____");
  $f90name = $name;
  print STDOUT "Cleaning $f90name ...\n";

  open(FILE_INPUT,"<_____temp08_____");
  open(FILEF90,">$f90name");

  $basename_obtained = basename($f90name,@suffixlist);

# open the input f90 file
  while($line = <FILE_INPUT>) {

# suppress trailing white spaces and carriage return
    $line =~ s/\s*$//;

# use new syntax of comparison operators, ignoring case in starting pattern (useful in case of mixed case)
    $line =~ s#\.le\.#<=#ogi;
    $line =~ s#\.ge\.#>=#ogi;
    $line =~ s#\.lt\.#<#ogi;
    $line =~ s#\.gt\.#>#ogi;
    $line =~ s#\.ne\.#/=#ogi;
    $line =~ s#\.eq\.#==#ogi;

# switch to lowercase for comparison operators
    $line =~ s#\.and\.#\.and\.#ogi;
    $line =~ s#\.or\.#\.or\.#ogi;
    $line =~ s#\.not\.#\.not\.#ogi;
    $line =~ s#\.eqv\.#\.eqv\.#ogi;
    $line =~ s#\.neqv\.#\.neqv\.#ogi;

    $line =~ s#\.true\.#\.true\.#ogi;
    $line =~ s#\.false\.#\.false\.#ogi;

# switch to Fortran2008 standard
    $line =~ s#call\s*getarg\(#call get_command_argument\(#ogi;

## DK DK reverting this because it breaks other compilers, since "\0" is not a true empty strings for other compilers.
## DK DK see https://github.com/SPECFEM/specfem3d/issues/818
# for null strings, make sure we fully conform to the standard otherwise the IBM xlf compiler refuses to compile
# do not do it for "Poisson''s ratio" in a Fortran print statement
#     if(index($line, "Poisson''s") == -1) {
#     if(index($line, "Young''s") == -1) {
#     if(index($line, "''''") == -1) {
#       $line =~ s#\"\"#\'\\0\'#ogi;
#       $line =~ s#\'\'#\'\\0\'#ogi;
#     }
#     }
#     }

# this is to revert the above fix for IBM xlf, which broke other compilers
    $line =~ s#\'\\0\'#\'\'#ogi;
    $line =~ s#\"\\0\"#\'\'#ogi;

#
# known problem: makes the changes also in constant strings, and not only
# in the code (which is dangerous, but really easier to program...)
#
# DK DK this could be dangerous if these words appear in strings or print statements
    $line =~ s#if\s*\(#if \(#ogi;
    $line =~ s#\)\s*then#\) then#ogi;
    $line =~ s#end\s*if#endif#ogi;
    $line =~ s#end\s*do#enddo#ogi;
    $line =~ s#else\s*if#else if#ogi;

# force lowercase keywords
    $line =~ s#subroutine#subroutine#ogi;
    $line =~ s#end\s*subroutine#end subroutine#ogi;
    $line =~ s#function#function#ogi;
    $line =~ s#end\s*function#end function#ogi;
    $line =~ s#continue#continue#ogi;
    $line =~ s#implicit none#implicit none#ogi;
    $line =~ s#implicit#implicit#ogi;
    $line =~ s#return#return#ogi;

    $line =~ s# go\s*to # goto #ogi;

    $line =~ s#use\s*::\s*mpi#use mpi#ogi;

    $line =~ s#,\s*only\s*:\s*#, only: #ogi;

    $line =~ s#print\s*\*#print \*#ogi;

    $line =~ s#NOISE_SOURCE_TIME_FUNCTION_TYPE#noise_source_time_function_type#ogi;

# do not move this before the above line in which we change the keyword "function"
    $line =~ s#use_ricker_time_function#USE_RICKER_TIME_FUNCTION#ogi;
    $line =~ s#print_source_time_function#PRINT_SOURCE_TIME_FUNCTION#ogi;
    $line =~ s#external_source_time_function#EXTERNAL_SOURCE_TIME_FUNCTION#ogi;
    $line =~ s#sourceTimeFunction#sourceTimeFunction#ogi;

    $line =~ s#external_stf#EXTERNAL_SOURCE_TIME_FUNCTION#ogi;
    $line =~ s#EXTERNAL_SOURCE_TIME_FUNCTION_filename#external_source_time_function_filename#ogi;
    $line =~ s#read_EXTERNAL_SOURCE_TIME_FUNCTION#read_external_source_time_function#ogi;
    $line =~ s#USE_MAP_function#USE_MAP_FUNCTION#ogi;

    $line =~ s#enddo_LOOP_IJK#ENDDO_LOOP_IJK#ogi;
    $line =~ s#OMP do#OMP DO#ogi;
    $line =~ s#OMP enddo#OMP ENDDO#ogi;

    $line =~ s#print\*#print \*#ogi;

    $line =~ s#spectral-elements#spectral elements#ogi;

    $line =~ s#gaussian#Gaussian#ogi;
    $line =~ s#hessian#Hessian#ogi;
    $line =~ s#cartesian#Cartesian#ogi;

# convert the GNU license from version 2+ to version 3 (we switched in June 2017)
    $line =~ s#either version 2 of the License#either version 3 of the License#ogi;

# for spaces around comparisons we exclude that file because it contains print statements that print
# XML file lines and thus it contains many < and > characters that must not be changed
    if($basename_obtained ne 'write_output_ASDF.f90' && $basename_obtained ne 'config.h.in') {

# make sure there is one white space on each side of comparison operators
      $line =~ s#\s*<\s*=\s*# <= #ogi;
      $line =~ s#\s*>\s*=\s*# >= #ogi;
      $line =~ s#\s*<\s*# < #ogi;
      $line =~ s#\s*/=\s*# /= #ogi;

# for these two we make sure the line is not a comment, because we have some comments that contain ===================== as paragraph delimiters
# and we also have Doxygen markers in comments that start with !>
      $linewithnospaceatall = $line;
      $linewithnospaceatall =~ s# ##ogi;
      my $first_letter = substr(($linewithnospaceatall),0,1);
      if($first_letter ne '!') {
        $line =~ s#\s*>\s*# > #ogi;
        $line =~ s#\s*==\s*# == #ogi;
      }

# restore operators that may have been split by the above introduction of white spaces
      $line =~ s#<\s*=#<=#ogi;
      $line =~ s#>\s*=#>=#ogi;
      $line =~ s#=\s*=#==#ogi;
      $line =~ s#/\s*=#/=#ogi;

# also restore bash file pipes that may appear in some print statements that save bash scripts to disk for future processing
      $line =~ s#>\s*&#>&#ogi;
      $line =~ s#<\s*&#<&#ogi;

# also restore xml-formatting strings '< and >'
      $line =~ s#\'\s*<\s*#\'<#ogi;
      $line =~ s#\s*>\s*\'#>\'#ogi;

# also restore -> and <- that may be used in comments; do this in comments only
# otherwise comparisons with negative numbers in source code may be affected (e.g. if (a < -100.) then...)
      if($first_letter eq '!') {
        $line =~ s#-\s*>#->#ogi;
        $line =~ s#<\s*-#<-#ogi;
      }

# for pointers
      $line =~ s#\s*=\s*>\s*# => #ogi;

    }

# for these ones we make sure that keyword is not the first of the sentence, otherwise we may break existing alignment
    $linewithnospaceatall = $line;
    $linewithnospaceatall =~ s# ##ogi;
    my $first_letter = substr(($linewithnospaceatall),0,1);
    if($first_letter ne '.') {
      $line =~ s#\s*\.and\.\s*# \.and\. #ogi;
      $line =~ s#\s*\.or\.\s*# \.or\. #ogi;
      $line =~ s#\s*\.not\.\s*# \.not\. #ogi;
      $line =~ s#\s*\.eqv\.\s*# \.eqv\. #ogi;
      $line =~ s#\s*\.neqv\.\s*# \.neqv\. #ogi;
    }

# suppress space between parenthesis and .not. (this can happen when testing logical operators)
    $line =~ s#\( \.not\. #\(\.not\. #ogi;

    $line =~ s#\)call#\) call#ogi;

# enforce upper case
    $line =~ s#CUSTOM_REAL#CUSTOM_REAL#ogi;

# do not use null strings, which are not part of the Fortran standard (and the IBM xlf compiler rejects them for instance)
    $line =~ s#print\s*\*\s*,\s*''#print \*#ogi;
    $line =~ s#write\s*\(\s*\*\s*,\s*\*\s*\)\s*''#print \*#ogi;
    $line =~ s#write\s*\(\s*IMAIN\s*,\s*\*\s*\)\s*''#write\(IMAIN,\*\)#ogi;
    $line =~ s#write\s*\(\s*IOUT\s*,\s*\*\s*\)\s*''#write\(IOUT,\*\)#ogi;

    $line =~ s#print\s*\*\s*,\s*""#print \*#ogi;
    $line =~ s#write\s*\(\s*\*\s*,\s*\*\s*\)\s*""#print \*#ogi;
    $line =~ s#write\s*\(\s*IMAIN\s*,\s*\*\s*\)\s*""#write\(IMAIN,\*\)#ogi;
    $line =~ s#write\s*\(\s*IOUT\s*,\s*\*\s*\)\s*""#write\(IOUT,\*\)#ogi;

# unit 6 means standard output, replace it with standard output symbol
    $line =~ s#write\s*\(\s*6\s*,\s*\*\s*\)#write\(\*,\*\)#ogi;
    $line =~ s#write\s*\(\s*6\s*,#write\(\*,#ogi;

# force space in , & at end of line
    $line =~ s#\s*,\s*&\s*$#, &#ogi;

# always use upper case for GLL when used as a word
    $line =~ s# gll # GLL #ogi;

    $line =~ s# mpi # MPI #ogi;

    $line =~ s# pml # PML #ogi;

# fix some typos I have found in the different codes, or non-US spelling.
# also switch to US spelling in order to have the same standard in all files.
    $line =~ s#regularisation#regularization#ogi;
    $line =~ s#optimisation#optimization#ogi;
    $line =~ s#analitical#analytical#ogi;
    $line =~ s#communIcation#communication#ogi;
    $line =~ s# in orfer # in order #ogi;
    $line =~ s# stepest # steepest #ogi;
    $line =~ s# stepest$# steepest#ogi;
    $line =~ s#aloow#allow#ogi;
    $line =~ s#neighbour#neighbor#ogi;
    $line =~ s#vecotr#vector#ogi;
    $line =~ s#computse#compute#ogi;
    $line =~ s#indicies#indices#ogi;
    $line =~ s#accordig#according#ogi;
    $line =~ s#paralell#parallel#ogi;
    $line =~ s#debbug#debug#ogi;
# do not suppress the white space here because it would then change "debugging" for instance
    $line =~ s#debugg #debug #ogi;
    $line =~ s#debugg$#debug#ogi;
    $line =~ s#familly#family#ogi;
    $line =~ s#warnning#warning#ogi;
    $line =~ s#elemement#element#ogi;
    $line =~ s#cartesion#Cartesian#ogi;
    $line =~ s#partiton#partition#ogi;
    $line =~ s#drection#direction#ogi;
    $line =~ s#seperation#separation#ogi;
    $line =~ s#inverision#inversion#ogi;
    $line =~ s#restauration#restoration#ogi;
    $line =~ s#restaure#restore#ogi;
    $line =~ s#memmory#memory#ogi;
    $line =~ s#convolution formation#convolution formulation#ogi;
    $line =~ s#fortran#Fortran#ogi;
    $line =~ s#adress#address#ogi;
    $line =~ s#gFortran#gfortran#ogi;

    $line =~ s# usefull # useful #ogi;
    $line =~ s# usefull$# useful#ogi;

# enforce upper case
    $line =~ s#MAX_neighborS#MAX_NEIGHBORS#ogi;

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

# When using this "find" command from Perl we need to use \\ instead of \ below otherwise Perl tries to interpret it.
# Purposely excluding Python files from this list, since Python can use tabs for indentation.
# Also make sure we exclude clean_listings_specfem.pl, i.e. this script, from the list, otherwise the script will clean itself and change the list of typos etc to the new one!
@objects = `find . -type f -name \\*Par_file\\* -print -o -name '.git' -prune -o -name 'm4' -prune -o -wholename './utils/clean_listings_specfem.pl' -prune -o -path './utils/ADJOINT_TOMOGRAPHY_TOOLS/flexwin' -prune -o -type f -regextype posix-extended -regex '.*\\.(bash|c|cpp|zsh|csh|cu|fh|f90|F90|h|h\\.in|fh\\.in|pl|tex|txt|sh)' -print`;

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

# convert the GNU license from version 2+ to version 3 (we switched in June 2017)
    $line =~ s#either version 2 of the License#either version 3 of the License#ogi;

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
# in the case of Python files, Makefiles and rules.mk files, do NOT expand the tabs into white spaces
#

# when using this "find" command from Perl we need to use \\ instead of \ below otherwise Perl tries to interpret it
@objects = `find . -type f -name \\*Par_file\\* -print -o -name '.git' -prune -o -name 'm4' -prune -o -path './utils/ADJOINT_TOMOGRAPHY_TOOLS/flexwin' -prune -o -type f -regextype posix-extended -regex '.*\\.(py)' -print -o -type f -regextype posix-extended -regex '.*Makefile*' -print -o -type f -regextype posix-extended -regex '.*Makefile\\.in' -print -o -type f -regextype posix-extended -regex '.*rules\\.mk' -print`;

foreach $name (@objects) {
  chop $name;
# create a copy of the input file
  system("cp $name _____temp08_____");
  $cname = $name;
  print STDOUT "Cleaning $cname ...\n";

  open(FILEDEPART,"<_____temp08_____");
  open(FILEC,">$cname");

# open the input C file
  while($line = <FILEDEPART>) {

# suppress trailing white spaces and carriage return
    $line =~ s/\s*$//;

# convert the GNU license from version 2+ to version 3 (we switched in June 2017)
    $line =~ s#either version 2 of the License#either version 3 of the License#ogi;

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

# when using this "find" command from Perl we need to use \\ instead of \ below otherwise Perl tries to interpret it
@objects = `find . -type f -name \\*Par_file\\* -print -o -name '.git' -prune -o -name 'm4' -prune -o -path './utils/ADJOINT_TOMOGRAPHY_TOOLS/flexwin' -prune -o -path './HDR_Vadim' -prune -o -type f -regextype posix-extended -regex '.*\\.(bash|c|cpp|zsh|csh|cu|fh|f90|F90|h|h\\.in|fh\\.in|pl|py|tex|txt|sh)' -print`;

system("rm -f _____temp08_____ _____temp09_____");

foreach $name (@objects) {
  chop $name;
  print STDOUT "Cleaning $name ...\n";

# this preserves file permissions of executable files
# the second "cp" is there just to be safe: if "iconv" is not installed on a system, then the last "cp", if executed, will copy the file identical to itself
  system("cp -p $name _____temp08_____ ; cp -p $name _____temp09_____ ; iconv -f utf-8 -t ASCII//TRANSLIT _____temp08_____ -o _____temp09_____ ; cp _____temp09_____ $name");

}

system("rm -f _____temp08_____ _____temp09_____");

