#! /usr/bin/perl
#
#    Copyright 2013, Tarje Nissen-Meyer, Alexandre Fournier, Martin van Driel
#                    Simon Stähler, Kasra Hosseini, Stefanie Hempel
#
#    This file is part of AxiSEM.
#    It is distributed from the webpage <http://www.axisem.info>
#
#    AxiSEM is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    AxiSEM is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with AxiSEM.  If not, see <http://www.gnu.org/licenses/>.
#
#    Generate a Makefile from the sources in the current directory.  The source
#    files may be in either C, FORTRAN 77, Fortran 90 or some combination of
#    these languages.  
#
#    Original version written by Michael Wester <wester@math.unm.edu> February 16, 1995
#    Cotopaxi (Consulting), Albuquerque, New Mexico
#
#    Modified by Martin van Driel, ETH Zürich and Simon Stähler, 
#    LMU München to fit the needs of Axisem. The compiler version is
#    now set in the file ../make_axisem.macros

open(MAKEFILE, "> Makefile");

print MAKEFILE "PROG = xmesh\n\n";

#
# Read header file with compiler names etc.
#
print MAKEFILE "include ../make_axisem.macros\n\n";

#
# Source listing
#
print MAKEFILE "SRCS =\t";
@srcs = <*.F90 *.f90 *.f *.F *.c>;
&PrintWords(8, 0, @srcs);
print MAKEFILE "\n\n";
#
# Object listing
#
print MAKEFILE "OBJS =\t";
@objs = @srcs;
foreach (@objs) { s/\.[^.]+$/.o/ };
&PrintWords(8, 0, @objs);
print MAKEFILE "\n\n";
#
# Define common macros
#
print MAKEFILE "ifeq (\$(strip \$(USE_NETCDF)),true)\n";
print MAKEFILE "   FFLAGS += -Dunc\n";
print MAKEFILE "   ifdef NETCDF_PATH\n";
print MAKEFILE "       LIBS = -L \$(strip \$(NETCDF_PATH))/lib -lnetcdff -Wl,-rpath,\$(strip \$(NETCDF_PATH))/lib\n";
print MAKEFILE "       INCLUDE = -I \$(strip \$(NETCDF_PATH))/include\n";
print MAKEFILE "   else\n";
print MAKEFILE "       LIBS = -lnetcdff\n";
print MAKEFILE "       INCLUDE = -I /usr/include\n";
print MAKEFILE "   endif\n";
print MAKEFILE "else\n";
print MAKEFILE "   LIBS = \n";
print MAKEFILE "   INCLUDE = \n";
print MAKEFILE "endif\n\n";


print MAKEFILE "\n\n";
print MAKEFILE "# cancel m2c implicit rule \n";
print MAKEFILE "%.o : %.mod \n ";
print MAKEFILE "\n\n";

#
# make
#
print MAKEFILE "all: \$(PROG)\n\n";
print MAKEFILE "\$(PROG): \$(OBJS)\n";
print MAKEFILE "\t\$(", &LanguageCompiler($ARGV[1], @srcs);
print MAKEFILE ") \$(LDFLAGS) -o \$@ \$(OBJS) \$(LIBS)\n\n";
#
# make clean
#
print MAKEFILE "clean:\n";
print MAKEFILE "\trm -f \$(PROG) \$(OBJS) *.M *.mod *.d *.il core \n\n";
#
# Make .f90 a valid suffix
#
print MAKEFILE ".SUFFIXES: \$(SUFFIXES) .f90 .F90\n\n";
#
# .f90 -> .o
#
print MAKEFILE ".f90.o:\n";
print MAKEFILE "\t\$(FC) \$(FFLAGS) -c \$(INCLUDE) \$<\n\n";
print MAKEFILE ".F90.o:\n";
print MAKEFILE "\t\$(FC) \$(FFLAGS) -c \$(INCLUDE) \$<\n\n";
#
# Dependency listings
#
&MakeDependsf90($ARGV[1]);
&MakeDepends("*.f *.F", '^\s*include\s+["\']([^"\']+)["\']');
&MakeDepends("*.c",     '^\s*#\s*include\s+["\']([^"\']+)["\']');

#
# &PrintWords(current output column, extra tab?, word list); --- print words
#    nicely
#
sub PrintWords {
   local($columns) = 78 - shift(@_);
   local($extratab) = shift(@_);
   local($wordlength);
   #
   print MAKEFILE @_[0];
   $columns -= length(shift(@_));
   foreach $word (@_) {
      $wordlength = length($word);
      if ($wordlength + 1 < $columns) {
         print MAKEFILE " $word";
         $columns -= $wordlength + 1;
         }
      else {
         #
         # Continue onto a new line
         #
         if ($extratab) {
            print MAKEFILE " \\\n\t\t$word";
            $columns = 62 - $wordlength;
            }
         else {
            print MAKEFILE " \\\n\t$word";
            $columns = 70 - $wordlength;
            }
         }
      }
   }

#
# &LanguageCompiler(compiler, sources); --- determine the correct language
#    compiler
#
sub LanguageCompiler {
   local($compiler) = &toLower(shift(@_));
   local(@srcs) = @_;
   #
   if (length($compiler) > 0) {
      CASE: {
         grep(/^$compiler$/, ("fc", "f77")) &&
            do { $compiler = "FC"; last CASE; };
         grep(/^$compiler$/, ("cc", "c"))   &&
            do { $compiler = "CC"; last CASE; };
         $compiler = "FC";
         }
      }
   else {
      CASE: {
         grep(/\.(f|F)90$/, @srcs)   && do { $compiler = "FC"; last CASE; };
         grep(/\.(f|F)$/, @srcs) && do { $compiler = "FC";  last CASE; };
         grep(/\.c$/, @srcs)     && do { $compiler = "CC";  last CASE; };
         $compiler = "???";
         }
      }
   $compiler;
   }

#
# &toLower(string); --- convert string into lower case
#
sub toLower {
   local($string) = @_[0];
   $string =~ tr/A-Z/a-z/;
   $string;
   }

#
# &uniq(sorted word list); --- remove adjacent duplicate words
#
sub uniq {
   local(@words);
   foreach $word (@_) {
      if ($word ne $words[$#words]) {
         push(@words, $word);
         }
      }
   @words;
   }

#
# &MakeDepends(language pattern, include file sed pattern); --- dependency
#    maker
#
sub MakeDepends {
   local(@incs);
   local($lang) = @_[0];
   local($pattern) = @_[1];
   #
   foreach $file (<${lang}>) {
      open(FILE, $file) || warn "Cannot open $file: $!\n";
      while (<FILE>) {
         /$pattern/i && push(@incs, $1);
         }
         $file =~ s/\.[^.]+$/.o/;
         print MAKEFILE "$file: ";
         &PrintWords(length($file) + 2, 0, @incs);
         print MAKEFILE " Makefile ../make_axisem.macros \n";
         undef @incs;
      }
   }

#
# &MakeDependsf90(f90 compiler); --- FORTRAN 90 dependency maker
#
sub MakeDependsf90 {
   local($compiler) = &toLower(@_[0]);
   local(@dependencies);
   local(%filename);
   local(@incs);
   local(@modules);
   local($objfile);
   #
   # Associate each module with the name of the file that contains it
   #
   foreach $file (<*.f90 *.F90>) {
      open(FILE, $file) || warn "Cannot open $file: $!\n";
      while (<FILE>) {
         /^\s*module\s+([^\s!]+)/i &&
            ($filename{&toLower($1)} = $file) =~ s/\.(f|F)90$/.o/;
         }
      }
   #
   # Print the dependencies of each file that has one or more include's or
   # references one or more modules
   #
   foreach $file (<*.f90 *.F90>) {
      open(FILE, $file);
      while (<FILE>) {
         /^\s*include\s+["\']([^"\']+)["\']/i && push(@incs,$1);
         /^\s*use\s+([^\s,!]+)/i && push(@modules, &toLower($1));
         }
      ($objfile = $file) =~ s/\.(f|F)90$/.o/;
         print MAKEFILE "$objfile: ";
         undef @dependencies;
         foreach $module (@modules) {
            push(@dependencies, $filename{$module});
            }
         @dependencies = &uniq(sort(@dependencies));
         &PrintWords(length($objfile) + 2, 0,
                     @dependencies, &uniq(sort(@incs)));
      print MAKEFILE " Makefile ../make_axisem.macros \n";
         undef @incs;
         undef @modules;
         #
         }
      }
   

print "\nCheck Makefile to make sure you're happy with it.\n\n";
system("perl -ni -e 'print unless /^kdtree2.o/' Makefile ");
