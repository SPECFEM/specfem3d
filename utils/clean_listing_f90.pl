#!/usr/bin/perl

#
#  Clean spaces, tabs and other non-standard or obsolete things in f90 files
#
#  Author : Dimitri Komatitsch, EPS - Harvard University, January 1998
#

#
# read and clean all f90 files in the current directory
#

      @objects = `ls *.f90 *.F90 *.h *.h.in */*.f90 */*.F90 */*.h */*.h.in */*/*.f90 */*/*.F90 */*/*.h */*/*.h.in`;

      foreach $name (@objects) {
            chop $name;
# change tabs to white spaces
            system("expand -2 < $name > _____tutu08_____");
            $f90name = $name;
            print STDOUT "Cleaning $f90name ...\n";

            open(FILE_INPUT,"<_____tutu08_____");
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

# change pause statements into stop statements
#     $line =~ s#pause#stop#ogi;

      $line =~ s#end do#enddo#ogi;
      $line =~ s#end if#endif#ogi;

      $line =~ s#elseif#else if#ogi;

      print FILEF90 "$line\n";

      }

            close(FILE_INPUT);
            close(FILEF90);

      }

            system("rm -f _____tutu08_____");

