#!/usr/bin/perl

#
# read and clean all Fortran files in the current directory and subdirectories
#

## DK DK June 2018: visualize (dry run) all the lines to which the other script will add checking of all allocate() statements in the code, calling exit_MPI_without_rank() in case of an error

## DK DK June 2018: three known minor issues for this script:
## - it makes replacements in all strings, even in print statements; thus if you have a print statement that contains the word "allocate()" (with parentheses) it will be changed :-(
## - since it works line by line, it will add a check even if there is an existing one already in the lines below the allocate() statement; this does not hurt, but can lead to duplicated lines or similar line
## - since it works line by line, it also has issues with allocate() statements that may extend over several lines; if so, it will likely generate something slightly incorrect,
#"      which you will have to fix manually when compiling the code for the first time (that should be easy to do)

# DK DK only do this in the "src" directory, otherwise independent programs in other directories such as "utils" will not have access to the "exit_MPI_without_rank()" subroutine

# when using this "find" command from Perl we need to use \\ instead of \ below otherwise Perl tries to interpret it
      @objects = `find 'src' -name '.git' -prune -o -name 'm4' -prune -o -type f -regextype posix-extended -regex '.*\\.(fh|f90|F90|h\\.in|fh\\.in)' -print`;

      foreach $name (@objects) {
            chop $name;
# change tabs to white spaces
            system("expand -2 < $name > _____temp08_____");
            $f90name = $name;
#####################            print STDOUT "Cleaning $f90name ...\n";

            open(FILE_INPUT,"<_____temp08_____");
###################            open(FILEF90,">$f90name");

# open the input f90 file
      while($line = <FILE_INPUT>) {

# suppress trailing white spaces and carriage return
      $line =~ s/\s*$//;

# clear the flag that says if we need to add an allocate statement check or not
      $need_to_add_an_allocate_statement_check = 0;

      $linewithnospaceatall = $line;
      $linewithnospaceatall =~ s# ##ogi;
      $first_letter = substr(($linewithnospaceatall),0,1);
# do not make replacements in comments
      if($first_letter ne '!') {

# test if the line contains an "allocate()" statement
        if (index($line, "allocate(") != -1) {
# and test that it is not a "deallocate()" statement
        if (index($line, "deallocate(") == -1) {

# suppress trailing white spaces, just in case we have added any in the above processing
          $line =~ s/\s*$//;

      print "$line\n";

        }
        }
      }

###################      print FILEF90 "$line\n";

      }

      close(FILE_INPUT);
###################      close(FILEF90);

      }

      system("rm -f _____temp08_____");

