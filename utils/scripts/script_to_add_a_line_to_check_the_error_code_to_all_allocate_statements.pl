#!/usr/bin/perl

#
# read and clean all Fortran files in the current directory and subdirectories
#

## DK DK June 2018: add checking of all allocate() statements in the code, calling exit_MPI_without_rank() in case of an error

## DK DK June 2018: three known minor issues for this script:
## - it makes replacements in all strings, even in print statements; thus if you have a print statement that contains the word "allocate()" (with parentheses) it will be changed :-(
## - since it works line by line, it will add a check even if there is an existing one already in the lines below the allocate() statement; this does not hurt, but can lead to duplicated lines or similar line
## - since it works line by line, it also has issues with allocate() statements that may extend over several lines; if so, it will likely generate something slightly incorrect,
#"      which you will have to fix manually when compiling the code for the first time (that should be easy to do)
## Because of these three known issues, I do not put this script in the automatic cleaning script of Buildbot (specfem3d/utils/clean_listings_specfem.pl)

# DK DK only do this in the "src" directory, otherwise independent programs in other directories such as "utils" will not have access to the "exit_MPI_without_rank()" subroutine

# when using this "find" command from Perl we need to use \\ instead of \ below otherwise Perl tries to interpret it
      @objects = `find 'src' -name '.git' -prune -o -name 'm4' -prune -o -type f -regextype posix-extended -regex '.*\\.(fh|f90|F90|h\\.in|fh\\.in)' -print`;

# create a counter to put in the error message printed if the allocate() fails, so that users can see which allocate() statement had a problem
      $counter_for_error_message = 0;

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

      $need_to_add_an_allocate_statement_check = 1;
      $counter_for_error_message = $counter_for_error_message + 1;

# count the number of white spaces until the beginning of the allocate() statement
# we will add the same number of white spaces in front of the if() test, to align it with it
      $number_of_white_spaces_before_allocate = index($line, "allocate(");
      if ($number_of_white_spaces_before_allocate < 0) {die "error when counting number of white spaces before allocate() statement";}

# if the stat=ier) exit code is not present at the end of the allocate() statement, add it by replacing the final parenthesis with ,stat=ier)
      if (index($line, "stat=ier") == -1) {
          $line =~ s/\)\s*$/\,stat=ier\)/;
         }

# suppress trailing white spaces, just in case we have added any in the above processing
          $line =~ s/\s*$//;

        }
        }
      }

      print FILEF90 "$line\n";
      if ($need_to_add_an_allocate_statement_check != 0) {
# put the right number of white spaces before the if() test in order to align it with the allocate() statement
        print FILEF90 ' ' x $number_of_white_spaces_before_allocate . "if (ier /= 0) call exit_MPI_without_rank('error allocating array $counter_for_error_message')\n";
        }

      }

      close(FILE_INPUT);
      close(FILEF90);

      }

      print "\n";
      print "A total of $counter_for_error_message error messages have been added to check all the allocate() statements of the code\n";
      print "\n";

      system("rm -f _____temp08_____");

