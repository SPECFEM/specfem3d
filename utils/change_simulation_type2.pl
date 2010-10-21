#!/usr/bin/perl

use Time::Local;
use Getopt::Std;
use POSIX;

sub Usage{
print STDERR <<END;

Usage:   change_simulation_type.pl  [-a|-f|-b|-F]
         Changes SIMULATION_TYPE in DATA/Par_file
         -a -- change type to run adjoint calculation(2)
         -f -- change type to run forward calculation(1)
         -b -- change type to run both simultaneously(3)
         -F -- change type to run forward calculation(1) with save_forward = .true.
END
exit(1);
}

@ARGV == 1 or Usage();
if(!getopts('abfF')) {die("check input arguments\n");}

open(IN,"DATA/Par_file");
@file_to_modify=<IN>;
close(IN);

foreach $file_to_modify (@file_to_modify){
  if($file_to_modify=~/SIMULATION_TYPE/){
    if(${opt_a}){
      print "Changed simulation_type to 2 in Par_file \n";
      $file_to_modify=~s/= 1/= 2/;
      $file_to_modify=~s/= 3/= 2/;
    }
    elsif(${opt_f}){
      print "Changed simulation_type to 1 in Par_file \n";
      $file_to_modify=~s/= 2/= 1/;
      $file_to_modify=~s/= 3/= 1/;
    }
    elsif(${opt_b}){
      print "Changed simulation_type to 3 in Par_file \n";
      $file_to_modify=~s/= 1/= 3/;
      $file_to_modify=~s/= 2/= 3/;
    }
    elsif(${opt_F}){
      print "Changed simulation_type to 1 and save_forward = .true. in Par_file \n";
      $file_to_modify=~s/= 2/= 1/;
      $file_to_modify=~s/= 3/= 1/;
    }
  }
  if ($file_to_modify=~/SAVE_FORWARD/) {
    if ($opt_F) { $file_to_modify=~s/= .*/= .true./; }
    else {$file_to_modify=~s/= .*/= .false./;}

  }
}

open(OUT,">DATA/Par_file") || die("No DATA/Par_file\n");
foreach $file_to_modify (@file_to_modify){
  print OUT "$file_to_modify";
}
close(OUT);
