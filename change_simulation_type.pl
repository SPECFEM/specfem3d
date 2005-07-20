#!/usr/bin/perl

use Time::Local;
use Getopt::Std;
use POSIX;

sub Usage{
print STDERR <<END;

Usage:   change_simulation_type.pl  [-a|-f|-b|-F]
         Changes SIMULATION_TYPE in OUTPUT_FILES/values_from_mesher.h
         -a -- change type to run adjoint calculation(2)
         -f -- change type to run forward calculation(1)
         -b -- change type to run both simultaneously(3)
         -F -- change type to run forward calculation(1) with save_forward = .true.
END
exit(1);
}

@ARGV == 1 or Usage();
if(!getopts('abfF')) {die("check input arguments\n");}

open(IN,"OUTPUT_FILES/values_from_mesher.h");
@vfm=<IN>;
close(IN);

foreach $vfm (@vfm){
  if($vfm=~/SIMULATION_TYPE/){
    if(${opt_a}){
      print "Changed simulation_type to 2 in values_from_mesher.h \n";
      $vfm=~s/= 1/= 2/;
      $vfm=~s/= 3/= 2/;
    }
    elsif(${opt_f}){
      print "Changed simulation_type to 1 in values_from_mesher.h \n";
      $vfm=~s/= 2/= 1/;
      $vfm=~s/= 3/= 1/;
    }
    elsif(${opt_b}){
      print "Changed simulation_type to 3 in values_from_mesher.h \n";
      $vfm=~s/= 1/= 3/;
      $vfm=~s/= 2/= 3/;
    }
    elsif(${opt_F}){
      print "Changed simulation_type to 1 and save_forward = .true. in values_from_mesher.h \n";
      $vfm=~s/= 2/= 1/;
      $vfm=~s/= 3/= 1/;
    }
  }
  if ($vfm=~/SAVE_FORWARD/) {
    if ($opt_F) { $vfm=~s/= .*/= .true./; }
    else {$vfm=~s/= .*/= .false./;}

  }
}

open(OUT,">OUTPUT_FILES/values_from_mesher.h") || die("No OUTPUT_FILES/values_from_mesher.h\n");
foreach $vfm (@vfm){
  print OUT "$vfm";
}
close(OUT);
