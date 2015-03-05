#!/usr/bin/perl

use Time::Local;
use Getopt::Std;
use POSIX;

sub Usage{
print STDERR <<END;

Usage:  cut_seis.pl -d SynDir data_files
  where
        -d SynDir -- directory of synthetic seismogram, like SYN/CMTSOLUTION200904060132A
        data_files -- name of data files to be processed (SAC format), the data like sta.net.comp.sac.bandpass
END
exit(1);
}

@ARGV > 1 or Usage();

if (!getopts('f:cd:l')) {die(" check input arguments\n");}

if ($opt_d and not -d $opt_d) {die("No such directory as $opt_d\n");}


foreach $file_dat (@ARGV) {
  if (! -f $file_dat) {die("No such file: $file_dat\n");}

  ($datname) = split(" ",`basename $file_dat`);
  ($sta,$net,$comp,undef,$ext)=split(/\./,$datname);
  $cmp=substr($comp,2,3);

  $file_syn="$opt_d/$sta.$net.LH$cmp.sem.sac.$ext";

  if ( -f $file_syn ) {
    print "Cutting file $file_syn\n";
    open(SAC,"|sac > /dev/null");
        print SAC "r $file_dat $file_syn \n";
        print SAC "cut off \n";
        print SAC "eval to cutbeg ( (max &1,b &2,b ) ) \n";
        print SAC "eval to cutend ( (min &1,e &2,e ) ) \n";
        print SAC "cut %cutbeg %cutend \n";
        print SAC "cut on \n";
        print SAC "r \n";
        print SAC "ch SCALE 1 \n";
        print SAC "w $file_dat $file_syn  \nquit\n";
    close(SAC);
  }
}
