#!/usr/bin/perl

#
#  Script to change the version number in f90 codes
#
#  Author : Dimitri Komatitsch, EPS - Harvard University, May 1998
#

#
# read all f90 and F90 (and *.h) files in the current directory
# f90 files are supposed to have the extension "*.f90" or "*.F90" or "*.h"
#

#
# known bug : does the changes also in constant strings, and not only
# in the code (which is dangerous, but really easier to program...)
#

#
# usage: ./update_headers_change_word_f90.pl 
#             run in directory root SPECFEM3D/
#


@objects = `ls src/*/*.f90 src/*/*.F90 src/*/*.h.in src/*/*.h src/*/*.c src/*/*.cu`;

$nlines_total = 0;
$nlines_noblank = 0;
$nlines_nocomment = 0;

foreach $name (@objects) {
  chop $name;

  
  # change tabs to white spaces
  if( 1 == 1 ){
    system("expand -2 < $name > _____tutu01_____");
    $f90name = $name;
    print STDOUT "Changing word in file $name ...\n";

    open(FILEF77,"<_____tutu01_____");
    open(FILEF90,">$f90name");

    # open the source file
    while($line = <FILEF77>) {
      chop $line;

      # suppress trailing white spaces and carriage return
      $line =~ s/\s*$//;

      # change the version number and copyright information
      #    $line =~ s#\(c\) California Institute of Technology and University of Pau, October 2007#\(c\) California Institute of Technology and University of Pau, November 2007#og;
      #    $line =~ s#rmass_sigma#rmass_time_integral_of_sigma#og;

      # write the modified line to the output file
      print FILEF90 "$line\n";

    }

    close(FILEF77);
    close(FILEF90);
  }

  # line count
  if( 1 == 0 ){
    print STDOUT "file $name : \n";
    
    # counts all lines in file
    $l = `wc -l $name | awk '{print \$1}'`;
    chomp $l;
    print "  lines = $l  \n";

    # without blank lines
    $lb = 0;
    # without comments
    $lc = 0;
    system("expand -2 < $name > _____tutu01_____");
    open(FILEF77,"<_____tutu01_____");
    # open the source file
    while($line = <FILEF77>) {
      chop $line;      
      chomp $line;
      # remove whitespace at start
      $line =~ s/^\s+//;
      # remove whitespace at end
      $line =~ s/\s+$//;
      # write the modified line to the output file
      if( $line ne ""){ $lb = $lb + 1; }
      if( ($line ne "") && (substr($line,0,1) ne "!") && (substr($line,0,1) ne "/") ){ $lc = $lc + 1; }
    }
    close(FILEF77);    
    print "  lines (no blank) = $lb \n";
    print "  lines (no comment) = $lc  \n";

    # summations
    $nlines_total = $nlines_total + $l;
    $nlines_noblank = $nlines_noblank + $lb;
    $nlines_nocomment = $nlines_nocomment + $lc;
    print "  total = $nlines_total \n";
    print "  total (no blank) = $nlines_noblank \n";
    print "  total (no comment) = $nlines_nocomment \n";
  }  
}

#line count output
if( 1 == 0 ){
  print "\ntotal number of lines: \n";
  print "  lines = $nlines_total \n";
  print "  lines (no blank) = $nlines_noblank \n";
  print "  lines (no comment) = $nlines_nocomment \n\n";
}

system("rm -f _____tutu01_____");

