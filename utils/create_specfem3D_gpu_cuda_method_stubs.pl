#!/usr/bin/perl

#
#  Script to extract the function declarations in cuda files
#
#
# usage: ./ceate_specfem3D_gpu_cuda_method_stubs.pl 
#             run in directory root SPECFEM3D/
#

$outfile = "src/cuda/specfem3D_gpu_cuda_method_stubs.c";


open(IOUT,"> _____temp_tutu_____");

$header = <<END;
/* 
!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  2 . 1
!               ---------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Princeton University, USA and University of Pau / CNRS / INRIA
! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
!                            April 2011
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================
*/

#include <stdio.h>
#include <stdlib.h>

#include "config.h"

typedef float realw;

END


$warning = <<END;
 fprintf(stderr,"ERROR: GPU_MODE enabled without GPU/CUDA Support. To enable GPU support, reconfigure with --with-cuda flag.\\n");
 exit(1);
END

print IOUT "$header \n";

$success = 0;

@objects = `ls src/cuda/*.cu`;

foreach $name (@objects) {  
  chop $name;
  print "extracting word in file $name ...\n";

  print IOUT "\n//\n// $name\n//\n\n";  
  
  # change tabs to white spaces
  system("expand -2 < $name > _____temp_tutu01_____");  
  open(IIN,"<_____temp_tutu01_____");

  
  # open the source file
  $success = 1;
  $do_extract = 0;
  while($line = <IIN>) {
    chop $line;
    
    # suppress trailing white spaces and carriage return
    $line =~ s/\s*$//;
    
    # change the version number and copyright information
    #    $line =~ s#\(c\) California Institute of Technology and University of Pau, October 2007#\(c\) California Institute of Technology and University of Pau, November 2007#og;
    #    $line =~ s#rmass_sigma#rmass_time_integral_of_sigma#og;
    
    if($line =~ /extern "C"/){
      # new function declaration starts  
      #print "$line\n";
      if( $line =~/FC_FUNC/ ){ 
        # function declaration on same line as extern, ask for line skip
        print "problem: please add a line break after extern 'C' here:";
        print "$line\n";
        $success = 0;
        close(IIN);  
        exit;
      }
      $do_extract = 1;
      next;          
    }
    
    # extract section
    if($do_extract == 1 ){
      # function declaration
      if($line =~ /{/){
        # function declaration ends
        if( $line =~ /INITIALIZE_CUDA_DEVICE/ ){
          # adds warning
          print IOUT "$line \n$warning\} \n\n";
        }else{
          print IOUT "$line\} \n\n";
        }
        $do_extract = 0;
      }else{
        # write line to the output file
        print IOUT "$line\n";  
      }
      next;
    }
  }
  close(IIN);  

  if( $success == 0 ){ exit; }
}

close(IOUT);
system("rm -f _____temp_tutu01_____");

# creates new stubs file if successful
if( $success == 1 ){
  print "\n\nsuccessfully extracted declarations \n\n";
  system("cp -p $outfile $outfile.bak");
  system("cp -p _____temp_tutu_____ $outfile");
  print "created new: $outfile \n";
}
system("rm -f _____temp_tutu_____");


