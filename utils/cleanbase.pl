#!/usr/bin/perl -w

# this script cleans the whole scratch disk under /scratch/user/
# Qinya Liu, May 2007, Caltech

if (@ARGV != 1) {die("cleanbase.pl machinefile\n");}
$mymachine = $ARGV[0];
if (not -f $mymachine) {die("check if $mymachine is a file or not\n");}

`shmux -M50 -Sall -c "rm -rf /scratch/$ENV{USER}/*; mkdir /scratch/$ENV{USER}/DATABASES_MPI;" - < $mymachine `;

