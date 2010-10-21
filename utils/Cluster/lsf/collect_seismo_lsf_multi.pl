#!/usr/bin/perl -w

# this script collects the seismograms from the scratch sub-directory given by Par_file
# Qinya Liu, May 2007, Caltech

if (@ARGV != 2) {die("collect_seismo_lsf_multi.pl lsf_machinefile Par_file\n");}

$machine = $ARGV[0];
$par_file = $ARGV[1];

open(FILE,"$machine") or die("Error opening file $machine\n");
($junk1) = <FILE>;
close(FILE);
(@junk2) = split(" ",$junk1);
for($i=0;$i<@junk2/2;$i++) {
$junk[$i] = $junk2[2*$i];}
  

for($i=0;$i<@junk;$i++) {
  ($node) = split(" ",$junk[$i]);
  push(@nodes,$node);
}

# now get the LOCAL_PATH
open(FILE3,"<$par_file") or die ("Fatal Error openning file $par_file\n");
while (<FILE3>) {
   if ($_ =~ /^LOCAL_PATH/) {
	chop;	
	@vals = split("=", $_);
	$mpidir = $vals[1];
	$mpidir =~ s/^\s+//;
	$mpidir =~ s/\s+$//;
	close(FILE3);
	last;
   }
}

open(FILE2,">procs.list");
for ($i=0;$i<@nodes;$i++) {print FILE2 "$nodes[$i]\n";}
close(FILE2);

print "@nodes\n";

foreach $node (@nodes) {
    system("scp $node:$mpidir/*sem* .");
    print "$node\n";}

#`shmux -M50 -Sall -c "rm -f $mpidir/*sem?" - < $machine`;

