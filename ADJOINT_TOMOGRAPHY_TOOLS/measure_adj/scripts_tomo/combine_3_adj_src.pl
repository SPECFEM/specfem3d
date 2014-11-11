#!/usr/bin/perl -w

#--------------------------------------------------------------------
# combine_3_adj_src.pl
# Qinya Liu, Carl Tape, 03-Oct-2008
#
# This script combines three sets of adjoint sources (for example, three
# different bandpassed ranges), and outputs the new adjoint sources,
# along with a new STATIONS_ADJOINT file to be used.
#
# To remove the cc label on files:
#   mv_files.pl "*.cc.adj" "*.adj"
#
# CALLED BY: combine_3_adj_src_all.pl
#
#--------------------------------------------------------------------

if (@ARGV != 7) {die("Usage: combine_3_adj_src.pl DIR_1 DIR_2 DIR_3 DIR_NEW type1(cc,ik,mtm) type2 type3\n");}

$dir1 = $ARGV[0];
$dir2 = $ARGV[1];
$dir3 = $ARGV[2];

$dir_new = $ARGV[3];
$type1 = $ARGV[4];
$type2 = $ARGV[5];
$type3 = $ARGV[6];

$sta1 = "$dir1/STATIONS_ADJOINT";
$sta2 = "$dir2/STATIONS_ADJOINT";
$sta3 = "$dir3/STATIONS_ADJOINT";
$sta_new = "$dir_new/STATIONS_ADJOINT";
if (not -f $sta1) {die("Check if $sta1 exist or not\n");}
if (not -f $sta2) {die("Check if $sta2 exist or not\n");}
if (not -f $sta3) {die("Check if $sta3 exist or not\n");}
if (not -d $dir_new) {system("mkdir $dir_new");}

#--------------

# open STATIONS file 1 -- this becomes the initial STATIONS list
open(FILE,"$sta1");
(@sta) = <FILE>;    # first line is the number of stations
close(FILE);
for ($i=1; $i<@sta; $i++){
  ($sta,$nt,$lat,$lon,$dep,$bur) = split(" ",$sta[$i]);
  $adj{$sta}{num} = 1;
  $adj{$sta}{net} = "$nt";
  $adj{$sta}{info} = "$nt $lat $lon $dep $bur";
  @files = glob("$dir1/$sta.$nt.BH?.*adj");
  if (@files == 0) {print("Check if $dir1/$sta.$nt.BH?.*adj exist or not\n");}
}

# open STATIONS file 2 -- update current STATIONS list with new stations
open(FILE,"$sta2");
(@sta) = <FILE>;    # first line is the number of stations
close(FILE);
for ($i=1; $i<@sta; $i++){
  ($sta,$nt,$lat,$lon,$dep,$bur) = split(" ",$sta[$i]);
  if (not defined $adj{$sta}{num}) {
    $adj{$sta}{num} = 1;
    $adj{$sta}{net} = "$nt";
    $adj{$sta}{info} = "$nt $lat $lon $dep $bur";}
  else {$adj{$sta}{num} ++ ; 
        # NOTE: this crashes if HAST.TA and HAST.BK are included
	if ($nt ne $adj{$sta}{net}) {die("Check if network name same for $sta\n");}
  }
  @files = glob("$dir2/$sta.$nt.BH?.*adj");
  if (@files == 0) {print("Check if $dir2/$sta.$nt.BH?.*adj exist or not\n");}
}

# open STATIONS file 3 -- update current STATIONS list with new stations
open(FILE,"$sta3");
(@sta) = <FILE>;    # first line is the number of stations
close(FILE);
for ($i=1; $i<@sta; $i++){
  ($sta,$nt,$lat,$lon,$dep,$bur) = split(" ",$sta[$i]);
  if (not defined $adj{$sta}{num}) {
    $adj{$sta}{num} = 1;
    $adj{$sta}{net} = "$nt";
    $adj{$sta}{info} = "$nt $lat $lon $dep $bur";}
  else {$adj{$sta}{num} ++ ; 
        # NOTE: this crashes if HAST.TA and HAST.BK are included
	if ($nt ne $adj{$sta}{net}) {die("Check if network name same for $sta\n");}
  }
  @files = glob("$dir3/$sta.$nt.BH?.*adj");
  if (@files == 0) {print("Check if $dir3/$sta.$nt.BH?.*adj exist or not\n");}
}

#--------------
# make new STATIONS_ADJOINT file, and sum the adjoint sources

open(FILE,">$dir_new/sta.tmp"); $nsta = 0;
foreach $sta (keys %adj) {
#foreach $sta ("TCF") {
  $net = $adj{$sta}{net}; $nsta ++ ;
  $nmatch = $adj{$sta}{num};
  print "$nsta, $sta, $net, $nmatch\n";
  print FILE "$sta $adj{$sta}{info} \n";

  foreach $comp ("BHE","BHN","BHZ") {

    $snc = "$sta.$net.$comp";
    $afile1 = "$dir1/$snc.$type1.adj";
    $afile2 = "$dir2/$snc.$type2.adj";
    $afile3 = "$dir3/$snc.$type3.adj";

    if ($nmatch == 3) {		# sum the three adjoint sources
      system("paste $afile1 $afile2 $afile3 | awk '{print \$1, \$2+\$4+\$6}' > $dir_new/$snc.adj");
    
    } elsif ($nmatch == 2) {	# sum the two adjoint sources
      if (-f $afile1 && -f $afile2) {
        system("paste $afile1 $afile2 | awk '{print \$1, \$2+\$4}' > $dir_new/$snc.adj");
      } elsif (-f $afile1 && -f $afile3) {
	system("paste $afile1 $afile3 | awk '{print \$1, \$2+\$4}' > $dir_new/$snc.adj");
      } elsif (-f $afile2 && -f $afile3) {
	system("paste $afile2 $afile3 | awk '{print \$1, \$2+\$4}' > $dir_new/$snc.adj");
      }
    
    } elsif ($nmatch == 1) {    # copy over the lone-station adjoint sources
      if (-f $afile1) {system("cp -f $afile1  $dir_new/$snc.adj >& /dev/null");}
      if (-f $afile2) {system("cp -f $afile2  $dir_new/$snc.adj >& /dev/null");}
      if (-f $afile3) {system("cp -f $afile3  $dir_new/$snc.adj >& /dev/null");}

    } else {			# error
      die("combine_3_adj_src.pl: NO MATCHES\n");
    }
  }
}
close(FILE);

print "\n number of stations in new STATIONS file should be $nsta\n";
system(" echo $nsta > $sta_new; cat $dir_new/sta.tmp >> $sta_new");

#----------------------------
