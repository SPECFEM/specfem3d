#!/usr/bin/perl -w

#--------------------------------------------------------------------
# combine_2_adj_src.pl
# Qinya Liu, 05-May-2007
#
# This script combines two sets of adjoint sources (for example, P wave time windows
# and surface wave time windows), and outputs the new adjoint sources, along with
# a new STATIONS file to be used.
#
# To remove the cc label on files:
#   mv_files.pl "*.cc.adj" "*.adj"
#
# CALLED BY: combine_2_adj_src_all.pl
#
#--------------------------------------------------------------------

if (@ARGV != 5) {die("Usage: combine_2_adj_src.pl DIR_1 DIR_2 DIR_NEW type1(cc,ik,mtm) type2(cc,ik,mtm)\n");}

$dir1 = $ARGV[0];
$dir2 = $ARGV[1];

$dir_new = $ARGV[2];
$type1 = $ARGV[3];
$type2 = $ARGV[4];

$sta1 = "$dir1/STATIONS_ADJOINT";
$sta2 = "$dir2/STATIONS_ADJOINT";
$sta_new = "$dir_new/STATIONS_ADJOINT";
if (not -f $sta1 or not -f $sta2) {die("Check if $sta1 and $sta2 exist or not\n");}
if (not -d $dir_new) {system("mkdir $dir_new");}

# open STATIONS file 1
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

# open STATIONS file 2
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


open(FILE,">$dir_new/sta.tmp"); $nsta = 0;
foreach $sta (keys %adj) {
#foreach $sta ("TCF") {
  $net = $adj{$sta}{net}; $nsta ++ ;
  print "$nsta, $sta, $net\n";
  print FILE "$sta $adj{$sta}{info} \n";
  if ($adj{$sta}{num} > 1) {
    # add up two adjoint sources
    foreach $comp ("BHE","BHN","BHZ") {
      system("paste $dir1/$sta.$net.$comp.$type1.adj $dir2/$sta.$net.$comp.$type2.adj | awk '{print \$1, \$2+\$4}' > $dir_new/$sta.$net.$comp.adj");
#      print "paste $dir1/$sta.$net.$comp.$type1.adj $dir2/$sta.$net.$comp.$type2.adj | awk '{print \$1, \$2+\$4}' > $dir_new/$sta.$net.$comp.adj\n";
    }
  }else{
    # copy over the adjoint source
    foreach $comp ("BHE","BHN","BHZ") {
       system("cp -f $dir1/$sta.$net.$comp.$type1.adj  $dir_new/$sta.$net.$comp.adj >& /dev/null");
       system("cp -f $dir2/$sta.$net.$comp.$type2.adj  $dir_new/$sta.$net.$comp.adj >& /dev/null");
     }
  }
}
close(FILE);

print "\n number of stations in new STATIONS file should be $nsta\n";
system(" echo $nsta > $sta_new; cat $dir_new/sta.tmp >> $sta_new");

#----------------------------
