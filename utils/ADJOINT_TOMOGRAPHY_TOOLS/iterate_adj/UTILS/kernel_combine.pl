#!/usr/bin/perl -w

#==========================================================
#
#  kernel_combine.pl
#  Carl Tape
#  19-June-2008
#
# EXAMPLE:
#    ~/UTILS/kernel_combine.pl 5 7 mu_all_kernel_set
#
# combines pdf files for a kernel plot
#==========================================================

if (@ARGV < 3) {die("Usage: kernel_combine.pl imodel_min imodel_max suffix\n")}
($imodel_min,$imodel_max,$ksuffix) = @ARGV;

# base directory
$dir0 = "/net/sierra/raid1/carltape/results/KERNELS";
$smodel = sprintf("m%2.2i",$imodel_max);
$kdirmax = "$dir0/kernel_${smodel}";
if (not -e $kdirmax) {die("check if kdirmax $kdirmax exist or not\n")}

# list of event IDs
#$file_eids = "/net/sierra/raid1/carltape/results/EID_LISTS/kernels_use_m05";
#if (not -f $file_eids) {die("\n check if $file_eids exists\n")}
#open(IN,$file_eids); @eids = <IN>; $numk = @eids;

@kfiles = glob("$kdirmax/*${ksuffix}.pdf");
$numk = @kfiles;
print "$numk kernel files in $kdirmax\n";

# final output file
$outfile = "$kdirmax/${ksuffix}_ALL.pdf";

#======================
# concatenate the PDF files into one document

# concatenate pdf files
@pdcat = "/home/carltape/bin/pdcat -r";
$p = 1;

for ($k = 1; $k <= $numk; $k++) {

   $kfilemax = $kfiles[$k-1]; chomp($kfilemax);
   $ktag = `basename $kfilemax`; chomp($ktag);
   print "$kfilemax -- $ktag --\n";

   for ($m = $imodel_min; $m <= $imodel_max; $m++) {
       $smod = sprintf("m%2.2i",$m);
       $kfilem = "$dir0/kernel_${smod}/$ktag";
       if (not -f $kfilem) {
          print "check if kfilem $kfilem exist or not\n";
       } else {
          $pdcat[$p] = $kfilem; $p = $p+1;
       }
   }
}

#======================

#`rm $outfile`;   # remove the output file if it exists
$pdcat[$p] = "$outfile";

print "\n @pdcat \n";
`@pdcat`;

print "\n ";

#=================================================================
