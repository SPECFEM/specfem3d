#!/usr/bin/perl -w

open(FILE,"stalist.TK");
@info = <FILE>;
close(FILE);
%stalist = ();

open(FILE,">new_stalist.TK");
foreach $line (@info) {
  ($sta,@else) = split(" ",$line);
  if ($else[0] !~/^\d/) {
    ($sta,$net) = split(" ",$line);}
  else {$net = "CI";}
  $stalist{$net}.=" $sta";
  $sta =~tr/[a-z]/[A-Z]/;
  printf FILE ("%-5s%-6s@else\n",$sta,$net);
}
close(FILE);

foreach $keyword (keys %stalist) {
  print "$keyword -- $stalist{$keyword} \n";
}

