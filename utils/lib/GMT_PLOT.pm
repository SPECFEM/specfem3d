package GMT_PLOT;

# Qinya Liu, Caltech, May 2007

use warnings;
use Exporter;
use lib '/opt/seismo-util/lib/perl';
use CMT_TOOLS;

@ISA = ('Exporter');
@EXPORT = ('plot_psbasemap','plot_pscoast','plot_psxy','plot_psxy_file',
	   'plot_pstext','plot_pstext_file','plot_grdimage','plot_psmeca',
	   'plot_psmeca_raw','plot_pssac2','plot_pssac2_raw',
	   'shift_xy');


# default paper orientation : landscape
$paper_orient = "";




# plot basemap
# default -JM -R -B
#

sub plot_psbasemap{
  my ($fh,$out,$opt) = @_;
  if (@_ != 3) {die("plot_psbasemap filehandle outfile options\n");}
  my ($J,$R,$B,$W,$O,$options,$arrow,@options,@others);
  @options = split(" ",$opt);
  $J = "-JM"; $R = "-R"; $B = "-B";
  foreach $options (@options) {
    if ($options =~ /-J/) {$J = $options;}
    elsif ($options =~ /-R/) {$R = $options;}
    elsif ($options =~ /-B/) {$B = $options;}
    elsif ($options =~ /-O/ or $options =~ /-K/ ) {$O .= "$options ";}
    else {push(@others,"$options ");}
  }
  if (not defined $O) {$O = "-K -O";}
  if ($O =~ /-K/ and $O !~ /-O/) {$arrow = ">";} else {$arrow = ">>";}
  print $fh "psbasemap $J $R $B @others $paper_orient $O -V   $arrow   $out\n"
    or die ("Error printing psbasemap\n");
}


# plot coastline
# default : -JM -R -W
# others : -Na -Dh -S255

sub plot_pscoast {

  my ($fh,$out,$opt) = @_;
  if (@_ != 3) {die("plot_pscoast filehandle outfile options\n");}
  my ($J,$R,$W,$O,$options,$arrow,@options,@others);
  @options = split(" ",$opt);
  $J = "-JM"; $R = "-R"; $W = "-W1";
  foreach $options (@options) {
    if ($options =~ /-J/) {$J = $options;}
    elsif ($options =~ /-R/) {$R = $options;}
    elsif ($options =~ /-O/ or $options =~ /-K/ ) {$O .= "$options ";}
    elsif ($options =~ /-W/) {$W = $options;}
    else {push(@others,"$options ");}
  }
  if (not defined $O) {$O = "-K -O";}
  if ($O =~ /-K/ and $O !~ /-O/) {$arrow = ">";} else {$arrow = ">>";}
  print $fh "pscoast $J $R $W @others $paper_orient $O -V   $arrow   $out\n"
    or die ("Error printing pscoast\n");
}



# plot psxy
# default -JM -R
# others -Sa0.2

sub plot_psxy {

  my ($fh,$out,$opt,$input) = @_;
  if (@_ != 4 ) {die("plot_psxy filehandle outfile options inputs \n");}
  my ($J,$R,$B,$W,$O,$options,$arrow,@options,@others);
  @options = split(" ",$opt);

  $J = "-JM"; $R = "-R";
  foreach $options (@options) {
    if ($options =~ /-J/) {$J = $options;}
    elsif ($options =~ /-R/) {$R = $options;}
    elsif ($options =~ /-O/ or $options =~ /-K/ ) {$O .= "$options ";}
    else {push(@others,"$options ");}
  }
  if (not defined $O) {$O = "-K -O";}
  if ($O =~ /-K/ and $O !~ /-O/) {$arrow = ">";} else {$arrow = ">>";}

  if ($input !~ /\n/ and -e $input) {
    print $fh "psxy $input $J $R @others $paper_orient $O -V   $arrow   $out\n";} else {
    if ($input =~ /^\s*$/) {$input = "";}
    else {$input = "$input\n";}
    print $fh "psxy $J $R @others $paper_orient $O -V   <<EOF  $arrow   $out\n${input}EOF\n" or die("Error printing psxy\n");}

}





sub plot_psxy_file {

  my($fh,$out,$opt,$input,$mask) = @_;
  if (@_ != 5) {die("plot_psxy_file filehandle outfile options inputs [mask]\n");}
  my ($J,$R,$B,$W,$O,$options,$arrow,@options,@others,@masks);
  @options = split(" ",$opt);
  if (defined $mask) {
    (@masks) = split(" ",$mask);
    if (@masks != 2) {die("psxy handle 2 entries\n");}}
  $J = "-JM"; $R = "-R";
  foreach $options (@options) {
    if ($options =~ /-J/) {$J = $options;}
    elsif ($options =~ /-R/) {$R = $options;}
    elsif ($options =~ /-O/ or $options =~ /-K/ ) {$O .= "$options ";}
    else {push(@others,"$options ");}
  }
  if (not defined $O) {$O = "-K -O";}
  if ($O =~ /-K/ and $O !~ /-O/) {$arrow = ">";} else {$arrow = ">>";}

  if (-f $input) {

    print $fh "awk '{print \$$masks[0],\$$masks[1]}' $input | psxy $J $R @others $paper_orient $O -V   $arrow   $out\n";
  } else {die("Check $input file is correct or not\n");}


}




# pstext
# necessary: -R -J
# optional :

sub plot_pstext {

  my ($fh,$out,$opt,$input) = @_;
  if (@_ != 4) {die("plot_pstext filehandle outfile options inputs\n");}
  my ($J,$R,$B,$W,$O,$options,$arrow,$dx,$dy,@options,@others,@masks);

  @options = split(" ",$opt);
  $J = "-JM"; $R = "-R";
  foreach $options (@options) {
    if ($options =~ /-J/) {$J = $options;}
    elsif ($options =~ /-R/) {$R = $options;}
    elsif ($options =~ /-O/ or $options =~ /-K/ ) {$O .= "$options ";}
    else {push(@others,"$options ");}
  }
  if (not defined $O) {$O = "-K -O";}
  if ($O =~ /-K/ and $O !~ /-O/) {$arrow = ">";} else {$arrow = ">>";}

  if ($input !~ /\n/ and -e $input) {
    print $fh "pstext $input $J $R @others $paper_orient $O -V   $arrow   $out\n";
  }else {
    print $fh "pstext $J $R @others $paper_orient $O -V   <<EOF  $arrow   $out\n$input \nEOF\n" or die("Error printing pstext\n");}

}


sub plot_pstext_file {

  my ($fh,$out,$opt,$input,$mask,$dxdy,$size) = @_;
  if (@_ != 6 and @_ != 7)
    {die("plot_pstext_file filehandle outfile options inputs $mask $dxdy $size\n");}
  my ($J,$R,$B,$W,$O,$options,$arrow,$dx,$dy,@options,@others,@masks);

  @options = split(" ",$opt);
  if (defined $mask) {
    (@masks) = split(" ",$mask);
    if (@masks != 3) {die("psxy handle 3 entries\n");}
    ($dx,$dy) = split(" ",$dxdy);
    if (not defined $dx or not defined $dy) {$dx = 0; $dy = 0;}
    if (not $size) {$size = "10 0 4 BL";}
  } else {die("Define $mask first\n");}

  $J = "-JM"; $R = "-R";
  foreach $options (@options) {
    if ($options =~ /-J/) {$J = $options;}
    elsif ($options =~ /-R/) {$R = $options;}
    elsif ($options =~ /-O/ or $options =~ /-K/ ) {$O .= "$options ";}
    else {push(@others,"$options ");}
  }
  if (not defined $O) {$O = "-K -O";}
  if ($O =~ /-K/ and $O !~ /-O/) {$arrow = ">";} else {$arrow = ">>";}

  if ($input !~ /\n/ and -e $input) {
    print $fh "awk '{print \$$masks[0]+$dx, \$$masks[1]+$dy, \"$size\", \$$masks[2]}' $input | pstext $J $R @others $paper_orient $O -V   $arrow   $out\n";}
  else {die("Check $input file is correct or not\n");}

}

# plot grdimage
# necessary: grdfile intfile
# optional: intfile

sub plot_grdimage {

  my ($fh,$out,$opt,$grd,$cpt) = @_;
  if (@_ != 5) {die("plot_grdimage filehandle outfile options grdfile cptfile\n");}
  my ($J,$R,$W,$O,$options,$arrow,@options,@others);
  @options = split(" ",$opt);
  $J = "-JM"; $R = "-R";
  foreach $options (@options) {
    if ($options =~ /-J/) {$J = $options;}
    elsif ($options =~ /-R/) {$R = $options;}
   elsif ($options =~ /-O/ or $options =~ /-K/ ) {$O .= "$options ";}
    else {push(@others,"$options ");}
  }
  if (not defined $O) {$O = "-K -O";}
  if ($O =~ /-K/ and $O !~ /-O/) {$arrow = ">";} else {$arrow = ">>";}

  print $fh "grdimage $grd $J $R -C$cpt @others $paper_orient $O -V   $arrow   $out\n"  or die ("Error printing grdimage\n");

}


# plot psmeca
# defaults: -JM -R -L
# others: -Sa/m -sdx/dy(shift the position to plot beach ball)

sub plot_psmeca{

  my ($fh,$out,$opt,@cmt_files) = @_;
  if (@_ < 3) {die("plot_psmeca filehandle outfile options [cmtfiles]\n");}
  my ($J,$R,$S,$O,$L,$options,$arrow,$moment,$mw,$temp_file,$input);
  my (@options,@others,@M,@sdr);
  @options = split(" ",$opt);
  $J = "-JM"; $R = "-R"; $L = "-L2"; my($add_name) = 0;
  foreach $options (@options) {
    if ($options =~ /-J/) {$J = $options;}
    elsif ($options =~ /-R/) {$R = $options;}
    elsif ($options =~ /-O/ or $options =~ /-K/ ) {$O .= "$options ";}
    elsif ($options =~ /-L/) {$L = $options;}
    elsif ($options =~ /-S/) {$S = $options;}
    elsif ($options =~ /-s/) {
      $s = $options; $s =~s/-s//;($dx,$dy) = split(/\//,$s);}
    elsif ($options =~ /-A/) {$add_name = 1;}
    elsif ($options !~ /^-/) {$input = $options;}
    else {push(@others,"$options ");}
  }
  if (not defined $O) {$O = "-K -O";}
  if ($O =~ /-K/ and $O !~ /-O/) {$arrow = ">";} else {$arrow = ">>";}
  if (not $dx or not $dy) {$dx = 0; $dy = 0;}

  if (not defined $input and not @cmt_files)
    {die("Define either input files or cmtsolution files\n");}
  if (not defined $S) {die("No -S is defined for psmeca\n");}

  if (not defined $input and @cmt_files) {
    print $fh "psmeca $J $R $S $L @others $paper_orient $O -V <<EOF $arrow $out \n";
    foreach $cmt_file (@cmt_files) {
      (@M[0..5],$elat,$elon,$edep) = get_cmt_short($cmt_file);
      $elat = $elat + $dy; $elon = $elon + $dx;
      if ($S =~ /-Sm/ or $S =~/-Sz/) {
	foreach $M (@M) {$M = $M /10;}
	print $fh "$elon $elat $edep @M[0..5] 1 0 0  ";
        if ($add_name) {print $fh " $cmt_file \n";} else {print $fh "\n";}
      } elsif ($S =~ /-Sa/) {
	($mw,$elat,$elon,$edep,@sdr[0..2]) = cmt2dc($cmt_file);
	if ($sdr[2] > 180) {$sdr[2] = $sdr[2] - 360;}
	print $fh "$elon $elat $edep @sdr[0..2] $mw 1 0 0 ";
        if ($add_name) {print $fh " $cmt_file \n";} else {print $fh "\n";}
      } else {die("plot_psmeca could only handle -Sm, -Sz and -Sa options\n");}
    }
    print $fh "EOF\n";
  }else {
    print $fh "psmeca $J $R $S $L @others $paper_orient $O -V $input $arrow $out \n" or die ("Error printing psmeca\n");}

}


# # specify the position of psmeca plots
sub plot_psmeca_raw{

  my ($fh,$out,$opt,$loc,$cmt_file) = @_;
  if (@_ != 5) {die("plot_psmeca_raw filehandle outfile options location cmtfile\n");}
  my ($J,$R,$S,$O,$L,$options,$arrow,$moment,$mw,$temp_file,$input);
  my (@options,@others,@M,@sdr);
  @options = split(" ",$opt);
  $J = "-JM"; $R = "-R"; $L = "-L2";
  ($xloc,$yloc) = split(" ",$loc);
  if (not defined $xloc) {$xloc = 0; $yloc = 0;}
  foreach $options (@options) {
    if ($options =~ /-J/) {$J = $options;}
    elsif ($options =~ /-R/) {$R = $options;}
    elsif ($options =~ /-O/ or $options =~ /-K/ ) {$O .= "$options ";}
    elsif ($options =~ /-L/) {$L = $options;}
    elsif ($options =~ /-S/) {$S = $options;}
    elsif ($options !~ /^-/) {die("no input file needed\n");}
    else {push(@others,"$options ");}
  }
  if (not defined $O) {$O = "-K -O";}
  if ($O =~ /-K/ and $O !~ /-O/) {$arrow = ">";} else {$arrow = ">>";}

  (@M[0..5],$elat,$elon,$edep) = get_cmt_short($cmt_file);
  if ($S =~ /-Sm/ or $S =~/-Sz/) {
    foreach $M (@M) {$M = $M /10;}
    $tmp =  "$xloc $yloc $edep @M[0..5] 1 0 0  \n";
  }elsif ($S =~ /-Sa/) {
    ($mw,$elat,$elon,$edep,@sdr[0..2]) = cmt2dc($cmt_file);
    if ($sdr[2] > 180) {$sdr[2] = $sdr[2] - 360;}
    $tmp = "$xloc $yloc $edep @sdr[0..2] $mw 1 0 0\n";
  } else {die("plot_psmeca could only handle -Sm, -Sz and -Sa options\n");}


  print $fh "psmeca $J $R $S $L @others $paper_orient $O -V <<EOF $arrow $out \n${tmp}EOF\n" or die ("Error printing psmeca\n");

}


# pssac2
# default : -JM -R -E -L -M
# others : -D<dx/dy> -S<shift_trace_by_seconds> -W<pen_attributes> -v(vertical)

sub plot_pssac2 {

  my ($fh,$out,$opt,$files) = @_;
  if (@_ < 4) {die("plot_pssac2 filehandle outfile options sacfiles\n");}
  my ($J,$R,$O,$options,$arrow,$E,$L,$M);
  my(@options,@others);

  @options = split(" ",$opt);
  $J = "-JM"; $R = "-R"; $E = "-Ekt-3"; $L = "-L100"; $M = "-M1e-5/0";
  foreach $options (@options) {
    if ($options =~ /^-J/) {$J = $options;}
    elsif ($options =~ /^-R/) {$R = $options;}
    elsif ($options =~ /^-O/ or $options =~ /-K/ ) {$O .= "$options ";}
    elsif ($options =~ /^-E/) {$E = $options;}
    elsif ($options =~ /^-L/) {$L = $options;}
    elsif ($options =~ /^-M/) {$M = $options;}
    else {push(@others,"$options ");}
  }
  if (not defined $O) {$O = "-K -O";}
  if ($O =~ /-K/ and $O !~ /-O/) {$arrow = ">";} else {$arrow = ">>";}
  print $fh "pssac2 $J $R $E $L $M @others $paper_orient $O -V $files $arrow  $out\n"
    or die ("Error printing pssac2\n");
}


# specifying pssac2 by specifying the location of plotted files
sub plot_pssac2_raw {

  my ($fh,$out,$opt,$file) = @_;
  if (@_ < 4) {die("plot_pssac2_raw filehandle outfile options sacfiles\n");}
  my ($J,$R,$O,$options,$arrow,$E,$L,$M);
  my(@options,@others);

  @options = split(" ",$opt);
  $J = "-JM"; $R = "-R"; $E = "-Ekt-3"; $L = "-L100"; $M = "-M1e-5/0";
  foreach $options (@options) {
    if ($options =~ /-J/) {$J = $options;}
    elsif ($options =~ /-R/) {$R = $options;}
    elsif ($options =~ /-O/ or $options =~ /-K/ ) {$O .= "$options ";}
    elsif ($options =~ /-E/) {$E = $options;}
    elsif ($options =~ /-L/) {$L = $options;}
    elsif ($options =~ /-M/) {$M = $options;}
    else {push(@others,"$options ");}
  }
  if (not defined $O) {$O = "-K -O";}
  if ($O =~ /-K/ and $O !~ /-O/) {$arrow = ">";} else {$arrow = ">>";}
  if ($file !~ /\s/ and -e $file) {
    print $fh "pssac2 $J $R $E $L $M @others $paper_orient $O -V < $file $arrow  $out\n"  or die ("Error printing pssac2\n");}
  else {
    print $fh "pssac2 $J $R $E $L $M @others $paper_orient $O -V <<EOF $arrow  $out\n$file\nEOF\n"  or die ("Error printing pssac2\n");
  }

}

## shift the location of -X -Y for subplot
# ($jx,$jy) = shift_xy("1 1","3 5","7 6",\@xy,"");
sub shift_xy {

  my($xy,$nxy,$dxy,$rxx,$dxy1) = @_;
  if (@_ < 4) {die("Check shift_xy(xy,nxy,dxy,rxx,[dxy1])\n");}

  my($x0,$y0,$dx,$dy,$dx1,$dy1,$jx,$jy,$ddx,$ddy,$tmp1,$tmp2);
  ($x0,$y0) = split(" ",$xy); if (not ($x0 and $y0)) {die("specify -X/Y\n");}
  ($nx,$ny) = split(" ",$nxy); if (not ($nx and $ny)) {die(" -Nx/Ny\n");}
  ($dx,$dy) = split(" ",$dxy);
  if (not ($dx and $dy)) {$dx = 7; $dy = 7; }
  if (defined $dxy1) {($dx1,$dy1) = split(" ",$dxy1);}
  if (not ($dx1 and $dy1)) {$dx1 = 1; $dy1 = 1;}

  $jx = $dx/$nx * $dx1; $jy = $dy/$ny * $dy1;
  $ddx = $dx/$nx; $ddy = $dy/$ny;

  for ($j=0;$j<$ny;$j++) { # rows
    for ($i=0;$i<$nx;$i++) { # columns
      $tmp1 = $x0 + $ddx * $i; $tmp2 = $y0 + $ddy * ($ny-1-$j);
      ${$rxx}[$j][$i] = "$tmp1/$tmp2";
    }
  }
  return($jx,$jy);
}
1;
