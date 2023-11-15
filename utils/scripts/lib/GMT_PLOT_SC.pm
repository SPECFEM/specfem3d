package GMT_PLOT_SC;

# Qinya Liu, Caltech, May 2007

use warnings;
use Exporter;
use lib '/opt/seismo-util/lib/perl';
use GMT_PLOT;
use SACLST;

@ISA = ('Exporter');
@EXPORT = ('plot_sc_faults','plot_sc_topo','plot_sc_all','plot_sc_pssac2_map');

$datalib_dir = "/opt/seismo-util/data/datalib_SC";
if (not -d $datalib_dir) {die("Check the datalib path\n");}

# plot faults
sub plot_sc_faults{

  if (@_ != 3) {die("plot_sc_faults filehandle,outfile options");}
  my($fh,$out,$opt) = @_;
  my($options,@options,$input);
  (@options) = split(" ",$opt);
  foreach $options (@options) {
    if ($options !~ /-/) {$input = $options;}
    else {push(@others,"$options");}
  }
  if (not $input or not -f $input) {$input = "${datalib_dir}/jennings.xy";}
  plot_psxy($fh,$out,"@others -M",$input);
}


# plot topo
sub plot_sc_topo{

  if (@_ != 3) {die("plot_sc_faults filehandle outfile options");}
  my($fh,$out,$opt) = @_;
  my($options,@options,$grd,$int,$cpt,@others);
  (@options) = split(" ",$opt);
  $grd = "${datalib_dir}/topo_cal.grd";
  $cpt = "${datalib_dir}/topo_cal_color.cpt";
  $int = "${datalib_dir}/topo_cal.int";
  foreach $options (@options) {
    if ($options !~ /-/)
      {$grd = $options;  if (not -f $grd) {die("Check $grd file\n");}}
    if ($options =~ /-I/)
      {$int = $options; $int =~s/^-.//;
       if (not -f $int) {die("Check $int file\n");}}
    if ($options =~ /-C/)
      {$cpt = $options; $cpt =~s/^-.//;
       if (not -f $cpt) {die("Check $cpt file\n");}}
    else {push(@others,"$options");}
  }

  plot_grdimage($fh,$out,"@others -I$int",$grd,$cpt);

}


# plot all
# -W gives the line attributes of coastline
# -w gives the line attributes of fault
# notice for labels,axis,and titles in -B, do not use 'space' or hyphen -,
# because it confuses the perl splitting and matching
sub plot_sc_all{

  if (@_ != 4) {die("plot_sc_all filehandle outfile options[JRBOICSW] mask\n");}
  my($fh,$out,$opt,$mask) = @_;
  my ($J,$R,$B,$O,$int,$cpt,$grd,$sea_color,$fault_color);
  my(@options,$options,@others);

  (@masks) = split(" ",$mask);
  if (@masks == 3) {
    if ($masks[0]) {$plot_topo = 1;}
    if ($masks[1]) {$plot_coast = 1;}
    if ($masks[2]) {$plot_faults = 1;}
  }else {$plot_topo = 1;$plot_coast = 1;$plot_faults = 1;}

  $J = "-JM6i";
  $R = "-R-120.3/-114.7/32.2/36.8";
  $B = "-B2/2";
  $O = "";
  $int = "";
  $cpt = "";
  $grd = "";
  $sea_color = "-S255";
  $coast_color = "-W4";
  $fault_color = "-W2.5";
  (@options) = split(" ",$opt);
  foreach $options (@options) {
    if ($options =~ /-J/) {$J = $options;}
    elsif ($options =~ /-R/) {$R = $options;}
    elsif ($options =~ /-B/) {$B = $options;}
    elsif ($options =~ /-O/ or $options =~ /-K/ ) {$O .= "$options ";}
    elsif ($options =~ /-I/) {$int = $options;}
    elsif ($options =~ /-C/) {$cpt = $options;}
    elsif ($options =~ /-S/) {$sea_color = $options;}
    elsif ($options =~ /-W/) {$coast_color = $options;}
    elsif ($options =~ /-w/) {$options =~s/-w/-W/; $fault_color = $options;}
    elsif ($options =~ /-g/) {$cpt = "-C/opt/seismo-util/data/datalib_SC/topo_cal_gray.cpt";}
    elsif ($options =~ /-/)  {push(@others,$options);}
    elsif ($options !~ /-/)  {$grd = $options;}
  }
  plot_psbasemap($fh,$out,"$J $R $B $O @others");
  if ($plot_topo) {plot_sc_topo($fh,$out,"$int $cpt $grd");}
  if ($plot_coast) {plot_pscoast($fh,$out,"-Na -Dh $coast_color $sea_color");}
  if ($plot_faults) {plot_sc_faults($fh,$out,"$fault_color");}
}


# default : -JM -R -Wdata -wsyn 
# options that has to be set -a[o,p,s,r,l]
# optional : -Ct1/t2 -C0/100
# this also checks if data and syn files are matching exactly
# in number and station location which may result in great difference
# in the result of plotting. To make sure everything is exact, we
# will cut syn and plot syn all according to data files instead of the
# sac headers.


sub plot_sc_pssac2_map {
  my($fh,$out,$opt,$data_files,$syn_files) = @_;
  my(@options,@data_files,@syn_files,$plot_syn,$J,$R,$O,$Wdata,$Wsyn,$E,$C,$S);
  my($options,$align,$cut1,$cut2,$cut_by_header,$h1,$h2,$h11);
  my($input,$stla,$stlo,$kstnm,$kcmpnm,$stla1,$stlo1,$kstnm1,$kcmpnm1);
  my($i,$data_file,$syn_file,@data_header_values,@syn_header_values);
  my($arrow,@others);

  (@options) = split(" ",$opt);
  (@data_files) = split(" ",$data_files);
  (@syn_files) = split(" ",$syn_files);
  if (@data_files == 0 and @syn_files == 0) {die("No data or syn files to plot\n");}
  if (@syn_files > 0) {
    $plot_syn = 1; if (@data_files != @syn_files)
      {die("Number of data and syn files need to match exactly\n"); }}
  else {$plot_syn=0;}

  $J = "-JM"; $R = "-R"; $Wdata = "-W21"; $Wsyn = "-W20/255/0/0"; $E="";$C="";
  $S = "";
  foreach $options (@options) {
    if ($options =~ /-J/) {$J = $options;}
    elsif ($options =~ /-R/) {$R = $options;}
    elsif ($options =~ /-O/ or $options =~ /-K/ ) {$O .= "$options ";}
    elsif ($options =~ /-W/) {$Wdata = $options;}
    elsif ($options =~ /-w/) {$options =~s/^-w/-W/; $Wsyn = $options;}
    elsif ($options =~ /-a/) {
      $align = $options; $align =~s/^-a//;
      if($align !~/^(o|p|s|r|l|\d*|\d*.\d*)$/) 
	{die "Align on o,p,s,r,l or number, not: $align \n";}
      if($align =~/o/){$E = "-Ent-3";}
      elsif($align =~/p/){$E = "-Ent1";}
      elsif($align =~/s/){$E = "-Ent2";}
      elsif($align =~/r/){$E = "-En3.8";}
      elsif($align =~/l/){$E = "-En4.4";}
      else{$E = "-En$align";}
    }elsif ($options =~ /-C/) {$C = $options;}
    elsif ($options =~ /-l/) {$l = $options;}
    elsif ($options =~ /-S/) {$S = $options;}
    elsif ($options =~ /-/) {push(@others,$options);}
  }
  if (not defined $O) {$O = "-K -O";}
  if ($O =~ /-K/ and $O !~ /-O/) {$arrow = ">";} else {$arrow = ">>";}

  $eps = 0.01;
  @data_header_values = sac_files_get_values
      ("stlo stla kstnm kcmpnm",@data_files);
  if ($plot_syn) {@syn_header_values = sac_files_get_values
		      ("stlo stla kstnm kcmpnm",@syn_files);}
  for ($i = 0; $i < @data_files; $i ++) {
    $data_file = $data_files[$i];
    ($stlo,$stla,$kstnm,$kcmpnm) = split(" ",$data_header_values[$i]);
    if ($plot_syn) {
      $syn_file = $syn_files[$i];
      ($stlo1,$stla1,$kstnm1,$kcmpnm1) = split(" ",$syn_header_values[$i]);
      if ($kstnm ne $kstnm1 or $kcmpnm ne $kcmpnm1) {die("Check station name and component name: $kstnm,$kcmpnm; $kstnm1,$kcmpnm1\n");}
      if (abs($stla-$stla1)>$eps or  abs($stlo-$stlo1)>$eps) {
	print("Check the WRONG station location for : $kstnm,$kcmpnm --- $stlo,$stla; $stlo1,$stla1\n"); return;}
      if ($stla ne $stla1 or $stlo ne $stlo1) {
	print "Check difference for $kstnm,$kcmpnm --- $stlo,$stla; $stlo1,$stla1\n";}
      $input .= "$syn_file $stlo $stla\n";
    }
   }
  plot_pssac2($fh,$out,"$J $R $Wdata $E $C $l @others",$data_files) or die "pssac2 failed";

  if ($plot_syn) {
    chomp($input);
    plot_pssac2_raw($fh,$out,"$J $R $Wsyn $E $C $S @others",$input) or die "pssac2 failed";}

}


1;
