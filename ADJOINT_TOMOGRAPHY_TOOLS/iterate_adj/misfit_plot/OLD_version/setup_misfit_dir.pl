#!/usr/bin/perl -w

#-----------------------------------
# Carl Tape, 03-Sept-2008
# setup_misfit_dir.pl
#
# This script xxx
#
# EXAMPLE: setup_misfit_dir.pl m00 m07 6/40
#
#-----------------------------------

if (@ARGV < 3) {die("Usage: setup_misfit_dir.pl smodel_min smodel_max Tmin/Tmax\n")}
($minm,$maxm,$Ts) = @ARGV;

$pwd = $ENV{PWD};
print "\n$pwd\n";

# labels for bandpass-filtered data and synthetics
($Tmin,$Tmax) = split("/",$Ts);
$sTmin = sprintf("T%3.3i",$Tmin);
$sTmax = sprintf("T%3.3i",$Tmax);
$Ttag = "${sTmin}_${sTmax}";

#-------------------------------------------
# USER INPUT

# directory containing all CMTSOLUTION files
$dir_src = "/net/sierra/raid1/carltape/results/SOURCES/socal_09";
$dir_cmt = "${dir_src}/CMT_files_pre_inverted";
if (not -e $dir_src) {die("check if dir_src $dir_src exist or not\n")}
if (not -e $dir_cmt) {die("check if dir_cmt $dir_cmt exist or not\n")}

# list of event IDs to use
$file_eid0 = "/net/sierra/raid1/carltape/results/EID_LISTS/syn_run_m07";
if (not -f ${file_eid0}) {die("check if file_eid ${file_eid0} exist or not\n")}

# FULL stations file
$file_stations0 = "/net/denali/home2/carltape/gmt/stations/seismic/Matlab_output/STATIONS_CALIFORNIA_TOMO_INNER_specfem";
if (not -f ${file_stations0}) {die("check if file_stations ${file_stations0} exist or not\n")}

# data and synthetics directories
$dir_data0  = "/net/sierra/raid1/carltape/socal/socal_3D/DATA/FINAL";
$dir_syn01 = "/net/sierra/raid1/carltape/socal/socal_3D/SYN/model_${minm}";
$dir_syn02 = "/net/sierra/raid1/carltape/socal/socal_3D/SYN/model_${maxm}";
if (not -e ${dir_data0}) {die("check if dir_data ${dir_data0} exist or not\n")}
if (not -e ${dir_syn01}) {die("check if dir_syn ${dir_syn01} exist or not\n")}
if (not -e ${dir_syn02}) {die("check if dir_syn ${dir_syn02} exist or not\n")}

# directory containing the output for windows, measurements, adjoint sources, etc
$dir_run = "/net/sierra/raid1/carltape/socal/socal_3D/RUNS";
if (not -e ${dir_run}) {die("check if dir_run ${dir_run} exist or not\n")}

#-------------------------------------------

# make the base output directory if it does not exist
$dir0 = "OUTPUT_MISFIT";
`mkdir -p $dir0`;
$dir1 = "$dir0/EVENTS";
`mkdir -p $dir1`;

# link various files
$file_eid = "$dir0/eid_list";
`ln -s ${file_eid0} ${file_eid}`;
$file_stations = "$dir0/STATIONS";
`ln -s ${file_stations0} ${file_stations}`;

# open EID list
open(IN,"${file_eid}"); @eids = <IN>; close(IN);
$nevent = @eids;
print "\n $nevent events in the list";

# open STATIONS file
open(IN,"${file_stations}"); @slines = <IN>; close(IN);
$nrec = @slines - 1;
print "\n $nrec stations in the list";

#-------------------------------------------

if (0==1) {
  for ($ievent = 1; $ievent <= $nevent; $ievent++) {
    $eid = $eids[$ievent-1]; chomp($eid);
    print "\n $ievent : $eid";
  }
  die("testing");
}

#=============================================
# MAKE DIRECTORIES

#$imin = 1; $imax = $nevent;
#$imin = 1; $imax = 100;
$imin = 142; $imax = $imin;

$isetup1 = 0;
if ($isetup1 == 1) {

  for ($ievent = $imin; $ievent <= $imax; $ievent++) {
    $eid = $eids[$ievent-1]; chomp($eid);
    print "\n $ievent : $eid";

    # make base directory for each event
    $dir2 = "$dir1/$eid";
    `mkdir -p $dir2`;

    # link CMTSOLUTION file
    $cmtfile0 = "${dir_cmt}/CMTSOLUTION_${eid}";
    $cmtfile = "$dir2/CMTSOLUTION";
    `ln -s $cmtfile0 $cmtfile`;

    # link STATIONS file from the windowing/measurement codes (different for each period range)
    # --> use the max model to show the stations used
    $stafile0 = "${dir_run}/$eid/$maxm/MEASURE_${sTmin}/ADJOINT_${sTmin}/STATIONS_ADJOINT";
    if (not -f $stafile0) {die("check if stafile0 $stafile0 exist or not\n")}
    $stafile = "$dir2/STATIONS_ADJOINT";
    `ln -s $stafile0 $stafile`;

    # make a directory for each station in the STATIONS list
    for ($irec = 1; $irec <= $nrec; $irec++) {
      ($sta,$net,$slat,$slon,undef,undef) = split(" ",$slines[$irec]);
      $dir_rec = "$dir2/$sta.$net";
      print "\n $irec out of $nrec -- $dir_rec";
      `mkdir -p $dir_rec`;
    }
  }

}

#=============================================
# COPY IN DATA AND SYNTHETICS
# NOTE: For now we are copying in all available records,
# but really we only want BHZ,BHR,BHT (or HHZ,HHR,HHT).

$isetup2 = 0;

if ($isetup2 == 1) {
  for ($ievent = $imin; $ievent <= $imax; $ievent++) {

    $eid = $eids[$ievent-1]; chomp($eid);
    print "\n $ievent : $eid";

    # data and synthetics directories for each event
    $dir_data  = "${dir_data0}/$eid/PROCESSED/PROCESSED_${Ttag}";
    $dir_syn1  = "${dir_syn01}/$eid/PROCESSED/PROCESSED_${Ttag}";
    $dir_syn2  = "${dir_syn02}/$eid/PROCESSED/PROCESSED_${Ttag}";
    if (not -e ${dir_data}) {die("check if dir_data ${dir_data} exist or not\n");}
    if (not -e ${dir_syn1}) {die("check if dir_syn ${dir_syn1} exist or not\n");}
    if (not -e ${dir_syn2}) {die("check if dir_syn ${dir_syn2} exist or not\n");}

    $dir2 = "$dir1/$eid";

    for ($irec = 1; $irec <= $nrec; $irec++) {
    #for ($irec = 128; $irec <= 128; $irec++) {

      cd_directory($pwd);
      ($sta,$net,$slat,$slon,undef,undef) = split(" ",$slines[$irec]);
      $station = "$sta.$net";
      $dir_rec = "${dir2}/${sta}.${net}";
      print "\n $irec out of $nrec -- ${dir_rec}";

      # link plotting script
      $plotfile0 = "${pwd}/plot_misfit.pl";
      $plotfile = "${dir_rec}/plot_misfit.pl";
      `rm $plotfile`;
      `ln -s $plotfile0 $plotfile`;

      $dtags = "${dir_data}/*$net.$sta*"; @dfiles = glob($dtags); `cp @dfiles ${dir_rec}`;

      $stags1 = "${dir_syn1}/$sta.$net*"; @sfiles1 = glob($stags1); #`cp @sfiles1 ${dir_rec}`;
      cd_directory($dir_rec);
      `cp @sfiles1 .`;
      `mv_files.pl -x "*semd.sac.d.${Ttag}" "*semd.sac.d.${Ttag}.$minm"`;

      $stags2 = "${dir_syn2}/$sta.$net*"; @sfiles2 = glob($stags2); #`cp @sfiles2 ${dir_rec}`;
      `cp @sfiles2 .`;
      `mv_files.pl -x "*semd.sac.d.${Ttag}" "*semd.sac.d.${Ttag}.$maxm"`;
    }
  }

}

#=============================================
# NOW RUN MATLAB CODE misfit_gmt_run.m

#=============================================
# GENERATE GMT FIGURES

$isetup3 = 0;

if ($isetup3 == 1) {
  for ($ievent = $imin; $ievent <= $imax; $ievent++) {

    $eid = $eids[$ievent-1]; chomp($eid);
    print "\n $ievent : $eid";
    $dir2 = "$dir1/$eid";

    for ($irec = 1; $irec <= $nrec; $irec++) {
    #for ($irec = 104; $irec <= 104; $irec++) {

      cd_directory($pwd);
      ($sta,$net,$slat,$slon,undef,undef) = split(" ",$slines[$irec]);
      $station = "$sta.$net";
      $dir_rec = "${dir2}/${sta}.${net}";
      print "\n $irec out of $nrec -- ${dir_rec}";

      cd_directory($dir_rec);
      @datfiles = glob("*.dat");
      $ndat = @datfiles;

      if($ndat > 0) {
        `plot_misfit.pl $eid $sta/$net $Ts`;
      }
    }
  }

}

#=============================================
# CONCATENATE PDF FIGURES

$isetup4 = 0;

if ($isetup4 == 1) {

  for ($ievent = $imin; $ievent <= $imax; $ievent++) {

    $eid = $eids[$ievent-1]; chomp($eid);
    print "\n $ievent : $eid";
    $dir2 = "$dir1/$eid";

    for ($irec = 1; $irec <= $nrec; $irec++) {

      ($sta,$net,$slat,$slon,undef,undef) = split(" ",$slines[$irec]);
      $station = "$sta.$net";
      $dir_rec = "${dir2}/${sta}.${net}";
      print "\n $irec out of $nrec -- ${dir_rec}";

      @pdffiles = glob("${dir_rec}/*${Ttag}.pdf");
      $nfile = @pdffiles;
      #$pdffile = "${dir_rec}/${sta}_${net}_${Ttag}.pdf";
      if($nfile == 1) {
         $pdffile = $pdffiles[0]; chomp($pdffile);
         `cp $pdffile $dir2`;
      }
    }

    # concatenate into a single pdf
    $outfile = "$dir1/${eid}_${Ttag}_waveform_az.pdf";
    `rm $outfile`;
    @pdfall = glob("${dir2}/*${Ttag}.pdf");
    print "\n /home/carltape/bin/pdcat -r @pdfall $outfile\n";
    `/home/carltape/bin/pdcat -r @pdfall $outfile`;
    `rm @pdfall`;
  }
}

#=============================================
# DELETE FILES

$iclean = 0;

if ($iclean == 1) {

  for ($ievent = $imin; $ievent <= $imax; $ievent++) {

    $eid = $eids[$ievent-1]; chomp($eid);
    print "\n $ievent : $eid";
    $dir2 = "$dir1/$eid";

    for ($irec = 1; $irec <= $nrec; $irec++) {
      #for ($irec = 128; $irec <= 128; $irec++) {

      ($sta,$net,$slat,$slon,undef,undef) = split(" ",$slines[$irec]);
      $station = "$sta.$net";
      $dir_rec = "${dir2}/${sta}.${net}";
      print "\n $irec out of $nrec -- ${dir_rec}";

      `rm ${dir_rec}/*pdf`; `rm ${dir_rec}/*ps`;
    }
  }
}

#================================================
print "\n done with setup_misfit_dir.pl\n\n";
#================================================

sub cd_directory {
    my($new) = @_;
    my $old = `pwd`;
    chomp($old);
    check_directory($new);
    #print "$prog: cd $new\n";
    chdir $new;
    return($old);
}

sub check_directory {
    if(! -e $_[0] ) {
        print "Directory not found: $_[0]\n";
        exit(-1);
    }
}

#=================================================================
