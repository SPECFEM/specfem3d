#!/usr/bin/perl -w

$ivert = 0;
$irun = 1;
$stip = sprintf("%2.2i",$irun);

#$ftag = "vert_${stip}_xc";
#$ftag = "horz_xc_vb_m16_m00";
$ftag = "horz_xc_vs_m01_m00";

@epsfiles = glob("${ftag}*.ps");

 $neps = @epsfiles;
 for ($k = 1; $k <= $neps; $k = $k+1) {
   $efile = $epsfiles[$k-1];
   `ps2eps -f $efile`;
 }

if ($irun==1 && $ivert==1) {
  `ps2ps vert_01_xc_HEC_CI_282_9968977_vs_m16_m00_003.ps temp.ps ; ps2eps temp.ps ; mv temp.eps vert_01_xc_HEC_CI_282_9968977_vs_m16_m00_003.eps`;
  `ps2ps vert_01_xc_LCG_CI_122_10215753_vs_m16_m00_021.ps temp.ps ; ps2eps temp.ps ; mv temp.eps vert_01_xc_LCG_CI_122_10215753_vs_m16_m00_021.eps`;
  `ps2ps vert_01_xc_LGU_CI_353_10097009_vs_m16_m00_067.ps temp.ps ; ps2eps temp.ps ; mv temp.eps vert_01_xc_LGU_CI_353_10097009_vs_m16_m00_067.eps`;
  `ps2ps vert_01_xc_BZN_AZ_306_14138080_vs_m16_m00_013.ps temp.ps ; ps2eps temp.ps ; mv temp.eps vert_01_xc_BZN_AZ_306_14138080_vs_m16_m00_013.eps`;
  `ps2ps vert_01_xc_EDW2_CI_134_14236768_vs_m16_m00_011.ps temp.ps ; ps2eps temp.ps ; mv temp.eps vert_01_xc_EDW2_CI_134_14236768_vs_m16_m00_011.eps`;
  `ps2ps vert_01_xc_ISA_CI_160_14383980_vs_m16_m00_049.ps temp.ps ; ps2eps temp.ps ; mv temp.eps vert_01_xc_ISA_CI_160_14383980_vs_m16_m00_049.eps`;
  `ps2ps vert_01_xc_HEC_CI_296_14418600_vs_m16_m00_007.ps temp.ps ; ps2eps temp.ps ; mv temp.eps vert_01_xc_HEC_CI_296_14418600_vs_m16_m00_007.eps`;
  `ps2ps vert_01_xc_ISA_CI_113_14418600_vs_m16_m00_044.ps temp.ps ; ps2eps temp.ps ; mv temp.eps vert_01_xc_ISA_CI_113_14418600_vs_m16_m00_044.eps`;
  `ps2ps vert_01_xc_DPP_CI_314_9753485_vs_m16_m00_079.ps temp.ps ; ps2eps temp.ps ; mv temp.eps vert_01_xc_DPP_CI_314_9753485_vs_m16_m00_079.eps`;
}

if ($irun==1 && $ivert==0) {
  `ps2ps ${ftag}_001.ps temp.ps ; ps2eps temp.ps ; mv temp.eps ${ftag}_001.eps`;
  `ps2ps ${ftag}_009.ps temp.ps ; ps2eps temp.ps ; mv temp.eps ${ftag}_009.eps`;
  `ps2ps ${ftag}_016.ps temp.ps ; ps2eps temp.ps ; mv temp.eps ${ftag}_016.eps`;
}

#==================================================
