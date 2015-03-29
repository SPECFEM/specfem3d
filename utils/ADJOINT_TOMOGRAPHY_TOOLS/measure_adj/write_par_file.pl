#!/usr/bin/perl -w
#
#==========================================================
#  write_par_file.pl
#
#  This script writes a parameters file to be used in the measurement code.
#  See the reading-in of this file in measure_adj.f90.
#
#  INPUT:
#     tstart/dt/ntime
#     imeas
#     channel (BH or LH)
#     TSHORT/TLONG
#     RUN_BANDPASS/DISPLAY_DETAILS/OUTPUT_MEASUREMENT_FILES/COMPUTE_ADJOINT_SOURCE
#     par1: TSHIFT_MIN/TSHIFT_MAX/DLNA_MIN/DLNA_MAX/CC_MIN
#     par2: ERROR_TYPE/DT_SIGMA_MIN/DLNA_SIGMA_MIN
#     par3: ITAPER/WTR/NPI/DT_FAC/ERR_FAC/DT_MAX_SCALE/NCYCLE_IN_WINDOW
#
#  EXAMPLE:
#    write_par_file.pl -0.585/0.011/18200 1 BH 6/30 0/1/1/1 -4.5/4.5/-1.5/1.5/0.69 1/1.0/0.5 1/0.02/2.5/2.0/2.5/3.5/1.5
#
#
#==========================================================

if (@ARGV < 8) {die("Usage: write_par_file.pl tstart/dt/npts imeas chan Tshort/Tlong
          RUN_BANDPASS/DISPLAY_DETAILS/OUTPUT_MEASUREMENT_FILES/COMPUTE_ADJONIT_SOURCE
          tshift_min/tshift_max/dlnA_min/dlnA_max/cc_min  error_type/dt_sigma_min/dlnA_sigma_min
          itaper/wtr/npi/dt_fac/err_fac/dt_max_scale/ncycle_in_window
             **** writes MEASUREMENT.PAR with given parameters *** \n\n");}
($tvec,$imeas,$chan,$Ts,$iparbools,$par1,$par2,$par3) = @ARGV;

# extract variables
($tstart,$dt,$ntime) = split("/",$tvec);
($TSHORT,$TLONG) = split("/",$Ts);
($ibool1,$ibool2,$ibool3,$ibool4) = split("/",$iparbools);
if($ibool1==1) {$RUN_BANDPASS = ".true."} else {$RUN_BANDPASS = ".false."}
if($ibool2==1) {$DISPLAY_DETAILS = ".true."} else {$DISPLAY_DETAILS = ".false."}
if($ibool3==1) {$OUTPUT_MEASUREMENT_FILES = ".true."} else {$OUTPUT_MEASUREMENT_FILES = ".false."}
if($ibool4==1) {$COMPUTE_ADJOINT_SOURCE = ".true."} else {$COMPUTE_ADJOINT_SOURCE = ".false."}
($TSHIFT_MIN,$TSHIFT_MAX,$DLNA_MIN,$DLNA_MAX,$CC_MIN) = split("/",$par1);
($ERROR_TYPE,$DT_SIGMA_MIN,$DLNA_SIGMA_MIN) = split("/",$par2);
($ITAPER,$WTR,$NPI,$DT_FAC,$ERR_FAC,$DT_MAX_SCALE,$NCYCLE_IN_WINDOW) = split("/",$par3);

# comments for the PAR_FILE
#$line00 = "output directory";  # cannot have a string here
$line01 = "tstart, DT, npts: time vector for simulations";
$line02 = "imeas (1-8; see manual)";
$line03 = "channel: BH or LH";
$line04 = "TLONG and TSHORT: band-pass periods for records";
$line05 = "RUN_BANDPASS: use band-pass on records";
$line06 = "DISPLAY_DETAILS";
$line07 = "OUTPUT_MEASUREMENT_FILES";
$line08 = "COMPUTE_ADJOINT_SOURCE";

$line09 = "TSHIFT_MIN; TSHIFT_MAX";
$line10 = "DLNA_MIN; DLNA_MAX";
$line11 = "CC_MIN";
$line12 = "ERROR_TYPE -- 0 none; 1 CC, MT-CC; 2 MT-jack-knife";
$line13 = "DT_SIGMA_MIN";
$line14 = "DLNA_SIGMA_MIN";

$line15 = "ITAPER -- taper type: 1 multi-taper; 2 cosine; 3 boxcar";
$line16 = "WTR, NPI (ntaper = 2*NPI)";
$line17 = "DT_FAC";
$line18 = "ERR_FAC";
$line19 = "DT_MAX_SCALE";
$line20 = "NCYCLE_IN_WINDOW";

#=============================================
$ofile = "./MEASUREMENT.PAR";
`\\rm $ofile`;   # remove file if it exists

print "\nWriting to $ofile ...\n";
open(OUT,">$ofile");

print OUT sprintf("%8.3f%7.4f%8i  # %s\n",$tstart,$dt,$ntime,$line01);
print OUT sprintf("%23i  # %s\n",$imeas,$line02);
print OUT sprintf("%23s  # %s\n",$chan,$line03);
print OUT sprintf("%12.3f%11.3f  # %s\n",$TLONG,$TSHORT,$line04);
print OUT sprintf("%23s  # %s\n",$RUN_BANDPASS,$line05);
print OUT sprintf("%23s  # %s\n",$DISPLAY_DETAILS,$line06);
print OUT sprintf("%23s  # %s\n",$OUTPUT_MEASUREMENT_FILES,$line07);
print OUT sprintf("%23s  # %s\n",$COMPUTE_ADJOINT_SOURCE,$line08);

print OUT sprintf("%12.4f%11.4f  # %s\n",$TSHIFT_MIN,$TSHIFT_MAX,$line09);
print OUT sprintf("%12.4f%11.4f  # %s\n",$DLNA_MIN,$DLNA_MAX,$line10);
print OUT sprintf("%23.3f  # %s\n",$CC_MIN,$line11);
print OUT sprintf("%23i  # %s\n",$ERROR_TYPE,$line12);
print OUT sprintf("%23.3f  # %s\n",$DT_SIGMA_MIN,$line13);
print OUT sprintf("%23.3f  # %s\n",$DLNA_SIGMA_MIN,$line14);

print OUT sprintf("%23i  # %s\n",$ITAPER,$line15);
print OUT sprintf("%17.3f%6.2f  # %s\n",$WTR,$NPI,$line16);
print OUT sprintf("%23.3f  # %s\n",$DT_FAC,$line17);
print OUT sprintf("%23.3f  # %s\n",$ERR_FAC,$line18);
print OUT sprintf("%23.3f  # %s\n",$DT_MAX_SCALE,$line19);
print OUT sprintf("%23.3f  # %s\n",$NCYCLE_IN_WINDOW,$line20);

#print OUT sprintf("%s\n",$odir);
#print OUT sprintf("%23i  # %s\n",$ITAPER,$line01);
#print OUT sprintf("%17.3f%6.2f  # %s\n",$WTR,$NPI,$line02);
#print OUT sprintf("%23i  # %s\n",$imeas,$line03);
#print OUT sprintf("%23s  # %s\n",$RUN_BANDPASS,$line04);
#print OUT sprintf("%12.3f%11.3f  # %s\n",$TLONG,$TSHORT,$line05);
#print OUT sprintf("%8.3f%7.4f%8i  # %s\n",$tstart,$dt,$ntime,$line06);
#print OUT sprintf("%23s  # %s\n",$DISPLAY_DETAILS,$line07);
#print OUT sprintf("%23s  # %s\n",$OUTPUT_MEASUREMENT_FILES,$line08);
#print OUT sprintf("%23i  # %s\n",$INCLUDE_ERROR,$line09);
#print OUT sprintf("%23.3f  # %s\n",$DT_FAC,$line10);
#print OUT sprintf("%23.3f  # %s\n",$ERR_FAC,$line11);
#print OUT sprintf("%23.3f  # %s\n",$DT_MAX_SCALE,$line12);
#print OUT sprintf("%23.3f  # %s\n",$NCYCLE_IN_WINDOW,$line13);
#print OUT sprintf("%12.4f%11.4f  # %s\n",$BEFORE_QUALITY,$AFTER_QUALITY,$line14);
#print OUT sprintf("%12.4f%11.4f  # %s\n",$BEFORE_TSHIFT,$AFTER_TSHIFT,$line15);
#print OUT sprintf("%12.4f%11.4f  # %s\n",$DT_SIGMA_MIN,$DLNA_SIGMA_MIN,$line16);
close(OUT);

#=================================================================
