package DELAZ5;
use strict;
use warnings;
use Math::Trig;
BEGIN {
    use Exporter();
    our($VERSION, @ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
    $VERSION     = 1.00;
    @ISA         = qw(Exporter);
    @EXPORT      = qw(&delaz5);
    %EXPORT_TAGS = ( );
    @EXPORT_OK   = ( );
}
our @EXPORT_OK;

# culled from Hiroos subroutines
sub delaz5 {
    my($thei, $alei, $thsi, $alsi, $i) = @_;
#	This is an old Fortran subroutine to compute the
#       distance and azimuths between two points.
#	Point A:  Latitude   THEI,  Longitude  ALEI
#	Point B:  Latitude   THSI,  Longitude  ALSI
# 	If the coordinates are in geographical coordinates in deg,  then I=0
#        If the coordinates are in geocentric in radia, then I=1
#	These are the input parameters.
#	Outputs are:
#		DELT=distance in radian
#  		DELTDG=distance in degree
#		DELTKM=distance in km
#               AZES=azimuth of Point B as viewed from Point A (radian)
#		AZESDEG=azimuth of Point B as viewed from Point A (radian) (degree)          
#		AZSE=azimuth of Point A as viewed from Point B (radian)
#		AZSEDEG=azimuth of Point A as viewed from Point B (radian) (degree)
    my($h, $delt, $deltkm, $azes, $azse);
    my($azsedg, $deltdg,$azesdg);
    my($the, $ale, $ths, $als, $aaa);
    my($c, $ak, $d, $e, $g);
    my($cp, $akp, $dp ,$ep, $ap, $bp, $gp, $hp);
    my($c1, $c2, $c3, $c4, $c5, $c6);
    
#   if(i) 50, 50, 51
    if($i <= 0.0) {
#   if  coordinates are geograph deg i=0
#   if coordinates are geocent radian  i=1
	$the=1.745329252e-2*$thei; # Line 50
	$ale=1.745329252e-2*$alei;
	$ths=1.745329252e-2*$thsi;
	$als=1.745329252e-2*$alsi;
	$aaa=0.9931177*tan($the);
	$the=atan($aaa);
	$aaa=0.9931177*tan($ths);
	$ths=atan($aaa);
#      go to 32
    } else {
	$the=$thei; # Line 51
	$ale=$alei;
	$ths=$thsi;
	$als=$alsi;
    }
#    32 continue
    $c= sin($the);
    $ak=-cos($the);
    $d=sin($ale);
    $e= -cos($ale);
    $a= $ak*$e;
    $b= -$ak*$d;
    $g=-$c*$e;
    $h=$c*$d;
    $cp=sin($ths);
    $akp=-cos($ths);
    $dp=sin($als);
    $ep = -cos($als);
    $ap = $akp*$ep;
    $bp=-$akp*$dp;
    $gp=-$cp*$ep;
    $hp=$cp*$dp;
    $c1=$a*$ap+$b*$bp+$c*$cp;
    
#      if( c1-0.94 )  30, 31, 31
    if($c1-0.94 < 0.0) {
#    30	if(c1+0.94) 28, 28, 29 
	if($c1+0.94 <= 0.0) {
	    # this used to be lower in the code delaz5.f
	    $c1=($a+$ap)**2+($b+$bp)**2+($c+$cp)**2; # Line 28
	    $c1 = sqrt($c1 );
	    $c1= $c1/2.0;
	    $delt = acos($c1);
	    $delt = 2.0*$delt;
#	    go to 33
	} else {
	    $delt=acos($c1); #Line 29
	}
    } else {
	# this used to be lower in the code delaz5.f
	$c1=($a-$ap)**2+($b-$bp)**2+($c-$cp)**2; # Line 31
	$c1= sqrt($c1);
	$c1=$c1/2.0;
	$delt = asin($c1);
	$delt= 2.0*$delt;
	#go to 33
    }
    $deltkm=6371.0*$delt; # Line 33
    $c3 = ($ap-$d)**2+($bp-$e)**2+$cp**2-2.0;
    $c4 = ($ap-$g)**2+($bp-$h)**2+($cp-$ak)**2-2.0;
    $c5 = ($a-$dp)**2+ ($b-$ep)**2+$c**2-2.0;
    $c6 = ($a-$gp)**2+($b-$hp)**2+($c-$akp)**2-2.0;
    $deltdg = 57.29577951*$delt;
    $azes = atan2($c3, $c4 );
#    if ( $azes ) 80, 81, 81
    if($azes < 0.0) {
	$azes = 6.283185308+ $azes; # Line 80
    }
    $azse = atan2( $c5, $c6 ); # Line 81

#    if ( azse ) 70, 71 , 71
    if($azse < 0.0) {
	$azse=6.283185308+$azse; #line 70
    }
    $azesdg=57.29577951*$azes; # Line 71
    $azsedg=57.29577951*$azse;
    return($delt, $deltdg, $deltkm, $azes, $azesdg, $azse, $azsedg);
}

1;
