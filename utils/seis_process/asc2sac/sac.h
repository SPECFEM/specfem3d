
#ifndef __SAC_H__
#define __SAC_H__
/*===========================================================================*/
/* SEED reader     |               sac.h                   |     header file */
/*===========================================================================*/
/*
  Name:   sac.h
  Purpose:  structure for header of a SAC (Seismic Analysis Code) data file
  Usage:    #include "sac.h"
  Input:    not applicable
  Output:   not applicable
  Warnings: not applicable
  Errors:   not applicable
  Called by:  output_data
  Calls to: none
  Algorithm:  not applicable
  Notes:    Key to comment flags describing each field:
        Column 1:
          R = required by SAC
            = (blank) optional
        Column 2:
          A = settable from a priori knowledge
          D = available in data
          F = available in or derivable from SEED fixed data header
          T = available in SEED header tables
            = (blank) not directly available from SEED data, header
              tables, or elsewhere
  Problems: none known
  References: O'Neill, D. (1987).  IRIS Interim Data Distribution Format
    (SAC ASCII), Version 1.0 (12 November 1987).  Incorporated
    Research Institutions for Seismology, 1616 North Fort Myer
    Drive, Suite 1440, Arlington, Virginia 22209.  11 pp.
    Tull, J. (1987).  SAC User's Manual, Version 10.2, October 7,
      1987.  Lawrence Livermore National Laboratory, L-205,
      Livermore, California 94550.  ??? pp.
  Language: C, hopefully ANSI standard
  Author: Dennis O'Neill
  Revisions:07/15/88  Dennis O'Neill  Initial preliminary release 0.9
    11/21/88  Dennis O'Neill  Production release 1.0
    09/19/89  Dennis O'Neill  corrected length of char strings
*/
/** Information for a sac file except data.
    The order of the variables are in the same order as
    would be for a sac file.
    R = required by SAC<p>
    A = settable from a priori knowledge<p>
    D = available in data<p>
    F = available in or derivable from SEED fixed data header<p>
    T = available in SEED header tables<p>
      = (blank) not directly available from SEED data, header tables, or elsewhere<p>

*/

struct sac
{
  /** RF time increment, sec    .
      @param Required
      @param SEED Available in or derivable from SEED fixed data header
  */
  float delta;
  /**    minimum amplitude      .    */  float  depmin;
  /**    maximum amplitude      .    */  float  depmax;
  /**    amplitude scale factor .    */  float  scale;
  /**    observed time inc      .    */  float  odelta;
  /** RD initial value, ampl.   .
      @param Required
      @param Data available in data
  */
  float b;
  /** RD final value, amplitude .
      @param Required
      @param Data available in data
  */
  float e;
  /**    event start, sec > 0   .    */  float  o;
  /**    1st arrival time       .    */  float  a;
  /**    internal use           .    */  float  internal1;
  /**    user-defined time pick .    */  float  t0;
  /**    user-defined time pick .    */  float  t1;
  /**    user-defined time pick .    */  float  t2;
  /**    user-defined time pick .    */  float  t3;
  /**    user-defined time pick .    */  float  t4;
  /**    user-defined time pick .    */  float  t5;
  /**    user-defined time pick .    */  float  t6;
  /**    user-defined time pick .    */  float  t7;
  /**    user-defined time pick .    */  float  t8;
  /**    user-defined time pick .    */  float  t9;
  /**    event end, sec > 0     .    */  float  f;
  /**    instrument respnse parm.    */  float  resp0;
  /**    instrument respnse parm.    */  float  resp1;
  /**    instrument respnse parm.    */  float  resp2;
  /**    instrument respnse parm.    */  float  resp3;
  /**    instrument respnse parm.    */  float  resp4;
  /**    instrument respnse parm.    */  float  resp5;
  /**    instrument respnse parm.    */  float  resp6;
  /**    instrument respnse parm.    */  float  resp7;
  /**    instrument respnse parm.    */  float  resp8;
  /**    instrument respnse parm.    */  float  resp9;
  /**  T station latititude     .
       @param Table available in SEED header tables
   */
  float stla;
  /**  T station longitude      .
       @param Table available in SEED header tables
  */
  float stlo;
  /**  T station elevation, m   .
       @param Table available in SEED header tables
  */
  float stel;
  /**  T station depth, m       .
       @param Table available in SEED header tables
  */
  float stdp;
  /**    event latitude         .    */  float  evla;
  /**    event longitude        .    */  float  evlo;
  /**    event elevation        .    */  float  evel;
  /**    event depth            .    */  float  evdp;
  /**    reserved for future use.    */  float  unused1;
  /**    available to user      .    */  float  user0;
  /**    available to user      .    */  float  user1;
  /**    available to user      .    */  float  user2;
  /**    available to user      .    */  float  user3;
  /**    available to user      .    */  float  user4;
  /**    available to user      .    */  float  user5;
  /**    available to user      .    */  float  user6;
  /**    available to user      .    */  float  user7;
  /**    available to user      .    */  float  user8;
  /**    available to user      .    */  float  user9;
  /**    stn-event distance, km .    */  float  dist;
  /**    event-stn azimuth      .    */  float  az;
  /**    stn-event azimuth      .    */  float  baz;
  /**    stn-event dist, degrees.    */  float  gcarc;
  /**    internal use           .    */  float  internal2;
  /**    internal use           .    */  float  internal3;
  /**    mean value, amplitude  .    */  float  depmen;
  /**  T component azimuth      .
       @param Table available in SEED header tables
  */
  float cmpaz;
  /**  T component inclination  .
       @param Table available in SEED header tables
  */
  float cmpinc;
  /**    reserved for future use.    */  float  unused2;
  /**    reserved for future use.    */  float  unused3;
  /**    reserved for future use.    */  float  unused4;
  /**    reserved for future use.    */  float  unused5;
  /**    reserved for future use.    */  float  unused6;
  /**    reserved for future use.    */  float  unused7;
  /**    reserved for future use.    */  float  unused8;
  /**    reserved for future use.    */  float  unused9;
  /**    reserved for future use.    */  float  unused10;
  /**    reserved for future use.    */  float  unused11;
  /**    reserved for future use.    */  float  unused12;
  /**  F zero time of file, yr  .
       @param SEED Available in or derivable from SEED fixed data header
  */
  long  nzyear;
  /**  F zero time of file, day .
       @param SEED Available in or derivable from SEED fixed data header
  */
  long  nzjday;
  /**  F zero time of file, hr  .
       @param SEED Available in or derivable from SEED fixed data header
  */
  long  nzhour;
  /**  F zero time of file, min .
       @param SEED Available in or derivable from SEED fixed data header
  */
  long  nzmin;
  /**  F zero time of file, sec .
       @param SEED Available in or derivable from SEED fixed data header
  */
  long  nzsec;
  /**  F zero time of file, msec.
       @param SEED Available in or derivable from SEED fixed data header
  */
  long  nzmsec;
  /**    internal use           .    */  long internal4;
  /**    internal use           .    */  long internal5;
  /**    internal use           .    */  long internal6;
  /** RF number of samples      .
      @param Required
      @param SEED Available in or derivable from SEED fixed data header
  */
  long  npts;
  /**    internal use           .    */  long internal7;
  /**    internal use           .    */  long internal8;
  /**    reserved for future use.    */  long unused13;
  /**    reserved for future use.    */  long unused14;
  /**    reserved for future use.    */  long unused15;
  /** RA type of file           .
      @param Required
      @param Prior settable from a priori knowledge
   */
  long  iftype;
  /**    type of amplitude      .    */  long idep;
  /**    zero time equivalence  .    */  long iztype;
  /**    reserved for future use.    */  long unused16;
  /**    recording instrument   .    */  long iinst;
  /**    stn geographic region  .    */  long istreg;
  /**    event geographic region.    */  long ievreg;
  /**    event type             .    */  long ievtyp;
  /**    quality of data        .    */  long iqual;
  /**    synthetic data flag    .    */  long isynth;
  /**    reserved for future use.    */  long unused17;
  /**    reserved for future use.    */  long unused18;
  /**    reserved for future use.    */  long unused19;
  /**    reserved for future use.    */  long unused20;
  /**    reserved for future use.    */  long unused21;
  /**    reserved for future use.    */  long unused22;
  /**    reserved for future use.    */  long unused23;
  /**    reserved for future use.    */  long unused24;
  /**    reserved for future use.    */  long unused25;
  /**    reserved for future use.    */  long unused26;
  /** RA data-evenly-spaced flag.
      @param Required
      @param Prior settable from a priori knowledge
  */
  long  leven;
  /**    station polarity flag  .    */  long lpspol;
  /**    overwrite permission   .    */  long lovrok;
  /**    calc distance, azimuth .    */  long lcalda;
  /**    reserved for future use.    */  long unused27;
  /**  F station name           .
       @param SEED Available in or derivable from SEED fixed data header
   */
  char  kstnm[8];
  /**    event name             .    */  char kevnm[16];
  /**    man-made event name    .    */  char khole[8];
  /**    event origin time id   .    */  char ko[8];
  /**    1st arrival time ident .    */  char ka[8];
  /**    time pick 0 ident      .    */  char kt0[8];
  /**    time pick 1 ident      .    */  char kt1[8];
  /**    time pick 2 ident      .    */  char kt2[8];
  /**    time pick 3 ident      .    */  char kt3[8];
  /**    time pick 4 ident      .    */  char kt4[8];
  /**    time pick 5 ident      .    */  char kt5[8];
  /**    time pick 6 ident      .    */  char kt6[8];
  /**    time pick 7 ident      .    */  char kt7[8];
  /**    time pick 8 ident      .    */  char kt8[8];
  /**    time pick 9 ident      .    */  char kt9[8];
  /**    end of event ident     .    */  char kf[8];
  /**    available to user      .    */  char kuser0[8];
  /**    available to user      .    */  char kuser1[8];
  /**    available to user      .    */  char kuser2[8];
  /**  F component name         .
       @param SEED Available in or derivable from SEED fixed data header
  */
  char  kcmpnm[8];
  /**    network name           .    */  char knetwk[8];
  /**    date data read         .    */  char kdatrd[8];
  /**    instrument name        .    */  char kinst[8];
};

#ifdef SAC_NULL
/** a SAC structure containing all null values.
    @see sac
*/
static struct sac sac_null = {
  -12345., -12345., -12345., -12345., -12345.,
  -12345., -12345., -12345., -12345., -12345.,
  -12345., -12345., -12345., -12345., -12345.,
  -12345., -12345., -12345., -12345., -12345.,
  -12345., -12345., -12345., -12345., -12345.,
  -12345., -12345., -12345., -12345., -12345.,
  -12345., -12345., -12345., -12345., -12345.,
  -12345., -12345., -12345., -12345., -12345.,
  -12345., -12345., -12345., -12345., -12345.,
  -12345., -12345., -12345., -12345., -12345.,
  -12345., -12345., -12345., -12345., -12345.,
  -12345., -12345., -12345., -12345., -12345.,
  -12345., -12345., -12345., -12345., -12345.,
  -12345., -12345., -12345., -12345., -12345.,
  -12345, -12345, -12345, -12345, -12345,
  -12345, -12345, -12345, -12345, -12345,
  -12345, -12345, -12345, -12345, -12345,
  -12345, -12345, -12345, -12345, -12345,
  -12345, -12345, -12345, -12345, -12345,
  -12345, -12345, -12345, -12345, -12345,
  -12345, -12345, -12345, -12345, -12345,
  -12345, -12345, -12345, -12345, -12345,
  { '-','1','2','3','4','5',' ',' ' },
  { '-','1','2','3','4','5',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ' },
  { '-','1','2','3','4','5',' ',' ' }, { '-','1','2','3','4','5',' ',' ' },
  { '-','1','2','3','4','5',' ',' ' }, { '-','1','2','3','4','5',' ',' ' },
  { '-','1','2','3','4','5',' ',' ' }, { '-','1','2','3','4','5',' ',' ' },
  { '-','1','2','3','4','5',' ',' ' }, { '-','1','2','3','4','5',' ',' ' },
  { '-','1','2','3','4','5',' ',' ' }, { '-','1','2','3','4','5',' ',' ' },
  { '-','1','2','3','4','5',' ',' ' }, { '-','1','2','3','4','5',' ',' ' },
  { '-','1','2','3','4','5',' ',' ' }, { '-','1','2','3','4','5',' ',' ' },
  { '-','1','2','3','4','5',' ',' ' }, { '-','1','2','3','4','5',' ',' ' },
  { '-','1','2','3','4','5',' ',' ' }, { '-','1','2','3','4','5',' ',' ' },
  { '-','1','2','3','4','5',' ',' ' }, { '-','1','2','3','4','5',' ',' ' },
  { '-','1','2','3','4','5',' ',' ' }
};
#endif

#define FLOAT_NULL   -12345.0
#define LONG_NULL    -12345
#define CHAR_NULL_8  "-12345  "
#define CHAR_NULL_16 "-12345          "

/* definitions of constants for SAC enumerated data values */
/* undocumented == couldn't find a definition for it (denio, 07/15/88) */
#define IREAL   0 /* undocumented              */
#define ITIME   1 /* file: time series data    */
#define IRLIM   2 /* file: real&imag spectrum  */
#define IAMPH   3 /* file: ampl&phas spectrum  */
#define IXY     4 /* file: gen'l x vs y data   */
#define IUNKN   5 /* x data: unknown type      */
/* zero time: unknown        */
/* event type: unknown       */
#define IDISP   6 /* x data: displacement (nm) */
#define IVEL    7 /* x data: velocity (nm/sec) */
#define IACC    8 /* x data: accel (cm/sec/sec)*/
#define IB      9 /* zero time: start of file  */
#define IDAY   10 /* zero time: 0000 of GMT day*/
#define IO     11 /* zero time: event origin   */
#define IA     12 /* zero time: 1st arrival    */
#define IT0    13 /* zero time: user timepick 0*/
#define IT1    14 /* zero time: user timepick 1*/
#define IT2    15 /* zero time: user timepick 2*/
#define IT3    16 /* zero time: user timepick 3*/
#define IT4    17 /* zero time: user timepick 4*/
#define IT5    18 /* zero time: user timepick 5*/
#define IT6    19 /* zero time: user timepick 6*/
#define IT7    20 /* zero time: user timepick 7*/
#define IT8    21 /* zero time: user timepick 8*/
#define IT9    22 /* zero time: user timepick 9*/
#define IRADNV 23 /* undocumented              */
#define ITANNV 24 /* undocumented              */
#define IRADEV 25 /* undocumented              */
#define ITANEV 26 /* undocumented              */
#define INORTH 27 /* undocumented              */
#define IEAST  28 /* undocumented              */
#define IHORZA 29 /* undocumented              */
#define IDOWN  30 /* undocumented              */
#define IUP    31 /* undocumented              */
#define ILLLBB 32 /* undocumented              */
#define IWWSN1 33 /* undocumented              */
#define IWWSN2 34 /* undocumented              */
#define IHGLP  35 /* undocumented              */
#define ISRO   36 /* undocumented              */
#define INUCL  37 /* event type: nuclear shot  */
#define IPREN  38 /* event type: nuke pre-shot */
#define IPOSTN 39 /* event type: nuke post-shot*/
#define IQUAKE 40 /* event type: earthquake    */
#define IPREQ  41 /* event type: foreshock     */
#define IPOSTQ 42 /* event type: aftershock    */
#define ICHEM  43 /* event type: chemical expl */
#define IOTHER 44 /* event type: other source  */ /* data quality: other problm*/
#define IGOOD  45 /* data quality: good        */
#define IGLCH  46 /* data quality: has glitches*/
#define IDROP  47 /* data quality: has dropouts*/
#define ILOWSN 48 /* data quality: low s/n     */
#define IRLDTA 49 /* data is real data         */
#define IVOLTS 50 /* file: velocity (volts)    */
#define INIV51 51 /* undocumented              */
#define INIV52 52 /* undocumented              */
#define INIV53 53 /* undocumented              */
#define INIV54 54 /* undocumented              */
#define INIV55 55 /* undocumented              */
#define INIV56 56 /* undocumented              */
#define INIV57 57 /* undocumented              */
#define INIV58 58 /* undocumented              */
#define INIV59 59 /* undocumented              */
#define INIV60 60 /* undocumented              */

/* True/false definitions */
#ifndef TRUE
#define FALSE 0
#define TRUE !FALSE
#endif

/* Format strings for writing headers for SAC ASCII files */
#define FCS "%15.7f%15.7f%15.7f%15.7f%15.7f\n"  /* for floats */
#define ICS "%10d%10d%10d%10d%10d\n"    /* for integers */
#define CCS1 "%-8.8s%-8.8s%-8.8s\n"   /* for strings */
#define CCS2 "%-8.8s%-16.16s\n"     /* for strings */

/* Construct for writing data to SAC ASCII files */
/*
 * for (i = 0; i < data_hdr->nsamples; i++)
 * {
 * fprintf (outfile, "%15.7f", *seismic_data_ptr++);
 * if ( ( ((i+1)%5) == 0) && (i > 0) ) fprintf (outfile, "\n");
 * }
 *
 */

#include <stdio.h>

FILE *sacopen(char *infile, char * mode);
int  sacclose(FILE *fp);
struct sac * sachdrread(FILE *fp);
float *sacdatamalloc(int npts);
void sacfree(float *data_ptr);
int sacdataread(FILE *fp, float * data_ptr, int npts, char *infile);
int sacdatawrite(FILE *fpout, float *data_ptr, int npts);
int sachdrwrite(FILE *fpout, struct sac *sac_hdr);
struct sac * sachdrmalloc();

#endif /* __SAC_H__ */



