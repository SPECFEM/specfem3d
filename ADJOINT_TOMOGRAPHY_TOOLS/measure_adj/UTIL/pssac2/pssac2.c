/*--------------------------------------------------------------------
 *    The GMT-system:	@(#)pssac.c	2.13  2/2/95
 *
 *    Copyright (c) 1991 by P. Wessel and W. H. F. Smith
 *    See README file for copying and redistribution conditions.
 *--------------------------------------------------------------------
 * pssac will read sacfiles and plot them
 * PostScript code is written to stdout.
 *
 * Revision:	Lupei Zhu from psxy
 * Modifications: Upgraded to new version of GMT (Brian Savage)
 *                Small updates
 *                Added ability to plot seismograms on maps
 *  2004 August 3   - Added Body wave scaling using JB Tables
 *  2006 December 5 - Will only use the absolute value of the distance/gcarc
 *                    Useful for negative distances
 */

#include <stdio.h>
#include <string.h>
#include <fcntl.h>
#include <math.h>

/* GMT and Postscript */
#include <gmt.h>
#include <pslib.h>

/* Local include files */
#include "sac.h"
#include "jb.h"

#define PI M_PI
#define LINE_SIZE 10512

#define GMT_VERSION_MAJOR  4
#define GMT_VERSION_MINOR  3

#define PSSAC "pssac2"

void GMT_history(int argc, char *argv[]);

static void
usage(void) {
  fprintf (stderr,"%s: Plot seismograms\n", PSSAC);
  fprintf(stderr,"Usage: %s <infiles> -C<cut0/cut1> -J<params> -R<w/e/s/n>\n", PSSAC);
  fprintf(stderr,"      [-B<tickinfo>] [-E(k|d|a|n|b)(t(0-9,-5(b),-3(o),-2(a))|vel)] [-K] [-V]\n");
  fprintf(stderr,"      [-Msize/alpha] [-Gr/g/b/c/t1/t2] [-O] [-P] [-Sshift] [-U[<label>]]\n");
  fprintf(stderr,"      [-W<pen>] [-#<ncopies>] [-X<x_shift>] [-Y<y_shift>] [-r]\n");
  fprintf(stderr,"\n");
  fprintf(stderr,"      sacfiles may be given on either the command line or\n");
  fprintf(stderr,"      if none are given they will be expected using standard in (stdin)\n");
  fprintf(stderr,"      Using stdin indicate (sacfile,x,y,pen)\n");
  fprintf(stderr,"          if x and y are given, it will override the position\n");
  fprintf(stderr,"               taken from the sacfile\n");
  fprintf(stderr,"          if pen is given the pen on the command line will be overridden\n");
  fprintf(stderr,"\n");
  fprintf(stderr,"      All options are in the same meaning as standard GMT EXCEPT:\n");
  fprintf(stderr,"      -C <cut_begin/cut_end>\n");
  fprintf(stderr,"              windowing data\n");
  fprintf(stderr,"              -Ct[0-9]/t[0-9] will use the t[0-9] header values\n");
  fprintf(stderr,"              or -Ct[0-9]/cut_end or -Ccut_begin/t[0-9]\n");
  fprintf(stderr,"      -E option determines\n");
  fprintf(stderr,"        (1) profile type:\n");
  fprintf(stderr,"	  a: azimuth profile\n");
  fprintf(stderr,"	  b: back-azimuth profile\n");
  fprintf(stderr,"	  d: distance profile with x in degrees\n");
  fprintf(stderr,"	  k: distance profile with x in kilometers\n");
  fprintf(stderr,"	  n: traces are numbered from 1 to N in y-axis\n");
  fprintf(stderr,"        (2) time alignment:\n");
  fprintf(stderr," 	  tn: align up with time mark tn in SAC head\n");
  fprintf(stderr,"	        n= -5(b), -3(o), -2(a), 0-9 (t0-t9)\n");
  fprintf(stderr,"	  vel: use reduced velocity\n");
  fprintf(stderr,"              Examples: -Ek8.0 -Eat-2 -Ed4.4 -Ent-3\n");
  fprintf(stderr,"      -M vertical scaling in sacfile_unit/MEASURE_UNIT = size<required> \n");
  fprintf(stderr,"          size: each trace will normalized to size (in MEASURE_UNIT)\n");
  fprintf(stderr,"              scale =  (1/size) * [data(max) - data(min)]\n");
  fprintf(stderr,"          size/alpha: plot absolute amplitude multiplied by (1/size)*r^alpha\n");
  fprintf(stderr,"              where r is the distance range in km across surface\n");
  fprintf(stderr,"              specifying alpha = 0.0 will give absolute amplitudes\n");
  fprintf(stderr,"              scale = (1/size) * r^alpha\n");
  fprintf(stderr,"          size/s: plot absolute amplitude multiplied by (1/size)*sqrt(sin(gcarc))\n");
  fprintf(stderr,"              where gcarc is the distance in degrees.\n");
  fprintf(stderr,"              scale = (1/size) * sqrt(sin(gcarc))\n");
  fprintf(stderr,"      -Gr/g/b/zero_line/t0/t1 positive phase painting\n");
  fprintf(stderr,"              paints the positive portion of the trace in color r/g/b\n");
  fprintf(stderr,"              from the zero_line to the top and between times t0 and t1\n");
  fprintf(stderr,"      -gr/g/b/zero_line/t0/t1 negative phase painting\n");
  fprintf(stderr,"              paints the negative portion of the trace in color r/g/b\n");
  fprintf(stderr,"              from the zero_line to the bottom and between times t0 and t1\n");
  fprintf(stderr,"      -N turn clipping off, data will be plotted outside the bounds\n");
  fprintf(stderr,"              clipping is on by default in X/Y plots and maps\n");
  fprintf(stderr,"      -S <seconds> shift the trace by seconds\n");
  fprintf(stderr,"              u: if next character after S is a u\n");
  fprintf(stderr,"                 then shift using the user variable\n");
  fprintf(stderr,"                 i.e. -Su9 will use 'user9' for shifting\n");
  fprintf(stderr,"                 positive shifts move the plot backwards in time\n");
  fprintf(stderr,"      -L <seconds per MEASURE_UNIT> while poltting on maps <required for maps>\n");
  fprintf(stderr,"              If your seismograms look choppy and pixelated\n");
  fprintf(stderr,"              Check the value of DOTS_PR_INCH in gmtdefaults\n");
  fprintf(stderr,"              and increase the value using gmtset\n");
  fprintf(stderr,"      -l <x/y/length/bar_length/font_size>\n");
  fprintf(stderr,"              plot a timescale in seconds per MEASURE_UNIT at (x,y)\n");
  fprintf(stderr,"              Only works while plotting on maps\n");
  fprintf(stderr,"      -D <dx/dy> shifts position of seismogram in MEASURE_UNIT\n");
  fprintf(stderr,"              Only works while plotting on maps\n");
  fprintf(stderr,"      -W <pen> the pen may be specified on the command line\n");
  fprintf(stderr,"              or it may be specified with individual files using stdin\n");
  fprintf(stderr,"      -n Distances and Great Circle Arcs may be negative, don't complain\n");
  fprintf(stderr,"      -r remove the mean value\n");
  fprintf(stderr,"      -s byte swap the data before plotting\n");  
  fprintf(stderr,"      -v plot the traces vertically\n");
  fprintf(stderr,"      -V Verbose (use more for more information)\n");
  exit(-1);
}

float 
get_sac_t_value(SACHEAD *h, int t) {
  float tval;
  float *fpt;
  fpt = (float *) h;
  tval = *(fpt + 10  + t);
  return(tval);
}

float 
get_sac_u_value(SACHEAD *h, int t) {
  float tval;
  float *fpt;
  fpt = (float *) h;
  tval = *(fpt + 40  + t);
  return(tval);
}

int
main (int argc, char **argv) {
  int 	i, ii, n, fno, nrdc;
  int	plot_n, pn, ierr;
  float	*yval;
  float *fpt;
  float t0;
  float y0, yy0;
  double *xx, *yy;
  double *xp ,*yp;
  double *xxx, *yyy;
  double *xxp, *yyp;
  double cam, cam0;
  double neg_cam0;
  double tt0, tt1;
  double neg_tt0, neg_tt1;
  char *c;
  char	line[LINE_SIZE];
  char sac_file[100];
  char penpen[20];
  char prf;
  char rdc;
  double cut0, cut1;
  int cut0t, cut1t;
  int user_variable;
  double t;
  double reduce_vel;
  double ysc;
  double alpha;
  double ty;

  /* For the actual line, used in GMT_plot_line */
  struct GMT_PEN pen;
  /* For phase painting, used in GMT_fill */
  struct GMT_FILL fill_pos, fill_neg;
  int *plot_pen;
  
  SACHEAD h;
  int rgb[3];
  float dt;
  float realt;
  float realdt;
  double secs_per_measure;

  /* TimeScale Variables */
  double timescale_x;
  double timescale_y;
  double timescale_len;
  double timescale_bar;
  int    timescale_font_size;
  char   timescale_txt[LINE_SIZE];

  double dx     = 0.0;
  double dy     = 0.0;
  double west   = 0.0;
  double east   = 0.0;
  double south  = 0.0;
  double north  = 0.0;
  double size   = 1.0;
  double delay  = 0.0;
  double yscale = 1.0;
  float  x0     = 0.0; 

  int red = 0;
  int blu = 0;
  int grn = 0;
  int n_files = 0;


  JBTable *jb;

  BOOLEAN negative_distances       = FALSE;
  BOOLEAN error                    = FALSE;
  BOOLEAN norm                     = FALSE;
  BOOLEAN reduce                   = FALSE;
  BOOLEAN window_cut               = FALSE;
  BOOLEAN phase_paint              = FALSE;
  BOOLEAN neg_phase_paint          = FALSE;
  BOOLEAN rmean                    = FALSE;
  BOOLEAN sin_scaling              = FALSE;
  BOOLEAN body_wave_scaling        = FALSE;
  BOOLEAN positioning              = FALSE;
  BOOLEAN vertical_trace           = FALSE;
  BOOLEAN plot_timescale           = FALSE;
  BOOLEAN option_M_specified       = FALSE;
  BOOLEAN option_L_specified       = FALSE;
  BOOLEAN clipping_on              = TRUE;
  BOOLEAN window_cut_use_headers_0 = FALSE;
  BOOLEAN window_cut_use_headers_1 = FALSE;
  BOOLEAN user_defined_shifts      = FALSE;
  int byte_swap                    = 0;
  int num_plotted                  = 0;
  int verbose_level                = 0;

  GMT_begin (argc,argv);
  GMT_init_pen (&pen,1);
  gmtdefs.verbose = FALSE;
  /* Check and interpret the command line arguments */
  for (i=1; !error && i < argc; i++) {
    if (argv[i][0] == '-') {
      switch(argv[i][1]) {
      /* Common parameters */
      case 'B': case 'J': case 'K': case 'O': case 'P':
      case 'R': case 'U': case 'X': case 'x': case 'Y':
      case 'y': case '#': case '\0':
	error += GMT_get_common_args (argv[i],&west,&east,&south,&north);
	break;
	
      /* Supplemental parameters,C,E,M,S,W,L */
      case 'M':		/* Multiple line segments,-Msize */
	if(strstr(&argv[i][2],"/s") == NULL &&
	   strstr(&argv[i][2],"/b") == NULL) { /* -Msize/alpha */
	  ierr=sscanf(&argv[i][2],"%lf/%lf",&size,&alpha);
	  if(ierr == 2) {
	    norm = TRUE;
	    size = 1/size;
	  } 
	} else if(strstr(&argv[i][2],"/b") != NULL) {  /* -Msize/b */
	  ierr=sscanf(&argv[i][2],"%lf",&size);
	  body_wave_scaling = TRUE;
	  size = 1/size;
	  jb = jbtable_load();
	} else if(strstr(&argv[i][2],"/s") != NULL) { /* -Msize/s */
	  ierr=sscanf(&argv[i][2],"%lf",&size);
	  sin_scaling = TRUE;
	  size = 1/size;
	} else {
	  fprintf(stderr, "%s: Error parsing -M option (%s)\n", PSSAC, argv[i]);
	  exit(-1);
	}
	option_M_specified = TRUE;
	break;
	break;
      case 'N':
	clipping_on = FALSE;
	break;
      case 'S':		/* Sdelay */
	if(argv[i][2] == 'u') { /* Align on User Variable */
	  sscanf(&argv[i][3],"%d", &user_variable);
	  user_defined_shifts = TRUE;
	} else {
	  sscanf(&argv[i][2],"%lf",&delay);
	}
	break;
      case 'G':		/* paint positive */
	sscanf(&argv[i][2],"%d/%d/%d/%lf/%lf/%lf",&red,&grn,&blu,&cam0,&tt0,&tt1);
	GMT_init_fill(&fill_pos, red, grn, blu);
	phase_paint = TRUE;
	break;
      case 'g':		/* paint positive */
	sscanf(&argv[i][2],"%d/%d/%d/%lf/%lf/%lf",&red,&grn,&blu,&neg_cam0,&neg_tt0,&neg_tt1);
	GMT_init_fill(&fill_neg, red, grn, blu);
	neg_phase_paint = TRUE;
	break;
      case 'C':		/* Ccut0/cut1 */
	c = &argv[i][2];
	if(c[0] == 't') {
	  window_cut_use_headers_0 = TRUE;
	  c++;
	}
	sscanf(c,"%lf",&cut0);
	c = strstr(&argv[i][2],"/");
	c++;
	if(c[0] == 't') {
	  window_cut_use_headers_1 = TRUE;
	  c++;
	}
	sscanf(c,"%lf",&cut1);
	window_cut = TRUE;
	if(window_cut_use_headers_0 == TRUE) {
	  cut0t = (int) cut0;
	  cut0 = 0.0;
	}
	if(window_cut_use_headers_1 == TRUE) {
	  cut1t = (int) cut1;
	  cut1 = 0.0;
	}
	break;
      case 'E':
	reduce=TRUE;
	prf=argv[i][2];
	rdc=argv[i][3];
	if (rdc == 't') {
	  if(sscanf(&argv[i][4],"%d",&nrdc) != 1) {
	    fprintf(stderr, "%s: error reading in SAC time variable\n", PSSAC);
	    exit(-1);
	  }
	} else {
	  if(sscanf(&argv[i][3],"%lf",&reduce_vel) != 1) {
	    fprintf(stderr, "%s: error reading in reduction velocity\n", PSSAC);
	    exit(-1);
	  }
	}
	break;
      case 'W':		/* Set line attributes */
	GMT_getpen (&argv[i][2],&pen);
	break;
      case 'n':
	negative_distances = TRUE;
	break;
      case 'D':		/* Set the position shifting */
	if(sscanf(&argv[i][2],"%lf/%lf",&dx, &dy) < 2) {
	  fprintf(stderr, "%s: error reading in position shifting variables\n", PSSAC);
	  exit(-1);
	}
	break;
      case 'l':		/* Set the timescale position */
	if(sscanf(&argv[i][2],"%lf/%lf/%lf/%lf/%d",
		  &timescale_x,
		  &timescale_y,
		  &timescale_len,
		  &timescale_bar,
		  &timescale_font_size) < 5) {
	  fprintf(stderr, "%s: error reading in timescale position\n", PSSAC);
	  exit(-1);
	}
	plot_timescale = TRUE;
	break;
      case 'L':		/* Set seconds per MEASURE_UNIT for map projection */
	option_L_specified = TRUE;
	sscanf(&argv[i][2],"%lf",&secs_per_measure);
	break;
      case 'r':
	rmean=TRUE;
	break;
      case 's':
	byte_swap = 1;
	break;
      case 'v':		/* plot trace vertically */
	vertical_trace = TRUE;
	break;
      case 'V':
	verbose_level = verbose_level + 1;
	break;
      default:		/* Options not recognized */
	error=TRUE;
	break;
      }
    }
    else n_files++;
  }
  if (argc == 1 || error) {	/* Display usage */
    usage();
  }
  if(option_M_specified == FALSE) {
    fprintf(stderr, "%s: GMT SYNTAX ERROR: Must specify -M\n", PSSAC);
    exit(-1);
  }
  if(option_L_specified == FALSE && project_info.projection >= 10) {
    fprintf(stderr, "%s: GMT SYNTAX ERROR: Must specify -L for plotting maps\n", PSSAC);
    exit(-1);
  }

#if ( GMT_VERSION_MAJOR < 4 )
  GMT_put_history(argc,argv); /* VERISON less and 4.0.0 */
#else
#if ( GMT_VERSION_MINOR < 2 ) 
  GMT_put_history(argc,argv); /* VERSION 4.1.4 and Less */
#else
  GMT_history(argc,argv); /* VERSION 4.2.0 and Greater */
#endif /* GMT_VERSION_MINOR < 4 */
#endif /* GMT_VERSION_MAJOR < 4 */

  GMT_map_setup (west,east,south,north);

#if ( GMT_VERSION_MAJOR < 4 )
  ps_plotinit (NULL,
	       gmtdefs.overlay,
	       gmtdefs.page_orientation,
	       gmtdefs.x_origin,
	       gmtdefs.y_origin,
	       gmtdefs.global_x_scale,
	       gmtdefs.global_y_scale,
	       gmtdefs.n_copies,
	       gmtdefs.dpi,
	       gmtdefs.measure_unit,
	       gmtdefs.paper_width,
	       gmtdefs.page_rgb,
	       GMT_epsinfo (argv[0]));
#else
  GMT_plotinit(argc, argv);
#endif
  GMT_echo_command (argc,argv);
  rgb[0] = rgb[1] = rgb[2] = -1;
  if(clipping_on == TRUE) {
    GMT_map_clip_on (rgb,3);
  }
  /* set pen attribute */ 
  GMT_setpen(&pen);

  yy0=0;
  fno=0;
  if(verbose_level > 1) {
    fprintf(stderr, "%s: Ready to recieve input \n", PSSAC);
  }

  while( (n_files && ++fno<argc) || (n_files==0 && fgets(line, LINE_SIZE, stdin))) {
    /* get sac file name */
    if (n_files) {	/* from command line */
      if(verbose_level > 1) {
	fprintf(stderr, "%s: Using commnd line \n", PSSAC);
      }
      if (argv[fno][0] == '-') continue;
      else strcpy(sac_file,argv[fno]);
    } else {		/* from standard in */
      if(verbose_level >= 1) {
	fprintf(stderr, "%s: Using stdin \n", PSSAC);
      }
      i=sscanf(line, "%s %f %f %s",sac_file,&x0,&y0,penpen);
      if (i>2) positioning = TRUE;
      if (i>3) {
	/* set pen attribute */
	GMT_getpen(penpen,&pen);
	GMT_setpen(&pen);
      }
      if(verbose_level > 2) {
	fprintf(stderr, "%s: Reading %s \n", PSSAC, sac_file);
      }
    }
    
    /* read in sac files */
    if(verbose_level > 1) {
      fprintf(stderr, "%s: Reading in %s ...\n", PSSAC, sac_file);
    }
    if( (yval = read_sac_swap(sac_file, &h, byte_swap)) == NULL ) {
      fprintf(stderr, "%s: Error reading in %s\n", PSSAC, sac_file);
      yy0++;
      continue;
    }
    fpt = (float *) &h;
    if(window_cut_use_headers_0 == TRUE) {
      cut0 = get_sac_t_value(&h, (int)cut0t);
    }
    if(window_cut_use_headers_1 == TRUE) {    
      cut1 = get_sac_t_value(&h, (int)cut1t);
    }
    
    /* determine time alignment and y-axis 
       location of the trace from profile type*/
    t=t0=h.b+delay;
    if(user_defined_shifts == TRUE) {
      if(get_sac_u_value(&h, (int)user_variable) != -12345.0) {
	t = t - get_sac_u_value(&h, (int)user_variable);
	t0 = t;
      }
    }
    if(verbose_level > 2) {
      fprintf(stderr, "%s: Before Time Retardation\n", PSSAC);
      fprintf(stderr, "%s: \tTime:    %f\n", PSSAC, t);
      fprintf(stderr, "%s: \tT0:      %f\n", PSSAC, t0);
      fprintf(stderr, "%s: \tb-value: %f\n", PSSAC, h.b);
      fprintf(stderr, "%s: \tdelay:   %f\n", PSSAC, delay);
      if(user_defined_shifts == TRUE) {
	fprintf(stderr, "%s: \tdelay u: %f\n", PSSAC, get_sac_u_value(&h, (int)user_variable));
      }
    }
    if (reduce) {
       switch (prf) {
       case 'd':
          yy0=h.gcarc;
	  if(verbose_level > 1) {
	    fprintf(stderr, "%s:\tGreat Circle Arc: %e\n", PSSAC, yy0);
	  }
	  if(h.gcarc < 0.0) {
	    fprintf(stderr, "%s: Great Circle Arc less than zero (%s) %e\n", PSSAC, sac_file, yy0);
	  }
          break;
       case 'k':
	  yy0=h.dist;
	  if(verbose_level > 1) {
	    fprintf(stderr, "%s:\tDistance: %e\n", PSSAC, yy0);
	  }
	  if(h.dist < 0.0 && !negative_distances) {
	    fprintf(stderr, "%s: Distance less than zero (%s) %e\n", PSSAC, sac_file, yy0);
	  }
          break;
       case 'a':
          yy0=h.az;
	  if(verbose_level > 1) {
	    fprintf(stderr, "%s:\tAzimuth: %e\n", PSSAC, yy0);
	  }
	  if(h.az < 0.0) {
	    fprintf(stderr, "%s: Azimuth less than zero (%s) %e\n", PSSAC, sac_file, yy0);
	  }
	  break;
       case 'b':
	  yy0=h.baz;
	  if(verbose_level > 1) {
	    fprintf(stderr, "%s:\tBack-Azimuth: %e\n", PSSAC, yy0);
	  }
	  if(h.baz < 0.0) {
	    fprintf(stderr, "%s: Back-Azimuth less than zero (%s) %e\n", PSSAC, sac_file, yy0);
	  }
	  break;
       case 'n':
	  if(verbose_level > 1) {
	    fprintf(stderr, "%s:\tNumber: %e\n", PSSAC, yy0);
	  }
          yy0++;
          break;
       default:
          fprintf(stderr,"%s: wrong choise of profile type (d|k|a|n)\n", PSSAC);
          exit(3);
       }
       if (rdc == 't') {
	 t = t0 - *(fpt + 10 + nrdc);
	 if(window_cut_use_headers_0 == TRUE) {
	   cut0 = cut0 - *(fpt + 10 + nrdc);
	 }
	 if(window_cut_use_headers_1 == TRUE) {
	   cut1 = cut1 - *(fpt + 10 + nrdc);
	 }
       }
       else {
	 t = t0 - fabs(h.dist/reduce_vel);
	 if(window_cut_use_headers_0 == TRUE) {
	   cut0 = cut0 - fabs(h.dist/reduce_vel);
	 }
	 if(window_cut_use_headers_1 == TRUE) {
	   cut1 = cut1 - fabs(h.dist/reduce_vel);
	 }
       }
    }
    if (! positioning ) y0 = yy0;

    xx=( double *) GMT_memory(CNULL,h.npts,sizeof(double),"pssac2");
    yy=( double *) GMT_memory(CNULL,h.npts,sizeof(double),"pssac2");
    if(verbose_level > 2) {
      fprintf(stderr, "%s: After Time Retardation\n", PSSAC);
      fprintf(stderr, "%s: \tt:     %f\n", PSSAC, t);
      fprintf(stderr, "%s: \tx0:    %f\n", PSSAC, x0);
    }
    /* Set the dt depending on the projections type */
    realdt = h.delta;
    realt  = t;

    /* At this point we need to determine 
       dt -- x spacing     in degrees (if projecting) or MEASURE_UNIT (if plotting)
       t  -- initial point in degrees (if projecting) or MEASURE_UNIT (if plotting)
       which depends on projection 
    */
    if(project_info.projection < 10) { 
      /* Linear, Log10, Pow */
      /* Tested with -JX  but nothing else */
      dt = h.delta;
      t += x0;
    } else { /* Other projections -- see gmt_project.h and MAP_PROJECTIONS */
      /* If the location was not specified, 
	 use the sac header stlo and stla 
      */
      if(positioning == FALSE) {
	x0 = h.stlo;
	y0 = h.stla;
      }
      /* Determine dt (in MEASURE_UNIT) */
      dt = (h.delta) / secs_per_measure;
      /* Convert original lat/lon point (x0,y0) to 
	 MEASURE_UNIT (t,ty=>y0) for plotting */
      if(vertical_trace == TRUE) {
	GMT_geo_to_xy(x0,y0,&ty,&t);
	y0 = t;
	t  += dy;
	y0 += dx;
      } else {
	GMT_geo_to_xy(x0,y0,&t,&ty);
	y0 = ty;
	t  += dx;
	y0 += dy;
      }
      /* Get position from t(position) and
	 realt (timing) and
	 secs_per_measure (time scaling) */
      t = t + (realt / secs_per_measure);
      if(verbose_level > 2) {
	fprintf(stderr, "%s: Paramters for map plotting \n", PSSAC);
	fprintf(stderr, "%s: \tsecs_per_measure: %f\n", PSSAC, secs_per_measure);
	fprintf(stderr, "%s: \tmeasure_unit: %s\n", PSSAC,
		GMT_unit_names[gmtdefs.measure_unit]);
      }
    }
    if(verbose_level > 2) {
      fprintf(stderr, "%s: After dt Determination\n", PSSAC);
      fprintf(stderr, "%s: \tt:          %f\n", PSSAC, t);
      fprintf(stderr, "%s: \tdt:         %f\n", PSSAC, dt);
      fprintf(stderr, "%s: \tprojection: %d\n", PSSAC, project_info.projection);
      fprintf(stderr, "%s: \tx0:         %f\n", PSSAC, x0);
      fprintf(stderr, "%s: \trealt:      %f\n", PSSAC, realt);
      fprintf(stderr, "%s: \trealdt:     %f\n", PSSAC, realdt);
      fprintf(stderr, "%s: \tcut0:       %f\n", PSSAC, cut0);
      fprintf(stderr, "%s: \tcut1:       %f\n", PSSAC, cut1);
      
    }
    for(n=0,i=0;i<h.npts;i++) {
      if( (window_cut && realt>=cut0 && realt<=cut1) || ( !window_cut) ){
	xx[n]=t;
	yy[n]=yval[i];
	n++;
      }
      t += dt;
      realt += realdt;
    }
    if(n==0) {
       free((char *)xx);
       free((char *)yy);
       if(verbose_level > 1) {
	 fprintf(stderr, "%s: Skipping %s (no points)\n", PSSAC, sac_file);
       }
       continue;
    }

    if(option_M_specified == TRUE){
      if(verbose_level > 1) {
	if(sin_scaling){
	  fprintf(stderr, "%s: Multi-Segment Scaled, trace multiplied by %lf *sqrt(sin(gcarc))\n",PSSAC, size);
	}
	else if(norm){
	  fprintf(stderr, "%s: Multi-Segment Scaled, trace multiplied by %lf *dist[km]^%lf\n",PSSAC,size,alpha);
	}
	else if(body_wave_scaling) {
	  fprintf(stderr, "%s: Multi-Segment Scaled, trace multiplied by %lf Geometrical Spreading Factor\n", PSSAC, size);
	}
	else{
	  fprintf(stderr, "%s: Multi-Segment Scaled, trace normalized to %lf inches\n", PSSAC, size);
	}
      }
      h.depmax=-1.e20;h.depmin=1.e20;h.depmen=0.;
      for(i=0;i<n;i++){
         h.depmax=h.depmax > yy[i] ? h.depmax : yy[i];
	 h.depmin=h.depmin < yy[i] ? h.depmin : yy[i];
	 h.depmen += yy[i];
      }
      h.depmen /= ( float ) n;
      if(!rmean) h.depmen=0.0;

      if(project_info.projection < 10) { /* X-Y-ish plots */
	ysc=size*fabs((north-south)/project_info.pars[1]);
	ysc=size*fabs((north-south)/GMT_map_height);
      } else { /* Map Projections */
	ysc = size;
      }
      if(norm) {
	if(h.dist < 0.0 && !negative_distances) {
	  fprintf(stderr, "%s: Warning -- Distance in km is less than zero (%s)\n", PSSAC, sac_file);
	}
	yscale=pow(fabs(h.dist),alpha)*ysc;
      }
      else if(sin_scaling) {
	if(h.gcarc < 0.0 && !negative_distances) {
	  fprintf(stderr, "%s: Warning -- Distance in degrees is less than zero (%s)\n", PSSAC, sac_file);
	}
	yscale=ysc*pow(sin(fabs(h.gcarc)*2*PI/360),0.5);
      }
      else if(body_wave_scaling) {
	if(fabs(h.gcarc) < 30.0 || fabs(h.gcarc) > 90.0) {
	  fprintf(stderr, "%s: Warning -- Distances < 30.0 or > 90.0 degrees are ill-behaved\n", PSSAC);
	}
	if(h.evdp < 0.0 || h.evdp > 800) {
	  fprintf(stderr, "%s: File: %s Depth: %f\n", PSSAC, sac_file, h.evdp);
	  fprintf(stderr, "%s: Warning -- Depths < 0.0 or > 800.0 km are ill-behaved\n", PSSAC);	  
	}
	yscale = ysc * gfact3(jb, fabs(h.gcarc), h.evdp);
      }
      else {
	yscale=ysc/fabs(h.depmax - h.depmin);
      }

      for(i=0;i<n;i++) yy[i]=(double) (yy[i]-h.depmen)*yscale + y0;
      if(verbose_level > 2) {
	fprintf(stderr, "%s: Vertical Scale Determination\n", PSSAC);
	fprintf(stderr, "%s: \tdepmin:  %f\n", PSSAC, h.depmin);
	fprintf(stderr, "%s: \tdepmax:  %f\n", PSSAC, h.depmax);
	fprintf(stderr, "%s: \tysc:     %f\n", PSSAC, ysc);
	fprintf(stderr, "%s: \tpars[0]: %f\n", PSSAC, project_info.pars[0]);
	fprintf(stderr, "%s: \tpars[1]: %f\n", PSSAC, project_info.pars[1]);
	fprintf(stderr, "%s: \tnorth:   %f\n", PSSAC, north);
	fprintf(stderr, "%s: \tsouth:   %f\n", PSSAC, south);
	fprintf(stderr, "%s: \twidth:   %f\n", PSSAC, GMT_map_width);
	fprintf(stderr, "%s: \theight:  %f\n", PSSAC, GMT_map_height);
	fprintf(stderr, "%s: \tyscale:  %f\n", PSSAC, yscale);
	fprintf(stderr, "%s: \ty0:      %f\n", PSSAC, y0);
      }
    }

    if (vertical_trace) {
       if(verbose_level > 1) {
	 fprintf(stderr, "%s: Vertical Trace\n", PSSAC);
       }
       xp = (double *) GMT_memory(CNULL, n, sizeof (double), PSSAC);
       memcpy((char *)xp, (char *)yy, n*sizeof(double));
       memcpy((char *)yy, (char *)xx, n*sizeof(double));
       memcpy((char *)xx, (char *)xp, n*sizeof(double));
       free(xp);
    }
    
    if(verbose_level > 3) {
      fprintf(stderr, "%s: projecting %d points\n", PSSAC, n);
      for(i = 0; i < n; i++) {
	fprintf(stderr, "%s: (pt) => %d: (x,y) => (%f, %f)\n", PSSAC, i,xx[i],yy[i]);
      }
    }
    if(project_info.projection < 10) {
      plot_n = GMT_geo_to_xy_line (xx,yy,n);
      if(verbose_level > 2) {
	fprintf(stderr, "%s: creating %d points\n", PSSAC, plot_n);
      }
      if(plot_n == 0) {
	fprintf(stderr, "%s: No points to plot -- %s\n", PSSAC, sac_file);
      }
      xp=(double *) GMT_memory (CNULL,plot_n,sizeof (double),PSSAC);
      yp=(double *) GMT_memory (CNULL,plot_n,sizeof (double),PSSAC);
      plot_pen = (int *) GMT_memory(CNULL, plot_n, sizeof(int), PSSAC);
      memcpy ((char *)xp,(char *)GMT_x_plot,plot_n * sizeof (double));
      memcpy ((char *)yp,(char *)GMT_y_plot,plot_n * sizeof (double));
      memcpy ((char *)plot_pen, (char *) GMT_pen, plot_n * sizeof(int));
    } else {
      plot_n = n;
      xp=(double *) GMT_memory (CNULL,plot_n,sizeof (double),PSSAC);
      yp=(double *) GMT_memory (CNULL,plot_n,sizeof (double),PSSAC);
      plot_pen = (int *) GMT_memory(CNULL, plot_n, sizeof(int), PSSAC);
      memcpy ((char *)xp,(char *)xx,plot_n * sizeof (double));
      memcpy ((char *)yp,(char *)yy,plot_n * sizeof (double));
      for(i = 0; i < plot_n; i++) {
	plot_pen[i] = 2;
      }
    }
 
    if (phase_paint) {
       if(verbose_level > 1) {
	 fprintf(stderr, "%s: Painting Phase\n", PSSAC);
       }
       xxx=(double *) GMT_memory(CNULL,n+2,sizeof (double),PSSAC);
       yyy=(double *) GMT_memory(CNULL,n+2,sizeof (double),PSSAC);
          cam = y0+cam0*yscale;
          for(i=1;i<n && xx[i]<tt1;i++) {
             ii = 0;
             if (yy[i]>cam && xx[i]>tt0) {
   	        yyy[ii] = cam;
   	        if (yy[i-1] < cam)
		  xxx[ii] = xx[i-1]+(cam-yy[i-1])*(xx[i]-xx[i-1])/(yy[i]-yy[i-1]);
		if (xxx[ii]<tt0) {
		  xxx[ii] = tt0;
		  yyy[ii] = yy[i-1]+(tt0-xx[i-1])*(yy[i]-yy[i-1])/(xx[i]-xx[i-1]);
		}
   	        else {
		  xxx[ii] = tt0;
		  if (xx[i-1]>tt0)
		    xxx[ii] = xx[i-1];
		  ii++;
		  xxx[ii] = xxx[ii-1];
		  yyy[ii] = yy[i-1]+(xxx[ii-1]-xx[i-1])*(yy[i]-yy[i-1])/(xx[i]-xx[i-1]);
   	        }
   	        ii++;
   	        while(yy[i]>cam && xx[i]<tt1 && i<n-1) {
		  xxx[ii] = xx[i];
		  yyy[ii] = yy[i];
		  ii++;
		  i++;
   	        }
   	        if (yy[i]>cam) {
		  xxx[ii] = tt1;
		  if (xx[i]<tt1)
		    xxx[ii] = xx[i];
		  yyy[ii] = yy[i-1]+(xxx[ii]-xx[i-1])*(yy[i]-yy[i-1])/(xx[i]-xx[i-1]);
		  ii++;
		  xxx[ii] = xxx[ii-1];
   	        } else {
		  xxx[ii] = xx[i-1] + (cam-yy[i-1])*(xx[i]-xx[i-1])/(yy[i]-yy[i-1]);
		}
		yyy[ii] = cam;
   	        ii++;
                if ((pn = GMT_geo_to_xy_line(xxx,yyy,ii)) < 3 ) continue;
                xxp=(double *) GMT_memory (CNULL,pn,sizeof (double),PSSAC);
                yyp=(double *) GMT_memory (CNULL,pn,sizeof (double),PSSAC);
                memcpy ((char *)xxp,(char *)GMT_x_plot,pn*sizeof (double));
                memcpy ((char *)yyp,(char *)GMT_y_plot,pn*sizeof (double));
		GMT_fill(xxp, yyp, pn, &fill_pos, FALSE);
		free((char *)xxp);
		free((char *)yyp);
             }
          }
       free((char *)xxx);
       free((char *)yyy);
    }
    if (neg_phase_paint) {
       if(verbose_level > 1) {
	 fprintf(stderr, "%s: Painting Phase Negative\n", PSSAC);
       }
       xxx=(double *) GMT_memory(CNULL,n+2,sizeof (double),PSSAC);
       yyy=(double *) GMT_memory(CNULL,n+2,sizeof (double),PSSAC);
       cam = y0+neg_cam0*yscale;
       /* Loop over all points until neg_tt1 */
       for(i=1;i<n && xx[i]<neg_tt1;i++) {
	 ii = 0;
	 /* Find a point greater than the initial time and greater that the upper limit */
	 if (yy[i] < cam && xx[i] > neg_tt0) {
	   yyy[ii] = cam;
	   /* If the previous y point is greater than the upper limit, linear interpolate */
	   if (yy[i-1] > cam)
	     xxx[ii] = xx[i-1]+(cam-yy[i-1])*(xx[i]-xx[i-1])/(yy[i]-yy[i-1]);
	   /* If the previous x point is less than neg_tt0, find the assicaiated y value */
	   if (xxx[ii] < neg_tt0) {
	     xxx[ii] = neg_tt0;
	     yyy[ii] = yy[i-1]+(neg_tt0-xx[i-1])*(yy[i]-yy[i-1])/(xx[i]-xx[i-1]);
	   }
	   else {
	     xxx[ii] = neg_tt0;
	     if (xx[i-1] > neg_tt0)
	       xxx[ii] = xx[i-1];
	     ii++;
	     xxx[ii] = xxx[ii-1];
	     yyy[ii] = yy[i-1]+(xxx[ii-1]-xx[i-1])*(yy[i]-yy[i-1])/(xx[i]-xx[i-1]);
	   }
	   ii++;
	   /* Continue until we find another y point above the upper limit
	      or Reach neg_tt1 in x or we reach the end of i	   */
	   while(yy[i] < cam && xx[i] < neg_tt1 && i < n-1) {
	     xxx[ii] = xx[i];
	     yyy[ii] = yy[i];
	     ii++;
	     i++;
	   }
	   
	   if (yy[i] < cam) {
	     xxx[ii] = neg_tt1;
	     if (xx[i] < neg_tt1)
	       xxx[ii] = xx[i];
	     yyy[ii] = yy[i-1]+(xxx[ii]-xx[i-1])*(yy[i]-yy[i-1])/(xx[i]-xx[i-1]);
	     ii++;
	     xxx[ii] = xxx[ii-1];
	   } else {
	     xxx[ii] = xx[i-1] + (cam-yy[i-1])*(xx[i]-xx[i-1])/(yy[i]-yy[i-1]);
	   }
	   /* Close the polygon */
	   yyy[ii] = cam;
	   ii++;
	   if ((pn = GMT_geo_to_xy_line(xxx,yyy,ii)) < 3 ) continue;
	   xxp=(double *) GMT_memory (CNULL,pn,sizeof (double),PSSAC);
	   yyp=(double *) GMT_memory (CNULL,pn,sizeof (double),PSSAC);
	   memcpy ((char *)xxp,(char *)GMT_x_plot,pn*sizeof (double));
	   memcpy ((char *)yyp,(char *)GMT_y_plot,pn*sizeof (double));
	   GMT_fill(xxp, yyp, pn, &fill_neg, FALSE);
	   free((char *)xxp);
	   free((char *)yyp);
	 } /* Finished painting this section of negative values */
       } /* Finished with the entire trace */
       free((char *)xxx);
       free((char *)yyy);
    }
    
    if(verbose_level > 2) {
      fprintf(stderr, "%s: plotting %d points\n", PSSAC, plot_n);
    }
    if(verbose_level > 3) {
      fprintf(stderr, "%s: projecting %d points\n", PSSAC, n);
      for(i = 0; i < n; i++) {
	fprintf(stderr, "%s: (pt) => %d: (x,y) => (%f, %f)\n", PSSAC, i,xp[i],yp[i]);
      }
    }
    GMT_plot_line(xp, yp, plot_pen, plot_n);
    num_plotted = num_plotted + 1;

    free((char *)xp);
    free((char *)yp);
    free((char *)yval);
    free((char *)xx);
    free((char *)yy);
  }
  if(plot_timescale) {
    ps_plot(timescale_x - timescale_len/2, timescale_y - timescale_bar/2, 3);
    ps_plot(timescale_x - timescale_len/2, timescale_y + timescale_bar/2, 2);
    ps_plot(timescale_x - timescale_len/2, timescale_y, 2);
    ps_plot(timescale_x + timescale_len/2, timescale_y, 2);
    ps_plot(timescale_x + timescale_len/2, timescale_y - timescale_bar/2, 2);
    ps_plot(timescale_x + timescale_len/2, timescale_y + timescale_bar/2, 2);
    sprintf(timescale_txt, "%.0f secs.", secs_per_measure * timescale_len);
#if ( GMT_VERSION_MAJOR < 4 )
    GMT_text3d(timescale_x, timescale_y - timescale_bar/2, project_info.z_level,
	       timescale_font_size, gmtdefs.anot_font, timescale_txt,
	       0.0, 10, 0);
#else
    GMT_text3D(timescale_x, timescale_y - timescale_bar/2, project_info.z_level,
	       timescale_font_size, gmtdefs.annot_font[0], timescale_txt,
	       0.0, 10, 0);
#endif
  }
  if(clipping_on == TRUE) {
    GMT_map_clip_off ();
  }
  
  GMT_init_pen(&pen,1);
  GMT_setpen(&pen);
  if (frame_info.plot) {
    GMT_map_basemap();
  }
  if(verbose_level > 0) {
    fprintf(stderr, "%s: plotted %d seismograms\n", PSSAC, num_plotted);
  }
#if (GMT_VERSION_MAJOR < 4)
  ps_plotend (gmtdefs.last_page);
  GMT_end (argc,argv);
#else 
  GMT_plotend();
#endif
  return(1);
}

