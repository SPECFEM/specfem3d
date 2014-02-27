#define TEST_WGRAPH

#include <sys/time.h>

#include "module.h"
#include "common.h"
#include "parser.h"
#include "graph.h"
#include "wgraph.h"
#include "wgraph_part_st.h"

/*********************/
/*                   */
/* The main routine. */
/*                   */
/*********************/

int
main (

int                 argc,
char *              argv[])
{
  FILE *              fileptr;
  FILE *              stream;
  Wgraph              wgrfdat;
  Strat *             straptr;
  double              clktime[2];
  double              runtime;
  struct timespec     tp;
  Gnum                partnbr;
  Gnum                partval;
  Gnum                vertnum;
  Gnum                comploadsum;
  Gnum                comploadavg;
  Gnum                comploadmax;
  Gnum                comploadmin;
  int                 i;
  char                c;
  const Gnum * restrict vlbltax;

  errorProg (argv[0]);

  if (argc != 5) {
    errorPrint ("main: usage is \"%s <nparts> <input source graph file> <output mapping file> <strategy>\"\n", argv[0]);
    exit       (1);
  }

  fprintf (stderr, "pid %d\n", getpid ());
  printf ("Waiting for key press...\n");
  scanf ("%c", &c);

  if (graphInit (&wgrfdat.s) != 0) {   /* Initialize source graph */
    errorPrint ("main: cannot initialize graph");
    return     (1);
  }

  if ((partnbr = atoi (argv[1])) < 1) {
    errorPrint ("main: invalid number of parts (\"%s\")", argv[1]);
    return     (1);
  }

  fileptr = NULL;

  if ((fileptr = fopen (argv[2], "r")) == NULL) {
    errorPrint ("main: cannot open graph file");
    return     (1);
  }

  if (graphLoad (&wgrfdat.s, fileptr, -1, 0) != 0) {
    errorPrint ("main: cannot load graph file");
    return     (1);
  }

  if (fileptr != NULL)
    fclose (fileptr);

  wgraphInit (&wgrfdat, &wgrfdat.s, partnbr);
/* wgraphZero (&wgrfdat); TODO to clean? */

  if ((straptr = stratInit (&wgraphpartststratab, argv[4])) == NULL) {
    errorPrint ("main: cannot initialize strategy");
    return     (1);
  }

  /* Init clock. */
  clktime[0] = 
  clktime[1] = 
  runtime    = 0;
 
  /* Start clock. */
  clock_gettime (CLOCK_REALTIME, &tp);            /* Elapsed time */
  clktime[0] = ((double) tp.tv_sec + (double) tp.tv_nsec * 1.0e-9L);
  wgraphPartSt (&wgrfdat, straptr);

  /* Stop clock. */
  clock_gettime (CLOCK_REALTIME, &tp);            /* Elapsed time */
  clktime[1] = ((double) tp.tv_sec + (double) tp.tv_nsec * 1.0e-9L);
  runtime += (clktime[1] - clktime[0]);

  comploadsum =
  comploadavg =
  comploadmax = 0;
  comploadmin = GNUMMAX;
  for (i = 0; i < partnbr; i ++) {
    comploadsum += wgrfdat.compload[i];
    if (wgrfdat.compload[i] > comploadmax)
      comploadmax = wgrfdat.compload[i];
    if (wgrfdat.compload[i] < comploadmin)
      comploadmin = wgrfdat.compload[i];
  }
  comploadavg = comploadsum / wgrfdat.partnbr;
  for (partval = 0; partval < wgrfdat.partnbr; partval ++) /* for each part */
    printf("\033[0;33mcompload[%d] %d %d\033[0m\n", partval, wgrfdat.compload[partval], wgrfdat.compsize[partval]);
  printf ("\033[0;32mCompLoad(max/avg)=%g,\t(min/avg)=%g\n",
	  ((double) comploadmax / (double) comploadavg),
	  ((double) comploadmin / (double) comploadavg));
  printf ("\033[0;32mFronLoad=%ld\n",
	  (long) wgrfdat.fronload);
  printf ("\033[0;32mTime=%f\n", runtime);
  stream = NULL;

  if ((stream = fopen (argv[3], "w")) == NULL) {
    errorPrint ("main: cannot open graph file");
    return     (1);
  }
  vlbltax = (wgrfdat.s.vlbltax != NULL) ? (wgrfdat.s.vlbltax) : NULL;

  if (fprintf (stream, "%ld\n", (long) wgrfdat.s.vertnbr) == EOF) {
    errorPrint ("main: bad output (1)");
    return     (1);
  }

  for (vertnum = wgrfdat.s.baseval; vertnum < wgrfdat.s.vertnnd; vertnum ++) {
    if (fprintf (stream, "%ld\t%ld\n",
                 (long) ((vlbltax != NULL) ? vlbltax[vertnum] : vertnum),
                 (long) wgrfdat.parttax[vertnum]) == EOF) {
      errorPrint ("main: bad output (2)");
      return     (1);
    }
  }

  if (stream != NULL)
    fclose (stream);

  wgrfdat.parttax = NULL;                      /* Do not free parttax as grouped with frontab */
  wgraphExit (&wgrfdat);
  stratExit  (straptr);
  exit       (0);
}
