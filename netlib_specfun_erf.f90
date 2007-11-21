
  subroutine calerf(ARG,RESULT,JINT)

!------------------------------------------------------------------
!
! This routine can be freely obtained from Netlib
! at http://www.netlib.org/specfun/erf
!
! Most Netlib software packages have no restrictions on their use
! but Netlib recommends that you check with the authors to be sure.
! See http://www.netlib.org/misc/faq.html#2.3 for details.
!
!------------------------------------------------------------------
!
!   This packet evaluates erf(x) for a real argument x.
!   It contains one FUNCTION type subprogram: ERF,
!   and one SUBROUTINE type subprogram, CALERF.  The calling
!   statements for the primary entries are:
!
!                   Y = ERF(X)
!
!   The routine  CALERF  is intended for internal packet use only,
!   all computations within the packet being concentrated in this
!   routine.  The function subprograms invoke  CALERF  with the
!   statement
!
!          call CALERF(ARG,RESULT,JINT)
!
!   where the parameter usage is as follows
!
!      Function                     Parameters for CALERF
!       call              ARG                  Result          JINT
!
!     ERF(ARG)      ANY REAL ARGUMENT         ERF(ARG)          0
!
!   The main computation evaluates near-minimax approximations
!   from "Rational Chebyshev approximations for the error function"
!   by William J. Cody, Math. Comp., 1969, PP. 631-638.  This
!   transportable program uses rational functions that theoretically
!   approximate  erf(x)  and  erfc(x)  to at least 18 significant
!   decimal digits.  The accuracy achieved depends on the arithmetic
!   system, the compiler, the intrinsic functions, and proper
!   selection of the machine-dependent constants.
!
!*******************************************************************
!*******************************************************************
!
! Explanation of machine-dependent constants
!
!   XMIN   = the smallest positive floating-point number.
!   XINF   = the largest positive finite floating-point number.
!   XNEG   = the largest negative argument acceptable to ERFCX;
!            the negative of the solution to the equation
!            2*exp(x*x) = XINF.
!   XSMALL = argument below which erf(x) may be represented by
!            2*x/sqrt(pi)  and above which  x*x  will not underflow.
!            A conservative value is the largest machine number X
!            such that   1.0 + X = 1.0   to machine precision.
!   XBIG   = largest argument acceptable to ERFC;  solution to
!            the equation:  W(x) * (1-0.5/x**2) = XMIN,  where
!            W(x) = exp(-x*x)/[x*sqrt(pi)].
!   XHUGE  = argument above which  1.0 - 1/(2*x*x) = 1.0  to
!            machine precision.  A conservative value is
!            1/[2*sqrt(XSMALL)]
!   XMAX   = largest acceptable argument to ERFCX; the minimum
!            of XINF and 1/[sqrt(pi)*XMIN].
!
!   Approximate IEEE double precision values are defined below.
!
!*******************************************************************
!*******************************************************************
!
! Error returns
!
!  The program returns  ERFC = 0      for  ARG >= XBIG;
!
!  Author: William J. Cody
!          Mathematics and Computer Science Division
!          Argonne National Laboratory
!          Argonne, IL 60439, USA
!
!  Latest modification: March 19, 1990
!
!  Converted to Fortran90 and slightly modified by
!  Dimitri Komatitsch, University of Pau, France, November 2007.
!
!------------------------------------------------------------------

  implicit none

  integer I,JINT
  double precision A,ARG,B,C,D,DEL,FOUR,HALF,P,ONE,Q,RESULT,SIXTEEN,SQRPI, &
       TWO,THRESHOLD,X,XBIG,XDEN,XHUGE,XINF,XMAX,XNEG,XNUM,XSMALL, &
       Y,YSQ,ZERO
  dimension A(5),B(4),C(9),D(8),P(6),Q(5)

!------------------------------------------------------------------
!  Mathematical constants
!------------------------------------------------------------------
  data FOUR,ONE,HALF,TWO,ZERO/4.0D0,1.0D0,0.5D0,2.0D0,0.0D0/, &
       SQRPI/5.6418958354775628695D-1/,THRESHOLD/0.46875D0/, &
       SIXTEEN/16.0D0/

!------------------------------------------------------------------
!  Machine-dependent constants
!------------------------------------------------------------------
  data XINF,XNEG,XSMALL/1.79D308,-26.628D0,1.11D-16/, &
       XBIG,XHUGE,XMAX/26.543D0,6.71D7,2.53D307/

!------------------------------------------------------------------
!  Coefficients for approximation to  erf  in first interval
!------------------------------------------------------------------
  data A/3.16112374387056560D00,1.13864154151050156D02, &
         3.77485237685302021D02,3.20937758913846947D03, &
         1.85777706184603153D-1/
  data B/2.36012909523441209D01,2.44024637934444173D02, &
         1.28261652607737228D03,2.84423683343917062D03/

!------------------------------------------------------------------
!  Coefficients for approximation to  erfc  in second interval
!------------------------------------------------------------------
  data C/5.64188496988670089D-1,8.88314979438837594D0, &
         6.61191906371416295D01,2.98635138197400131D02, &
         8.81952221241769090D02,1.71204761263407058D03, &
         2.05107837782607147D03,1.23033935479799725D03, &
         2.15311535474403846D-8/
  data D/1.57449261107098347D01,1.17693950891312499D02, &
         5.37181101862009858D02,1.62138957456669019D03, &
         3.29079923573345963D03,4.36261909014324716D03, &
         3.43936767414372164D03,1.23033935480374942D03/

!------------------------------------------------------------------
!  Coefficients for approximation to  erfc  in third interval
!------------------------------------------------------------------
  data P/3.05326634961232344D-1,3.60344899949804439D-1, &
         1.25781726111229246D-1,1.60837851487422766D-2, &
         6.58749161529837803D-4,1.63153871373020978D-2/
  data Q/2.56852019228982242D00,1.87295284992346047D00, &
         5.27905102951428412D-1,6.05183413124413191D-2, &
         2.33520497626869185D-3/

  X = ARG
  Y = ABS(X)
  if (Y <= THRESHOLD) then

!------------------------------------------------------------------
!  Evaluate  erf  for  |X| <= 0.46875
!------------------------------------------------------------------
      YSQ = ZERO
      if (Y > XSMALL) YSQ = Y * Y
      XNUM = A(5)*YSQ
      XDEN = YSQ

      do I = 1, 3
         XNUM = (XNUM + A(I)) * YSQ
         XDEN = (XDEN + B(I)) * YSQ
      enddo

      RESULT = X * (XNUM + A(4)) / (XDEN + B(4))
      if (JINT  /=  0) RESULT = ONE - RESULT
      if (JINT  ==  2) RESULT = EXP(YSQ) * RESULT
      goto 800

!------------------------------------------------------------------
!  Evaluate  erfc  for 0.46875 <= |X| <= 4.0
!------------------------------------------------------------------
   else if (Y <= FOUR) then
      XNUM = C(9)*Y
      XDEN = Y

      do I = 1, 7
         XNUM = (XNUM + C(I)) * Y
         XDEN = (XDEN + D(I)) * Y
      enddo

      RESULT = (XNUM + C(8)) / (XDEN + D(8))
      if (JINT  /=  2) then
         YSQ = AINT(Y*SIXTEEN)/SIXTEEN
         DEL = (Y-YSQ)*(Y+YSQ)
         RESULT = EXP(-YSQ*YSQ) * EXP(-DEL) * RESULT
      endif

!------------------------------------------------------------------
!  Evaluate  erfc  for |X| > 4.0
!------------------------------------------------------------------
   else
      RESULT = ZERO
      if (Y >= XBIG) then
         if (JINT /= 2 .OR. Y >= XMAX) goto 300
         if (Y >= XHUGE) then
            RESULT = SQRPI / Y
            goto 300
         endif
      endif
      YSQ = ONE / (Y * Y)
      XNUM = P(6)*YSQ
      XDEN = YSQ

      do I = 1, 4
         XNUM = (XNUM + P(I)) * YSQ
         XDEN = (XDEN + Q(I)) * YSQ
      enddo

      RESULT = YSQ *(XNUM + P(5)) / (XDEN + Q(5))
      RESULT = (SQRPI -  RESULT) / Y
      if (JINT /= 2) then
         YSQ = AINT(Y*SIXTEEN)/SIXTEEN
         DEL = (Y-YSQ)*(Y+YSQ)
         RESULT = EXP(-YSQ*YSQ) * EXP(-DEL) * RESULT
      endif
  endif

!------------------------------------------------------------------
!  Fix up for negative argument, erf, etc.
!------------------------------------------------------------------
  300 if (JINT == 0) then
      RESULT = (HALF - RESULT) + HALF
      if (X < ZERO) RESULT = -RESULT
   else if (JINT == 1) then
      if (X < ZERO) RESULT = TWO - RESULT
   else
      if (X < ZERO) then
         if (X < XNEG) then
               RESULT = XINF
            else
               YSQ = AINT(X*SIXTEEN)/SIXTEEN
               DEL = (X-YSQ)*(X+YSQ)
               Y = EXP(YSQ*YSQ) * EXP(DEL)
               RESULT = (Y+Y) - RESULT
         endif
      endif
  endif

  800 return

  end subroutine calerf

!--------------------------------------------------------------------

  double precision function netlib_specfun_erf(X)

! This subprogram computes approximate values for erf(x).
!   (see comments heading CALERF).
!
!   Author/date: William J. Cody, January 8, 1985

  implicit none

  integer JINT
  double precision X, RESULT

  JINT = 0
  call calerf(X,RESULT,JINT)
  netlib_specfun_erf = RESULT

  end function netlib_specfun_erf

!
! Subject: RE: Can one freely use and redistribute Fortran routines "specfun" from Netlib?
! From: Jack Dongarra
! Date: Wed, 21 Nov 2007 10:33:45 -0500
! To: Rusty Lusk, Dimitri Komatitsch
!
! Yes the code can freely be used and incorporated into other software. You
! should of course acknowledge the use of the software.
!
! Hope this helps,
!
! Jack Dongarra
!
! **********************************************************************
! Prof. Jack Dongarra; Innovative Computing Laboratory; EECS Department;
! 1122 Volunteer Blvd; University of Tennessee; Knoxville TN 37996-3450;
! +1-865-974-8295; http://www.cs.utk.edu/~dongarra/
!
! -----Original Message-----
! From: Rusty Lusk
! Sent: Wednesday, November 21, 2007 10:29 AM
! To: Dimitri Komatitsch
! Cc: Jack Dongarra
! Subject: Re: Can one freely use and redistribute Fortran routines "specfun"
! from Netlib?
!
! Netlib is managed at the University of Tennesee, not Argonne at this
! point. I have copied Jack Dongarra on this reply; he should be able
! to answer questions about licensing issues for code from Netlib.
!
! Regards,
! Rusty
!
! On Nov 21, 2007, at 8:36 AM, Dimitri Komatitsch wrote:
!
! >
! > Dear Sir,
! >
! > Can one freely use and redistribute Fortran routines "specfun" from
! > Netlib http://netlib2.cs.utk.edu/specfun/
! > which were written back in 1985-1990 by William J. Cody
! > from the Mathematics and Computer Science Division at Argonne?
! >
! > We use one of these routines (the error function, erf())
! > in one of our source codes, which we would like to
! > release as open source under GPL v2+, and we therefore
! > wonder if we could include that erf() routine in the
! > package in a separate file (of course saying in a comment in the
! > header that it comes from Netlib and was written by William J. Cody from
! > Argonne).
! >
! > Thank you,
! > Best regards,
! >
! > Dimitri Komatitsch.
! >
! > --
! > Dimitri Komatitsch - dimitri.komatitsch aT univ-pau.fr
! > Professor, University of Pau, Institut universitaire de France
! > and INRIA Magique3D, France   http://www.univ-pau.fr/~dkomati1
! >
