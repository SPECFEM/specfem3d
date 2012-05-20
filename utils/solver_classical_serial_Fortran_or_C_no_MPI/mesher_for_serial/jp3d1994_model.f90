!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  4 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory, California Institute of Technology, USA
!             and University of Pau / CNRS / INRIA, France
! (c) California Institute of Technology and University of Pau / CNRS / INRIA
!                            February 2008
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================


!=====================================================================
!
!         Last Time Modified by Min Chen, Caltech, 03/14/2008
!
!          Program ----- veljp3d.f -----
!
!       This program is used to calculate 3-D P-wave velocity
!    distribution beneath the Japan Islands which is obtained
!    by a simultaneous inversion of arrival time data from local,
!    regional and teleseismic events.  For details, see "Deep
!    structure of the Japan subduction zone as derived from local,
!    regional, and teleseismic events" by Zhao, Hasegawa & Kanamori,
!    JGR, 99, 22313-22329, 1994.
!
!       The meaningful range of this model is as follows:
!        latitude : 32 - 45 N
!        longitude: 130-145 E
!        depth    : 0  - 500 km
!
!                            Dapeng Zhao
!                            Dept. of Earth & Planet. Sci
!                            Washington University
!                            St. Louis, MO 63130
!                            U.S.A.
!                            dapeng@izu.wustl.edu
!=========================================================================
subroutine read_iso3d_dpzhao_model(JP3DM_V)

  implicit none

  include "constants.h"
! jp3d_model_variables
  type jp3d_model_variables
    sequence
! vmod3d
  integer :: NPA
  integer :: NRA
  integer :: NHA
  integer :: NPB
  integer :: NRB
  integer :: NHB
  double precision :: PNA(MPA)
  double precision :: RNA(MRA)
  double precision :: HNA(MHA)
  double precision :: PNB(MPB)
  double precision :: RNB(MRB)
  double precision :: HNB(MHB)
  double precision :: VELAP(MPA,MRA,MHA)
  double precision :: VELBP(MPB,MRB,MHB)
! discon
  double precision :: PN(51)
  double precision :: RRN(63)
  double precision :: DEPA(51,63)
  double precision :: DEPB(51,63)
  double precision :: DEPC(51,63)
! locate
  integer :: IPLOCA(MKA)
  integer :: IRLOCA(MKA)
  integer :: IHLOCA(MKA)
  integer :: IPLOCB(MKB)
  integer :: IRLOCB(MKB)
  integer :: IHLOCB(MKB)
  double precision :: PLA
  double precision :: RLA
  double precision :: HLA
  double precision :: PLB
  double precision :: RLB
  double precision :: HLB
! weight
  integer :: IP
  integer :: JP
  integer :: KP
  integer :: IP1
  integer :: JP1
  integer :: KP1
  double precision :: WV(8)
! prhfd
  double precision :: P
  double precision :: R
  double precision :: H
  double precision :: PF
  double precision :: RF
  double precision :: HF
  double precision :: PF1
  double precision :: RF1
  double precision :: HF1
  double precision :: PD
  double precision :: RD
  double precision :: HD
! jpmodv
  double precision :: VP(29)
  double precision :: VS(29)
  double precision :: RA(29)
  double precision :: DEPJ(29)
  end type jp3d_model_variables

  type (jp3d_model_variables) JP3DM_V
! jp3d_model_variables

      OPEN(2,FILE="DATA/Zhao_JP_model/m3d1341")
      OPEN(3,FILE="DATA/Zhao_JP_model/datadis")

      CALL INPUTJP(JP3DM_V)
      CALL INPUT1(JP3DM_V)
      CALL INPUT2(JP3DM_V)

end subroutine read_iso3d_dpzhao_model
!==========================================================================
subroutine iso3d_dpzhao_model(radius,theta,phi,vp,vs,dvp,dvs,rho,found_crust,JP3DM_V)
  implicit none

  include "constants.h"

! jp3d_model_variables
  type jp3d_model_variables
    sequence
! vmod3d
  integer :: NPA
  integer :: NRA
  integer :: NHA
  integer :: NPB
  integer :: NRB
  integer :: NHB
  double precision :: PNA(MPA)
  double precision :: RNA(MRA)
  double precision :: HNA(MHA)
  double precision :: PNB(MPB)
  double precision :: RNB(MRB)
  double precision :: HNB(MHB)
  double precision :: VELAP(MPA,MRA,MHA)
  double precision :: VELBP(MPB,MRB,MHB)
! discon
  double precision :: PN(51)
  double precision :: RRN(63)
  double precision :: DEPA(51,63)
  double precision :: DEPB(51,63)
  double precision :: DEPC(51,63)
! locate
  integer :: IPLOCA(MKA)
  integer :: IRLOCA(MKA)
  integer :: IHLOCA(MKA)
  integer :: IPLOCB(MKB)
  integer :: IRLOCB(MKB)
  integer :: IHLOCB(MKB)
  double precision :: PLA
  double precision :: RLA
  double precision :: HLA
  double precision :: PLB
  double precision :: RLB
  double precision :: HLB
! weight
  integer :: IP
  integer :: JP
  integer :: KP
  integer :: IP1
  integer :: JP1
  integer :: KP1
  double precision :: WV(8)
! prhfd
  double precision :: P
  double precision :: R
  double precision :: H
  double precision :: PF
  double precision :: RF
  double precision :: HF
  double precision :: PF1
  double precision :: RF1
  double precision :: HF1
  double precision :: PD
  double precision :: RD
  double precision :: HD
! jpmodv
  double precision :: VP(29)
  double precision :: VS(29)
  double precision :: RA(29)
  double precision :: DEPJ(29)
  end type jp3d_model_variables

  type (jp3d_model_variables) JP3DM_V
! jp3d_model_variables

  logical found_crust
  double precision :: radius,theta,phi,vp,vs,dvs,dvp,rho
  double precision :: PE,RE,HE,H1,H2,H3,scaleval
  integer :: LAY


  found_crust = .false.

  PE = theta
  RE = phi
  HE = (ONE - radius)*R_EARTH_KM
!  calculate depths of the Conrad, the Moho and
!  the plate boundary beneath the location (PHI,RAM)
  CALL HLAY(PE,RE,H1,1,JP3DM_V)
  CALL HLAY(PE,RE,H2,2,JP3DM_V)
  CALL HLAY(PE,RE,H3,3,JP3DM_V)
!   when LAY = 1, the focus is in the upper crust;
!   when LAY = 2, the focus is in the lower crust;
!   when LAY = 3, the focus is in the mantle wedge;
!   when LAY = 4, the focus is beneath the plate boundary.
  IF(HE.LE.H1)                   THEN
     LAY = 1
     found_crust = .true.
  ELSE IF(HE.GT.H1.AND.HE.LE.H2) THEN
     LAY = 2
     found_crust = .true.
  ELSE IF(HE.GT.H2.AND.HE.LE.H3) THEN
     LAY = 3
  ELSE
     LAY = 4
  END IF
  CALL VEL1D(HE,vp,LAY,1,JP3DM_V)
  CALL VEL1D(HE,vs,LAY,2,JP3DM_V)

  CALL VEL3(PE,RE,HE,dvp,LAY,JP3DM_V)
  dvp = 0.01d0*dvp
  dvs = 1.5d0*dvp
  vp = vp*(1.0d0+dvp)
  vs = vs*(1.0d0+dvs)

! determine rho
  if(LAY .eq. 1) then
     rho=2.6
  endif
  if(LAY .eq. 2) then
     rho=2.9
  endif
  if(LAY .GT. 2) then
     rho=3.3+(vs-4.4)*0.66667
  endif
! non-dimensionalize
! time scaling (s^{-1}) is done with scaleval
  scaleval=dsqrt(PI*GRAV*RHOAV)
  rho=rho*1000.0d0/RHOAV
  vp=vp*1000.0d0/(R_EARTH*scaleval)
  vs=vs*1000.0d0/(R_EARTH*scaleval)
END subroutine iso3d_dpzhao_model

!---------------------------------------------------------------

      SUBROUTINE INPUT1(JP3DM_V)
   implicit none

   include "constants.h"
! jp3d_model_variables
  type jp3d_model_variables
    sequence
! vmod3d
  integer :: NPA
  integer :: NRA
  integer :: NHA
  integer :: NPB
  integer :: NRB
  integer :: NHB
  double precision :: PNA(MPA)
  double precision :: RNA(MRA)
  double precision :: HNA(MHA)
  double precision :: PNB(MPB)
  double precision :: RNB(MRB)
  double precision :: HNB(MHB)
  double precision :: VELAP(MPA,MRA,MHA)
  double precision :: VELBP(MPB,MRB,MHB)
! discon
  double precision :: PN(51)
  double precision :: RRN(63)
  double precision :: DEPA(51,63)
  double precision :: DEPB(51,63)
  double precision :: DEPC(51,63)
! locate
  integer :: IPLOCA(MKA)
  integer :: IRLOCA(MKA)
  integer :: IHLOCA(MKA)
  integer :: IPLOCB(MKB)
  integer :: IRLOCB(MKB)
  integer :: IHLOCB(MKB)
  double precision :: PLA
  double precision :: RLA
  double precision :: HLA
  double precision :: PLB
  double precision :: RLB
  double precision :: HLB
! weight
  integer :: IP
  integer :: JP
  integer :: KP
  integer :: IP1
  integer :: JP1
  integer :: KP1
  double precision :: WV(8)
! prhfd
  double precision :: P
  double precision :: R
  double precision :: H
  double precision :: PF
  double precision :: RF
  double precision :: HF
  double precision :: PF1
  double precision :: RF1
  double precision :: HF1
  double precision :: PD
  double precision :: RD
  double precision :: HD
! jpmodv
  double precision :: VP(29)
  double precision :: VS(29)
  double precision :: RA(29)
  double precision :: DEPJ(29)
  end type jp3d_model_variables

  type (jp3d_model_variables) JP3DM_V
! jp3d_model_variables

100     FORMAT(3I3)
      READ(2,100)  JP3DM_V%NPA,JP3DM_V%NRA,JP3DM_V%NHA
      CALL PUT1(JP3DM_V%NPA,JP3DM_V%NRA,JP3DM_V%NHA,JP3DM_V%PNA,JP3DM_V%RNA,JP3DM_V%HNA,JP3DM_V%VELAP)
      READ(2,100)  JP3DM_V%NPB,JP3DM_V%NRB,JP3DM_V%NHB
      CALL PUT1(JP3DM_V%NPB,JP3DM_V%NRB,JP3DM_V%NHB,JP3DM_V%PNB,JP3DM_V%RNB,JP3DM_V%HNB,JP3DM_V%VELBP)
      CALL BLDMAP(JP3DM_V)
      RETURN
    END SUBROUTINE INPUT1

      SUBROUTINE PUT1(NPX,NRX,NHX,PNX,RNX,HNX,VELXP)
      integer :: NPX,NRX,NHX,K,I,J
      double precision ::  VELXP(NPX,NRX,NHX), &
                PNX(NPX),RNX(NRX),HNX(NHX)
      READ(2,110) (PNX(I),I=1,NPX)
      READ(2,110) (RNX(I),I=1,NRX)
      READ(2,120) (HNX(I),I=1,NHX)
      DO K = 1,NHX
         DO I = 1,NPX
            READ(2,140) (VELXP(I,J,K),J=1,NRX)
110         FORMAT(6(9F7.2/))
120         FORMAT(3(8F7.2/))
140         FORMAT(4(14F5.2/))
         enddo
      enddo
    END SUBROUTINE PUT1

      SUBROUTINE INPUT2(JP3DM_V)
  implicit none

  include "constants.h"

! jp3d_model_variables
  type jp3d_model_variables
    sequence
! vmod3d
  integer :: NPA
  integer :: NRA
  integer :: NHA
  integer :: NPB
  integer :: NRB
  integer :: NHB
  double precision :: PNA(MPA)
  double precision :: RNA(MRA)
  double precision :: HNA(MHA)
  double precision :: PNB(MPB)
  double precision :: RNB(MRB)
  double precision :: HNB(MHB)
  double precision :: VELAP(MPA,MRA,MHA)
  double precision :: VELBP(MPB,MRB,MHB)
! discon
  double precision :: PN(51)
  double precision :: RRN(63)
  double precision :: DEPA(51,63)
  double precision :: DEPB(51,63)
  double precision :: DEPC(51,63)
! locate
  integer :: IPLOCA(MKA)
  integer :: IRLOCA(MKA)
  integer :: IHLOCA(MKA)
  integer :: IPLOCB(MKB)
  integer :: IRLOCB(MKB)
  integer :: IHLOCB(MKB)
  double precision :: PLA
  double precision :: RLA
  double precision :: HLA
  double precision :: PLB
  double precision :: RLB
  double precision :: HLB
! weight
  integer :: IP
  integer :: JP
  integer :: KP
  integer :: IP1
  integer :: JP1
  integer :: KP1
  double precision :: WV(8)
! prhfd
  double precision :: P
  double precision :: R
  double precision :: H
  double precision :: PF
  double precision :: RF
  double precision :: HF
  double precision :: PF1
  double precision :: RF1
  double precision :: HF1
  double precision :: PD
  double precision :: RD
  double precision :: HD
! jpmodv
  double precision :: VP(29)
  double precision :: VS(29)
  double precision :: RA(29)
  double precision :: DEPJ(29)
  end type jp3d_model_variables

  type (jp3d_model_variables) JP3DM_V
! jp3d_model_variables

      integer :: NP,NNR,I,J
      READ(3,100)  NP,NNR
      READ(3,110) (JP3DM_V%PN(I),I=1,NP)
      READ(3,120) (JP3DM_V%RRN(I),I=1,NNR)
      DO 1  I = NP,1,-1
      READ(3,130) (JP3DM_V%DEPA(I,J),J=1,NNR)
1     CONTINUE
      DO 2  I = NP,1,-1
      READ(3,130) (JP3DM_V%DEPB(I,J),J=1,NNR)
2     CONTINUE
      DO 3  I = NP,1,-1
      READ(3,130) (JP3DM_V%DEPC(I,J),J=1,NNR)
3     CONTINUE
100   FORMAT(2I6)
110   FORMAT(5(10F7.2/),F7.2)
120   FORMAT(6(10F7.2/),3F7.2)
130   FORMAT(6(10F7.1/),3F7.1)
      RETURN
      END

      SUBROUTINE BLDMAP(JP3DM_V)
  implicit none

  include "constants.h"
! jp3d_model_variables
  type jp3d_model_variables
    sequence
! vmod3d
  integer :: NPA
  integer :: NRA
  integer :: NHA
  integer :: NPB
  integer :: NRB
  integer :: NHB
  double precision :: PNA(MPA)
  double precision :: RNA(MRA)
  double precision :: HNA(MHA)
  double precision :: PNB(MPB)
  double precision :: RNB(MRB)
  double precision :: HNB(MHB)
  double precision :: VELAP(MPA,MRA,MHA)
  double precision :: VELBP(MPB,MRB,MHB)
! discon
  double precision :: PN(51)
  double precision :: RRN(63)
  double precision :: DEPA(51,63)
  double precision :: DEPB(51,63)
  double precision :: DEPC(51,63)
! locate
  integer :: IPLOCA(MKA)
  integer :: IRLOCA(MKA)
  integer :: IHLOCA(MKA)
  integer :: IPLOCB(MKB)
  integer :: IRLOCB(MKB)
  integer :: IHLOCB(MKB)
  double precision :: PLA
  double precision :: RLA
  double precision :: HLA
  double precision :: PLB
  double precision :: RLB
  double precision :: HLB
! weight
  integer :: IP
  integer :: JP
  integer :: KP
  integer :: IP1
  integer :: JP1
  integer :: KP1
  double precision :: WV(8)
! prhfd
  double precision :: P
  double precision :: R
  double precision :: H
  double precision :: PF
  double precision :: RF
  double precision :: HF
  double precision :: PF1
  double precision :: RF1
  double precision :: HF1
  double precision :: PD
  double precision :: RD
  double precision :: HD
! jpmodv
  double precision :: VP(29)
  double precision :: VS(29)
  double precision :: RA(29)
  double precision :: DEPJ(29)
  end type jp3d_model_variables

  type (jp3d_model_variables) JP3DM_V
! jp3d_model_variables

      CALL LOCX(JP3DM_V%PNA,JP3DM_V%RNA,JP3DM_V%HNA,JP3DM_V%NPA,JP3DM_V%NRA,JP3DM_V%NHA,MKA, &
           JP3DM_V%PLA,JP3DM_V%RLA,JP3DM_V%HLA,JP3DM_V%IPLOCA,JP3DM_V%IRLOCA,JP3DM_V%IHLOCA)
      CALL LOCX(JP3DM_V%PNB,JP3DM_V%RNB,JP3DM_V%HNB,JP3DM_V%NPB,JP3DM_V%NRB,JP3DM_V%NHB,MKB, &
           JP3DM_V%PLB,JP3DM_V%RLB,JP3DM_V%HLB,JP3DM_V%IPLOCB,JP3DM_V%IRLOCB,JP3DM_V%IHLOCB)
      RETURN
      END

      SUBROUTINE LOCX(PNX,RNX,HNX,NPX,NRX,NHX,MKX, &
                 PLX,RLX,HLX,IPLOCX,IRLOCX,IHLOCX)
     integer ::  NPX,NRX,NHX,MKX,IPLOCX(MKX),IRLOCX(MKX),IHLOCX(MKX)
     integer ::  IPMAX,IP,IP1,IRMAX,IR,IR1,IH1,IH,IHMAX,I
     double precision :: PNX(NPX),RNX(NRX),HNX(NHX)
     double precision :: PLX,RLX,HLX,PNOW,RNOW,HNOW
      PLX      = 1.0-PNX(1)*100.0
      IPMAX    = IDNINT(PNX(NPX)*100.0+PLX)
      IP       = 1
      DO 10 I  = 1,IPMAX
      IP1      = IP+1
      PNOW     = (FLOAT(I)-PLX)/100.0
      IF(PNOW.GE.PNX(IP1))   IP = IP1
      IPLOCX(I)= IP
10    CONTINUE
      RLX      = 1.0-RNX(1)*100.0
      IRMAX    = IDNINT(RNX(NRX)*100.0+RLX)
      IR       = 1
      DO 20 I  = 1,IRMAX
      IR1      = IR+1
      RNOW     = (FLOAT(I)-RLX)/100.0
      IF(RNOW.GE.RNX(IR1))   IR = IR1
      IRLOCX(I)= IR
20    CONTINUE
      HLX      = 1.0-HNX(1)
      IHMAX    = IDNINT(HNX(NHX)+HLX)
      IH       = 1
      DO 30 I  = 1,IHMAX
      IH1      = IH+1
      HNOW     = FLOAT(I)-HLX
      IF(HNOW.GE.HNX(IH1))   IH = IH1
      IHLOCX(I)= IH
30    CONTINUE
      RETURN
      END

      SUBROUTINE VEL3(PE,RE,HE,V,LAY,JP3DM_V)
        implicit none

        include "constants.h"
! jp3d_model_variables
  type jp3d_model_variables
    sequence
! vmod3d
  integer :: NPA
  integer :: NRA
  integer :: NHA
  integer :: NPB
  integer :: NRB
  integer :: NHB
  double precision :: PNA(MPA)
  double precision :: RNA(MRA)
  double precision :: HNA(MHA)
  double precision :: PNB(MPB)
  double precision :: RNB(MRB)
  double precision :: HNB(MHB)
  double precision :: VELAP(MPA,MRA,MHA)
  double precision :: VELBP(MPB,MRB,MHB)
! discon
  double precision :: PN(51)
  double precision :: RRN(63)
  double precision :: DEPA(51,63)
  double precision :: DEPB(51,63)
  double precision :: DEPC(51,63)
! locate
  integer :: IPLOCA(MKA)
  integer :: IRLOCA(MKA)
  integer :: IHLOCA(MKA)
  integer :: IPLOCB(MKB)
  integer :: IRLOCB(MKB)
  integer :: IHLOCB(MKB)
  double precision :: PLA
  double precision :: RLA
  double precision :: HLA
  double precision :: PLB
  double precision :: RLB
  double precision :: HLB
! weight
  integer :: IP
  integer :: JP
  integer :: KP
  integer :: IP1
  integer :: JP1
  integer :: KP1
  double precision :: WV(8)
! prhfd
  double precision :: P
  double precision :: R
  double precision :: H
  double precision :: PF
  double precision :: RF
  double precision :: HF
  double precision :: PF1
  double precision :: RF1
  double precision :: HF1
  double precision :: PD
  double precision :: RD
  double precision :: HD
! jpmodv
  double precision :: VP(29)
  double precision :: VS(29)
  double precision :: RA(29)
  double precision :: DEPJ(29)
  end type jp3d_model_variables

  type (jp3d_model_variables) JP3DM_V
! jp3d_model_variables

       double precision :: PE,RE,HE,V

       integer :: LAY

        JP3DM_V%P     = 90.0-PE/DEGREES_TO_RADIANS
        JP3DM_V%R     = RE/DEGREES_TO_RADIANS
        JP3DM_V%H     = HE
        IF(LAY.LE.3)       THEN
           CALL PRHF(JP3DM_V%IPLOCA,JP3DM_V%IRLOCA,JP3DM_V%IHLOCA,JP3DM_V%PLA,JP3DM_V%RLA,JP3DM_V%HLA, &
                JP3DM_V%PNA,JP3DM_V%RNA,JP3DM_V%HNA,MPA,MRA,MHA,MKA,JP3DM_V)
        ELSE IF(LAY.EQ.4)  THEN
           CALL PRHF(JP3DM_V%IPLOCB,JP3DM_V%IRLOCB,JP3DM_V%IHLOCB,JP3DM_V%PLB,JP3DM_V%RLB,JP3DM_V%HLB, &
                JP3DM_V%PNB,JP3DM_V%RNB,JP3DM_V%HNB,MPB,MRB,MHB,MKB,JP3DM_V)
        ELSE
        END IF
        JP3DM_V%WV(1) = JP3DM_V%PF1*JP3DM_V%RF1*JP3DM_V%HF1
        JP3DM_V%WV(2) = JP3DM_V%PF*JP3DM_V%RF1*JP3DM_V%HF1
        JP3DM_V%WV(3) = JP3DM_V%PF1*JP3DM_V%RF*JP3DM_V%HF1
        JP3DM_V%WV(4) = JP3DM_V%PF*JP3DM_V%RF*JP3DM_V%HF1
        JP3DM_V%WV(5) = JP3DM_V%PF1*JP3DM_V%RF1*JP3DM_V%HF
        JP3DM_V%WV(6) = JP3DM_V%PF*JP3DM_V%RF1*JP3DM_V%HF
        JP3DM_V%WV(7) = JP3DM_V%PF1*JP3DM_V%RF*JP3DM_V%HF
        JP3DM_V%WV(8) = JP3DM_V%PF*JP3DM_V%RF*JP3DM_V%HF
        !   calculate velocity
        IF(LAY.LE.3)      THEN
           CALL VABPS(MPA,MRA,MHA,JP3DM_V%VELAP,V,JP3DM_V)
        ELSE IF(LAY.EQ.4) THEN
           CALL VABPS(MPB,MRB,MHB,JP3DM_V%VELBP,V,JP3DM_V)
        ELSE
        END IF

        RETURN
      END SUBROUTINE VEL3

      SUBROUTINE VABPS(MP,MR,MH,V,VEL,JP3DM_V)
  implicit none

  include "constants.h"


! jp3d_model_variables
  type jp3d_model_variables
    sequence
! vmod3d
  integer :: NPA
  integer :: NRA
  integer :: NHA
  integer :: NPB
  integer :: NRB
  integer :: NHB
  double precision :: PNA(MPA)
  double precision :: RNA(MRA)
  double precision :: HNA(MHA)
  double precision :: PNB(MPB)
  double precision :: RNB(MRB)
  double precision :: HNB(MHB)
  double precision :: VELAP(MPA,MRA,MHA)
  double precision :: VELBP(MPB,MRB,MHB)
! discon
  double precision :: PN(51)
  double precision :: RRN(63)
  double precision :: DEPA(51,63)
  double precision :: DEPB(51,63)
  double precision :: DEPC(51,63)
! locate
  integer :: IPLOCA(MKA)
  integer :: IRLOCA(MKA)
  integer :: IHLOCA(MKA)
  integer :: IPLOCB(MKB)
  integer :: IRLOCB(MKB)
  integer :: IHLOCB(MKB)
  double precision :: PLA
  double precision :: RLA
  double precision :: HLA
  double precision :: PLB
  double precision :: RLB
  double precision :: HLB
! weight
  integer :: IP
  integer :: JP
  integer :: KP
  integer :: IP1
  integer :: JP1
  integer :: KP1
  double precision :: WV(8)
! prhfd
  double precision :: P
  double precision :: R
  double precision :: H
  double precision :: PF
  double precision :: RF
  double precision :: HF
  double precision :: PF1
  double precision :: RF1
  double precision :: HF1
  double precision :: PD
  double precision :: RD
  double precision :: HD
! jpmodv
  double precision :: VP(29)
  double precision :: VS(29)
  double precision :: RA(29)
  double precision :: DEPJ(29)
  end type jp3d_model_variables

  type (jp3d_model_variables) JP3DM_V
! jp3d_model_variables
      double precision :: VEL
      integer :: MP,MR,MH
      double precision :: V(MP,MR,MH)
      VEL = JP3DM_V%WV(1)*V(JP3DM_V%IP,JP3DM_V%JP,JP3DM_V%KP)  + JP3DM_V%WV(2)*V(JP3DM_V%IP1,JP3DM_V%JP,JP3DM_V%KP) &
          + JP3DM_V%WV(3)*V(JP3DM_V%IP,JP3DM_V%JP1,JP3DM_V%KP) + JP3DM_V%WV(4)*V(JP3DM_V%IP1,JP3DM_V%JP1,JP3DM_V%KP) &
          + JP3DM_V%WV(5)*V(JP3DM_V%IP,JP3DM_V%JP,JP3DM_V%KP1) + JP3DM_V%WV(6)*V(JP3DM_V%IP1,JP3DM_V%JP,JP3DM_V%KP1) &
          + JP3DM_V%WV(7)*V(JP3DM_V%IP,JP3DM_V%JP1,JP3DM_V%KP1)+ JP3DM_V%WV(8)*V(JP3DM_V%IP1,JP3DM_V%JP1,JP3DM_V%KP1)
      RETURN
      END

      SUBROUTINE INTMAP(R,IRLOC,NNR,RL,IR)
      integer :: NNR,IRLOC(NNR),IS,IR
      double precision :: R,RL
      IS      = IDNINT(R+RL)
      IR      = IRLOC(IS)
      RETURN
      END

      SUBROUTINE PRHF(IPLOCX,IRLOCX,IHLOCX,PLX,RLX,HLX, &
                      PNX,RNX,HNX,MPX,MRX,MHX,MKX,JP3DM_V)
  implicit none

  include "constants.h"

! jp3d_model_variables
  type jp3d_model_variables
    sequence
! vmod3d
  integer :: NPA
  integer :: NRA
  integer :: NHA
  integer :: NPB
  integer :: NRB
  integer :: NHB
  double precision :: PNA(MPA)
  double precision :: RNA(MRA)
  double precision :: HNA(MHA)
  double precision :: PNB(MPB)
  double precision :: RNB(MRB)
  double precision :: HNB(MHB)
  double precision :: VELAP(MPA,MRA,MHA)
  double precision :: VELBP(MPB,MRB,MHB)
! discon
  double precision :: PN(51)
  double precision :: RRN(63)
  double precision :: DEPA(51,63)
  double precision :: DEPB(51,63)
  double precision :: DEPC(51,63)
! locate
  integer :: IPLOCA(MKA)
  integer :: IRLOCA(MKA)
  integer :: IHLOCA(MKA)
  integer :: IPLOCB(MKB)
  integer :: IRLOCB(MKB)
  integer :: IHLOCB(MKB)
  double precision :: PLA
  double precision :: RLA
  double precision :: HLA
  double precision :: PLB
  double precision :: RLB
  double precision :: HLB
! weight
  integer :: IP
  integer :: JP
  integer :: KP
  integer :: IP1
  integer :: JP1
  integer :: KP1
  double precision :: WV(8)
! prhfd
  double precision :: P
  double precision :: R
  double precision :: H
  double precision :: PF
  double precision :: RF
  double precision :: HF
  double precision :: PF1
  double precision :: RF1
  double precision :: HF1
  double precision :: PD
  double precision :: RD
  double precision :: HD
! jpmodv
  double precision :: VP(29)
  double precision :: VS(29)
  double precision :: RA(29)
  double precision :: DEPJ(29)
  end type jp3d_model_variables

  type (jp3d_model_variables) JP3DM_V
! jp3d_model_variables

        integer :: MPX,MRX,MHX,MKX
        integer ::  IPLOCX(MKX),IRLOCX(MKX),IHLOCX(MKX)
        double precision :: PNX(MPX),RNX(MRX),HNX(MHX)
        double precision :: PLX,RLX,HLX
      CALL LIMIT(PNX(1),PNX(MPX),JP3DM_V%P)
      CALL LIMIT(RNX(1),RNX(MRX),JP3DM_V%R)
      CALL LIMIT(HNX(1),HNX(MHX),JP3DM_V%H)
      CALL INTMAP(JP3DM_V%P*100.0,IPLOCX,MKX,PLX,JP3DM_V%IP)
      CALL INTMAP(JP3DM_V%R*100.0,IRLOCX,MKX,RLX,JP3DM_V%JP)
      CALL INTMAP(JP3DM_V%H,IHLOCX,MKX,HLX,JP3DM_V%KP)
      JP3DM_V%IP1   = JP3DM_V%IP+1
      JP3DM_V%JP1   = JP3DM_V%JP+1
      JP3DM_V%KP1   = JP3DM_V%KP+1
      JP3DM_V%PD    = PNX(JP3DM_V%IP1)-PNX(JP3DM_V%IP)
      JP3DM_V%RD    = RNX(JP3DM_V%JP1)-RNX(JP3DM_V%JP)
      JP3DM_V%HD    = HNX(JP3DM_V%KP1)-HNX(JP3DM_V%KP)
      JP3DM_V%PF    = (JP3DM_V%P-PNX(JP3DM_V%IP))/JP3DM_V%PD
      JP3DM_V%RF    = (JP3DM_V%R-RNX(JP3DM_V%JP))/JP3DM_V%RD
      JP3DM_V%HF    = (JP3DM_V%H-HNX(JP3DM_V%KP))/JP3DM_V%HD
      JP3DM_V%PF1   = 1.0-JP3DM_V%PF
      JP3DM_V%RF1   = 1.0-JP3DM_V%RF
      JP3DM_V%HF1   = 1.0-JP3DM_V%HF
      RETURN
      END

      SUBROUTINE HLAY(PE,RE,HE,IJK,JP3DM_V)
  implicit none

  include "constants.h"
! jp3d_model_variables
  type jp3d_model_variables
    sequence
! vmod3d
  integer :: NPA
  integer :: NRA
  integer :: NHA
  integer :: NPB
  integer :: NRB
  integer :: NHB
  double precision :: PNA(MPA)
  double precision :: RNA(MRA)
  double precision :: HNA(MHA)
  double precision :: PNB(MPB)
  double precision :: RNB(MRB)
  double precision :: HNB(MHB)
  double precision :: VELAP(MPA,MRA,MHA)
  double precision :: VELBP(MPB,MRB,MHB)
! discon
  double precision :: PN(51)
  double precision :: RRN(63)
  double precision :: DEPA(51,63)
  double precision :: DEPB(51,63)
  double precision :: DEPC(51,63)
! locate
  integer :: IPLOCA(MKA)
  integer :: IRLOCA(MKA)
  integer :: IHLOCA(MKA)
  integer :: IPLOCB(MKB)
  integer :: IRLOCB(MKB)
  integer :: IHLOCB(MKB)
  double precision :: PLA
  double precision :: RLA
  double precision :: HLA
  double precision :: PLB
  double precision :: RLB
  double precision :: HLB
! weight
  integer :: IP
  integer :: JP
  integer :: KP
  integer :: IP1
  integer :: JP1
  integer :: KP1
  double precision :: WV(8)
! prhfd
  double precision :: P
  double precision :: R
  double precision :: H
  double precision :: PF
  double precision :: RF
  double precision :: HF
  double precision :: PF1
  double precision :: RF1
  double precision :: HF1
  double precision :: PD
  double precision :: RD
  double precision :: HD
! jpmodv
  double precision :: VP(29)
  double precision :: VS(29)
  double precision :: RA(29)
  double precision :: DEPJ(29)
  end type jp3d_model_variables

  type (jp3d_model_variables) JP3DM_V
! jp3d_model_variables
        double precision :: PE,RE,HE,WV1,WV2,WV3,WV4,P,R,PF,RF,PF1,RF1
        integer :: IJK,J,J1,I,I1
        P = 90.0-PE/DEGREES_TO_RADIANS
        R = RE/DEGREES_TO_RADIANS
        CALL LIMIT(JP3DM_V%PN(1),JP3DM_V%PN(51),P)
        CALL LIMIT(JP3DM_V%RRN(1),JP3DM_V%RRN(63),R)
        DO 1 I = 1,50
           I1     = I+1
           IF(P.GE.JP3DM_V%PN(I).AND.P.LT.JP3DM_V%PN(I1)) GO TO 11
1          CONTINUE
11         CONTINUE
           DO 2 J = 1,62
              J1     = J+1
              IF(R.GE.JP3DM_V%RRN(J).AND.R.LT.JP3DM_V%RRN(J1)) GO TO 22
2             CONTINUE
22            CONTINUE
              PF    = (P-JP3DM_V%PN(I))/(JP3DM_V%PN(I1)-JP3DM_V%PN(I))
              RF    = (R-JP3DM_V%RRN(J))/(JP3DM_V%RRN(J1)-JP3DM_V%RRN(J))
              PF1   = 1.0-PF
              RF1   = 1.0-RF
              WV1   = PF1*RF1
              WV2   = PF*RF1
              WV3   = PF1*RF
              WV4   = PF*RF
              IF(IJK.EQ.1)       THEN
                 HE  = WV1*JP3DM_V%DEPA(I,J)  + WV2*JP3DM_V%DEPA(I1,J) &
                      + WV3*JP3DM_V%DEPA(I,J1) + WV4*JP3DM_V%DEPA(I1,J1)
              ELSE IF(IJK.EQ.2)  THEN
                 HE  = WV1*JP3DM_V%DEPB(I,J)  + WV2*JP3DM_V%DEPB(I1,J) &
                      + WV3*JP3DM_V%DEPB(I,J1) + WV4*JP3DM_V%DEPB(I1,J1)
              ELSE IF(IJK.EQ.3)  THEN
                 HE  = WV1*JP3DM_V%DEPC(I,J)  + WV2*JP3DM_V%DEPC(I1,J) &
                      + WV3*JP3DM_V%DEPC(I,J1) + WV4*JP3DM_V%DEPC(I1,J1)
              ELSE
              END IF
              RETURN
            END SUBROUTINE HLAY

      SUBROUTINE LIMIT(C1,C2,C)
      double precision :: A1,A2,C1,C2,C
      A1    = dmin1(C1,C2)
      A2    = dmax1(C1,C2)
      IF(C.LT.A1)   C = A1
      IF(C.GT.A2)   C = A2
    END SUBROUTINE LIMIT

      SUBROUTINE VEL1D(HE,V,LAY,IPS,JP3DM_V)
  implicit none

  include "constants.h"
! jp3d_model_variables
  type jp3d_model_variables
    sequence
! vmod3d
  integer :: NPA
  integer :: NRA
  integer :: NHA
  integer :: NPB
  integer :: NRB
  integer :: NHB
  double precision :: PNA(MPA)
  double precision :: RNA(MRA)
  double precision :: HNA(MHA)
  double precision :: PNB(MPB)
  double precision :: RNB(MRB)
  double precision :: HNB(MHB)
  double precision :: VELAP(MPA,MRA,MHA)
  double precision :: VELBP(MPB,MRB,MHB)
! discon
  double precision :: PN(51)
  double precision :: RRN(63)
  double precision :: DEPA(51,63)
  double precision :: DEPB(51,63)
  double precision :: DEPC(51,63)
! locate
  integer :: IPLOCA(MKA)
  integer :: IRLOCA(MKA)
  integer :: IHLOCA(MKA)
  integer :: IPLOCB(MKB)
  integer :: IRLOCB(MKB)
  integer :: IHLOCB(MKB)
  double precision :: PLA
  double precision :: RLA
  double precision :: HLA
  double precision :: PLB
  double precision :: RLB
  double precision :: HLB
! weight
  integer :: IP
  integer :: JP
  integer :: KP
  integer :: IP1
  integer :: JP1
  integer :: KP1
  double precision :: WV(8)
! prhfd
  double precision :: P
  double precision :: R
  double precision :: H
  double precision :: PF
  double precision :: RF
  double precision :: HF
  double precision :: PF1
  double precision :: RF1
  double precision :: HF1
  double precision :: PD
  double precision :: RD
  double precision :: HD
! jpmodv
  double precision :: VP(29)
  double precision :: VS(29)
  double precision :: RA(29)
  double precision :: DEPJ(29)
  end type jp3d_model_variables

  type (jp3d_model_variables) JP3DM_V
! jp3d_model_variables

      integer :: IPS,LAY
      double precision :: HE,V,VM,HM
      IF(LAY.EQ.1)      THEN
        V    = 6.0
        IF(IPS.EQ.2)    V = 3.5
      ELSE IF(LAY.EQ.2) THEN
        V    = 6.7
        IF(IPS.EQ.2)    V = 3.8
      ELSE IF(LAY.GE.3) THEN
        HM   = 40.0
        IF(HE.LT.HM)    THEN
          CALL JPMODEL(IPS,HM,VM,JP3DM_V)
          V  = VM-(HM-HE)*0.003
        ELSE
          CALL JPMODEL(IPS,HE,V,JP3DM_V)
        END IF
      ELSE
      END IF
      RETURN
      END

      SUBROUTINE INPUTJP(JP3DM_V)
  implicit none

  include "constants.h"
! jp3d_model_variables
  type jp3d_model_variables
    sequence
! vmod3d
  integer :: NPA
  integer :: NRA
  integer :: NHA
  integer :: NPB
  integer :: NRB
  integer :: NHB
  double precision :: PNA(MPA)
  double precision :: RNA(MRA)
  double precision :: HNA(MHA)
  double precision :: PNB(MPB)
  double precision :: RNB(MRB)
  double precision :: HNB(MHB)
  double precision :: VELAP(MPA,MRA,MHA)
  double precision :: VELBP(MPB,MRB,MHB)
! discon
  double precision :: PN(51)
  double precision :: RRN(63)
  double precision :: DEPA(51,63)
  double precision :: DEPB(51,63)
  double precision :: DEPC(51,63)
! locate
  integer :: IPLOCA(MKA)
  integer :: IRLOCA(MKA)
  integer :: IHLOCA(MKA)
  integer :: IPLOCB(MKB)
  integer :: IRLOCB(MKB)
  integer :: IHLOCB(MKB)
  double precision :: PLA
  double precision :: RLA
  double precision :: HLA
  double precision :: PLB
  double precision :: RLB
  double precision :: HLB
! weight
  integer :: IP
  integer :: JP
  integer :: KP
  integer :: IP1
  integer :: JP1
  integer :: KP1
  double precision :: WV(8)
! prhfd
  double precision :: P
  double precision :: R
  double precision :: H
  double precision :: PF
  double precision :: RF
  double precision :: HF
  double precision :: PF1
  double precision :: RF1
  double precision :: HF1
  double precision :: PD
  double precision :: RD
  double precision :: HD
! jpmodv
  double precision :: VP(29)
  double precision :: VS(29)
  double precision :: RA(29)
  double precision :: DEPJ(29)
  end type jp3d_model_variables

  type (jp3d_model_variables) JP3DM_V
! jp3d_model_variables
      double precision :: VP1(29),VS1(29),RA1(29)
      integer :: L
      DATA VP1/7.75, 7.94, 8.13, 8.33, 8.54, 8.75, 8.97, &
               9.50, 9.91,10.26,10.55,10.99,11.29,11.50, &
              11.67,11.85,12.03,12.20,12.37,12.54,12.71, &
              12.87,13.02,13.16,13.32,13.46,13.60,13.64,13.64/
      DATA VS1/4.353,4.444,4.539,4.638,4.741,4.850,4.962, &
               5.227,5.463,5.670,5.850,6.125,6.295,6.395, &
               6.483,6.564,6.637,6.706,6.770,6.833,6.893, &
               6.953,7.012,7.074,7.137,7.199,7.258,7.314,7.304/
      DATA RA1/1.00,0.99,0.98,0.97,0.96,0.95,0.94,0.93, &
               0.92,0.91,0.90,0.88,0.86,0.84,0.82,0.80, &
               0.78,0.76,0.74,0.72,0.70,0.68,0.66,0.64, &
               0.62,0.60,0.58,0.56,0.55/
      DO 1 L  = 1,29
      JP3DM_V%VP(L)   = VP1(L)
      JP3DM_V%VS(L)   = VS1(L)
      JP3DM_V%RA(L)   = RA1(L)
      JP3DM_V%DEPJ(L) = 40.0+6325.59*(1.0-RA1(L))
1     CONTINUE
      RETURN
      END

      SUBROUTINE JPMODEL(IPS,H,V,JP3DM_V)
  implicit none

  include "constants.h"
! jp3d_model_variables
  type jp3d_model_variables
    sequence
! vmod3d
  integer :: NPA
  integer :: NRA
  integer :: NHA
  integer :: NPB
  integer :: NRB
  integer :: NHB
  double precision :: PNA(MPA)
  double precision :: RNA(MRA)
  double precision :: HNA(MHA)
  double precision :: PNB(MPB)
  double precision :: RNB(MRB)
  double precision :: HNB(MHB)
  double precision :: VELAP(MPA,MRA,MHA)
  double precision :: VELBP(MPB,MRB,MHB)
! discon
  double precision :: PN(51)
  double precision :: RRN(63)
  double precision :: DEPA(51,63)
  double precision :: DEPB(51,63)
  double precision :: DEPC(51,63)
! locate
  integer :: IPLOCA(MKA)
  integer :: IRLOCA(MKA)
  integer :: IHLOCA(MKA)
  integer :: IPLOCB(MKB)
  integer :: IRLOCB(MKB)
  integer :: IHLOCB(MKB)
  double precision :: PLA
  double precision :: RLA
  double precision :: HLA
  double precision :: PLB
  double precision :: RLB
  double precision :: HLB
! weight
  integer :: IP
  integer :: JP
  integer :: KP
  integer :: IP1
  integer :: JP1
  integer :: KP1
  double precision :: WV(8)
! prhfd
  double precision :: P
  double precision :: R
  double precision :: H
  double precision :: PF
  double precision :: RF
  double precision :: HF
  double precision :: PF1
  double precision :: RF1
  double precision :: HF1
  double precision :: PD
  double precision :: RD
  double precision :: HD
! jpmodv
  double precision :: VP(29)
  double precision :: VS(29)
  double precision :: RA(29)
  double precision :: DEPJ(29)
  end type jp3d_model_variables

  type (jp3d_model_variables) JP3DM_V
! jp3d_model_variables
      integer :: IPS,K,K1
      double precision :: H1,H2,H12,H,V
      DO 2 K = 1,28
      K1     = K+1
      H1     = JP3DM_V%DEPJ(K)
      H2     = JP3DM_V%DEPJ(K1)
      IF(H.GE.H1.AND.H.LT.H2) GO TO 3
2     CONTINUE
3     CONTINUE
      H12    = (H-H1)/(H2-H1)
      IF(IPS.EQ.1)  THEN
         V   = (JP3DM_V%VP(K1)-JP3DM_V%VP(K))*H12+JP3DM_V%VP(K)
      ELSE
         V   = (JP3DM_V%VS(K1)-JP3DM_V%VS(K))*H12+JP3DM_V%VS(K)
      END IF
      RETURN
      END


