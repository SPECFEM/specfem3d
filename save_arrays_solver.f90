!=====================================================================
!
!          S p e c f e m 3 D  B a s i n  V e r s i o n  1 . 4
!          --------------------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!         (c) California Institute of Technology September 2006
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================

  subroutine save_arrays_solver(flag_sediments,not_fully_in_bedrock,rho_vp,rho_vs,prname,xixstore,xiystore,xizstore, &
            etaxstore,etaystore,etazstore, &
            gammaxstore,gammaystore,gammazstore,jacobianstore, &
            xstore,ystore,zstore,kappastore,mustore, &
            ANISOTROPY, &
            c11store,c12store,c13store,c14store,c15store,c16store, &
            c22store,c23store,c24store,c25store,c26store,c33store,c34store,c35store,c36store, &
            c44store,c45store,c46store,c55store,c56store,c66store, &
            ibool,idoubling,rmass,rmass_ocean_load,npointot_oceans, &
            ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top, &
            nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax, &
            normal_xmin,normal_xmax,normal_ymin,normal_ymax,normal_bottom,normal_top, &
            jacobian2D_xmin,jacobian2D_xmax,jacobian2D_ymin,jacobian2D_ymax, &
            jacobian2D_bottom,jacobian2D_top, &
            iMPIcut_xi,iMPIcut_eta,nspec,nglob, &
            NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP,OCEANS)

  implicit none

  include "constants.h"

  integer nspec,nglob
  integer NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP
  integer npointot_oceans

  logical OCEANS
  logical ANISOTROPY

! arrays with jacobian matrix
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec) :: &
    xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore, &
    gammaxstore,gammaystore,gammazstore,jacobianstore

! arrays with mesh parameters
  double precision xstore(NGLLX,NGLLY,NGLLZ,nspec)
  double precision ystore(NGLLX,NGLLY,NGLLZ,nspec)
  double precision zstore(NGLLX,NGLLY,NGLLZ,nspec)

  real(kind=CUSTOM_REAL) kappastore(NGLLX,NGLLY,NGLLZ,nspec)
  real(kind=CUSTOM_REAL) mustore(NGLLX,NGLLY,NGLLZ,nspec)

  real(kind=CUSTOM_REAL) c11store(NGLLX,NGLLY,NGLLZ,nspec)
  real(kind=CUSTOM_REAL) c12store(NGLLX,NGLLY,NGLLZ,nspec)
  real(kind=CUSTOM_REAL) c13store(NGLLX,NGLLY,NGLLZ,nspec)
  real(kind=CUSTOM_REAL) c14store(NGLLX,NGLLY,NGLLZ,nspec)
  real(kind=CUSTOM_REAL) c15store(NGLLX,NGLLY,NGLLZ,nspec)
  real(kind=CUSTOM_REAL) c16store(NGLLX,NGLLY,NGLLZ,nspec)
  real(kind=CUSTOM_REAL) c22store(NGLLX,NGLLY,NGLLZ,nspec)
  real(kind=CUSTOM_REAL) c23store(NGLLX,NGLLY,NGLLZ,nspec)
  real(kind=CUSTOM_REAL) c24store(NGLLX,NGLLY,NGLLZ,nspec)
  real(kind=CUSTOM_REAL) c25store(NGLLX,NGLLY,NGLLZ,nspec)
  real(kind=CUSTOM_REAL) c26store(NGLLX,NGLLY,NGLLZ,nspec)
  real(kind=CUSTOM_REAL) c33store(NGLLX,NGLLY,NGLLZ,nspec)
  real(kind=CUSTOM_REAL) c34store(NGLLX,NGLLY,NGLLZ,nspec)
  real(kind=CUSTOM_REAL) c35store(NGLLX,NGLLY,NGLLZ,nspec)
  real(kind=CUSTOM_REAL) c36store(NGLLX,NGLLY,NGLLZ,nspec)
  real(kind=CUSTOM_REAL) c44store(NGLLX,NGLLY,NGLLZ,nspec)
  real(kind=CUSTOM_REAL) c45store(NGLLX,NGLLY,NGLLZ,nspec)
  real(kind=CUSTOM_REAL) c46store(NGLLX,NGLLY,NGLLZ,nspec)
  real(kind=CUSTOM_REAL) c55store(NGLLX,NGLLY,NGLLZ,nspec)
  real(kind=CUSTOM_REAL) c56store(NGLLX,NGLLY,NGLLZ,nspec)
  real(kind=CUSTOM_REAL) c66store(NGLLX,NGLLY,NGLLZ,nspec)

! Stacey
  real(kind=CUSTOM_REAL) rho_vp(NGLLX,NGLLY,NGLLZ,nspec)
  real(kind=CUSTOM_REAL) rho_vs(NGLLX,NGLLY,NGLLZ,nspec)

! flag indicating whether point is in the sediments
  logical flag_sediments(NGLLX,NGLLY,NGLLZ,nspec)
  logical not_fully_in_bedrock(nspec)

  integer ibool(NGLLX,NGLLY,NGLLZ,nspec)

! doubling mesh flag
  integer idoubling(nspec)

! mass matrix
  real(kind=CUSTOM_REAL) rmass(nglob)

! additional ocean load mass matrix
  real(kind=CUSTOM_REAL) rmass_ocean_load(npointot_oceans)

! boundary parameters locator
  integer ibelm_xmin(NSPEC2DMAX_XMIN_XMAX),ibelm_xmax(NSPEC2DMAX_XMIN_XMAX)
  integer ibelm_ymin(NSPEC2DMAX_YMIN_YMAX),ibelm_ymax(NSPEC2DMAX_YMIN_YMAX)
  integer ibelm_bottom(NSPEC2D_BOTTOM),ibelm_top(NSPEC2D_TOP)

! normals
  real(kind=CUSTOM_REAL) normal_xmin(NDIM,NGLLY,NGLLZ,NSPEC2DMAX_XMIN_XMAX)
  real(kind=CUSTOM_REAL) normal_xmax(NDIM,NGLLY,NGLLZ,NSPEC2DMAX_XMIN_XMAX)
  real(kind=CUSTOM_REAL) normal_ymin(NDIM,NGLLX,NGLLZ,NSPEC2DMAX_YMIN_YMAX)
  real(kind=CUSTOM_REAL) normal_ymax(NDIM,NGLLX,NGLLZ,NSPEC2DMAX_YMIN_YMAX)
  real(kind=CUSTOM_REAL) normal_bottom(NDIM,NGLLX,NGLLY,NSPEC2D_BOTTOM)
  real(kind=CUSTOM_REAL) normal_top(NDIM,NGLLX,NGLLY,NSPEC2D_TOP)

! jacobian on 2D edges
  real(kind=CUSTOM_REAL) jacobian2D_xmin(NGLLY,NGLLZ,NSPEC2DMAX_XMIN_XMAX)
  real(kind=CUSTOM_REAL) jacobian2D_xmax(NGLLY,NGLLZ,NSPEC2DMAX_XMIN_XMAX)
  real(kind=CUSTOM_REAL) jacobian2D_ymin(NGLLX,NGLLZ,NSPEC2DMAX_YMIN_YMAX)
  real(kind=CUSTOM_REAL) jacobian2D_ymax(NGLLX,NGLLZ,NSPEC2DMAX_YMIN_YMAX)
  real(kind=CUSTOM_REAL) jacobian2D_bottom(NGLLX,NGLLY,NSPEC2D_BOTTOM)
  real(kind=CUSTOM_REAL) jacobian2D_top(NGLLX,NGLLY,NSPEC2D_TOP)

! number of elements on the boundaries
  integer nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax

! MPI cut-planes parameters along xi and along eta
  logical iMPIcut_xi(2,nspec),iMPIcut_eta(2,nspec)

  integer i,j,k,ispec,iglob

! processor identification
  character(len=150) prname

! xix
  open(unit=27,file=prname(1:len_trim(prname))//'xix.bin',status='unknown',form='unformatted')
  write(27) xixstore
  close(27)

! xiy
  open(unit=27,file=prname(1:len_trim(prname))//'xiy.bin',status='unknown',form='unformatted')
  write(27) xiystore
  close(27)

! xiz
  open(unit=27,file=prname(1:len_trim(prname))//'xiz.bin',status='unknown',form='unformatted')
  write(27) xizstore
  close(27)

! etax
  open(unit=27,file=prname(1:len_trim(prname))//'etax.bin',status='unknown',form='unformatted')
  write(27) etaxstore
  close(27)

! etay
  open(unit=27,file=prname(1:len_trim(prname))//'etay.bin',status='unknown',form='unformatted')
  write(27) etaystore
  close(27)

! etaz
  open(unit=27,file=prname(1:len_trim(prname))//'etaz.bin',status='unknown',form='unformatted')
  write(27) etazstore
  close(27)

! gammax
  open(unit=27,file=prname(1:len_trim(prname))//'gammax.bin',status='unknown',form='unformatted')
  write(27) gammaxstore
  close(27)

! gammay
  open(unit=27,file=prname(1:len_trim(prname))//'gammay.bin',status='unknown',form='unformatted')
  write(27) gammaystore
  close(27)

! gammaz
  open(unit=27,file=prname(1:len_trim(prname))//'gammaz.bin',status='unknown',form='unformatted')
  write(27) gammazstore
  close(27)

! jacobian
  open(unit=27,file=prname(1:len_trim(prname))//'jacobian.bin',status='unknown',form='unformatted')
  write(27) jacobianstore
  close(27)

! flag_sediments
  open(unit=27,file=prname(1:len_trim(prname))//'flag_sediments.bin',status='unknown',form='unformatted')
  write(27) flag_sediments
  close(27)

! not_fully_in_bedrock
  open(unit=27,file=prname(1:len_trim(prname))//'not_fully_in_bedrock.bin',status='unknown',form='unformatted')
  write(27) not_fully_in_bedrock
  close(27)

! rho_vs
! Stacey
! rho_vp
  open(unit=27,file=prname(1:len_trim(prname))//'rho_vp.bin',status='unknown',form='unformatted')
  write(27) rho_vp
  close(27)

! rho_vs
  open(unit=27,file=prname(1:len_trim(prname))//'rho_vs.bin',status='unknown',form='unformatted')
  write(27) rho_vs
  close(27)

!!$! vp (for checking the mesh and model)
!!$  open(unit=27,file=prname(1:len_trim(prname))//'vp.bin',status='unknown',form='unformatted')
!!$  write(27) (FOUR_THIRDS * mustore + kappastore) / rho_vp
!!$  close(27)
!!$
!!$! vs (for checking the mesh and model)
!!$  open(unit=27,file=prname(1:len_trim(prname))//'vs.bin',status='unknown',form='unformatted')
!!$  write(27) mustore / rho_vs
!!$  close(27)

! kappa
  open(unit=27,file=prname(1:len_trim(prname))//'kappa.bin',status='unknown',form='unformatted')
  write(27) kappastore
  close(27)

! mu
  open(unit=27,file=prname(1:len_trim(prname))//'mu.bin',status='unknown',form='unformatted')
  write(27) mustore
  close(27)

! ibool
  open(unit=27,file=prname(1:len_trim(prname))//'ibool.bin',status='unknown',form='unformatted')
  write(27) ibool
  close(27)

! doubling
  open(unit=27,file=prname(1:len_trim(prname))//'idoubling.bin',status='unknown',form='unformatted')
  write(27) idoubling
  close(27)

! mass matrix
  open(unit=27,file=prname(1:len_trim(prname))//'rmass.bin',status='unknown',form='unformatted')
  write(27) rmass
  close(27)

! For anisotropy
  if(ANISOTROPY) then
     ! c11
     open(unit=27,file=prname(1:len_trim(prname))//'c11.bin',status='unknown',form='unformatted')
     write(27) c11store
     close(27)

     ! c12
     open(unit=27,file=prname(1:len_trim(prname))//'c12.bin',status='unknown',form='unformatted')
     write(27) c12store
     close(27)

     ! c13
     open(unit=27,file=prname(1:len_trim(prname))//'c13.bin',status='unknown',form='unformatted')
     write(27) c13store
     close(27)

     ! c14
     open(unit=27,file=prname(1:len_trim(prname))//'c14.bin',status='unknown',form='unformatted')
     write(27) c14store
     close(27)

     ! c15
     open(unit=27,file=prname(1:len_trim(prname))//'c15.bin',status='unknown',form='unformatted')
     write(27) c15store
     close(27)

     ! c16
     open(unit=27,file=prname(1:len_trim(prname))//'c16.bin',status='unknown',form='unformatted')
     write(27) c16store
     close(27)

     ! c22
     open(unit=27,file=prname(1:len_trim(prname))//'c22.bin',status='unknown',form='unformatted')
     write(27) c22store
     close(27)

     ! c23
     open(unit=27,file=prname(1:len_trim(prname))//'c23.bin',status='unknown',form='unformatted')
     write(27) c23store
     close(27)

     ! c24
     open(unit=27,file=prname(1:len_trim(prname))//'c24.bin',status='unknown',form='unformatted')
     write(27) c24store
     close(27)

     ! c25
     open(unit=27,file=prname(1:len_trim(prname))//'c25.bin',status='unknown',form='unformatted')
     write(27) c25store
     close(27)

     ! c26
     open(unit=27,file=prname(1:len_trim(prname))//'c26.bin',status='unknown',form='unformatted')
     write(27) c26store
     close(27)

     ! c33
     open(unit=27,file=prname(1:len_trim(prname))//'c33.bin',status='unknown',form='unformatted')
     write(27) c33store
     close(27)

     ! c34
     open(unit=27,file=prname(1:len_trim(prname))//'c34.bin',status='unknown',form='unformatted')
     write(27) c34store
     close(27)

     ! c35
     open(unit=27,file=prname(1:len_trim(prname))//'c35.bin',status='unknown',form='unformatted')
     write(27) c35store
     close(27)

     ! c36
     open(unit=27,file=prname(1:len_trim(prname))//'c36.bin',status='unknown',form='unformatted')
     write(27) c36store
     close(27)

     ! c44
     open(unit=27,file=prname(1:len_trim(prname))//'c44.bin',status='unknown',form='unformatted')
     write(27) c44store
     close(27)

     ! c45
     open(unit=27,file=prname(1:len_trim(prname))//'c45.bin',status='unknown',form='unformatted')
     write(27) c45store
     close(27)

     ! c46
     open(unit=27,file=prname(1:len_trim(prname))//'c46.bin',status='unknown',form='unformatted')
     write(27) c46store
     close(27)

     ! c55
     open(unit=27,file=prname(1:len_trim(prname))//'c55.bin',status='unknown',form='unformatted')
     write(27) c55store
     close(27)

     ! c56
     open(unit=27,file=prname(1:len_trim(prname))//'c56.bin',status='unknown',form='unformatted')
     write(27) c56store
     close(27)

     ! c66
     open(unit=27,file=prname(1:len_trim(prname))//'c66.bin',status='unknown',form='unformatted')
     write(27) c66store
     close(27)

  endif

! additional ocean load mass matrix if oceans
  if(OCEANS) then
    open(unit=27,file=prname(1:len_trim(prname))//'rmass_ocean_load.bin',status='unknown',form='unformatted')
    write(27) rmass_ocean_load
    close(27)
  endif

! boundary parameters
  open(unit=27,file=prname(1:len_trim(prname))//'ibelm.bin',status='unknown',form='unformatted')
  write(27) ibelm_xmin
  write(27) ibelm_xmax
  write(27) ibelm_ymin
  write(27) ibelm_ymax
  write(27) ibelm_bottom
  write(27) ibelm_top
  close(27)

  open(unit=27,file=prname(1:len_trim(prname))//'normal.bin',status='unknown',form='unformatted')
  write(27) normal_xmin
  write(27) normal_xmax
  write(27) normal_ymin
  write(27) normal_ymax
  write(27) normal_bottom
  write(27) normal_top
  close(27)

  open(unit=27,file=prname(1:len_trim(prname))//'jacobian2D.bin',status='unknown',form='unformatted')
  write(27) jacobian2D_xmin
  write(27) jacobian2D_xmax
  write(27) jacobian2D_ymin
  write(27) jacobian2D_ymax
  write(27) jacobian2D_bottom
  write(27) jacobian2D_top
  close(27)

  open(unit=27,file=prname(1:len_trim(prname))//'nspec2D.bin',status='unknown',form='unformatted')
  write(27) nspec2D_xmin
  write(27) nspec2D_xmax
  write(27) nspec2D_ymin
  write(27) nspec2D_ymax
  close(27)

! MPI cut-planes parameters along xi and along eta
  open(unit=27,file=prname(1:len_trim(prname))//'iMPIcut_xi.bin',status='unknown',form='unformatted')
  write(27) iMPIcut_xi
  close(27)

  open(unit=27,file=prname(1:len_trim(prname))//'iMPIcut_eta.bin',status='unknown',form='unformatted')
  write(27) iMPIcut_eta
  close(27)

! mesh arrays used in the solver to locate source and receivers
! use rmass for temporary storage to perform conversion, since already saved

!--- x coordinate
  rmass(:) = 0._CUSTOM_REAL
  do ispec = 1,nspec
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          iglob = ibool(i,j,k,ispec)
! distinguish between single and double precision for reals
          if(CUSTOM_REAL == SIZE_REAL) then
            rmass(iglob) = sngl(xstore(i,j,k,ispec))
          else
            rmass(iglob) = xstore(i,j,k,ispec)
          endif
        enddo
      enddo
    enddo
  enddo
  open(unit=27,file=prname(1:len_trim(prname))//'x.bin',status='unknown',form='unformatted')
  write(27) rmass
  close(27)

!--- y coordinate
  rmass(:) = 0._CUSTOM_REAL
  do ispec = 1,nspec
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          iglob = ibool(i,j,k,ispec)
! distinguish between single and double precision for reals
          if(CUSTOM_REAL == SIZE_REAL) then
            rmass(iglob) = sngl(ystore(i,j,k,ispec))
          else
            rmass(iglob) = ystore(i,j,k,ispec)
          endif
        enddo
      enddo
    enddo
  enddo
  open(unit=27,file=prname(1:len_trim(prname))//'y.bin',status='unknown',form='unformatted')
  write(27) rmass
  close(27)

!--- z coordinate
  rmass(:) = 0._CUSTOM_REAL
  do ispec = 1,nspec
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          iglob = ibool(i,j,k,ispec)
! distinguish between single and double precision for reals
          if(CUSTOM_REAL == SIZE_REAL) then
            rmass(iglob) = sngl(zstore(i,j,k,ispec))
          else
            rmass(iglob) = zstore(i,j,k,ispec)
          endif
        enddo
      enddo
    enddo
  enddo
  open(unit=27,file=prname(1:len_trim(prname))//'z.bin',status='unknown',form='unformatted')
  write(27) rmass
  close(27)

  end subroutine save_arrays_solver

