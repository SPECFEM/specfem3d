!
!    Copyright 2013, Tarje Nissen-Meyer, Alexandre Fournier, Martin van Driel
!                    Simon Stahler, Kasra Hosseini, Stefanie Hempel
!
!    This file is part of AxiSEM.
!    It is distributed from the webpage < http://www.axisem.info>
!
!    AxiSEM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    AxiSEM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with AxiSEM.  If not, see < http://www.gnu.org/licenses/>.
!

!=========================================================================================
module meshgen

  use global_parameters, only: sp, dp
  implicit none

  public :: generate_skeleton
  private

  real(kind=dp)   , dimension(:), allocatable :: aspect_ratio

  !REGION BY REGION PARAMETERS
  integer, dimension(:,:), allocatable          :: lnodeso   ! OUTER SHELL
  character(len=6), dimension(:), allocatable   :: eltypeo
  logical, dimension(:), allocatable            :: coarsingo

  integer                                       :: nelo
  real(kind=dp)   , dimension(:), allocatable   :: so,zo

  integer, dimension(:,:), allocatable          :: lnodesi   ! INNER SHELL
  character(len=6), dimension(:), allocatable   :: eltypei
  logical, dimension(:), allocatable            :: coarsingi
  integer                                       :: neli
  real(kind=dp)   , dimension(:), allocatable   :: si,zi

  integer, dimension(:,:), allocatable          :: lnodessq  ! INNER SQUARE
  character(len=6), dimension(:), allocatable   :: eltypesq
  integer                                       :: nelsq
  real(kind=dp)   , dimension(:), allocatable   :: ssq,zsq

  integer, dimension(:,:), allocatable          :: lnodesbuf ! BUFFER SHELL
  character(len=6), dimension(:), allocatable   :: eltypebuf
  integer                                       :: nelbuf
  real(kind=dp)   , dimension(:), allocatable   :: sbuf,zbuf

contains

!-----------------------------------------------------------------------------------------
subroutine generate_skeleton

  use data_grid
  use data_diag
  use data_mesh
  use data_numbering
  use numbering

  ! OUTER SHELL GENERATION
  !  define reference grid
  if (dump_mesh_info_screen) then
    write(*,*)'generating reference spherical grid....';call flush(6)
  endif
  call def_reference_spherical_grid_discont

  ! Aspect ratio
  write(*,*)'estimate aspect ratio...';call flush(6)
  call estimate_aspect_ratio

  ! final outer shell w/ potential coarsening w/ depth
  write(*,*)'define spherical shell...';call flush(6)
  call define_spherical_shell
  deallocate(crd_grds,crd_grdc)
  deallocate(z_unif,s_unif)

  ! CENTRAL REGION GENERATION (inner shell + central square + buffer layer)
  write(*,*)'define central region...';call flush(6)
  call define_central_region

  ! GATHER INFORMATION FROM DIFFERENT REGIONS IN GLOBAL ARRAYS
  write(*,*)'gather skeleton...';call flush(6)
  call gather_skeleton

  ! must add a flag there
  if (southern) then
   write(*,*)'generate southern hemisphere...';call flush(6)
   call generate_southern_hemisphere
  else
   call donot_generate_southern_hemisphere
  endif

  write(*,*)
  write(*,"(10x,' THE TOTAL NUMBER OF ELEMENTS IS ',i10 )") neltot
  write(*,*)

  call generate_serendipity(npointot, neltot, sg, zg)

end subroutine generate_skeleton
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine generate_serendipity(npoin, nel, sg, zg)
! 08/09/2004: generates mesh database compatible w/ former version of mesher
!            (element topology defined by 8 control points instead of 4)

  use numbering
  use data_grid
  use data_diag

  integer, intent(in)                           :: npoin, nel
  real(kind=dp)   , intent(in)                  :: sg(4,nel),zg(4,nel)
  real(kind=dp)   , dimension(:,:), allocatable :: sg2,zg2
  integer, dimension(:,:), allocatable          :: lnods
  integer                                       :: iel, inode, ipt,npoin2

  ! Numbering arrays
  integer, dimension(:), allocatable :: ipt_iglob,el_iglob,inode_iglob
  integer                            :: nglob_serend
  integer, dimension(:), allocatable :: iglob_serend,loc_serend
  logical, dimension(:), allocatable :: ifseg_serend

  npoin2 = 8*npoin/4

  if (dump_mesh_info_screen) write(*,*) ' NPOIN 2 is ', npoin2,nel*8
  allocate(sg2(8,nel),zg2(8,nel))

  ! The even coordinate indices below are linearly interpolated, i.e. DO NOT
  ! represent the correct location for any spheroidally shaped element,
  ! only for the linear shapes at the center (the only case for which
  ! these serendipity nodes are actually needed)!

  do iel = 1, nel
     sg2(1,iel) = sg(1,iel)
     zg2(1,iel) = zg(1,iel)
     sg2(3,iel) = sg(2,iel)
     zg2(3,iel) = zg(2,iel)
     sg2(5,iel) = sg(3,iel)
     zg2(5,iel) = zg(3,iel)
     sg2(7,iel) = sg(4,iel)
     zg2(7,iel) = zg(4,iel)

     sg2(2,iel) = .5d0 * ( sg2(1,iel) + sg2(3,iel) )
     zg2(2,iel) = .5d0 * ( zg2(1,iel) + zg2(3,iel) )
     sg2(4,iel) = .5d0 * ( sg2(3,iel) + sg2(5,iel) )
     zg2(4,iel) = .5d0 * ( zg2(3,iel) + zg2(5,iel) )
     sg2(6,iel) = .5d0 * ( sg2(5,iel) + sg2(7,iel) )
     zg2(6,iel) = .5d0 * ( zg2(5,iel) + zg2(7,iel) )
     sg2(8,iel) = .5d0 * ( sg2(7,iel) + sg2(1,iel) )
     zg2(8,iel) = .5d0 * ( zg2(7,iel) + zg2(1,iel) )


     ! TNM: make sure axial points are equal to zero
     if (sg2(1,iel) <= 0.1d0*dabs(sg2(2,iel)-sg2(1,iel))) sg2(1,iel) = 0.d0

     if (sg2(7,iel) <= 0.1d0*dabs(sg2(6,iel)-sg2(7,iel))) sg2(7,iel) = 0.d0

     if (sg2(8,iel) <= 0.1d0*dabs(sg2(6,iel)-sg2(7,iel))) sg2(8,iel) = 0.d0
  enddo

  ! write out meshes: entire domain, central region, crust, coarsening level
  if (dump_mesh_info_files) then
    call write_serendipity_meshes(nel,sg2,zg2)
  endif

  allocate(iglob_serend(npoin2))
  iglob_serend(:) = 0
  allocate(loc_serend(npoin2))
  loc_serend(:) = 0
  allocate(ifseg_serend(npoin2))

  if (dump_mesh_info_screen) write(*,*) 'CALLING GLOBAL NUMBERING'

  call get_global(nel, sg2, zg2, iglob_serend, loc_serend, ifseg_serend, &
                  nglob_serend, npoin2, 8)
  if (dump_mesh_info_screen) write(*,*) 'NGLOB SERENDIPITY IS ' , nglob_serend


  allocate (ipt_iglob(nglob_serend),el_iglob(nglob_serend),inode_iglob(nglob_serend))

  do iel = 1, nel
     do inode = 1, 8
        ipt = (iel-1)*8 + inode
        ipt_iglob(iglob_serend(ipt))   = ipt
        el_iglob(iglob_serend(ipt))    = iel
        inode_iglob(iglob_serend(ipt)) = inode
     enddo
  enddo

  allocate(lnods(8,nel))
  do iel = 1, nel
     do inode = 1, 8
        ipt = (iel-1) * 8  + inode
        lnods(inode,iel) = iglob_serend(ipt)
     enddo
  enddo

  deallocate(inode_iglob,el_iglob)
  deallocate(ifseg_serend,loc_serend,iglob_serend)
  deallocate(sg2,zg2)

end subroutine generate_serendipity
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine write_serendipity_meshes(nel,sg2,zg2)
  use data_bkgrdmodel
  use data_diag
  use data_grid, only: ri,router

  integer, intent(in)           :: nel
  real(kind=dp)   , intent(in)  :: sg2(8,nel),zg2(8,nel)
  integer                       :: iel

  write(*,*)'writing all elements....'
  open(unit=157,file=diagpath(1:lfdiag)//'/global_skel.dat')
  do iel = 1, nel
     write(157,*) sg2(1,iel), zg2(1,iel)
     write(157,*) sg2(2,iel), zg2(2,iel)
     write(157,*) sg2(3,iel), zg2(3,iel)
     write(157,*) sg2(4,iel), zg2(4,iel)
     write(157,*) sg2(5,iel), zg2(5,iel)
     write(157,*) sg2(6,iel), zg2(6,iel)
     write(157,*) sg2(7,iel), zg2(7,iel)
     write(157,*) sg2(8,iel), zg2(8,iel)
     write(157,*) sg2(1,iel), zg2(1,iel)
     write(157,*)
  enddo
  close(157)

  write(*,*)'writing regions of elements...'; call flush(6)
  open(unit=2157,file=diagpath(1:lfdiag)//'/foc_skel.dat')
  open(unit=3157,file=diagpath(1:lfdiag)//'/smcic_skel.dat')
  open(unit=1157,file=diagpath(1:lfdiag)//'/center_skel.dat')
  open(unit=1257,file=diagpath(1:lfdiag)//'/um_antipode.dat')
  open(unit=1357,file=diagpath(1:lfdiag)//'/midmantle.dat')
  open(unit=1457,file=diagpath(1:lfdiag)//'/uppermantle.dat')
  open(unit=1557,file=diagpath(1:lfdiag)//'/uppermantle_north.dat')

  do iel = 1, nel

     ! fluid core
     if ( sqrt(sg2(5,iel)**2+zg2(5,iel)**2) <= ri) then
        write(2157,*) sg2(1,iel), zg2(1,iel)
        write(2157,*) sg2(2,iel), zg2(2,iel)
        write(2157,*) sg2(3,iel), zg2(3,iel)
        write(2157,*) sg2(4,iel), zg2(4,iel)
        write(2157,*) sg2(5,iel), zg2(5,iel)
        write(2157,*) sg2(6,iel), zg2(6,iel)
        write(2157,*) sg2(7,iel), zg2(7,iel)
        write(2157,*) sg2(8,iel), zg2(8,iel)
        write(2157,*) sg2(1,iel), zg2(1,iel)
        write(2157,*)
     endif

     ! solid regions: mantle & inner core
     if ( sqrt(sg2(5,iel)**2+zg2(5,iel)**2) >= ri) then
         write(3157,*) sg2(1,iel), zg2(1,iel)
         write(3157,*) sg2(2,iel), zg2(2,iel)
         write(3157,*) sg2(3,iel), zg2(3,iel)
         write(3157,*) sg2(4,iel), zg2(4,iel)
         write(3157,*) sg2(5,iel), zg2(5,iel)
         write(3157,*) sg2(6,iel), zg2(6,iel)
         write(3157,*) sg2(7,iel), zg2(7,iel)
         write(3157,*) sg2(8,iel), zg2(8,iel)
         write(3157,*) sg2(1,iel), zg2(1,iel)
         write(3157,*)
     endif

     ! write only IC section
     if ( sqrt(sg2(5,iel)**2+zg2(5,iel)**2) <= discont(ndisc)/router) then
        write(1157,*) sg2(1,iel), zg2(1,iel)
        write(1157,*) sg2(2,iel), zg2(2,iel)
        write(1157,*) sg2(3,iel), zg2(3,iel)
        write(1157,*) sg2(4,iel), zg2(4,iel)
        write(1157,*) sg2(5,iel), zg2(5,iel)
        write(1157,*) sg2(6,iel), zg2(6,iel)
        write(1157,*) sg2(7,iel), zg2(7,iel)
        write(1157,*) sg2(8,iel), zg2(8,iel)
        write(1157,*) sg2(1,iel), zg2(1,iel)
        write(1157,*)
     endif

     ! plot only upper mantle near antipode
     if ( sqrt(sg2(5,iel)**2+zg2(5,iel)**2) >= 0.9d0 .and. &
         sg2(5,iel) <= 0.12 .and. zg2(5,iel) < 0.d0 ) then
        write(1257,*) sg2(1,iel), zg2(1,iel)
        write(1257,*) sg2(2,iel), zg2(2,iel)
        write(1257,*) sg2(3,iel), zg2(3,iel)
        write(1257,*) sg2(4,iel), zg2(4,iel)
        write(1257,*) sg2(5,iel), zg2(5,iel)
        write(1257,*) sg2(6,iel), zg2(6,iel)
        write(1257,*) sg2(7,iel), zg2(7,iel)
        write(1257,*) sg2(8,iel), zg2(8,iel)
        write(1257,*) sg2(1,iel), zg2(1,iel)
        write(1257,*)
     endif

     ! plot only upper mantle
     if ( sqrt(sg2(5,iel)**2+zg2(5,iel)**2) >= 0.9d0) then
        write(1457,*) sg2(1,iel), zg2(1,iel)
        write(1457,*) sg2(2,iel), zg2(2,iel)
        write(1457,*) sg2(3,iel), zg2(3,iel)
        write(1457,*) sg2(4,iel), zg2(4,iel)
        write(1457,*) sg2(5,iel), zg2(5,iel)
        write(1457,*) sg2(6,iel), zg2(6,iel)
        write(1457,*) sg2(7,iel), zg2(7,iel)
        write(1457,*) sg2(8,iel), zg2(8,iel)
        write(1457,*) sg2(1,iel), zg2(1,iel)
        write(1457,*)
     endif

     ! plot only upper mantle in north
      if ( sqrt(sg2(5,iel)**2+zg2(5,iel)**2) >= 0.95d0 .and. &
          sg2(5,iel) <= 1000./router .and. zg2(5,iel) > 0.d0  ) then
       write(1557,*) sg2(1,iel), zg2(1,iel)
       write(1557,*) sg2(2,iel), zg2(2,iel)
       write(1557,*) sg2(3,iel), zg2(3,iel)
       write(1557,*) sg2(4,iel), zg2(4,iel)
       write(1557,*) sg2(5,iel), zg2(5,iel)
       write(1557,*) sg2(6,iel), zg2(6,iel)
       write(1557,*) sg2(7,iel), zg2(7,iel)
       write(1557,*) sg2(8,iel), zg2(8,iel)
       write(1557,*) sg2(1,iel), zg2(1,iel)
       write(1557,*)
      endif

   enddo
   close(1557)
   close(1457)
   close(1357)
   close(1257)
   close(1157)
   close(3157)
   close(2157)

end subroutine write_serendipity_meshes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine def_reference_spherical_grid_discont
! ALEX 08/04/2004
! We make the assumption of a uniform rectangular spacing
! in the associated cylindrical grid, which gets
! mapped into a spherical shell grid of inner
! radius ri and outer radius ro.

  use data_bkgrdmodel
  use data_grid
  use data_diag
  use analytic_spheroid_mapping

  real(kind=dp)   , dimension(8,2)              :: crd_control_nodes(8,2)
  real(kind=dp)   , dimension(:), allocatable   :: dz
  integer :: npts
  integer :: iz

  ! FIRST DEFINE PARAMETERS FOR ASSOCIATED CYLINDRICAL/Cartesian GRID

  ns = ns_glob ! inherited from discont_meshing routine
  nz = nz_glob ! inherited from discont_meshing routine
  ri = rmin/router
  ro = 1.
  allocate(dz(1:nz))
  do iz = 1, nz
    dz(iz) = 2.d0 * dz_glob(nz-iz+1) / (router-rmin)
  enddo

  if (dump_mesh_info_screen) then
     write(*,*)  'ns,nz,ri,ro:', ns,nz,ri,ro
     write(*,*)  'dz', dz(:)
     write(*,*)  'SUM(dz):', SUM(dz)
  endif

  ! Total number of points (soon to be corners)  in the non-coarsened shell
  npts = (ns + 1) * (nz + 1)

  if (dump_mesh_info_screen) then
     write(*,*)'CHECK PARAMS 4 ri,ro,router,nz,ns:',ri,ro,router,nz,ns
     call flush(6)
  endif

  allocate(s_unif (npts),z_unif (npts))
  s_unif (:) = 0.d0
  z_unif (:) = 0.d0

  allocate(crd_grdc(1:ns+1,1:nz+1,2),crd_grds(1:ns+1,1:nz+1,2))
  crd_grdc(:,:,:) = 0.d0
  crd_grds(:,:,:) = 0.d0

  ! Define coordinates in the parent square
  call def_ref_cart_coordinates_discont(ns, nz, crd_grdc, dz)

  ! In order to use analytic mapping to get coordinates
  ! 1) Define control nodes to define shell geometry (Northern H)
  call def_control_nodes(crd_control_nodes, ri, ro)

  ! 2) Use the map_spheroid function
  call def_mapped_coordinates(ns,nz,crd_grds,crd_grdc,crd_control_nodes)

  ! At this stage we know the coordinates in the physical domain of the
  ! 8 control points that define a spectral element which belongs to the
  ! spherical-shelly part of the mesh
  ! 3) Use lexicographical ordering to create 1d arrays of coordinates
  call def_global_coordinates(npts,s_unif,z_unif,ns,nz,crd_grds)
  ! should probably add a possible geometrical diag output at this stage
  ! ...
end subroutine def_reference_spherical_grid_discont
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine def_global_coordinates(npts, s, z, ns1, nz1, crd)
  integer, intent(in)                                     :: npts, ns1, nz1
  real(kind=dp), dimension(1:npts), intent(out)           :: s, z
  real(kind=dp), dimension(1:ns1+1,1:nz1+1,2), intent(in) :: crd

  integer :: ipt, is, iz

  do iz = 1,nz1+1
     do is = 1, ns1+1
        ipt = uniform_nodenumber(is,iz,ns1)
        s(ipt) = crd(is,iz,1)
        z(ipt) = crd(is,iz,2)
     enddo
  enddo

end subroutine def_global_coordinates
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine def_mapped_coordinates(ns1, nz1, crds, crdc, crd_cont)
  use analytic_spheroid_mapping

  integer, intent(in) :: ns1, nz1

  real(kind=dp)   , dimension(1:ns1+1,1:nz1+1,2), intent(out) :: crds
  real(kind=dp)   , dimension(1:ns1+1,1:nz1+1,2), intent(in)  :: crdc
  real(kind=dp)   , dimension(8,2), intent(in)                :: crd_cont

  integer           :: is, iz
  real(kind=dp)     :: xi, eta

  do iz = 1,nz1+1
     do is = 1,ns1+1
        xi = crdc(is,iz,1)
        eta = crdc(is,iz,2)
        crds(is,iz,1) = map_spheroid(xi, eta, crd_cont, 1)
        crds(is,iz,2) = map_spheroid(xi, eta, crd_cont, 2)
     enddo
  enddo

end subroutine def_mapped_coordinates
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine def_control_nodes(crd, ri1, ro1)
  real(kind=dp)   , intent(in) :: ri1, ro1
  real(kind=dp)   , dimension(8,2), intent(out) :: crd
  !hemispherical case

  crd(1,1) = 0.d0
  crd(1,2) = ri1
  crd(2,1) = ri1*.5*dsqrt(2.d0)
  crd(2,2) = ri1*.5d0*dsqrt(2.d0)
  crd(3,1) = ri1
  crd(3,2) = 0.d0
  crd(4,1) = .5d0*(ri1+ro1)
  crd(4,2) = 0.d0
  crd(5,1) = ro1
  crd(5,2) = 0.d0
  crd(6,1) = ro1*.5*dsqrt(2.d0)
  crd(6,2) = ro1*.5d0*dsqrt(2.d0)
  crd(7,1) = 0.d0
  crd(7,2) = ro1
  crd(8,1) = 0.d0
  crd(8,2) = .5d0*(ro1+ri1)

end subroutine def_control_nodes
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine estimate_aspect_ratio
  use data_grid
  use data_diag
  real(kind=dp)     :: s1, z1, s2, z2, hr
  integer           :: iz

  allocate(aspect_ratio(nz))
  aspect_ratio(:) = 0.

  if (dump_mesh_info_files) open(unit=16,file=diagpath(1:lfdiag)//'/fort.16')

  do iz=1,nz
    s1 = .5d0*(crd_grds(1,iz,1)+crd_grds(1,iz+1,1))
    s2 = .5d0*(crd_grds(2,iz,1)+crd_grds(2,iz+1,1))
    z1 = .5d0*(crd_grds(1,iz,2)+crd_grds(1,iz+1,2))
    z2 = .5d0*(crd_grds(2,iz,2)+crd_grds(2,iz+1,2))
    hr = dsqrt( (s2-s1)**2 + (z2-z1)**2 ) ! WIDTH IS COMPUTED
    s1 = .5d0*(crd_grds(1,iz  ,1)+crd_grds(2,iz  ,1))
    s2 = .5d0*(crd_grds(1,iz+1,1)+crd_grds(2,iz+1,1))
    z1 = .5d0*(crd_grds(1,iz  ,2)+crd_grds(2,iz  ,2))
    z2 = .5d0*(crd_grds(1,iz+1,2)+crd_grds(2,iz+1,2))
    hr = hr / dsqrt( (s2-s1)**2 + (z2-z1)**2 ) ! HR = WIDTH / HEIGHT
    if (dump_mesh_info_files) write(16,*) iz, hr
    aspect_ratio(iz) = hr
  enddo

  if (dump_mesh_info_files) close(16)

end subroutine estimate_aspect_ratio
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine define_spherical_shell

  use data_bkgrdmodel
  use data_grid
  use data_diag
  use data_coarse
  implicit none

  integer :: nel
  integer :: inode, iel, ipt,ipto, ic
  ! This will have to be parameterised later on

  nc=nc_glob
  allocate (iclev(0:nc+1))
  iclev(0)=nz_glob+1
  do ic = 1, nc
     iclev(ic)=iclev_glob(ic)
  enddo
  iclev(nc+1)=1

  ! COMPATIBILITY CONDITIONS HAVE TO BE DEFINED HERE
  ! IS NS a multiple of 2**nc ? (an even multiple!)
  ! IS NZ compatible with the coarsening, etc.
  ! compute new number of elements

  call compute_nelclean(nel, nc, iclev, ns, nz)
  ! define elements topology
  allocate(lnodeso(4,nel))
  lnodeso(:,:) = 0
  ! array to characterize element geometry
  allocate(eltypeo(nel))
  ! call define_lnodes

  allocate(coarsingo(nel))
  coarsingo = .false.

  call define_lnodesclean(nel, lnodeso, eltypeo, coarsingo, nc, iclev, ns, nz)

  nelo  = nel
  ! Fill so and zo arrays
  allocate(so(4*nelo),zo(4*nelo))
  so(:) = 0.
  zo(:) = 0.

  do iel = 1, nel
     do inode = 1, 4
        ipt = lnodeso(inode,iel)
        ipto = (iel-1)*4 + inode
        so(ipto) = s_unif (ipt)
        zo(ipto) = z_unif (ipt)
     enddo
  enddo

  ! gnuplot dump
  if (dump_mesh_info_files) then
     open(unit=3,file=diagpath(1:lfdiag)//'/testcrd.dat',STATUS="UNKNOWN",POSITION="REWIND")
     do iel = 1, nel
        do inode = 1, 4
           ipt = lnodeso(inode,iel)
           ipto = (iel-1)*4 + inode
           write(3,*) so(ipto),zo(ipto)
        enddo
          ipt = lnodeso(1,iel)
          ipto = (iel-1)*4 + 1
          write(3,*) so(ipto),zo(ipto)
          write(3,*)
     enddo
  endif

end subroutine define_spherical_shell
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
integer function uniform_nodenumber(is, iz, ns, istart)
  !uses lexicographical order

  integer, intent(in)           :: is, iz, ns
  integer, intent(in), optional :: istart
  integer                       :: i0

  i0 = 0
  if (PRESENT(istart)) i0 = istart
  uniform_nodenumber = (iz-1)*(ns+1) + is + i0

end function uniform_nodenumber
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine compute_nelclean(nel, nc1, iclev1, ns1, nz1)
  use data_diag
  !returns the number of elements defining the new spherical shell grid.
  integer, intent(out) :: nel
  integer, intent(in) :: nc1
  integer, dimension(0:nc1+1) :: iclev1
  integer, intent(in) :: ns1, nz1
  integer :: iz,ic,icc,icold,nelregion(nc1+1),nelabove,nelregionsum


  nel = 0
  nelregion(:)=0

  icold = 0

  if ( nc1 == 0 ) then
    nel = ns1*nz1
    return
  endif


  icc =1
  do iz = nz1, iclev1(icc),-1
     if ( iz > iclev1(icc) ) nel = nel + ns1
     if ( iz == iclev1(icc) ) nel = nel + 3*ns1/2
  enddo
  icold = 1
  nelregion(1)=nel
  nelabove=nel
  do icc = 2,nc1
     ic = 2**(icc-1)
     do iz = iclev1(icc-1)-2, iclev1(icc),-1
        if ( iz > iclev1(icc) ) nel = nel + ns1/ic ! TEST
        if ( iz == iclev1(icc) ) nel = nel + 3*ns1/(2*ic)
     enddo
     nelregion(icc)=nel-nelabove
     nelabove=nel
     icold = ic
  enddo

  icc = nc1+1
  ic = 2**(icc-1)
  do iz = iclev1(icc-1)-2, 1,-1
     nel = nel + ns1/ic ! TEST
  enddo

  nelregion(icc)=nel-nelabove
  if (dump_mesh_info_screen) write(*,*) 'nel =',nel, 'instead of', nz1*ns1
  nelregionsum=sum(nelregion)

  if (dump_mesh_info_screen) then
     write(*,*) 'SUM # EL. ALL SPHERICAL REGIONS: ',nelregionsum
     do icc=1,nc1+1
        write(*,*)'# EL. in REGION ',icc,' :',nelregion(icc),' (', &
                   real(nelregion(icc))/real(nelregionsum)*100.,' percent)'
     enddo
  endif

  if (dump_mesh_info_files) then
   open(unit=9999,file=diagpath(1:lfdiag)//'/fort.9999')
    do icc=1,nc1+1
     write(9999,*)icc,nelregion(icc)
    enddo
   close(9999)
  endif

end subroutine compute_nelclean
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine define_lnodesclean(nel, lnodes, eltype, coarsingloc, nc, iclev, ns, nz)

  ! The topology of the spherical shell and the array characterizing
  ! the geometry of each element is defined here.

  integer, intent(in)                           :: nel
  integer, dimension(4,nel),intent(out)         :: lnodes
  character(len=6), dimension(nel), intent(out) :: eltype
  logical, dimension(nel), intent(out)          :: coarsingloc
  integer, intent(in)                           :: nc

  integer, dimension(0:nc+1),intent(in)         :: iclev
  integer, intent(in)                           :: ns, nz
  integer                                       :: iz, ic, icc, iel, is, icold

  iel = 0
  ! no coarsening case
  if (nc == 0) then
     do iz = nz,1,-1
        do is = 1, ns
           iel = iel + 1
           lnodes(1,iel) = uniform_nodenumber(is   ,iz,ns)
           lnodes(2,iel) = uniform_nodenumber(is+1,iz,ns)
           lnodes(3,iel) = uniform_nodenumber(is+1,iz+1,ns)
           lnodes(4,iel) = uniform_nodenumber(is   ,iz+1,ns)
           eltype(iel) = 'curved'
        enddo
     enddo
     return
  endif

  ! Going from top to bottom
  icc = 1
  ic = 1
  do iz = nz, iclev(icc),-ic
     if ( iz > iclev(icc) ) then
        do is = 1, ns, ic
           iel = iel + 1
           lnodes(1,iel) = uniform_nodenumber(is   ,iz,ns)
           lnodes(2,iel) = uniform_nodenumber(is+ic,iz,ns)
           lnodes(3,iel) = uniform_nodenumber(is+ic,iz+ic,ns)
           lnodes(4,iel) = uniform_nodenumber(is   ,iz+ic,ns)
           eltype(iel) = 'curved'
        enddo
     else if ( iz == iclev(icc) ) then
        !This takes care of the "conformal mortars"
        !This case here corresponds to the first remeshing/coarsening
        !when going down.
        do is = 1, ns, ic*4

           iel = iel + 1
           lnodes(1,iel) = uniform_nodenumber(is   ,iz,ns   )
           lnodes(2,iel) = uniform_nodenumber(is+ic,iz,ns   )
           lnodes(3,iel) = uniform_nodenumber(is+ic,iz+1,ns)
           lnodes(4,iel) = uniform_nodenumber(is   ,iz+1,ns)
           eltype(iel) = 'curved'
           coarsingloc(iel)=.true.

           iel = iel + 1
           lnodes(1,iel) = uniform_nodenumber(is     ,iz-1,ns )
           lnodes(2,iel) = uniform_nodenumber(is+2*ic,iz-1,ns )
           lnodes(3,iel) = uniform_nodenumber(is+  ic,iz,ns   )
           lnodes(4,iel) = uniform_nodenumber(is     ,iz,ns   )
           eltype(iel) = 'curved'
           coarsingloc(iel)=.true.

           iel = iel + 1
           lnodes(1,iel) = uniform_nodenumber(is  +ic,iz,ns   )
           lnodes(2,iel) = uniform_nodenumber(is+2*ic,iz-1,ns )
           lnodes(3,iel) = uniform_nodenumber(is+2*ic,iz+1,ns )
           lnodes(4,iel) = uniform_nodenumber(is  +ic,iz+1,ns )
           eltype(iel) = 'semino'
           coarsingloc(iel)=.true.

           iel = iel + 1
           lnodes(1,iel) = uniform_nodenumber(is+2*ic,iz- 1,ns)
           lnodes(2,iel) = uniform_nodenumber(is+3*ic,iz,ns   )
           lnodes(3,iel) = uniform_nodenumber(is+3*ic,iz+ 1,ns)
           lnodes(4,iel) = uniform_nodenumber(is+2*ic,iz+ 1,ns)
           eltype(iel) = 'semino'
           coarsingloc(iel)=.true.

           iel = iel + 1
           lnodes(1,iel) = uniform_nodenumber(is+2*ic,iz-1,ns )
           lnodes(2,iel) = uniform_nodenumber(is+4*ic,iz-1,ns )
           lnodes(3,iel) = uniform_nodenumber(is+4*ic,iz ,ns  )
           lnodes(4,iel) = uniform_nodenumber(is+3*ic,iz ,ns  )
           eltype(iel) = 'curved'
           coarsingloc(iel)=.true.

           iel = iel + 1
           lnodes(1,iel) = uniform_nodenumber(is+3*ic,iz   ,ns)
           lnodes(2,iel) = uniform_nodenumber(is+4*ic,iz   ,ns)
           lnodes(3,iel) = uniform_nodenumber(is+4*ic,iz+1 ,ns)
           lnodes(4,iel) = uniform_nodenumber(is+3*ic,iz+1 ,ns)
           eltype(iel) = 'curved'
           coarsingloc(iel)=.true.

        enddo
     endif
  enddo

  ! intermediate domains, bound by coarsening levels above & below
  icold = 1
  do icc = 2, nc
     ic = 2**(icc-1)
       do iz = iclev(icc-1)-2, iclev(icc),-1
          if ( iz > iclev(icc) ) then
             do is = 1, ns, ic
                 iel = iel + 1
                 lnodes(1,iel) = uniform_nodenumber(is   ,iz,ns)
                 lnodes(2,iel) = uniform_nodenumber(is+ic,iz,ns)
                 lnodes(3,iel) = uniform_nodenumber(is+ic,iz+1 ,ns)
                 lnodes(4,iel) = uniform_nodenumber(is   ,iz+1 ,ns)
                 eltype(iel) = 'curved'
             enddo
          else if ( iz == iclev(icc) ) then
             ! Remeshing
             do is = 1, ns, ic*4

                iel = iel + 1
                lnodes(1,iel) = uniform_nodenumber(is   ,iz   ,ns)
                lnodes(2,iel) = uniform_nodenumber(is+ic,iz   ,ns)
                lnodes(3,iel) = uniform_nodenumber(is+ic,iz+1 ,ns)
                lnodes(4,iel) = uniform_nodenumber(is   ,iz+1 ,ns)
                eltype(iel) = 'curved'
                coarsingloc(iel)=.true.

                iel = iel + 1
                lnodes(1,iel) = uniform_nodenumber(is     ,iz-1 ,ns)
                lnodes(2,iel) = uniform_nodenumber(is+2*ic,iz-1 ,ns)
                lnodes(3,iel) = uniform_nodenumber(is+  ic,iz   ,ns)
                lnodes(4,iel) = uniform_nodenumber(is     ,iz   ,ns)
                eltype(iel) = 'curved'
                coarsingloc(iel)=.true.

                iel = iel + 1
                lnodes(1,iel) = uniform_nodenumber(is  +ic,iz   ,ns)
                lnodes(2,iel) = uniform_nodenumber(is+2*ic,iz-1 ,ns)
                lnodes(3,iel) = uniform_nodenumber(is+2*ic,iz+1 ,ns)
                lnodes(4,iel) = uniform_nodenumber(is  +ic,iz+1 ,ns)
                eltype(iel) = 'semino'
                coarsingloc(iel)=.true.

                iel = iel + 1
                lnodes(1,iel) = uniform_nodenumber(is+2*ic,iz-1 ,ns)
                lnodes(2,iel) = uniform_nodenumber(is+3*ic,iz   ,ns)
                lnodes(3,iel) = uniform_nodenumber(is+3*ic,iz+1 ,ns)
                lnodes(4,iel) = uniform_nodenumber(is+2*ic,iz+1 ,ns)
                eltype(iel) = 'semino'
                coarsingloc(iel)=.true.

                iel = iel + 1
                lnodes(1,iel) = uniform_nodenumber(is+2*ic,iz-1 ,ns)
                lnodes(2,iel) = uniform_nodenumber(is+4*ic,iz-1 ,ns)
                lnodes(3,iel) = uniform_nodenumber(is+4*ic,iz   ,ns)
                lnodes(4,iel) = uniform_nodenumber(is+3*ic,iz   ,ns)
                eltype(iel) = 'curved'
                coarsingloc(iel)=.true.

                iel = iel + 1
                lnodes(1,iel) = uniform_nodenumber(is+3*ic,iz   ,ns)
                lnodes(2,iel) = uniform_nodenumber(is+4*ic,iz   ,ns)
                lnodes(3,iel) = uniform_nodenumber(is+4*ic,iz+1 ,ns)
                lnodes(4,iel) = uniform_nodenumber(is+3*ic,iz+1 ,ns)
                eltype(iel) = 'curved'
                coarsingloc(iel)=.true.

             enddo
          endif
     enddo
     icold = ic
  enddo
  ! Last series of layer after last remeshing
  icc = nc + 1
  ic = 2**(icc-1)
  do iz = iclev(icc-1)-2, 1,-1
     do is = 1, ns, ic
       iel = iel + 1
       lnodes(1,iel) = uniform_nodenumber(is   ,iz ,ns)
       lnodes(2,iel) = uniform_nodenumber(is+ic,iz ,ns)
       lnodes(3,iel) = uniform_nodenumber(is+ic,iz+1,ns)
       lnodes(4,iel) = uniform_nodenumber(is   ,iz+1,ns)
       eltype(iel) = 'curved'
     enddo
  enddo

end subroutine define_lnodesclean
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine define_central_region

  use data_bkgrdmodel
  use data_coarse
  use data_diag
  use data_grid
  use global_parameters
  use data_spec
  use data_mesh

  real(kind=dp)   , dimension(:,:,:), allocatable :: crd_cyl
  real(kind=dp)   , dimension(:,:,:), allocatable :: crd_cent
  real(kind=dp)   , dimension(8,2)                :: crd_control_nodes
  real(kind=dp)   , dimension(:), allocatable     :: scc,zcc
  integer           :: nr, nex, nzs, nzss
  real(kind=dp)     :: ric, roc
  integer           :: npts, nelc, nelctot, nptscc
  integer           :: ipt, inode, iel, is, iz, ipti, iptsq, iptbuf, ieltest
  real(kind=dp)     :: test

  integer           :: ix, iy, maxy
  real(kind=dp)     :: rad, angle, p, x, y

  integer :: ntheta_opt, ntheta_opt_buff, nn

  ! number of latitude slices at shell - central region boundary (termed ib)
  ns_ib = ns / 2**nc
  if (dump_mesh_info_screen) write(*,*) 'ns_ib is ', ns_ib

  ! define dimensions of central square (just one!)
  lsq = lsq_fac * ri

  ! define number of divisions in one direction for central square
  ndivs = ns_ib/2

  if (only_suggest_ntheta) then
     ntheta_opt_buff = -1
     write(*,*)
     write(*,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
     write(*,*) '   suggested number of theta slices for optimal mesh decomposition:'
     do nn=1, 10
         ntheta_opt = ndivs / (2 * nn)
         if (mod(ntheta_opt, 4) > 0) ntheta_opt = ntheta_opt + 4 - mod(ntheta_opt, 4)
         if (ntheta_opt > 4) then
            if (ntheta_opt /= ntheta_opt_buff) write(*,*) ntheta_opt
         else
            exit
         endif
         ntheta_opt_buff = ntheta_opt
     enddo
     write(*,*) '   1, 2 and 4 are always decomposed optimally'
     write(*,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
     write(*,*)
     write(*,*) 'ONLY_SUGGEST_NTHETA was set, hence stopping now. Set to false to actually generate a mesh!'
     call exit()
  endif

  ! compute number of necessary extra coarsenings nex
  nex = ns_ib/(2*ndivs)-1
  if (dump_mesh_info_screen) write(*,*) ' NEX IS ', nex

  ! compute number of radial levels from center of the earth to
  ! shell boundary (equidistant spacing)
  ! PREM CASE
  ! dr = radius(2) - radius(1) ; write(*,*) ' DR IS ' ,dr
  ! nzs = int((ri-lsq)/dr)
  ! ALEX GRE; af 2009, 4 years later, that seems sensible but I do not remember why I
  ! picked this value. Might be useful to keep what follows compatible with nzs=/1
  ! MvD 2013: but I think it's not
  nzs = 1

  ! check
  test = (dsqrt(2.d0)*lsq) / (.9d0 *(lsq + (ri-lsq)/dble(nzs)))
  if ( test >= 1. ) then
     write(*,*) ' STOP: failed test for generation of central region'
     write(*,*) ' test is ',test
     stop
  endif

  ! nzs = int(real(ndivs)*(ri-lsq)/lsq)
  ! nzs = 2
  nr  = ndivs + nzs  ! will come from xtoto (newdiscont_meshing)
  if (dump_mesh_info_screen) write(*,*) 'nr is ', nr

  ! roc = ri ; ric = roc - real((nzs-1))*(roc-lsq)/dble(nzs)
  ! Assuming nzs = 1
  roc = ri ; ric = roc*(dble(nr-1)/dble(nr))
  if (dump_mesh_info_screen) then
     write(*,*) 'RIC = ', RIC
     write(*,*) 'ROC = ', ROC
     write(*,*) 'nzs = ', nzs
  endif

  if ( nzs < 2 * nex + 1) then
   write(*,*) 'incompatibility in central square between nex and nr'
   stop
  endif
  ! The '+1' indicates that we want the  bottommost "spherical layer" to
  ! be distinguished from the others, as it will be connected to the central square.

  ! We define the (nzs-1) spherical layers of the central region
  ! as for the exterior spherical shell.
  ! we define, again, the size of the associated uniform grid.
  ! nzss is the number of purely spherical levels
  nzss = nzs - 1
  if (dump_mesh_info_screen) write(*,*) 'nzss = ', nzss

  ! MvD: this is never true for nzs = 1 as hardcoded above
  !      nzs /= 1 however crashes
  if (nzss /= 0) then

     npts = (nzss+1) * (ns_ib+1)

     !!
     !! TNM JULY 2009: Old method, redundant...
     !!
     allocate( crd_cyl(1:ns_ib+1,1:nzss+1,2)) ;  crd_cyl(:,:,:)=0.
     allocate(crd_cent(1:ns_ib+1,1:nzss+1,2)) ; crd_cent(:,:,:)=0.

     allocate(s_unif (1:npts), z_unif (1:npts))

     ! define reference cylindrical coordinates
     call def_ref_cart_coordinates(ns_ib,nzss,crd_cyl,.true.)
     call def_control_nodes(crd_control_nodes,ric,roc)

     ! 2) Use the map_spheroid function
     call def_mapped_coordinates(ns_ib,nzss,crd_cent,crd_cyl,crd_control_nodes)
     !
     ! 3) Use lexicographical ordering to create 1d arrays of coordinates
      call def_global_coordinates(npts,s_unif,z_unif,ns_ib,nzss,crd_cent)
     !
     ! Again, we define the levels at which the
     ! potential coarsenings will take place
     allocate(iclevc(0:nex+1)) ; iclevc(:) = 0
     iclevc(0)     = nzss + 1
     iclevc(1)     = nzss
     iclevc(nex+1) = 1
     ! compute new number of elements
     call compute_nelclean(nelc,nex,iclevc,ns_ib,nzss)
     write(*,*) 'nelc ', nelc , 'npts'
     ! define elements topology and type
     allocate(lnodesi(4,nelc)) ; lnodesi(:,:) = 0
     allocate(eltypei(nelc),coarsingi(nelc))
     ! call define_lnodesclean(nelc,lnodesi,eltypei,coarsingi,nex,iclevc,ns_ib,nzss)
     ! fill si and zi arrays
     neli = nelc
     allocate(si(4*neli),zi(4*neli)) ; si(:) = 0. ; zi(:) = 0.

     do iel = 1, nelc
        do inode = 1, 4
           ipt = lnodesi(inode,iel)
           ipti = (iel-1)*4 + inode
           si(ipti) = s_unif (ipt)
           zi(ipti) = z_unif (ipt)
        enddo
     enddo

     write(*,*) 'npts = ', npts

  else
     npts = 0
     nelc = 0
  endif

  ! define grid points for central region
  ! Points numbering and coordinates
  nptscc = (ndivs+1)**2

  allocate(scc(nptscc+npts),zcc(nptscc+npts))
  scc(:) = 0.d0
  zcc(:) = 0.d0

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !        NEW METHOD: |x|^p + |y|^p = r^p
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  allocate(s_arr(ndivs+1,ndivs+1),z_arr(ndivs+1,ndivs+1))
  s_arr = zero
  z_arr = zero

  if (dump_mesh_info_screen) write(*,*)'CENTR: ndivs,nr',ndivs,nr
  if (dump_mesh_info_screen) write(*,*)'CENTR: ri',ri

  if (dump_mesh_info_files) then
     open(unit=46,file=diagpath(1:lfdiag)//'/fort.46')
     open(unit=47,file=diagpath(1:lfdiag)//'/fort.47')
     open(unit=48,file=diagpath(1:lfdiag)//'/fort.48')
  endif

  do ix=1, 2*ndivs+1
     maxy=int(min(ix,ceiling(dble(ix)/2.)+1))
     if (dump_mesh_info_screen) write(*,*)'CENTR: ix,maxy',ix,maxy
     do iy=1, maxy
        rad = dble(ix) / dble(2*ndivs+1) * (ri-maxh_icb/router)
        ! linear
        p = dble(ix) / dble(2*ndivs)+1.

        if (iy > 1) then
          angle = tan( pi/2.d0 * ( 1 - dble(iy-1)/dble(ix) ) )
          y = rad / (angle**p + 1.d0)**(1.d0/p)
        else
          angle = 0.d0
          y = 0.d0;
        endif
        x = dble(( rad**p-abs(y)**p )**(1.d0/p))

        ! compute s,z coordinates in rotated frame and indices
        if (mod(ix,2) == 0) then

           !   below diagonal, s >= z
           s_arr(int(ix/2)+1,maxy-iy+1) = x + y
           z_arr(int(ix/2)+1,maxy-iy+1) = x - y
           !   above diagonal, s < z (switched indices, negative y)
           s_arr(maxy-iy+1,int(ix/2)+1) = x - y
           z_arr(maxy-iy+1,int(ix/2)+1) = x + y

           if (dump_mesh_info_files) then
              write(46,*)y,x
              write(47,*)x+y,x-y
              write(48,*)x-y,x+y
           endif

        endif
     enddo
  enddo
  if (dump_mesh_info_files) then
     close(46)
     close(47)
     close(48)
  endif

  ! earth center
  s_arr(1,1) = zero
  z_arr(1,1) = zero

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! TNM MAY 2007: Apply stretching for 45-deg pathological elements
  !
  ! This is a trial-and-error fix to the triangularly deformed elements
  ! along the 45 deg diagonal inside the central cube which lead to
  ! grid spacing of diagonal points (is,iz) vs. (is+1,iz+1) of about
  ! a factor of 2-4 lower than the predicted/expected value (for higher
  ! resolution grids at least). The solver used to blow up as soon as
  ! reaching the inner core if this is left out. One should be careful
  ! and examine this region for each new high resolution grid before
  ! running the solver to be sure it looks properly (i.e. now crossing
  ! element boundaries etc).
  if (dump_mesh_info_files) then
     open(unit=5559,file=diagpath(1:lfdiag)//'/fort.5559')
  endif

  ! loop along diagonal: is==iz
  do is=2, ndivs+1
     p = 0.5d0 * (s_arr(is,is) - s_arr(is-1,is-1)) * dble(is) / dble(ndivs)
     if (dump_mesh_info_files) write(5559,*)is,is,p,s_arr(is,is)
     s_arr(is,is) = s_arr(is,is) + p
     z_arr(is,is) = z_arr(is,is) + p
     if (dump_mesh_info_files) write(5559,*)is,is,p,s_arr(is,is)
  enddo

  if (dump_mesh_info_files) then
     close(5559)
  endif

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !@TODO add here shrinking of axial elements as in the sperical shell
  !      might not be necessary, but keeping this in a comment for now
  !do is=ndivs+1, 2, -1
  !   do iz=1, ndivs+1
  !      s_arr(is,iz) = s_arr(is,iz) - s_arr(2,iz) * 0.25 * (ndivs + 2 - is) / ndivs
  !   enddo
  !enddo
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! make sure max val is ri
  s_arr = s_arr / maxval(abs(s_arr)) * (ri - maxh_icb/router)
  z_arr = z_arr / maxval(abs(z_arr)) * (ri - maxh_icb/router)

  if (dump_mesh_info_files) then
     open(unit=44,file=diagpath(1:lfdiag)//'/fort.44')
     open(unit=45,file=diagpath(1:lfdiag)//'/fort.45')
  endif

  do is=1,ndivs+1
     do iz=1,ndivs+1
        ipt = uniform_nodenumber(is,iz,ndivs,npts)
        scc(ipt) = s_arr(is,iz)
        zcc(ipt) = z_arr(is,iz)
        if (dump_mesh_info_files) write(44,*) scc(ipt), zcc(ipt)
        if (dump_mesh_info_files) write(45,*) is,iz,s_arr(is,iz),z_arr(is,iz)
     enddo
  enddo

  if (dump_mesh_info_files) then
     close(44)
     close(45)
  endif

  ! TNM JULY 2009: Need this for gather_skeleton
  if (nelc /= 0) then
     do iel = 1, nelc
        write(*,*) 'af iel ', iel
        do inode = 1, 4
           ipt = lnodesi(inode,iel)
           ipti = (iel-1)*4 + inode
           si(ipti) = scc(ipt)
           zi(ipti) = zcc(ipt)
        enddo
     enddo
  endif

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !        END NEW METHOD: |x|^p + |y|^p = r^p
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! total number of elements
  nelctot = nelc + 2*ndivs + ndivs**2
  if (dump_mesh_info_screen) write(*,*) 'nelctot ', nelctot
  nelsq = ndivs**2
  allocate(lnodessq(4, nelsq))
  allocate(eltypesq(nelsq))
  iel = 0
  do iz = 1, ndivs
     do is = 1, ndivs
        iel = iel + 1
        lnodessq(1,iel) = uniform_nodenumber(is  ,iz  ,ndivs,npts)
        lnodessq(2,iel) = uniform_nodenumber(is+1,iz  ,ndivs,npts)
        lnodessq(3,iel) = uniform_nodenumber(is+1,iz+1,ndivs,npts)
        lnodessq(4,iel) = uniform_nodenumber(is  ,iz+1,ndivs,npts)
        eltypesq(iel) = 'linear'
     enddo
  enddo

  ! FILL ARRAY SSQ AND ZSQ
  allocate(ssq(4*nelsq),zsq(4*nelsq))
  ssq(:) = 0.
  zsq(:) = 0.
  do iel = 1, nelsq
     do inode = 1, 4
        ipt = lnodessq(inode,iel)
        iptsq = 4*(iel-1)+inode
        ssq(iptsq) = scc(ipt)
        zsq(iptsq) = zcc(ipt)
     enddo
  enddo

  ! gnuplot dump
  if (dump_mesh_info_files) then
     open(unit=3,file=diagpath(1:lfdiag)//'/testcrd.dat',STATUS="OLD",POSITION="APPEND")
     do iel = 1, ndivs**2
        do inode = 1, 4
           ipt = lnodessq(inode,iel)
           write(3,*) scc(ipt),zcc(ipt)
        enddo
        ipt = lnodessq(1,iel)
        write(3,*) scc(ipt),zcc(ipt)
        write(3,*)
     enddo
  endif

  ! THE REST IN THIS ROUTINE IS ABOUT THE buffer layer**************************
  nelbuf = 2*ndivs
  allocate(lnodesbuf(4,nelbuf)) ; lnodesbuf(:,:) = 0
  allocate(eltypebuf(nelbuf))
  iel = 0
  do is = 1, ndivs
     iel = iel + 1
     lnodesbuf(1,iel) = uniform_nodenumber(is,ndivs+1,ndivs,npts)
     lnodesbuf(2,iel) = uniform_nodenumber(is+1,ndivs+1,ndivs,npts)
     if ( neli > 0 ) then
        lnodesbuf(3,iel) = uniform_nodenumber(1+  (is)*2**nex,1,ns_ib)
        lnodesbuf(4,iel) = uniform_nodenumber(1+(is-1)*2**nex,1,ns_ib)
     else ! CONNECT WITH OUTER SHELL
        ieltest = nelo + is -2*ndivs
        lnodesbuf(3,iel) = (ieltest-1)*4 + 2
        lnodesbuf(4,iel) = (ieltest-1)*4 + 1
     endif
     eltypebuf(iel) = 'semino'
  enddo
  do iz = 1, ndivs
   iel = iel + 1
   !lnodesbuf(1,iel) = uniform_nodenumber(ndivs+1,iz,ndivs,npts)
   !lnodesbuf(2,iel) = uniform_nodenumber(ns_ib+1-(iz-1)*2**nex,1,ns_ib)
   !lnodesbuf(3,iel) = uniform_nodenumber(ns_ib+1-(iz)*2**nex,1,ns_ib)
   !lnodesbuf(4,iel) = uniform_nodenumber(ndivs+1,iz+1,ndivs,npts)
   lnodesbuf(4,iel) = uniform_nodenumber(ndivs+1,iz,ndivs,npts)
   if ( neli > 0 ) then
      lnodesbuf(1,iel) = uniform_nodenumber(ns_ib+1-(iz-1)*2**nex,1,ns_ib)
      lnodesbuf(2,iel) = uniform_nodenumber(ns_ib+1-(iz)*2**nex,1,ns_ib)
   else ! CONNECT WITH OUTER ELEMENT
      ieltest = nelo -( iz - 1 )
      lnodesbuf(1,iel) = (ieltest-1)*4 + 2
      lnodesbuf(2,iel) = (ieltest-1)*4 + 1
   endif
   lnodesbuf(3,iel) = uniform_nodenumber(ndivs+1,iz+1,ndivs,npts)
   eltypebuf(iel) = 'semiso'
  enddo

  ! fill array sbuf and zbuf
  allocate(sbuf(4*nelbuf),zbuf(4*nelbuf))
  sbuf(:)=0.d0
  zbuf(:)=0.d0
  do iel = 1, nelbuf

     if (neli > 0) then
        do inode = 1, 4
           ipt = lnodesbuf(inode,iel)
           iptbuf = 4*(iel-1) + inode
           sbuf(iptbuf) = scc(ipt)
           zbuf(iptbuf) = zcc(ipt)
        enddo
     else ! CONNECT WITH OUTER SHELL
        if ( iel < (nelbuf/2 + 1) ) then
           do inode = 1,2
              ipt = lnodesbuf(inode,iel)
              iptbuf = 4*(iel-1) + inode
              sbuf(iptbuf) = scc(ipt)
              zbuf(iptbuf) = zcc(ipt)
           enddo
           do inode = 3,4
              ipt = lnodesbuf(inode,iel)
              iptbuf = 4*(iel-1) + inode
              sbuf(iptbuf) = so(ipt)
              zbuf(iptbuf) = zo(ipt) !; write(*,*) zo(ipt), ipt
           enddo
        else
           do inode = 1,2
              ipt = lnodesbuf(inode,iel)
              iptbuf = 4*(iel-1) + inode
              sbuf(iptbuf) = so(ipt)
              zbuf(iptbuf) = zo(ipt) !; write(*,*) zo(ipt), ipt
           enddo
           do inode = 3,4
              ipt = lnodesbuf(inode,iel)
              iptbuf = 4*(iel-1) + inode
              sbuf(iptbuf) = scc(ipt)
              zbuf(iptbuf) = zcc(ipt)
           enddo
        endif
     endif
  enddo

  ! gnuplot dump
  if (dump_mesh_info_files) then
     open(unit=3,file=diagpath(1:lfdiag)//'/testcrd.dat',STATUS="OLD",POSITION="APPEND")
     do iel = 1, nelbuf
        do inode = 1, 4
           ipt = lnodesbuf(inode,iel)
           iptbuf = 4*(iel-1) + inode
           write(3,*) sbuf(iptbuf),zbuf(iptbuf)
        enddo
        ipt = lnodesbuf(1,iel)
        iptbuf = 4*(iel-1) + 1
        write(3,*) sbuf(iptbuf),zbuf(iptbuf)
        write(3,*)
     enddo
  endif

end subroutine define_central_region
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine def_ref_cart_coordinates(nst, nzt, crd, inner_shell)
  use data_grid
  integer, intent(in)   :: nst, nzt
  real(kind=dp)   , dimension(1:nst+1,1:nzt+1,2), intent(out) :: crd
  logical, optional     :: inner_shell

  integer               :: is, iz
  real(kind=dp)         :: ds, dz

  ! uniform grid in s and z
  if (nst /= 0 .and. nzt /= 0) then
     ds = 2./dble(nst)
     dz = 2./dble(nzt)
  else
     write(*,*)'nst and nzt are ZERO!'
     stop
  endif

  do iz = 1, nzt+1
     do is = 1, nst+1
        crd(is,iz,1) = -1. + dble(is-1) * ds
        crd(is,iz,2) = -1. + dble(iz-1) * dz
        if (PRESENT(inner_shell)) crd(is,iz,2) = -1. + dble(iz-1) * dz
     enddo
  enddo
end subroutine def_ref_cart_coordinates
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine def_ref_cart_coordinates_discont(nst, nzt, crd, dz)

  use data_grid, only: axisfac
  use data_bkgrdmodel, only: nc_init, nthetaslices
  use data_pdb, only: theta_max_proc, theta_min_proc

  integer, intent(in) :: nst, nzt
  real(kind=dp), dimension(1:nst+1,1:nzt+1,2), intent(out) :: crd
  real(kind=dp), dimension(1:nzt) :: dz

  integer           :: is, iz
  real(kind=dp)     :: ds1, ds2

  real(kind=dp)     :: pi2
  integer           :: iproc


  !Make axial elements a bit smaller, to avoid artifacts from axial integration scheme
  ds1 = 2.d0 / dble(nst) * axisfac
  ds2 = (2.d0 - 2**nc_init * ds1) / dble(nst - 2**nc_init)

  do is = 1, nst+1
     do iz = 1, nzt+1
        if (is <= 2**nc_init) then
           crd(is,iz,1) = -1.d0 + dble(is-1) * ds1
        else
           crd(is,iz,1) = -1.d0 + 2**nc_init * ds1 + dble(is-1 - 2**nc_init) * ds2
        endif
     enddo
     crd(is,1,2) = -1.d0
     do iz = 2, nzt+1
        crd(is,iz,2) = crd(is,iz-1,2) +  dz(iz-1)
     enddo
  enddo

  ! Create colatitude bounds array for outer shell
  ds1 = 1.d0 / dble(nst) * axisfac
  ds2 = (1.d0 - 2**nc_init * ds1) / dble(nst - 2**nc_init)

  write(*,*) ds1, ds2, nst

  pi2 = 2.d0 * dasin(1.d0)
  allocate(theta_min_proc(0:nthetaslices-1), theta_max_proc(0:nthetaslices-1))
  theta_min_proc(:) = 0.d0
  theta_max_proc(:) = 0.d0
  theta_max_proc(nthetaslices-1) = pi2

  do iproc = 0, nthetaslices-2
     theta_min_proc(iproc+1) = 0.5d0 * pi2 * (ds1 * 2**nc_init &
                       + ds2 * (nst * 2 / nthetaslices * (iproc + 1) - 2**nc_init))
     theta_max_proc(iproc)   = 0.5d0 * pi2 * (ds1 * 2**nc_init &
                       + ds2 * (nst * 2 / nthetaslices * (iproc + 1) - 2**nc_init))
  enddo

end subroutine def_ref_cart_coordinates_discont
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine gather_skeleton
! This routine defines global arrays to assemble skeleton.

  use data_mesh
  use data_grid, only: ndivs

  integer :: ipt, inode, iel, istart, is, iz

  write(*,*)
  write(*,"(10x,'SKELETON INFORMATIONS (NORTHERN HEMISPHERE ONLY)')")
  write(*,*)
  write(*,"(10x,'Number of elements in the outer shell:    ',i10)")  nelo
  write(*,"(10x,'Number of elements in the inner shell:    ',i10)")  neli
  write(*,"(10x,'Number of elements in the buffer layer:   ',i10)")  nelbuf
  write(*,"(10x,'Number of elements in the central square: ',i10)")  nelsq

  neltot = nelo + neli + nelsq + nelbuf
  write(*,*)
  write(*,"(10x,'Total num. of elements in northern skel.: ',i10)")  neltot
  write(*,*)
  write(*,*)

  npointot = 4 * neltot
  allocate(sg(npointot),zg(npointot))
  sg(:) = 0.d0
  zg(:) = 0.d0

  ! outer shell
  istart = 1
  if (allocated(sg)) sg(istart:4*nelo) = so(1:4*nelo)
  if (allocated(zg)) zg(istart:4*nelo) = zo(1:4*nelo)

  ! inner shell
  istart = 4*nelo + 1
  if (allocated(si)) sg(istart:istart+4*neli-1) = si(1:4*neli)
  if (allocated(zi)) zg(istart:istart+4*neli-1) = zi(1:4*neli)

  ! buffer
  istart = 4*(nelo+neli) + 1
  if (allocated(sbuf)) sg(istart:istart+4*nelbuf-1) = sbuf(1:4*nelbuf)
  if (allocated(zbuf)) zg(istart:istart+4*nelbuf-1) = zbuf(1:4*nelbuf)


  ! central square region
  !istart = 4*(nelo+neli) + 1
  istart = 4*(nelo+neli+nelbuf) + 1
  if (allocated(ssq)) sg(istart:istart+4*nelsq-1) = ssq(1:4*nelsq)
  if (allocated(zsq)) zg(istart:istart+4*nelsq-1) = zsq(1:4*nelsq)

  ! define a mapping array for domain decomposition: is,iz = => global iel
  allocate(central_is_iz_to_globiel(1:ndivs,1:ndivs))
  iel = 0
  do iz = 1, ndivs
     do is = 1, ndivs
       iel = iel + 1

       central_is_iz_to_globiel(is,iz) = nelo + neli + nelbuf + iel
    enddo
  enddo


  ! gather element type
  allocate(lnodesg(4,neltot))
  lnodesg(:,:) = 0
  allocate(eltypeg(neltot),coarsing(neltot))
  if (allocated(eltypeo)) eltypeg(1:nelo) = eltypeo(1:nelo)
  if (allocated(eltypei)) eltypeg(nelo+1:nelo+neli) = eltypei(1:neli)
  if (allocated(eltypebuf)) eltypeg(nelo+neli+1:nelo+neli+nelbuf) = eltypebuf(1:nelbuf)
  if (allocated(eltypeg)) eltypeg(nelo+neli+nelbuf+1:neltot) = eltypesq(1:nelsq)

  coarsing = .false.
  coarsing(1:nelo) = coarsingo(1:nelo)

  do iel = 1, neltot
     do inode = 1, 4
        ipt = (iel-1)*4 + inode
        lnodesg(inode,iel) = ipt
     enddo
  enddo

end subroutine gather_skeleton
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine generate_southern_hemisphere

  use data_mesh
  use data_diag

  integer                                     :: neltot2, npointot2, ipt, iel,inode
  real(kind=dp)   , dimension(:), allocatable :: sg2,zg2
  integer, dimension(:,:), allocatable        :: lnodesg2
  character(len=6), dimension(:), allocatable :: eltypeg2
  logical, dimension(:), allocatable          :: coarsing2

  npointot2 = 2 * npointot
  neltot2 = 2*neltot
  allocate(sg2(npointot2),zg2(npointot2)) ; sg2(:)=0.d0 ; zg2(:)=0.d0
  allocate(lnodesg2(4,neltot2)) ; lnodesg2(:,:) = 0
  allocate(eltypeg2(neltot2),coarsing2(neltot2))
  do ipt = 1, npointot
     sg2(ipt) = sg(ipt)
     zg2(ipt) = zg(ipt)
     sg2(ipt+npointot) =  sg(ipt)
     zg2(ipt+npointot) = -zg(ipt)
  enddo

  do iel = 1, neltot

     do inode = 1, 4
        ipt = 4*(iel-1) + inode
        lnodesg2(inode,iel) = ipt
        lnodesg2(5-inode,iel+neltot) = ipt + npointot
     enddo

     eltypeg2(iel) = eltypeg(iel)
     coarsing2(iel) = coarsing(iel)
     coarsing2(iel+neltot) = coarsing(iel)

     if (eltypeg(iel) == 'semiso') then
        eltypeg2(iel+neltot) = 'semino'
     else if (eltypeg(iel) == 'semino') then
        eltypeg2(iel+neltot) = 'semiso'
     else
        eltypeg2(iel+neltot) = eltypeg(iel)
     endif

  enddo
  ! copy back into original arrays
  deallocate(sg,zg,lnodesg,eltypeg,coarsing)
  npointot = npointot2
  neltot = neltot2

  allocate(sg(npointot),zg(npointot)) ; sg(:)=0.d0 ; zg(:)=0.d0
  allocate(lnodesg(4,neltot)) ; lnodesg(:,:) = 0
  allocate(eltypeg(neltot),coarsing(neltot))

  sg(:) = sg2(:)
  zg(:) = zg2(:)
  lnodesg(:,:) = lnodesg2(:,:)
  eltypeg(:) = eltypeg2(:)
  coarsing(:)=coarsing2(:)

  deallocate(sg2,zg2,lnodesg2,eltypeg2,coarsing2)

end subroutine generate_southern_hemisphere
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine donot_generate_southern_hemisphere

  use data_mesh
  use data_diag

  real(kind=dp)   , dimension(:), allocatable   :: sg2, zg2
  integer, dimension(:,:), allocatable          :: lnodesg2
  character(len=6), dimension(:), allocatable   :: eltypeg2
  integer :: neltot2
  integer :: npointot2
  integer :: ipt
  integer :: iel, inode

  npointot2 =  npointot
  neltot2 = neltot
  allocate(sg2(npointot2),zg2(npointot2))
  sg2(:) = 0.d0
  zg2(:) = 0.d0
  allocate(lnodesg2(4,neltot2))
  lnodesg2(:,:) = 0
  allocate(eltypeg2(neltot2))

  do ipt = 1, npointot
     sg2(ipt) = sg(ipt)
     zg2(ipt) = zg(ipt)
  enddo

  do iel = 1, neltot
     do inode = 1, 4
        ipt = 4*(iel-1) + inode
        lnodesg2(inode,iel) = ipt
        eltypeg2(iel) = eltypeg(iel)
     enddo
  enddo
  ! copy back into original arrays
  deallocate(sg,zg,lnodesg,eltypeg)
  npointot = npointot2
  neltot = neltot2

  allocate(sg(npointot),zg(npointot))
  sg(:) = 0.d0
  zg(:) = 0.d0
  allocate(lnodesg(4,neltot))
  lnodesg(:,:) = 0
  allocate(eltypeg(neltot))

  sg(:) = sg2(:)
  zg(:) = zg2(:)
  lnodesg(:,:) = lnodesg2(:,:)
  eltypeg(:) = eltypeg2(:)
  deallocate(sg2,zg2,lnodesg2,eltypeg2)

end subroutine donot_generate_southern_hemisphere
!-----------------------------------------------------------------------------------------

end module meshgen
!=========================================================================================
