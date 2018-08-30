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
module mesh_info

  use data_gllmesh
  use data_mesh
  use data_spec
  use data_bkgrdmodel
  use data_grid
  use data_diag
  use global_parameters
  use data_numbering

  implicit none
  public :: define_regions, def_fluid_regions, def_solid_regions
  public :: define_boundaries

  private
contains

!-----------------------------------------------------------------------------------------
subroutine define_regions
  ! In this routine we establish the number of elements that
  ! belong to each region of the radial earth model of interest
  ! by considering the position of its center of mass
  integer                        :: iel, i, ipol, jpol, count
  real(kind=dp)                  :: s1, z1, theta
  real(kind=dp)   , dimension(4) :: rtmp

  ! define center of mass for each element
  allocate(scom(neltot))
  scom(:) = 0.d0
  allocate(zcom(neltot))
  zcom(:) = 0.d0
  allocate(thetacom(neltot))
  thetacom(:) = 0.d0
  allocate(rcom(neltot))
  rcom(:) = 0.d0
  allocate(rmin_el(neltot))
  rmin_el(:) = 0.d0
  allocate(rmax_el(neltot))
  rmax_el(:) = 0.d0

  do iel = 1, neltot
     count = 0
     ! find rmin_el and rmax_el
     do jpol=0, npol, npol
        do ipol=0, npol, npol
           count = count + 1
           rtmp(count) = dsqrt(sgll(ipol,jpol,iel)**2 + zgll(ipol,jpol,iel)**2)
        enddo
     enddo
     rmin_el(iel) = minval(rtmp)
     rmax_el(iel) = maxval(rtmp)

     scom(iel) = sgll(npol/2,npol/2,iel)
     zcom(iel) = zgll(npol/2,npol/2,iel)

     rcom(iel) = dsqrt(scom(iel)**2 + zcom(iel)**2)

     s1 = scom(iel)
     z1 = zcom(iel)
     theta = dacos(z1 / dsqrt(z1**2 + s1**2))
     if (theta < 0) then
        write(*,*) 'theta < 0', theta, s1, z1
        stop
     else if (theta > pi) then
        write(*,*) 'theta > pi', theta, s1, z1
        stop
     endif
     thetacom(iel) = theta
  enddo

  allocate(region(neltot), nel_region(neltot))
  region(:) = 0
  nel_region(:) = 0

  do iel = 1, neltot
     call assign_region(region(iel), scom(iel), zcom(iel))
     if (region(iel) > 0 .and. region(iel) < ndisc + 2) then
        nel_region(region(iel)) = nel_region(region(iel)) + 1
     else
        write(*,*) ' problem with assigning region to element number ', iel
        write(*,*) scom(iel), zcom(iel)
        stop
     endif
  enddo

  if (dump_mesh_info_screen) then
     print *
     write(*,*)'NUMBER OF ELEMENTS IN EACH SUBREGION:'
     do i=1, ndisc-1
        write(*,10)nel_region(i),discont(i)/1000., &
                  discont(i+1)/1000.
        write(*,11)real(nel_region(i))/real(neltot)*100., &
                   (discont(i)**2-discont(i+1)**2)/discont(1)**2*100.
        print *
     enddo
     write(*,12)nel_region(ndisc),discont(ndisc)/1000.
     write(*,11)real(nel_region(ndisc))/real(neltot)*100., &
                discont(ndisc)**2/discont(1)**2*100.
     print *
     write(*,*)'Elements summed over all regions:',sum(nel_region)
     write(*,*)'Total # Elements:',neltot
  endif

  if (neltot /= sum(nel_region)) then
     write(*,*)'In subroutine define_regions: Sum of elements in all regions not equal to total number'
     write(*,*) neltot, sum(nel_region)
     stop
  endif

  print *

10 format(i8,' elements between ',f6.1,' and ',f6.1)
11 format(f7.1,' element percent ',f4.1,' area percent')
12 format(i8,' elements below ',f6.1)

end subroutine define_regions
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine assign_region(ireg, s, z)

  integer, intent(out)          :: ireg
  real(kind=dp)   , intent(in)  :: s, z
  real(kind=dp)                 :: r
  integer                       :: idisc

  r = router*dsqrt(s**2+z**2) ! s and z are nondimensionalized by outer rad.
  do idisc = 2, ndisc
   if ( (r < discont(idisc-1)) .and. (r > discont(idisc)) ) ireg = idisc - 1
  enddo
  if ( r < discont(ndisc) ) ireg = ndisc

end subroutine assign_region
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine def_fluid_regions

  integer :: ielem, ifluid
  logical :: test_fluid

  allocate(fluid(neltot))
  fluid(:) = .false.
  neltot_fluid = 0

  do ielem = 1, neltot
     call check_fluid(test_fluid,ielem)
     if (test_fluid) then
        fluid(ielem) = .true.
        neltot_fluid = neltot_fluid + 1
     endif
  enddo

  if (dump_mesh_info_screen) then
     write(*,*) ' FLUID ELEMENTS '
     write(*,*) ' NUMBER OF FLUID ELEMENTS '
     write(*,*) neltot_fluid
  endif

  allocate(ielem_fluid(neltot_fluid))
  allocate(inv_ielem_fluid(neltot)); inv_ielem_fluid(:)=0
  ifluid = 0
  do ielem = 1, neltot
     if (fluid(ielem)) then
        ifluid = ifluid + 1
        ielem_fluid(ifluid) = ielem
        inv_ielem_fluid(ielem) = ifluid
     endif
  enddo
end subroutine def_fluid_regions
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine def_solid_regions
  integer :: ielem, isolid
  logical :: test_solid

  allocate(solid(neltot))
  solid(:) = .false.
  neltot_solid = 0

  do ielem = 1, neltot
     call check_solid(test_solid,ielem)
     if (test_solid) then
        solid(ielem) = .true.
        neltot_solid = neltot_solid + 1
     endif
  enddo

  if (dump_mesh_info_screen) then
     write(*,*) ' solid ELEMENTS '
     write(*,*) ' NUMBER OF solid ELEMENTS '
     write(*,*) neltot_solid
  endif
  allocate(ielem_solid(neltot_solid))
  allocate(inv_ielem_solid(neltot)); inv_ielem_solid(:)=0
  ! ielem_solid: given solid el number, return global
  ! inv_ielem_solid: given global el number, return solid
  isolid = 0
  do ielem = 1, neltot
     if (solid(ielem)) then
        isolid = isolid + 1
        ielem_solid(isolid) = ielem
        inv_ielem_solid(ielem) = isolid
     endif
  enddo
end subroutine def_solid_regions
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine check_fluid(test_fluid,  iel)
!
  logical, intent(out)  :: test_fluid
  integer, intent(in)   :: iel

  test_fluid = .false.
  ! FOR PREM
  if (.not. solid_domain(region(iel))) test_fluid = .true.

end subroutine check_fluid
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine check_solid(test_solid,iel)

  logical, intent(out) :: test_solid
  integer, intent(in) :: iel
  test_solid = .false.
  ! FOR PREM
  if (solid_domain(region(iel))) test_solid = .true.

end subroutine check_solid
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine define_boundaries
! Defines number, size (number of elements) and allocates arrays
! pertaining to the solid-fluid boundaries.
! The actual arrays needed by the solver
! are computed in define_my_boundary_neighbor.

  integer           :: j,ipol,ibelem
  real(kind=dp)     :: dmax, rbound
  integer           :: nbelemmax
  real(kind=dp), allocatable :: bdry_radius(:)

  if (neltot_fluid > 0 .and. neltot_solid > 0 ) then
     ! TODO: This is only true, if the fluid is not layered!!!
     nbcnd = 2*nfluidregions  ! 1=CMB; 2=ICB

     if (.not. solid_domain(ndisc)) nbcnd = nbcnd - 1

     write(*,*) '..... the number of fluid boundaries is not general enough....'
     write(*,*) '.....should insert a test on whether the fluid is indeed completely embedded!'
  else if (neltot_solid == 0) then
     nbcnd = 0
  else if (neltot_fluid == 0 ) then
     nbcnd = 0
  endif

  ! Allocate memory for the number of boundary elements
  allocate(nbelem(nbcnd),bdry_radius(nbcnd))
  ! set number of boundary elements to zero
  nbelem(:) = 0 ! careful: is a global variable!!!

  write(*,*)'size idom_fluid:',size(idom_fluid),nbcnd

  if (have_fluid) then
     dmax = min_distance_nondim
     do j = 1, nbcnd

        if (mod(j,2) /= 0) then ! upper boundary of fluid region
           rbound = discont(idom_fluid(j))/router
           if (dump_mesh_info_screen) then
              write(*,*)'ABOVE SOLID-FLUID BOUNDARY:',j,idom_fluid(j), &
                                                     rbound*router/1000.
              call flush(6)
           endif
        else  ! lower boundary of fluid region
           rbound = discont(idom_fluid(j-1)+1)/router
           if (dump_mesh_info_screen) then
              write(*,*)'BELOW SOLID-FLUID BOUNDARY:',j,idom_fluid(j-1)+1, &
                                                      rbound*router/1000.
              call flush(6)
           endif
        endif

        call belem_count_new(dmax, j, rbound)
        bdry_radius(j) = rbound
     enddo
     nbelemmax = maxval(nbelem(:))
     allocate (belem(nbelemmax,nbcnd))
     do j = 1, nbcnd
        if (mod(j,2) /= 0) then ! upper boundary of fluid region
           rbound = discont(idom_fluid(j))/router
        else  ! lower boundary of fluid region
           rbound = discont(idom_fluid(j-1)+1)/router
        endif
        call belem_list_new(dmax,j,rbound)
     enddo

     if (dump_mesh_info_files) then
        open(unit=6968,file=diagpath(1:lfdiag)//'/fort.6968')
        do j=1,nbcnd
           do ibelem=1,nbelem(j)
              write(6968,*)j,ibelem,sqrt(sgll(npol/2,0,belem(ibelem,j))**2 + &
                     zgll(npol/2,0,belem(ibelem,j))**2)*router/1000., &
                     sqrt(sgll(npol/2,npol,belem(ibelem,j))**2 + &
                     zgll(npol/2,npol,belem(ibelem,j))**2)*router/1000.
           enddo
        enddo
        close(6968)
     endif

     write(*,*)'allocating boundary arrays....',nbelemmax; call flush(6)

     allocate(my_neighbor(nbelemmax,nbcnd))
     allocate(bdry_above_el(nbelemmax/2,nbcnd),bdry_below_el(nbelemmax/2,nbcnd))
     allocate(bdry_solid_el(nbelemmax/2,nbcnd),bdry_fluid_el(nbelemmax/2,nbcnd))
     allocate(bdry_s(0:npol,nbelemmax/2,nbcnd),bdry_z(0:npol,nbelemmax/2,nbcnd))
     allocate(bdry_jpol_solid(nbelemmax/2,nbcnd))
     allocate(bdry_jpol_fluid(nbelemmax/2,nbcnd))
     allocate(bdry_globnum_above(0:npol,nbelemmax/2,nbcnd))
     allocate(bdry_globnum_below(0:npol,nbelemmax/2,nbcnd))
     allocate(bdry_locnum_above(0:npol,nbelemmax/2,nbcnd))
     allocate(bdry_locnum_below(0:npol,nbelemmax/2,nbcnd))

     allocate(bdry_globnum_solid(0:npol,nbelemmax/2,nbcnd))
     allocate(bdry_globnum_fluid(0:npol,nbelemmax/2,nbcnd))
     allocate(bdry_locnum_solid(0:npol,nbelemmax/2,nbcnd))
     allocate(bdry_locnum_fluid(0:npol,nbelemmax/2,nbcnd))

     call define_my_boundary_neighbor

  else ! only solid
     return
  endif

  !=======OUTPUT BOUNDARY INFO AND ARRAYS===========================
  ! MvD: crashes for nbcnd == 1, but do we need this anyway?
  if (dump_mesh_info_files .and. have_solid) then
     open(unit=67660,file=diagpath(1:lfdiag)//"/bdry_info.dat")
     !write(67660,*)nbcnd
     do j=1,nbcnd
        write(67660,*)j,bdry_radius(j),nbelem(j)
     enddo
     close(67660)

     open(unit=67661,file=diagpath(1:lfdiag)//"/bdry_el_neighbors.dat")
     open(unit=67662,file=diagpath(1:lfdiag)//"/bdry_el_above_below.dat")
     open(unit=67666,file=diagpath(1:lfdiag)//"/bdry_el_solid_fluid.dat")
     open(unit= 7775,file=diagpath(1:lfdiag)//"/bdry_globcoord_solel.dat")
     open(unit= 7776,file=diagpath(1:lfdiag)//"/bdry_globcoord_fluel.dat")
     open(unit= 7777,file=diagpath(1:lfdiag)//"/bdry_globcoord_myel.dat")
     open(unit= 7778,file=diagpath(1:lfdiag)//"/bdry_globcoord_neighborel.dat")
     open(unit=67671,file=diagpath(1:lfdiag)//"/bdry_sz.dat")
     open(unit=67672,file=diagpath(1:lfdiag)//"/bdry_glob_above_below.dat")
     open(unit=67673,file=diagpath(1:lfdiag)//"/bdry_loc_above_below.dat")
     open(unit=67674,file=diagpath(1:lfdiag)//"/bdry_glob_solid_fluid.dat")
     open(unit=67675,file=diagpath(1:lfdiag)//"/bdry_loc_solid_fluid.dat")
     do j=1,nbcnd
        do ibelem=1,nbelem(j)
           write(67661,*)belem(ibelem,j),belem(my_neighbor(ibelem,j),j)
        enddo
        do ibelem=1,nbelem(j)/2
           write(67662,*)bdry_above_el(ibelem,j),bdry_below_el(ibelem,j)
           write(67666,*)bdry_solid_el(ibelem,j),bdry_fluid_el(ibelem,j)

           write(7775,*)sgll(npol/2,npol/2,ielem_solid(bdry_solid_el(ibelem,j))), &
                        zgll(npol/2,npol/2,ielem_solid(bdry_solid_el(ibelem,j))), &
                        ielem_solid(bdry_solid_el(ibelem,j)),bdry_solid_el(ibelem,j)
           write(7776,*)sgll(npol/2,npol/2,ielem_fluid(bdry_fluid_el(ibelem,j))), &
                        zgll(npol/2,npol/2,ielem_fluid(bdry_fluid_el(ibelem,j))), &
                        ielem_fluid(bdry_fluid_el(ibelem,j)),bdry_fluid_el(ibelem,j)

           write(7777,*)sgll(npol/2,npol/2,bdry_above_el(ibelem,j)), &
                        zgll(npol/2,npol/2,bdry_above_el(ibelem,j)),bdry_above_el(ibelem,j)
           write(7778,*)sgll(npol/2,npol/2,bdry_below_el(my_neighbor(ibelem,j),j)), &
                           zgll(npol/2,npol/2,bdry_below_el(my_neighbor(ibelem,j),j)), &
                           bdry_below_el(my_neighbor(ibelem,j),j)

           do ipol=0,npol
              write(67671,*)bdry_s(ipol,ibelem,j),bdry_z(ipol,ibelem,j)
              write(67672,*)bdry_globnum_above(ipol,ibelem,j),bdry_globnum_below(ipol,ibelem,j)
              write(67673,*)bdry_locnum_above(ipol,ibelem,j),bdry_locnum_below(ipol,ibelem,j)
              write(67674,*)bdry_globnum_solid(ipol,ibelem,j),bdry_globnum_fluid(ipol,ibelem,j)
              write(67675,*)bdry_locnum_solid(ipol,ibelem,j),bdry_locnum_fluid(ipol,ibelem,j)
           enddo !ipol
        enddo !ibelem

     enddo !j
     close(67661)
     close(67662)
     close(67666)
     close(67671)
     close(67672)
     close(67673)
     close(67674)
     close(67675)
     close(7775)
     close(7776)
     close(7777)
     close(7778)

     ! reference coordinates
     open(unit=67663,file=diagpath(1:lfdiag)//'/bdry_refcoord.dat')
     open(unit=67664,file=diagpath(1:lfdiag)//'/bdry_refglob.dat')
     open(unit=67665,file=diagpath(1:lfdiag)//'/bdry_refslob.dat')
     open(unit=67667,file=diagpath(1:lfdiag)//'/bdry_refsolflu_el.dat')
     open(unit=67668,file=diagpath(1:lfdiag)//'/bdry_refflob.dat')

     do ibelem = 1,neltot

        if ( abs(rcom(ibelem)-bdry_radius(1)) <&
                abs( zgll(0,0,belem(1,1) )-zgll(npol,npol,belem(1,1) ) ) ) then
           write(67667,*)inv_ielem_solid(ibelem),inv_ielem_fluid(ibelem)
           do ipol=0,npol
              if (zcom(ibelem) > 0.) then ! northern hemisphere
                 write(67663,*)sgll(ipol,0,ibelem),zgll(ipol,0,ibelem)

                 if (rcom(ibelem) > bdry_radius(1)) then ! above, solid
                    write(67664,*)ibelem,iglob( (ibelem-1)*(npol+1)**2 + ipol + 1 )
                    write(67665,*)sgll(ipol,0,ibelem),zgll(ipol,0,ibelem), &
                         iglob_solid( (inv_ielem_solid(ibelem)-1)*(npol+1)**2 + ipol + 1 )
                 endif
              else ! southern hemisphere

                 write(67663,*)sgll(ipol,npol,ibelem),zgll(ipol,npol,ibelem)

                 if (rcom(ibelem) > bdry_radius(1)) then ! above, solid
                    write(67664,*)ibelem,iglob( (ibelem-1)*(npol+1)**2 &
                         + npol*(npol+1) + ipol + 1 )
                    write(67665,*)sgll(ipol,npol,ibelem),zgll(ipol,npol,ibelem), &
                         iglob_solid( (inv_ielem_solid(ibelem)-1)*(npol+1)**2  &
                         + npol*(npol+1) + ipol + 1 )

                 else ! below, fluid
                    write(67664,*)ibelem,iglob( (ibelem-1)*(npol+1)**2 + ipol + 1 )
                 endif

              endif
           enddo
        endif
        if (abs(rcom(ibelem)-bdry_radius(2)) < abs( zgll(0,0,belem(1,2) )-&
                                           zgll(npol,npol,belem(1,2) ) ) ) then

           write(67667,*)inv_ielem_solid(ibelem),inv_ielem_fluid(ibelem)
           do ipol=0,npol
              if (zcom(ibelem) > 0.) then ! northern hemisphere
                 write(67663,*)sgll(ipol,0,ibelem),zgll(ipol,0,ibelem)

              if (rcom(ibelem) > bdry_radius(2)) then ! above, fluid
                 write(67664,*)ibelem,iglob( (ibelem-1)*(npol+1)**2 + ipol + 1)
                 ! TNM 2011: Some issue with lower bound of iglob_fluid.... we
                 ! render this unimportant at this point.
              else ! below, solid
                 write(67664,*)ibelem,iglob((ibelem-1)*(npol+1)**2+&
                               npol*(npol+1) + ipol + 1 )
                 write(67665,*)sgll(ipol,npol,ibelem),zgll(ipol,npol,ibelem), &
                      iglob_solid((inv_ielem_solid(ibelem)-1)*(npol+1)**2+ &
                      npol*(npol+1) + ipol + 1 )
              endif

              else ! southern hemisphere
                 write(67663,*)sgll(ipol,npol,ibelem),zgll(ipol,npol,ibelem)

              if (rcom(ibelem) > bdry_radius(2)) then ! above, fluid
                 write(67664,*)ibelem,iglob((ibelem-1)*(npol+1)**2+ &
                               npol*(npol+1) + ipol + 1 )
                 ! TNM 2011: Some issue with lower bound of iglob_fluid.... we
                 ! render this unimportant at this point.
                 ! Maybe not unimportant. But probably caused by some ill-posed
                 ! indexing in this loop, like neltot...
              else ! below, solid
                 write(67664,*)ibelem,iglob( (ibelem-1)*(npol+1)**2+ ipol + 1 )
                 write(67665,*)sgll(ipol,0,ibelem),zgll(ipol,0,ibelem), &
                      iglob_solid( (inv_ielem_solid(ibelem)-1)*&
                                    (npol+1)**2 + ipol + 1 )
              endif

              endif
           enddo
        endif
     enddo

     close(67663)
     close(67664)
     close(67665)
     close(67667)
     close(67668)
  endif


end subroutine define_boundaries
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine define_my_boundary_neighbor
! Defines & computes all arrays needed to identify the elements hugging
! all solid-fluid boundaries, their respective cross-boundary neighbors,
! their respective global-global numbers, their respective flobal/slobal
! numbers and s,z coordinates of all GLL points on the boundaries.
! This will be used in the database to store for the solver.
! Note that nbelem(j) is the total amount of ALL elements at boundary j,
! i.e. counting above AND below.

  integer :: ibelem,jbelem,ipol
  integer :: myel,herel
  integer :: j,abovecount,belowcount
  real(kind=dp)    :: mytheta, hertheta
  real(kind=dp)    :: tolerance,rbound
  real(kind=dp)    :: rup,thetaup,rdown,thetadown

  ! respective global numbers
  integer :: ipt_glob_ab,ipt_glob_be
  integer :: ipt_slob_ab,ipt_slob_be
  integer :: ipt_flob_ab,ipt_flob_be

  integer :: jbelemmin
  real(kind=dp)    :: distmin, dist

  logical :: foundone

  tolerance = min_distance_nondim
  if (dump_mesh_info_screen) &
       write(*,*) 'Tolerance to find boundary-hugging partner:', tolerance

  do j=1, nbcnd

     if (mod(j,2) /= 0) then ! upper boundary of fluid region
        rbound = discont(idom_fluid(j))/router
     else  ! lower boundary of fluid region
        rbound = discont(idom_fluid(j-1)+1)/router
     endif

     ! should add conditional statement about the chosen one
     ! as the min_(jneqi) r_{ij}

     do ibelem = 1, nbelem(j)
        foundone = .false.
        myel = belem(ibelem,j)
        mytheta = thetacom(myel)

        do jbelem = 1, nbelem(j)
           herel = belem(jbelem,j)
           hertheta = thetacom(herel)
           dist = dabs(hertheta-mytheta)*rbound
           if (dist < tolerance .and. myel /= herel ) then
               my_neighbor(ibelem,j) = jbelem
               foundone = .true.
               exit
           endif
        enddo

        if (.not. foundone) then
           distmin = rbound
           do  jbelem = 1, nbelem(j)
              herel = belem(jbelem,j)
              hertheta = thetacom(herel)
              dist = dabs(hertheta-mytheta)*rbound
              if (dist < distmin .and. myel /= herel ) then
                 jbelemmin = jbelem
                 distmin = dist
              endif
           enddo
           my_neighbor(ibelem,j) = jbelemmin
        endif
        call flush(6)
     enddo
  enddo

  if (dump_mesh_info_files) then
     open(unit=131313,file=diagpath(1:lfdiag)//'/bdry_elems.dat')
     open(unit=6565,file=diagpath(1:lfdiag)//'/bdry_partners.dat')
     open(unit=6464,file=diagpath(1:lfdiag)//'/bdry_coords_elems_me.dat')
     open(unit=6465,file=diagpath(1:lfdiag)//'/bdry_coords_elems_her.dat')
  endif
  ! Easy check ... will that last?
  do j = 1, nbcnd
        abovecount= 0; belowcount=0
      if (mod(j,2) /= 0) then ! upper boundary of fluid region
         rbound = discont(idom_fluid(j))/router
      else  ! lower boundary of fluid region
         rbound = discont(idom_fluid(j-1)+1)/router
      endif

      do ibelem = 1, nbelem(j)

         ! easy check
         if (dump_mesh_info_files) &
         write(6565,*)'Bdry partners:',j,ibelem,my_neighbor(ibelem,j), &
                                    my_neighbor(my_neighbor(ibelem,j),j)
         if ( ibelem /= my_neighbor(my_neighbor(ibelem,j),j) ) then
            write(*,*)'Problem in boundary-partner mapping:'
            write(*,*)'Me:',ibelem
            write(*,*)"My partner's partner:", &
                       my_neighbor(my_neighbor(ibelem,j),j)
            stop
         endif

         if (dump_mesh_info_files) &
         write(6464,16)scom(belem(ibelem,j))*router/1000., &
                       zcom(belem(ibelem,j))*router/1000.,ibelem

         if (dump_mesh_info_files) &
         write(6465,16)scom(belem(my_neighbor(ibelem,j),j))*router/1000., &
                       zcom(belem(my_neighbor(ibelem,j),j))*router/1000., &
                       my_neighbor(ibelem,j)

16       format(2(1pe12.4),i4)

         ! determine whether element is above boundary
         if (rcom(belem(ibelem,j)) > rbound ) then

             ! I am above
             rup=rmin_el(belem(ibelem,j))*router/1000.d0
             thetaup=thetacom(belem(ibelem,j))*180.d0/pi
             rdown=rmax_el(belem(my_neighbor(ibelem,j),j))*router/1000.d0
             thetadown=thetacom(belem(my_neighbor(ibelem,j),j))*180.d0/pi
             abovecount= abovecount+1

             !================================================================
             ! ARRAYS NEEDED BY THE SOLVER: above_el, below_el, sbdry, zbdry
             !================================================================

             ! global element number
             bdry_above_el(abovecount,j)=belem(ibelem,j)
             bdry_below_el(abovecount,j)=belem(my_neighbor(ibelem,j),j)

             ! solid/fluid element & jpol numbers: these arrays are the ones needed
             ! in the time stepping algorithm to map between boundary and
             ! solid/fluid coordinate systems
             if (solid(belem(ibelem,j))) then ! Solid above
                bdry_solid_el(abovecount,j)=inv_ielem_solid(belem(ibelem,j))
                bdry_fluid_el(abovecount,j)=&
                                  inv_ielem_fluid(belem(my_neighbor(ibelem,j),j))
                bdry_jpol_solid(abovecount,j)=0
                if (zcom(belem(ibelem,j)) < 0.d0) bdry_jpol_solid(abovecount,j)=npol
                bdry_jpol_fluid(abovecount,j)=npol
                if (zcom(belem(ibelem,j)) < 0.) bdry_jpol_fluid(abovecount,j)=0

             else !Fluid above
                bdry_solid_el(abovecount,j)=&
                                  inv_ielem_solid(belem(my_neighbor(ibelem,j),j))
                bdry_fluid_el(abovecount,j)=inv_ielem_fluid(belem(ibelem,j))

                bdry_jpol_solid(abovecount,j)=npol
                if (zcom(belem(ibelem,j)) < 0.) bdry_jpol_solid(abovecount,j)=0

                bdry_jpol_fluid(abovecount,j)=0
                if (zcom(belem(ibelem,j)) < 0.) bdry_jpol_fluid(abovecount,j)=npol

                ! very crude way to accomodate the case of having the buffer layer just
                ! below the ICB. ...i.e., the jpol indices are switched for 45 < theta <
                ! 135 deg
                if (eltypeg(belem(my_neighbor(ibelem,j),j)) /= 'curved') then
                   if (dump_mesh_info_screen) write(*,*)'ELTYPE:',j,abovecount, &
                          eltypeg(belem(my_neighbor(ibelem,j),j))
                   if ( scom(belem(my_neighbor(ibelem,j),j)) >&
                        abs(zcom(belem(my_neighbor(ibelem,j),j))) ) then
                      bdry_jpol_solid(abovecount,j)=abs(bdry_jpol_solid(abovecount,j)-npol)
                   endif
                endif

             endif

             ! boundary coordinates
             if (zcom(belem(ibelem,j)) > 0.) then ! northern hemisphere
               bdry_s(0:npol,abovecount,j)=sgll(0:npol,0,belem(ibelem,j))
               bdry_z(0:npol,abovecount,j)=zgll(0:npol,0,belem(ibelem,j))
             else ! southern hemisphere
               bdry_s(0:npol,abovecount,j)=sgll(0:npol,npol,belem(ibelem,j))
               bdry_z(0:npol,abovecount,j)=zgll(0:npol,npol,belem(ibelem,j))
             endif

             ! arrays inside this loop are actually not needed elsewhere....
             ! global numbering is not invoked on the boundary!
             do ipol=0,npol

                ! global number above
                ipt_glob_ab = (belem(ibelem,j)-1)*(npol+1)**2 + ipol + 1
                if (zcom(belem(ibelem,j)) < 0.) then  ! south
                     ipt_glob_ab = (belem(ibelem,j)-1)* &
                                   (npol+1)**2 + npol*(npol+1) + ipol + 1
                endif

                ! global number below
                ipt_glob_be = (belem(my_neighbor(ibelem,j),j)-1)* &
                              (npol+1)**2 + npol*(npol+1) + ipol + 1
                if (zcom(belem(ibelem,j)) < 0.) then  ! south
                     ipt_glob_be = (belem(my_neighbor(ibelem,j),j)-1)* &
                                   (npol+1)**2  + ipol + 1
                endif

                if (solid(belem(ibelem,j))) then ! Solid above
                  ! slobal number above
                  ipt_slob_ab=(inv_ielem_solid(belem(ibelem,j))-1)*(npol+1)**2+ipol+1

                   if (zcom(belem(ibelem,j)) < 0.) then  ! south
                        ipt_slob_ab = ( inv_ielem_solid(belem(ibelem,j)) -1 )* &
                                      (npol+1)**2 + npol*(npol+1) + ipol + 1
                   endif
                   ! flobal number below
                   ipt_flob_be=(inv_ielem_fluid(belem(my_neighbor(ibelem,j),j))-1)*&
                                 (npol+1)**2 + npol*(npol+1) + ipol + 1
                    if (zcom(belem(ibelem,j)) < 0.) then  ! south
                      ipt_flob_be = &
                           (inv_ielem_fluid(belem(my_neighbor(ibelem,j),j))-1)* &
                           (npol+1)**2 + ipol + 1
                   endif

                   ! Fluid above
                else ! take my neighbor below
                   ! slobal number below
                   ipt_slob_be = (inv_ielem_solid(belem(my_neighbor(ibelem,j),j))-1)* &
                        (npol+1)**2 + npol*(npol+1) + ipol+1
                   if (zcom(belem(ibelem,j)) < 0.) then  ! south
                        ipt_slob_be = (inv_ielem_solid(belem(my_neighbor(ibelem,j),j))-1) &
                        * (npol+1)**2 + ipol + 1
                   endif
                   ! flobal number above
                   ipt_flob_ab = (inv_ielem_fluid(belem(ibelem,j))-1)* &
                        (npol+1)**2 + ipol+1
                    if (zcom(belem(ibelem,j)) < 0.) then  ! south
                        ipt_flob_ab = (inv_ielem_fluid(belem(ibelem,j))-1)* &
                        (npol+1)**2 + npol*(npol+1) + ipol + 1
                   endif
                endif ! solid/fluid above

                ! global number arrays
                bdry_globnum_above(ipol,abovecount,j)=iglob(ipt_glob_ab)
                bdry_globnum_below(ipol,abovecount,j)=iglob(ipt_glob_be)

                ! slobal/flobal number arrays
                if (solid(belem(ibelem,j))) then ! solid above

                  ! global numbering: solid
                   bdry_globnum_solid(ipol,abovecount,j)=iglob(ipt_glob_ab)
                   bdry_globnum_fluid(ipol,abovecount,j)=iglob(ipt_glob_be)

                   ! local numbering (i.e., refering to fluid or solid)
                   bdry_locnum_above(ipol,abovecount,j)=iglob_solid(ipt_slob_ab)
                   bdry_locnum_below(ipol,abovecount,j)=iglob_fluid(ipt_flob_be)

                   ! local numbering: solid (above or below)
                   bdry_locnum_solid(ipol,abovecount,j)=iglob_solid(ipt_slob_ab)
                   bdry_locnum_fluid(ipol,abovecount,j)=iglob_fluid(ipt_flob_be)

                else ! fluid above

                   ! global numbering: fluid
                   bdry_globnum_fluid(ipol,abovecount,j)=iglob(ipt_glob_ab)
                   bdry_globnum_solid(ipol,abovecount,j)=iglob(ipt_glob_be)

                   ! local numbering (i.e., refering to fluid or solid)
                   bdry_locnum_above(ipol,abovecount,j)=iglob_fluid(ipt_flob_ab)
                   bdry_locnum_below(ipol,abovecount,j)=iglob_solid(ipt_slob_be)

                   ! local numbering: (above or below)
                   bdry_locnum_solid(ipol,abovecount,j)=iglob_solid(ipt_slob_be)
                   bdry_locnum_fluid(ipol,abovecount,j)=iglob_fluid(ipt_flob_ab)
                endif

             enddo ! ipol

             !================================================================
             ! END OF ARRAYS NEEDED BY THE SOLVER
             !================================================================

         else
         ! I am below (only for output check)
            rdown=rmax_el(belem(ibelem,j))*router/1000.d0
            thetadown=thetacom(belem(ibelem,j))*180.d0/pi
            rup=rmin_el(belem(my_neighbor(ibelem,j),j))*router/1000.d0
            thetaup=thetacom(belem(my_neighbor(ibelem,j),j))*180.d0/pi
         endif

         if (dump_mesh_info_files) write(131313,30)rup,rdown,thetaup,thetadown

         if ( ibelem /= my_neighbor(my_neighbor(ibelem,j),j) ) then
           write(*,*)'ibelem different from my neighbor)',ibelem,j,nbelem(j)
           stop
         endif

      enddo ! ibelem (boundary elements)
  enddo ! j (# boundaries)

  close(131313)
  close(6565)
  close(6464)
  close(6465)

  write(*,*)'Done with boundary neighbor search.'; call flush(6)

30 format(12(1pe12.5,2x))

end subroutine define_my_boundary_neighbor
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine belem_count_new(dmax, j, rbound)
! This routine determines the number of boundary elements
! for the j-th component of the velocity field, in the spheroidal
! container case.
!
  integer, intent(in)           :: j
  real(kind=dp)   , intent(in)  :: dmax, rbound
  integer                       :: ibd, ielem, ipol, jpol, itest
  real(kind=dp)                 :: s, z

  do ielem = 1, neltot
     ibd = 0
     do jpol = 0, npol, npol
        do ipol = 0, npol, npol
           s=sgll(ipol,jpol,ielem)
           z=zgll(ipol,jpol,ielem)
           call check_boundary(itest,s,z,dmax,rbound)
           ibd = max(itest,ibd)
        enddo
     enddo
     if (ibd == 1) nbelem(j) = nbelem(j) + 1
  enddo
  if (dump_mesh_info_screen) write(*,*)nbelem(j),'boundary elements for interface',j

end subroutine belem_count_new
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine belem_list_new(dmax, j, rbound)
! This routine records the indices for boundary elements
! for the j-th component of the velocity field, in the spheroidal
! container case.

  integer, intent(in)           :: j
  real(kind=dp)   , intent(in)  :: dmax, rbound
  integer                       :: ibd, ielem, ipol, jpol, itest, icount
  real(kind=dp)                 :: s, z

  icount = 0
  do ielem = 1, neltot
     ibd = 0
     do jpol = 0, npol, npol
        do ipol = 0, npol, npol
           s=sgll(ipol,jpol,ielem)
           z=zgll(ipol,jpol,ielem)
           call check_boundary(itest,s,z,dmax,rbound)
           ibd = max(itest,ibd)
        enddo
     enddo
     if (ibd == 1) then
        icount = icount + 1
        belem(icount,j) = ielem
     endif
  enddo
end subroutine belem_list_new
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine check_boundary(itest, s, z, dmax, rbound)

  integer, intent(out) :: itest
  real(kind=dp)   , intent(in) :: s,z,dmax,rbound

  itest = 0

  if (dabs(dsqrt(s**2+z**2)-rbound) < dmax) itest = 1

end subroutine check_boundary
!-----------------------------------------------------------------------------------------

end module mesh_info
!=========================================================================================
