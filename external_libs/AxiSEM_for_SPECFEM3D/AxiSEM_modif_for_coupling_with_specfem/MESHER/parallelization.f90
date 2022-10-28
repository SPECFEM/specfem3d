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
module parallelization

  use data_grid
  use data_gllmesh
  use data_mesh
  use data_spec
  use data_pdb
  use data_diag
  use data_bkgrdmodel

  implicit none

  public :: create_domain_decomposition

  private

contains

!-----------------------------------------------------------------------------------------
subroutine create_domain_decomposition
! nel:    number of glocal elements, i.e. total number of a processor's elements
! nelmax: maximal number of glocal elements
! neltot: global total number of elements
! procel:

  integer :: iproc, iel, nelmax, nelmax_fluid, nelmax_solid
  integer :: itheta, irad
  integer :: nelsolid_per_tsl, nelsolid_per_trsl
  integer :: nelfluid_per_tsl, nelfluid_per_trsl
  logical :: attributed(1:neltot)

  write(*,*)'     creating domain decomposition....'

  ! check if number of elements at ICB is multiple of nthetaslices
  call check_nproc(nthetaslices)
  nproc = nthetaslices * nradialslices

  attributed(:) = .false.

  allocate(nel(0:nproc-1))
  allocate(nel_fluid(0:nproc-1))
  allocate(nel_solid(0:nproc-1))

  ! This must have been always the case as long as the two test statements
  ! above where passed and neltot = neltot_solid + neltot_fluid
  if (mod(neltot,nproc) == 0 .and. mod(neltot_solid,nproc) == 0 &
        .and. mod(neltot_fluid,nproc) == 0) then
     nel(:) = neltot / nproc
     nel_fluid(:) = neltot_fluid / nproc
     nel_solid(:) = neltot_solid / nproc
  else
     ! exact load balancing not possible. This scheme attributes the same number
     ! of elements +/- 1 to all processors. The element numbers are the same in
     ! each thetaslice
     if (mod(neltot,nthetaslices) /= 0 .or. mod(neltot_solid,nthetaslices) /= 0 &
           .or. mod(neltot_fluid,nthetaslices) /= 0) then
        write(*,*) 'ERROR: programming error, should not end up here'
        write(*,*) '       - each theta slice should have same number of elements!'
        stop
     endif

     nelsolid_per_tsl = neltot_solid / nthetaslices
     nelfluid_per_tsl = neltot_fluid / nthetaslices

     nelsolid_per_trsl = nelsolid_per_tsl / nradialslices
     nelfluid_per_trsl = nelfluid_per_tsl / nradialslices

     nel_solid(:) = nelsolid_per_trsl
     nel_fluid(:) = nelfluid_per_trsl

     do itheta = 0, nthetaslices-1
        do irad = 0, nradialslices-1
           iproc = itheta * nradialslices + irad
           if (sum(nel_solid(itheta*nradialslices:(itheta+1)*nradialslices-1)) &
 == nelsolid_per_tsl) exit
           nel_solid(iproc) = nel_solid(iproc) + 1
        enddo
     enddo

     do itheta = 0, nthetaslices-1
        do irad = 0, nradialslices-1
           iproc = itheta * nradialslices + irad
           if (sum(nel_fluid(itheta*nradialslices:(itheta+1)*nradialslices-1)) &
 == nelfluid_per_tsl) exit
           nel_fluid(iproc) = nel_fluid(iproc) + 1
        enddo
     enddo

     nel(:) = nel_solid(:) + nel_fluid(:)
  endif


  nelmax = maxval(nel)
  nelmax_solid = maxval(nel_solid)
  nelmax_fluid = maxval(nel_fluid)

  if (dump_mesh_info_screen) then

     write(*,*)
     write(*,*)'************** VARIOUS NUMBERS OF ELEMENTS**********************'
     do iproc=0, nproc-1
        write(*,13) iproc, 'has tot/sol/flu number       =', &
                      nel(iproc), nel_solid(iproc), nel_fluid(iproc)
     enddo
     write(*,12) 'Maximal glocal total number nelmax       =', &
                  nelmax,nelmax_solid,nelmax_fluid
     write(*,*)
     write(*,14) 'Sum over total elements for all procs    =', sum(nel)
     write(*,14) 'Global, total number neltot              =', neltot
     write(*,14) 'Sum over solid elements for all procs    =', sum(nel_solid)
     write(*,14) 'Global, solid-domain number neltot_solid =', neltot_solid
     write(*,14) 'Sum over fluid elements for all procs    =', sum(nel_fluid)
     write(*,14) 'Global, fluid-domain number neltot_fluid =', neltot_fluid

13   format(i4,a37,3(i8))
12   format(a41,3(i8))
14   format(a41,i8)
     write(*,*)'*****************************************************************'
     write(*,*)

  endif

  allocate(procel(nelmax,0:nproc-1))

  allocate(procel_fluid(nelmax_fluid,0:nproc-1))
  allocate(procel_solid(nelmax_solid,0:nproc-1))
  allocate(el2proc(neltot))
  allocate(inv_procel(neltot))

  procel_fluid = -1
  procel_solid = -1
  el2proc      = -1
  inv_procel   = -1

  ! Decompose such that each processor owns a cake piece in the theta direction,
  ! i.e. same amount of elements in solid and fluid respectively.
  ! The inner cube is done such that each processor maximally has 2 neighbors.
  if (nradialslices == 1) then
     call domain_decomposition_theta(attributed, nproc)
  else
     call domain_decomposition_theta_r(attributed, nproc, nthetaslices, nradialslices, &
                                       nelmax, nelmax_fluid, nelmax_solid)
  endif

  ! write out procel arrays
  if (dump_mesh_info_files) then
     open(unit=666,file=diagpath(1:lfdiag)//'/inv_procel.dat')
     do iproc=0, nproc-1
       do iel=1, neltot
          if (el2proc(iel) == iproc) then
             write(666,*) iproc, iel, inv_procel(iel), &
                          procel(inv_procel(iel),iproc)
          endif
       enddo
     enddo
     close(666)

     open(unit=666,file=diagpath(1:lfdiag)//'/procel.dat')
     do iproc=0, nproc-1
       do iel=1, nel(iproc)
          write(666,*) iproc, iel, procel(iel,iproc), inv_procel(procel(iel,iproc))
       enddo
     enddo
     close(666)
  endif

  ! check that every element has been assigned
  if (any(.not. attributed)) then
     do iel = 1, neltot
         if (.not. attributed(iel) ) then
            write(*,*) ' NOT ATTRIBUTED '
            write(*,*) iel, thetacom(iel), solid(iel), fluid(iel)
         endif
     enddo
     stop
  endif

  if (minval(el2proc) == -1) then
     write(*,*) ' '
     write(*,*) 'Element(s) not assigned to any processor:', minloc(el2proc)
     stop
  endif

  if (dump_mesh_info_screen) then
     write(*,*)
     write(*,*) 'NUMBER OF ELEMENTS IN EACH SUBDOMAIN:'
     do iproc=0, nproc-1
        write(*,'("Proc ",i3, " has ",i8, " solid,",i6," fluid,",i9," total elements")') &
                iproc, nel_solid(iproc), nel_fluid(iproc), nel(iproc)
     enddo
     write(*,*)
     call flush(6)
  endif

  if (dump_mesh_vtk) call plot_dd_vtk

end subroutine create_domain_decomposition
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine plot_dd_vtk
  use test_bkgrdmodel, only: write_VTK_bin_scal_old, write_VTK_bin_scal

  real(kind=realkind), dimension(:), allocatable :: wel2proc
  real(kind=sp),   allocatable, dimension(:,:)   :: mesh1
  integer                                        :: iel
  character(len=200)                             :: fname

  integer                                        :: ct
  real(kind=sp), allocatable, dimension(:)       :: x, y, z

  ! write VTK with point data

  allocate(mesh1(neltot,2))
  do iel=1, neltot
      mesh1(iel,1) = real(sgll(npol/2,npol/2,iel))
      mesh1(iel,2) = real(zgll(npol/2,npol/2,iel))
  enddo

  fname = trim(diagpath)//'/mesh_domaindecomposition'
  allocate(wel2proc(1:neltot))
  wel2proc = real(el2proc)
  call write_VTK_bin_scal_old(wel2proc, mesh1, neltot, fname)
  deallocate(wel2proc)

  ! write VTK with cell data

  allocate(wel2proc(1:neltot*4))

  fname = trim(diagpath)//'/mesh_domaindecomposition_cell'

  allocate(x(neltot*4), y(neltot*4), z(neltot*4))
  z = 0.d0
  ct = 0

  do iel=1, neltot
      x(ct+1) = sgll(0,0,iel)
      x(ct+2) = sgll(npol,0,iel)
      x(ct+3) = sgll(npol,npol,iel)
      x(ct+4) = sgll(0,npol,iel)
      y(ct+1) = zgll(0,0,iel)
      y(ct+2) = zgll(npol,0,iel)
      y(ct+3) = zgll(npol,npol,iel)
      y(ct+4) = zgll(0,npol,iel)
      wel2proc(ct+1) = real(el2proc(iel))
      wel2proc(ct+2) = real(el2proc(iel))
      wel2proc(ct+3) = real(el2proc(iel))
      wel2proc(ct+4) = real(el2proc(iel))

      ct = ct + 4
  enddo

  call write_VTK_bin_scal(x, y, z, wel2proc, neltot, fname)

  deallocate(x, y, z)
  deallocate(wel2proc)
  deallocate(mesh1)

end subroutine plot_dd_vtk
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine domain_decomposition_theta(attributed, nprocl)
! nel:    number of glocal elements, i.e. total number of a processor's elements
! nelmax: maximal number of glocal elements
! neltot: global total number of elements
! procel:

  logical, intent(inout)    :: attributed(:)
  integer, intent(in)       :: nprocl

  integer                   :: iproc, iiproc, iel
  integer                   :: mycount
  integer, allocatable      :: central_count(:)
  real(kind=dp)             :: pi2

  pi2 = two * dasin(one)

  write(*,*)'     THETA-SLICING as domain decomposition....'

  allocate(central_count(0:nprocl-1))

  ! Create colatitude bounds array for outer shell
  ! theta_min_proc and theta_max_proc are now filled up in
  ! meshgen.f90:def_ref_cart_coordinates_discont

  ! **************** INNER CUBE **********************
  if (nprocl == 1 .or. nprocl == 2) then
      ! define quadratic functions to delineate processor boundaries.
      ! Works for nprocl = 1, 2
      call decompose_inner_cube_quadratic_fcts(central_count, attributed, nprocl, &
                                               procel_solid, procel_fluid)

  else if (nprocl >= 4 .and. (nprocl / 4) * 4 == nprocl) then
      ! newest version of inner core decomposition (nprocl needs to be multiple of 4)
      call decompose_inner_cube_opt(central_count, attributed, nprocl, &
                                    procel_solid, procel_fluid)
  else
      write(*,*)
      write(*,*)'PROBLEM: central cube decomposition not implemented for nprocl = ',nprocl
      write(*,*)'         nprocl should be in 1, 2 or a multiple of 4!'
      stop
  endif

  ! **************** END OF INNER CUBE****************

  do iproc = 0, nprocl -1

     if (solid_domain(ndisc)) then
        mycount = 0
     else
        mycount = central_count(iproc)
     endif

     do iel = 1, neltot

        ! add the extra requirement that element iel to be in appropriate theta slice
        if (fluid(iel) .and. .not. attributed(iel) .and. &
            (thetacom(iel) >= theta_min_proc(iproc)) .and. &
            (thetacom(iel) <= theta_max_proc(iproc)) ) then
            mycount = mycount + 1
            procel_fluid(mycount,iproc) = iel
            attributed(iel) = .true.
        endif

        if ( mycount == nel_fluid(iproc) ) exit
     enddo ! iel

     if (solid_domain(ndisc)) then
        mycount = central_count(iproc)
     else
        mycount = 0
     endif

     !  At this stage we have assigned nel_fluid fluid elements
     !  to each processor, stored in a procel_fluid(1:nel_fluid,iproc)
     !  Here we start the loop over solid elements and try to assign them to iproc

     do iel = 1, neltot
        if (.not. fluid(iel) .and. .not. attributed(iel) .and. &
             (thetacom(iel) >= theta_min_proc(iproc)) .and. &
             (thetacom(iel) <= theta_max_proc(iproc)) ) then
           mycount = mycount + 1
           procel_solid(mycount,iproc) = iel
           attributed(iel) = .true.
        endif

        if ( mycount == nel_solid(iproc) ) then
           if (dump_mesh_info_screen) then
              write(*,*) ' PROC ', iproc ,' has everybody it needs ', mycount, &
                         nel_solid(iproc)
              call flush(6)
           endif
           exit
        endif
     enddo

     if (mycount < nel_solid(iproc)) then
        write(*,*)
        write(*,*) 'Problem: not all solid elements attributed for proc', iproc, &
                    mycount, nel_solid(iproc)
        do iiproc=0, nprocl-1
           write(*,*) 'nel_solid(iproc), centralcount:', &
                     iiproc, nel_solid(iiproc), central_count(iiproc)
        enddo
        stop
     endif

     ! procel contains
     ! procel(1:nel_fluid) : the nel_fluid element numbers pertaining to iproc
     ! procel(nel_fluid+1:nel(iproc)) : the nel_solid solid element numbers
     ! belonging to iproc
     ! Element numbers are defined in a global sense (solid+fluid whole mesh)
     do iel = 1, nel(iproc)
        if (iel <= nel_fluid(iproc) ) then
           procel(iel,iproc) = procel_fluid(iel,iproc)
        else
           procel(iel,iproc) = procel_solid(iel-nel_fluid(iproc),iproc)
           if (procel(iel,iproc) <= 0) then
              write(*,*) 'PROCEL ZERO!', iproc, eltypeg(iel), iel, nel_fluid(iproc)
              write(*,*) '            ', thetacom(iel), pi/nprocl * iproc, &
                         pi/nprocl * (iproc + 1)
              stop
           endif
        endif
        el2proc(procel(iel,iproc)) = iproc
        inv_procel(procel(iel,iproc)) = iel
     enddo

  enddo !nprocl-1

  ! OUTPUT OF DOMAIN DECOMPOSITION

  if (dump_mesh_info_files) then

     write(*,*)'Writing out the domain decomposition...'; call flush(6)

     ! central dd only, but checking through whole grid
     open(unit=647,file=diagpath(1:lfdiag)//'/dd_central_sz_ielglob_iproc.dat')
     do iel=1,neltot
        if (eltypeg(iel) == 'linear') write(647,14)scom(iel),zcom(iel),iel, &
                                                 el2proc(iel)
     enddo
     close(647)
14    format(2(1pe13.3),2(i8))

     ! THIS ONE WILL BE KINDA LARGE (but still only one per element)
     open(unit=648,file=diagpath(1:lfdiag)//'/dd_sz_ielglob_iproc.dat')
     do iel=1,neltot
        write(648,14) scom(iel), zcom(iel), iel, el2proc(iel)
     enddo
     close(648)

  endif

end subroutine domain_decomposition_theta
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine domain_decomposition_theta_r(attributed, nprocl, nthetal, nrl, &
                                        nelmax, nelmax_fluid, nelmax_solid)
! nel:     number of glocal elements, i.e. total number of a processor's elements
! nelmax: maximal number of glocal elements
! neltot: global total number of elements
! procel:

  use sorting

  logical, intent(inout)    :: attributed(:)
  integer, intent(in)       :: nprocl, nthetal, nrl
  integer, intent(in)       :: nelmax, nelmax_fluid, nelmax_solid

  integer                   :: itheta, iitheta, iel, nel_fluid_theta
  integer                   :: irad, iproc
  integer                   :: mycount, nicb, ncmb
  integer                   :: iprocb(2), mycountb(2), j1, j2
  real(kind=dp)             :: deltatheta
  integer, allocatable      :: central_count(:)
  real(kind=dp)             :: pi2

  integer, allocatable          :: thetaslel(:,:), thetaslel_fluid(:,:), thetaslel_solid(:,:)
  integer, allocatable          :: inner_core_buf(:)
  real(kind=dp), allocatable    :: inner_core_r(:)
  integer, allocatable          :: el2thetaslel(:)

  allocate(thetaslel(nelmax*nrl,0:nthetal-1))
  allocate(thetaslel_fluid(nelmax_fluid*nrl,0:nthetal-1))
  allocate(thetaslel_solid(nelmax_solid*nrl,0:nthetal-1))
  allocate(el2thetaslel(neltot))

  thetaslel       = -1
  thetaslel_fluid = -1
  thetaslel_solid = -1
  el2thetaslel    = -1

  pi2 = two * dasin(one)

  write(*,*)'     THETA+RADIAL-SLICING as domain decomposition....'

  if (nbcnd > 2) then
     write(*,*) 'ERROR: radial slicing only implemented for a single fluid layer'
     write(*,*) '       workaround: set NRADIAL_SLICES to 1'
     stop
  endif

  allocate(central_count(0:nthetal-1))

  ! Create colatitude bounds array for outer shell
  ! theta_min_proc and theta_max_proc are now filled up in
  ! meshgen.f90:def_ref_cart_coordinates_discont

  ! **************** INNER CUBE **********************
  if (nthetal == 1 .or. nthetal == 2) then
      ! define quadratic functions to delineate processor boundaries.
      ! Works for nthetal = 1, 2
      call decompose_inner_cube_quadratic_fcts(central_count, attributed, nthetal, &
                                               thetaslel_solid, thetaslel_fluid)

  else if (nthetal >= 4 .and. (nthetal / 4) * 4 == nthetal) then
      ! newest version of inner core decomposition (nthetal needs to be multiple of 4)
      call decompose_inner_cube_opt(central_count, attributed, nthetal, &
                                    thetaslel_solid, thetaslel_fluid)
  else
      write(*,*)
      write(*,*) 'PROBLEM: central cube decomposition not implemented for nthetal = ', nthetal
      write(*,*) '         nthetal should be in 1, 2 or a multiple of 4!'
      stop
  endif
  ! **************** END OF INNER CUBE****************



  ! add the extra requirement that element iel to be in appropriate theta slice
  do itheta = 0, nthetal-1

     if (solid_domain(ndisc)) then
        mycount = 0
     else
        mycount = central_count(itheta)
     endif

     do iel = 1, neltot
        if (fluid(iel) .and. .not. attributed(iel) .and. &
               (thetacom(iel) >= theta_min_proc(itheta)) .and. &
               (thetacom(iel) <= theta_max_proc(itheta)) ) then
            mycount = mycount + 1
            thetaslel_fluid(mycount,itheta) = iel
            attributed(iel) = .true.
        endif
        if ( mycount == sum(nel_fluid(itheta:itheta+nrl-1)) ) then
           write(*,*) 'all fluid elements assigned'
           exit
        endif
     enddo ! iel

     if (neltot_solid > 0 ) then
        if (solid_domain(ndisc)) then
            mycount = central_count(itheta)
        else
            mycount = 0
        endif

        do iel = 1, neltot

           if (.not. fluid(iel) .and. .not. attributed(iel) .and. &
                 (thetacom(iel) >= theta_min_proc(itheta)) .and. &
                 (thetacom(iel) <= theta_max_proc(itheta)) ) then
              if (solid_domain(ndisc)) then
                 thetaslel_solid(sum(nel_solid(itheta:itheta+nrl-1)) &
                                 - mycount + central_count(itheta),itheta) = iel
              else
                 thetaslel_solid(sum(nel_solid(itheta:itheta+nrl-1)) &
                                 - mycount,itheta) = iel
              endif

              mycount = mycount + 1
              attributed(iel) = .true.
           endif

           if ( mycount == sum(nel_solid(itheta:itheta+nrl-1)) ) then
              if (dump_mesh_info_screen) then
                 write(*,*) ' THETASL ', itheta ,' has everybody it needs ', mycount, &
                            sum(nel_solid(itheta:itheta+nrl-1))
                 call flush(6)
              endif
              exit
           endif
        enddo
     else
        mycount = 0
     endif

     if (mycount < sum(nel_solid(itheta:itheta+nrl-1))) then
        write(*,*)
        write(*,*) 'Problem: not all solid elements attributed for thetaslice', itheta, &
                    mycount, sum(nel_solid(itheta:itheta+nrl-1))
        do iitheta=0, nthetal-1
           write(*,*) 'nel_solid(itheta), centralcount:', &
                     iitheta, sum(nel_solid(iitheta:itheta+nrl-1)), central_count(iitheta)
        enddo
        stop
     else if (mycount > sum(nel_solid(itheta:itheta+nrl-1))) then
        write(*,*)
        write(*,*) 'Problem: too many solid elements attributed for thetaslice', itheta, &
                    mycount, sum(nel_solid(itheta:itheta+nrl-1))
        do iitheta=0, nthetal-1
           write(*,*) 'nel_solid(itheta), centralcount:', &
                     iitheta, sum(nel_solid(iitheta:itheta+nrl-1)), central_count(iitheta)
        enddo
        stop
     endif

     ! thetaslel contains
     ! thetaslel(1:nel_fluid) : the nel_fluid element numbers pertaining to itheta
     ! thetaslel(nel_fluid+1:nel(itheta)) : the nel_solid solid element numbers
     ! belonging to itheta
     ! Element numbers are defined in a global sense (solid+fluid whole mesh)
     do iel = 1, sum(nel(itheta:itheta+nrl-1))
        if (iel <= sum(nel_fluid(itheta:itheta+nrl-1))) then
           thetaslel(iel,itheta) = thetaslel_fluid(iel,itheta)
        else
           thetaslel(iel,itheta) = thetaslel_solid(iel-sum(nel_fluid(itheta:itheta+nrl-1)),itheta)
           if (thetaslel(iel,itheta) <= 0) then
              write(*,*) 'PROCEL ZERO!', itheta, eltypeg(iel), iel, nel_fluid(itheta)
              write(*,*) '            ', thetacom(iel), pi/nthetal * itheta, &
                         pi/nthetal * (itheta + 1)
              write(*,*) iel, sum(nel(itheta:itheta+nrl-1))
              stop
           endif
        endif
        el2thetaslel(thetaslel(iel,itheta)) = itheta
     enddo

  enddo !nthetal-1

  if (any(.not. attributed)) then
     write(*,*) 'ERROR: not all elements assigned to a theta-slice'
     stop
  endif


  ! sort inner core elements according to radius
  ! Using same stupid choice as in inner core decomposition

  if (solid_domain(ndisc)) then
     do itheta = 0, nthetal-1
        allocate(inner_core_buf(central_count(itheta)))
        allocate(inner_core_r(central_count(itheta)))
        inner_core_buf(:) = thetaslel_solid(1:central_count(itheta),itheta)

        do iel = 1, central_count(itheta)
           inner_core_r(iel) = rcom(inner_core_buf(iel))
        enddo

        call mergesort_3(inner_core_r, il=inner_core_buf, p=4)

        thetaslel_solid(1:central_count(itheta),itheta) = inner_core_buf(:)

        deallocate(inner_core_buf)
        deallocate(inner_core_r)
     enddo
  else
     do itheta = 0, nthetal-1
        nel_fluid_theta = sum(nel_fluid(itheta:itheta+nrl-1))
        allocate(inner_core_buf(nel_fluid_theta))
        allocate(inner_core_r(nel_fluid_theta))
        inner_core_buf(:) = thetaslel_fluid(1:nel_fluid_theta,itheta)

        do iel = 1, nel_fluid_theta
           ! sort by radius, if radius is the same, theta makes the difference
           inner_core_r(iel) = rcom(inner_core_buf(iel)) + 1e-10 * thetacom(inner_core_buf(iel))
        enddo

        call mergesort_3(inner_core_r, il=inner_core_buf, p=4)

        thetaslel_fluid(1:nel_fluid_theta,itheta) = inner_core_buf(:)

        deallocate(inner_core_buf)
        deallocate(inner_core_r)
     enddo
  endif

  ! reset, as we have to touch each element again!
  attributed = .false.

  ! Now decomposition in radius
  do itheta = 0, nthetal-1

     ! SOLID domain first
     do irad = 0, nrl-1
        iproc = itheta * nrl + irad
        mycount = 1
        do iel = 1, neltot_solid / nthetal
           if (mycount == nel_solid(iproc) + 1) exit
           if (.not. attributed(thetaslel_solid(iel,itheta))) then
              procel_solid(mycount, iproc) = thetaslel_solid(iel,itheta)
              attributed(thetaslel_solid(iel,itheta)) = .true.
              mycount = mycount + 1
           endif
        enddo
     enddo

     ! fill solid mapping arrays (used in the fluid decomposition)
     do irad = 0, nrl-1
        iproc = itheta * nrl + irad

        do iel = nel_fluid(iproc)+1, nel(iproc)
           procel(iel,iproc) = procel_solid(iel-nel_fluid(iproc),iproc)
           if (procel(iel,iproc) <= 0) then
              write(*,*) 'PROCEL ZERO!', iproc, eltypeg(iel), iel, nel_fluid(iproc)
              write(*,*) '            ', thetacom(iel), pi/nprocl * iproc, &
                         pi/nprocl * (iproc + 1)
              write(*,*) iel, nel(iproc)
              stop
           endif
           el2proc(procel(iel,iproc)) = iproc
           inv_procel(procel(iel,iproc)) = iel
        enddo
     enddo

     ! FLUID
     ! special treatment of the processors with sol/flu boundary
     ! kind of hacky, but should work for most cases of earth like models
     iprocb = -1
     mycountb = 1

     if (nbcnd == 1) then
        nicb = 0

        ! take the upper most layer in the fluid (CMB) and account to the
        ! processor having the solid neighbor
        do iel = 1, nbelem(1)
           if (fluid(belem(iel,1)) .and. el2thetaslel(belem(iel,1)) == itheta &
                   .and. .not. attributed(belem(iel,1))) then
               iproc = el2proc(belem(my_neighbor(iel,1),1))
               if (iprocb(1) == iproc .or. iprocb(1) == -1) then
                  iprocb(1) = iproc
                  procel_fluid(mycountb(1), iproc) = belem(iel,1)
                  mycountb(1) = mycountb(1) + 1
               else if (iprocb(2) == iproc .or. iprocb(2) == -1) then
                  iprocb(2) = iproc
                  procel_fluid(mycountb(2), iproc) = belem(iel,1)
                  mycountb(2) = mycountb(2) + 1
               else
                  write(*,*) 'ERROR: more then two procs involved in s/f boundary of one theta slice'
               endif
               attributed(belem(iel,1)) = .true.
           endif
        enddo
        ncmb = mycountb(1) + mycountb(2) - 2 - nicb
     else if (nbcnd == 2) then

        ! take the lower most layer in the fluid (ICB) and account to the
        ! processor having the solid neighbor
        do iel = 1, nbelem(2)
           if (fluid(belem(iel,2)) .and. el2thetaslel(belem(iel,2)) == itheta &
                   .and. .not. attributed(belem(iel,2))) then
               iproc = el2proc(belem(my_neighbor(iel,2),2))
               if (iprocb(1) == iproc .or. iprocb(1) == -1) then
                  iprocb(1) = iproc
                  procel_fluid(mycountb(1), iproc) = belem(iel,2)
                  mycountb(1) = mycountb(1) + 1
               else if (iprocb(2) == iproc .or. iprocb(2) == -1) then
                  iprocb(2) = iproc
                  procel_fluid(mycountb(2), iproc) = belem(iel,2)
                  mycountb(2) = mycountb(2) + 1
               else
                  write(*,*) 'ERROR: more then two procs involved in s/f boundary of one theta slice'
               endif
               attributed(belem(iel,2)) = .true.
           endif
        enddo

        nicb = mycountb(1) + mycountb(2) - 2

        ! take the upper most layer in the fluid (CMB) and account to the
        ! processor having the solid neighbor
        do iel = 1, nbelem(1)
           if (fluid(belem(iel,1)) .and. el2thetaslel(belem(iel,1)) == itheta &
                   .and. .not. attributed(belem(iel,1))) then
               iproc = el2proc(belem(my_neighbor(iel,1),1))
               if (iprocb(1) == iproc .or. iprocb(1) == -1) then
                  iprocb(1) = iproc
                  procel_fluid(mycountb(1), iproc) = belem(iel,1)
                  mycountb(1) = mycountb(1) + 1
               else if (iprocb(2) == iproc .or. iprocb(2) == -1) then
                  iprocb(2) = iproc
                  procel_fluid(mycountb(2), iproc) = belem(iel,1)
                  mycountb(2) = mycountb(2) + 1
               else
                  write(*,*) 'ERROR: more then two procs involved in s/f boundary of one theta slice'
               endif
               attributed(belem(iel,1)) = .true.
           endif
        enddo
        ncmb = mycountb(1) + mycountb(2) - 2 - nicb
     else if (nbcnd == 0) then
        nicb = 0
        ncmb = 0
     else
        write(*,*) 'domain decomposition only implemented for 0 or 2 solid fluid boundaries'
        stop
     endif


     if (any(mycountb - 1 > nel_fluid(iproc))) then
        write(*,*) 'ERROR: more boundary elements than fluid elements. try less NRADIAL_SLICES'
        stop
     endif

     ! depending on whether a possible proc boundary on the sf boundary (on the
     ! solid side) is on CMB or ICB decide which processor should take stuff
     ! from below CMB / above ICB
     ! this choice leads to one processor having just one contiguous domain and
     ! the other having two domains, one of which is very small (just boundary
     ! elements)
     if (mycountb(2) - 1 > nicb) then
        j1 = 2
        j2 = 1
     else
        j1 = 1
        j2 = 2
     endif

     ! fill up j1 with stuff from below CMB
     if (iprocb(j1) /= -1) then
        do iel = 1, neltot_fluid / nthetal
           if (mycountb(j1) == nel_fluid(iproc) + 1) exit
           if (.not. attributed(thetaslel_fluid(iel,itheta))) then
              procel_fluid(mycountb(j1), iprocb(j1)) = thetaslel_fluid(iel,itheta)
              attributed(thetaslel_fluid(iel,itheta)) = .true.
              mycountb(j1) = mycountb(j1) + 1
           endif
        enddo
     endif

     ! fill up j2 with stuff from above ICB
     if (iprocb(j2) /= -1) then
        do iel = neltot_fluid / nthetal, 1, -1
           if (mycountb(j2) == nel_fluid(iproc) + 1) exit
           if (.not. attributed(thetaslel_fluid(iel,itheta))) then
              procel_fluid(mycountb(j2), iprocb(j2)) = thetaslel_fluid(iel,itheta)
              attributed(thetaslel_fluid(iel,itheta)) = .true.
              mycountb(j2) = mycountb(j2) + 1
           endif
        enddo
     endif

     ! fill up the bulk with the other processors
     do irad = 0, nrl-1
        iproc = itheta * nrl + irad
        if (any(iproc == iprocb)) cycle

        ! just go inwards radially
        mycount = 1
        do iel = 1, neltot_fluid / nthetal
           if (mycount == nel_fluid(iproc) + 1) exit
           if (.not. attributed(thetaslel_fluid(iel,itheta))) then
              procel_fluid(mycount, iproc) = thetaslel_fluid(iel,itheta)
              attributed(thetaslel_fluid(iel,itheta)) = .true.
              mycount = mycount + 1
           endif
        enddo
     enddo

     ! fill up fluid part of the mapping arrays
     do irad = 0, nrl-1
        iproc = itheta * nrl + irad
        do iel = 1, nel_fluid(iproc)
           procel(iel,iproc) = procel_fluid(iel,iproc)
           el2proc(procel(iel,iproc)) = iproc
           inv_procel(procel(iel,iproc)) = iel
        enddo
     enddo

  enddo ! itheta

end subroutine domain_decomposition_theta_r
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine decompose_inner_cube_quadratic_fcts(central_count, attributed, nthetal, &
                                               procel_solidl, procel_fluidl)

  integer, intent(in)       :: nthetal
  integer, intent(out)      :: central_count(0:nthetal-1)
  logical, intent(inout)    :: attributed(:)
  integer, intent(out)      :: procel_solidl(:,0:), procel_fluidl(:,0:)

  integer :: iproc, is, iz, nthetal2
  integer,allocatable :: proc_central(:,:),num_columns(:),upper_boundary_el(:)
  integer,allocatable :: num_columns_hi(:),num_columns_lo(:),num_el(:)
  integer,allocatable :: count_assi(:)

  if (dump_mesh_info_screen) then
     write(*,*)
     write(*,*)'<> < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > <>'
     write(*,*)'CENTRAL LINEAR DOMAIN: simple decomposition!'
     write(*,*)'ndivs,nthetal:',ndivs,nthetal
     write(*,*)' = => each processor should have el=',ndivs**2/nthetal*2
     write(*,*)'<> < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > <>'
     write(*,*)
  endif

  if (nthetal > 2) then
     write(*,*) 'simple mesh decomposition cannot handle more than 2 processors'
     write(*,*) 'use optimization scheme instead'
     stop
  endif

  nthetal2 = nthetal / 2 - 1
  if (nthetal == 1) nthetal2 = 0

  allocate(proc_central(1:ndivs,1:ndivs))
  proc_central(1:ndivs,1:ndivs)=-1

  allocate(num_columns(1:ndivs))
  allocate(upper_boundary_el(1:ndivs))
  allocate(num_columns_hi(1:ndivs),num_columns_lo(1:ndivs))
  allocate(num_el(0:nthetal2))

  proc_central(1:ndivs,1:ndivs) = 0
  num_el(0) = ndivs**2

  allocate(count_assi(0:nthetal2))
  count_assi = 0
  ! count respective processors elements
  num_el(0:nthetal2) = 0

  do is=1, ndivs
     do iz=1, ndivs
        do iproc=0, nthetal2
           if (proc_central(is,iz) == iproc) then
              num_el(iproc) = num_el(iproc) + 1
              count_assi(iproc) = count_assi(iproc) + 1
           endif
        enddo
        if (proc_central(is,iz) < 0 .or. proc_central(is,iz) > nthetal2) then
           write(*,*) 'Problem:', is, iz, 'has no processor!'
           stop
        endif
     enddo
  enddo

  do iproc=0, nthetal2
    if (count_assi(iproc) /= ndivs**2/(nthetal2+1)) then
       write(*,*) 'Problem: Not every element is assigned to processor', iproc
       write(*,*) 'Counted assigned els/total els:', count_assi(iproc), ndivs**2/(nthetal2+1)
       if (iproc < nthetal2) &
            write(*,*) 'els for other procs:', count_assi(iproc+1:nthetal2)
       stop
    endif
  enddo

  if (dump_mesh_info_screen) then
     write(*,*)
     do iproc=0, nthetal2
        write(*,12) iproc,num_el(iproc)
     enddo
12   format('Central cube: proc', i3, ' has', i6, ' elements')
  endif

  ! connect these processor dependencies to the global element numbering scheme
  central_count(0:nthetal-1) = 0
  do iz = 1, ndivs
   do is = 1, ndivs
      if (proc_central(is,iz) /= -1) then
         central_count(proc_central(is,iz)) = central_count(proc_central(is,iz)) + 1

         attributed(central_is_iz_to_globiel(is,iz)) = .true.
         attributed(central_is_iz_to_globiel(is,iz) + neltot / 2)=.true.
      else
         write(*,*) 'Unassigned element in the central cube!'
         write(*,*) 'is,iz:',is,iz
         stop
      endif

     if (solid_domain(ndisc)) then
        procel_solidl(central_count(proc_central(is,iz)),proc_central(is,iz)) = &
                                                central_is_iz_to_globiel(is,iz)
     else
        procel_fluidl(central_count(proc_central(is,iz)),proc_central(is,iz)) = &
                                                central_is_iz_to_globiel(is,iz)
     endif

     ! South: inverted copy
     if (nthetal > 1) then
        if (solid_domain(ndisc)) then
           procel_solidl(central_count(proc_central(is,iz)), &
                nthetal-1-proc_central(is,iz)) = &
                central_is_iz_to_globiel(is,iz) + neltot / 2
        else
           procel_fluidl(central_count(proc_central(is,iz)), &
                nthetal-1-proc_central(is,iz)) = &
                central_is_iz_to_globiel(is,iz) + neltot / 2
        endif
     endif
    enddo
  enddo

  ! South:
  if (nthetal > 1) then
     do iproc=0, nthetal / 2 - 1
        central_count(nthetal-iproc-1) = central_count(iproc)
     enddo
  endif

  ! special case one processor... still needs to count the south!
  if (nthetal == 1) then
    do iz = 1, ndivs
     do is = 1, ndivs
          central_count(proc_central(is,iz)) = central_count(proc_central(is,iz)) + 1
          if (solid_domain(ndisc)) then
             procel_solidl(central_count(proc_central(is,iz)),0) = &
                  central_is_iz_to_globiel(is,iz) + neltot / 2
          else
             procel_fluidl(central_count(proc_central(is,iz)),0) = &
                  central_is_iz_to_globiel(is,iz) + neltot / 2
          endif
      enddo
    enddo
  endif

  ! check if all central-cube elements are assigned
  ! MvD: this test assumes all linear elements are in the inner core!
  do is=1, neltot
     if (eltypeg(is) == 'linear' ) then
        if (.not. attributed(is)) then
           write(*,*)
           write(*,*) 'Problem: Central cube element not assigned!', is
           stop
        endif
     endif
  enddo

  ! write out the central cube decomposition
  if (dump_mesh_info_files) then
     open(unit=10008,file=diagpath(1:lfdiag)//'/central_locind_locsz_iproc.dat')
     do iz = 1, ndivs
        do is = 1, ndivs
           write(10008,13) is, iz, s_arr(is,iz), z_arr(is,iz), proc_central(is,iz)
        enddo
     enddo
     do iz = 1, ndivs
        do is = 1, ndivs
           write(10008,13)is, iz, s_arr(is,iz), -z_arr(is,iz), nthetal-1-proc_central(is,iz)
        enddo
     enddo
     close(10008)
13   format(2(i4),2(f9.3),i4)
  endif

  if (dump_mesh_info_screen) then
     write(*,*)
     write(*,*)'<> < > < > < > Finished the central domain decomposition! < > < > < > <>'
     write(*,*)
     call flush(6)
  endif

end subroutine decompose_inner_cube_quadratic_fcts
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine check_my_els(iproc, proc_central)

  integer, intent(in) :: iproc, proc_central(1:ndivs,1:ndivs)
  integer             :: procelcount, is, iz

  ! Check if proc has right amount of elements
  procelcount=0
  do is=1,ndivs
     do iz=1,ndivs
        if (proc_central(is,iz) == iproc) procelcount=procelcount+1
     enddo
  enddo
  if (procelcount /= 2*ndivs**2/nproc) then
     write(*,*)
     write(*,12)iproc,procelcount
     write(*,*)'Needed:',2*ndivs**2/nproc
     stop
  else
     if (dump_mesh_info_screen) then
      write(*,*)
      write(*,*)' >  >>', iproc,' has right amount of elements:',2*ndivs**2/nproc
      call flush(6)
     endif
  endif

12 format('Problem: Processor',i4,' has ',i6,' elements!')

end subroutine check_my_els
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine check_nproc(np)
  use data_coarse
  integer :: np,nb,itest,ritest
  nb = ns_ib ! Number of elements at the inner spherical boundary (defined while
             ! generating the skeleton)
  ritest = dble(2*nb)/dble(np)
  itest = 2*nb/np
  if (ritest /= dble(itest) ) then
     write(*,*) ritest,itest,2*nb,np
     write(*,*) 2*nb,np
     write(*,*) ' Number of processes not compliant with mesh topology '
     stop
  endif
end subroutine check_nproc
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine decompose_inner_cube_opt(central_count, attributed, nthetal, &
                                    procel_solidl, procel_fluidl)

  integer, intent(in)        :: nthetal
  integer, intent(out)       :: central_count(0:nthetal-1)
  logical, intent(inout)     :: attributed(:)
  integer, intent(out)       :: procel_solidl(:,0:), procel_fluidl(:,0:)

  integer                    :: nthetal2, nlinsteps, ndivsppx0, ip, &
                                ncorrections = 0, npart, is, iz, sign_buff, n, &
                                ids, idz
  real(kind=sp)              :: r1, r2, stepsize, dphi
  real(kind=sp), allocatable :: x0(:), x1(:), x2(:), x3(:), z0(:), z1(:), z2(:), &
                                z3(:), phi(:)
  integer, allocatable       :: proc(:,:), nelem(:)
  logical, allocatable       :: proc_iq_min(:,:), proc_iq_max(:,:), elems(:,:)
  logical                    :: exit_buff

  if (dump_mesh_info_screen) then
      write(*,*)
      write(*,*)'<> < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > <>'
      write(*,*)'CENTRAL LINEAR DOMAIN: decomposing using an '
      write(*,*)'      optimization scheme!'
      write(*,*)'ndivs,nthetal:',ndivs,nthetal
      write(*,*)' = => each processor should have el=',ndivs**2/nthetal*2
      write(*,*)'<> < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > < > <>'
      write(*,*)
  endif


  if (mod(2*ndivs**2,nthetal) /= 0) then
     write(*,*)'Central cube area not divisible into equal areas!'
     stop
  endif

  if (mod(4*ndivs,nthetal) /= 0) then
     write(*,*)'PROBLEM with number of central-region elements and nthetal:'
     write(*,*)'ndivs,nthetal:',ndivs,nthetal
     write(*,*)'ndivs (number of northern elements in one direction)'
     write(*,*)'needs to be multiple of nthetal/4...'
     stop
  endif

  ! area well defined?
  if (mod(2*ndivs**2,nthetal) /= 0) then
     write(*,*)'PROBLEM with number of central-region elements and nthetal:'
     write(*,*)'els,nthetal:',ndivs**2,nthetal
     write(*,*)'number of elements needs to be multiple of nthetal/2...'
     stop
  endif

  nthetal2 = nthetal / 2

  nlinsteps = 50

  ndivsppx0 = ndivs / 2 / nthetal2
  if (ndivsppx0 < 2) then
      ndivsppx0 = 2
  endif

  ! treat 4 procs as special case
  if (nthetal2 == 2) then
      ndivsppx0 = 2
  endif

  if (dump_mesh_info_screen) then
      write(*,*) 'number of axis elements per processor = ', ndivsppx0
  endif

  r1 = ndivsppx0 * (nthetal2 - .9)
  r2 = 0.85 * ndivs

  allocate(x0(0:nthetal2-2))
  allocate(x1(0:nthetal2-2))
  allocate(x2(0:nthetal2-2))
  allocate(x3(0:nthetal2-2))
  allocate(z0(0:nthetal2-2))
  allocate(z1(0:nthetal2-2))
  allocate(z2(0:nthetal2-2))
  allocate(z3(0:nthetal2-2))
  allocate(phi(0:nthetal2-2))

  x0 = 0.
  x3 = real(ndivs)
  do ip = nthetal2 / 2, nthetal2 - 2, 1
      x3(ip) = real(ndivs - (ip + 1 - (nthetal2) / 2) * ndivs / (nthetal2 / 2)) + .5
  enddo

  do ip = 0, nthetal2 - 2, 1
      z0(ip) = real(ndivsppx0 * (ip + 1)) + 0.1
      x1(ip) = ndivsppx0 * (nthetal2 / 2 + 1) + (nthetal2 / 2 - ip) * ndivsppx0 * 1.1
      z3(ip) = x3(nthetal2 - 2 - ip)
      phi(ip) = pi / 2. / real(nthetal2) * real(ip + 1)
      x2(ip) = (r2 + ndivs * 0.3 * sin(2. * phi(ip))**4) * cos(phi(ip))
      z2(ip) = (r2 + ndivs * 0.3 * sin(2. * phi(ip))**4) * sin(phi(ip))
  enddo

  z1 = z0

  if (dump_mesh_info_screen) then
      write(*,*) 'x0 = ', x0
      write(*,*) 'x1 = ', x1
      write(*,*) 'x2 = ', x2
      write(*,*) 'x3 = ', x3
      write(*,*) 'z0 = ', z0
      write(*,*) 'z1 = ', z1
      write(*,*) 'z2 = ', z2
      write(*,*) 'z3 = ', z3
  endif

  allocate(proc(0:ndivs-1,0:ndivs-1))
  allocate(proc_iq_min(0:ndivs-1,0:ndivs-1))
  allocate(proc_iq_max(0:ndivs-1,0:ndivs-1))
  allocate(elems(0:ndivs-1,0:ndivs-1))
  allocate(nelem(0:nthetal2-1))

  proc = -1
  nelem = 0
  elems = .false.

  npart = ndivs**2 / nthetal2

  if (dump_mesh_info_screen) then
      write(*,*) 'ntot  = ', ndivs**2
      write(*,*) 'npart = ', npart
  endif


  do ip = 0, nthetal2 - 1, 1
      proc_iq_min = .false.
      proc_iq_max = (proc == - 1)

      ! initialize all elements within two elements distance of the last
      ! processor to this proc, as otherwise they had two neighbors
      ! relevant at high frequencies and lots of procs only
      if (ip > 0) then
          do is = 0, ndivs - 1, 1
              do iz = 0, ndivs - 1, 1
                  if (proc(is, iz) == nthetal2 - ip) then
                      do ids = -2, 2, 1
                          do idz = -2, 2, 1
                              if ((is + ids < 0) .or. (iz + idz < 0) .or. &
                                      (is + ids > ndivs - 1) .or. (iz + idz > ndivs - 1)) then
                                  continue
                              else if (proc(is + ids,iz + idz) == -1) then
                                  proc_iq_min(is + ids, iz + idz) = .true.
                              endif
                          enddo
                      enddo
                  endif
              enddo
          enddo
      endif

      stepsize = pi / 2. / 18.
      sign_buff = 0

      ! optimize number of elements per processor

      if (ip < nthetal2 - 1) then
          do n = 1, nlinsteps, 1
              call nelem_under(ndivs, x0(ip), x1(ip), x2(ip), x3(ip), z0(ip), &
                  z1(ip), z2(ip), z3(ip), proc_iq_max, proc_iq_min, nelem(ip), &
                  elems)
              if (nelem(ip) == npart) then
                  exit
              else
                  if (npart > nelem(ip)) then
                      if (sign_buff == -1) then
                          stepsize = stepsize / 2.
                      endif
                      dphi = stepsize
                      sign_buff = 1
                  else
                      if (sign_buff == 1) then
                          stepsize = stepsize / 2.
                      endif
                      dphi = -1. * stepsize
                      sign_buff = -1
                  endif
                  phi(ip) = phi(ip) + dphi
                  x2(ip) = (r2 + ndivs * 0.3 * sin(2. * phi(ip))**4) * cos(phi(ip))
                  z2(ip) = (r2 + ndivs * 0.3 * sin(2. * phi(ip))**4) * sin(phi(ip))
                  if (dphi < 0) then
                      proc_iq_max = elems
                  else
                      proc_iq_min = elems
                  endif
              endif
          enddo

          ! if the optimization did not converge, take or remove as many from
          ! the elements in question as needed

          do while (nelem(ip) < npart)
              ncorrections = ncorrections + 1
              exit_buff = .false.
              do is = ndivs - 1, 0, -1
                  do iz = ndivs - 1, 0, -1
                      if (XOR(proc_iq_max(is,iz), proc_iq_min(is,iz))) then
                          proc_iq_min(is,iz) = .true.
                          proc_iq_max(is,iz) = .true.
                          elems(is,iz) = .true.
                          nelem(ip) = nelem(ip) + 1
                          exit_buff = .true.
                          exit
                      endif
                  enddo
                  if (exit_buff) then
                      exit
                  endif
              enddo
          enddo

          do while (nelem(ip) > npart)
              ncorrections = ncorrections + 1
              exit_buff = .false.
              do is = ndivs - 1, 0, -1
                  do iz = ndivs - 1, 0, -1
                      if (XOR(proc_iq_max(is,iz), proc_iq_min(is,iz))) then
                          proc_iq_min(is,iz) = .true.
                          proc_iq_max(is,iz) = .true.
                          elems(is,iz) = .false.
                          nelem(ip) = nelem(ip) - 1
                          exit_buff = .true.
                          exit
                      endif
                  enddo
                  if (exit_buff) then
                      exit
                  endif
              enddo
          enddo

          ! fill up the array

          do is = 0, ndivs - 1, 1
              do iz = 0, ndivs - 1, 1
                  if (proc(is,iz) == -1 .and. elems(is,iz)) then
                          proc(is,iz) = nthetal2 - ip - 1
                  endif
              enddo
          enddo

      ! for the last processor just take the remaining elements
      else
          do is = 0, ndivs - 1, 1
              do iz = 0, ndivs - 1, 1
                  if (proc(is,iz) == -1) then
                      nelem(ip) = nelem(ip) + 1
                      proc(is,iz) = nthetal2 - ip - 1
                  endif
              enddo
          enddo
      endif
  enddo

  if (dump_mesh_info_screen) then
      write(*,*) 'ncorrections = ', ncorrections
      write(*,*) 'sum   = ', sum(nelem)
      write(*,*) 'proc  = '
      write(*,*) 'x0 = ', x0
      write(*,*) 'x1 = ', x1
      write(*,*) 'x2 = ', x2
      write(*,*) 'x3 = ', x3
      write(*,*) 'z0 = ', z0
      write(*,*) 'z1 = ', z1
      write(*,*) 'z2 = ', z2
      write(*,*) 'z3 = ', z3

      do ip = 0 , nthetal2 - 1, 1
              write(*,12)ip,nelem(ip)
      enddo
      call ascii_print(ndivs, proc, 3)
  endif
12 format('Central cube: proc', i3, ' has', i6, ' elements')

  exit_buff = test_decomp(ndivs, proc, npart, nthetal2)


  ! connect these processor dependencies to the global element numbering scheme
  central_count(0:nthetal-1) = 0
  do iz = 1, ndivs
      do is = 1, ndivs
          central_count(proc(is-1,iz-1)) = &
              central_count(proc(is-1,iz-1)) + 1
          attributed(central_is_iz_to_globiel(is,iz)) = .true.
          attributed(central_is_iz_to_globiel(is,iz) + neltot / 2) = .true.

          if (solid_domain(ndisc)) then
              procel_solidl(central_count(proc(is-1,iz-1)), &
                      proc(is-1,iz-1)) = central_is_iz_to_globiel(is,iz)
          else
              procel_fluidl(central_count(proc(is-1,iz-1)), &
                      proc(is-1,iz-1)) = central_is_iz_to_globiel(is,iz)
          endif

          ! South: inverted copy
          if (nthetal > 1) then
              if (solid_domain(ndisc)) then
                  procel_solidl(central_count(proc(is-1,iz-1)), &
                      nthetal - 1 - proc(is-1,iz-1)) = &
                      central_is_iz_to_globiel(is,iz) + neltot / 2
              else
                  procel_fluidl(central_count(proc(is-1,iz-1)), &
                      nthetal - 1 - proc(is-1,iz-1)) = &
                      central_is_iz_to_globiel(is,iz) + neltot / 2
              endif
          endif
      enddo
  enddo

  ! South:
  if (nthetal > 1) then
      do ip = 0, nthetal / 2 - 1
          central_count(nthetal - ip - 1) = central_count(ip)
      enddo
  endif

  ! check if all central-cube elements are assigned
  do is = 1, neltot
      if (eltypeg(is) == 'linear' ) then
          if (.not. attributed(is)) then
              write(*,*)
              write(*,*) 'Problem: Central cube element not assigned!',is
              stop
          endif
      endif
  enddo

  ! write out the central cube decomposition

  if (dump_mesh_info_files) then
      open(unit=10008,file=diagpath(1:lfdiag)//'/central_locind_locsz_iproc.dat')
      do iz = 1, ndivs
          do is = 1, ndivs
              write(10008,13) is, iz, s_arr(is,iz), z_arr(is,iz), &
                  proc(is-1,iz-1)
          enddo
      enddo
      do iz = 1, ndivs
          do is = 1, ndivs
              write(10008,13) is, iz, s_arr(is,iz), -z_arr(is,iz), &
                  nthetal - 1 - proc(is-1,iz-1)
          enddo
      enddo
      close(10008)
  endif
13 format(2(i4),2(f9.3),i4)

  if (dump_mesh_info_screen) then
      write(*,*)
      write(*,*)'<> < > < > < > Finished the central domain decomposition! < > < > < > <>'
      write(*,*)
      call flush(6)
  endif

  deallocate(x0, x1, x2, x3, z0, z1, z2, z3, phi)
  deallocate(proc, proc_iq_min, proc_iq_max, elems, nelem)

end subroutine decompose_inner_cube_opt
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
logical function test_decomp(ndivs, proc, npart, nproc2)

  integer, intent(in)     :: ndivs, proc(0:ndivs-1,0:ndivs-1), &
                             nproc2, npart
  integer                 :: is, iz, idx, idz, ip, nelem(0:nproc2), &
                             neighbor_buff

  !test processor bounds
  do is = 0, ndivs - 1, 1
      do iz = 0, ndivs - 1, 1
          if (proc(is,iz) < 0 .or. proc(is,iz) > nproc2 - 1) then
              write(*,*) 'Problem: element (', is, ',', iz, ') has no processor'
              stop
          endif
      enddo
  enddo

  !test number of elements per processor
  nelem = 0
  do ip = 0, nproc2-1, 1
      do is = 0, ndivs - 1, 1
          do iz = 0, ndivs - 1, 1
              if (proc(is,iz) == ip) then
                  nelem(ip) = nelem(ip) + 1
              endif
          enddo
      enddo
      if (nelem(ip) /= npart) then
          write(*,*) 'Problem: number of elements of proc ', ip, ' is ', &
              nelem(ip), ' should be ', npart
          stop
      endif
  enddo

  !test number of neighbors from different domains
  do is = 0, ndivs - 1, 1
      do iz = 0, ndivs - 1, 1
          neighbor_buff = - 1
          do idx = -1, 1, 1
              do idz = -1, 1, 1
                  if ((is + idx < 0) .or. (iz + idz < 0) .or. &
                          (is + idx > ndivs - 1) .or. (iz + idz > ndivs - 1)) then
                      continue
                  else if (proc(is + idx,iz + idz) /= proc(is,iz)) then
                      if (neighbor_buff == -1) then
                          neighbor_buff = proc(is + idx,iz + idz)
                      else if (neighbor_buff /= proc(is + idx,iz + idz)) then
                          call ascii_print_markregion(ndivs, proc, is, iz)
                          write(*,*) 'Problem: element (', is, ',', iz, &
                                     ') has neighbors from two other regions!'
                          stop
                      endif
                  endif
              enddo
          enddo
      enddo
  enddo

  test_decomp = .true.

end function test_decomp
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure logical function below(x0, y0, x1, y1, xp, yp)

  real(kind=sp), intent(in) :: x0, y0, x1, y1, xp, yp
  real(kind=sp)             :: m, b

  m = (y1 - y0) / (x1 - x0)
  b = (x1*y0 - x0*y1) / (x1 - x0)

  if (yp < (m*xp + b)) then
      below = .true.
  else
      below = .false.
  endif
end function below
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure logical function rightof(x0, y0, x1, y1, xp, yp)

  real(kind=sp), intent(in) :: x0, y0, x1, y1, xp, yp
  real(kind=sp)             :: m, b

  m = (y1 - y0) / (x1 - x0)
  b = (x1*y0 - x0*y1) / (x1 - x0)

  if (xp > ((yp - b) / m)) then
      rightof = .true.
  else
      rightof = .false.
  endif
end function rightof
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine nelem_under(ndivs, x0, x1, x2, x3, z0, z1, z2, z3, proc_max, &
                            proc_min, nelem, proc)

  real(kind=sp), intent(in) :: x0, x1, x2, x3, z0, z1, z2, z3
  integer, intent(in)       :: ndivs
  logical, intent(in)       :: proc_max(0:ndivs-1,0:ndivs-1), &
                               proc_min(0:ndivs-1,0:ndivs-1)
  integer, intent(out)      :: nelem
  logical, intent(out)      :: proc(0:ndivs-1,0:ndivs-1)
  integer                   :: is, iz

  nelem = 0

  do is = 0, ndivs - 1, 1
      do iz = 0, ndivs - 1, 1
          if (proc_min(is,iz)) then
              nelem = nelem + 1
          endif
      enddo
  enddo

  proc = proc_min

  do is = 1, ndivs, 1
      do iz = 1, ndivs, 1
          if (XOR(proc_max(is-1,iz-1), proc_min(is-1,iz-1))) then
              ! first segment
              !if (is < x1 .and. below(x0, z0, x1, z1, is, iz)):
              if (below(x0, z0, x1, z1, real(is), real(iz))) then
                  nelem = nelem + 1
                  proc(is-1,iz-1) = .true.
              ! second segment
              else if (x2 > x1 .and. is < x2 .and. is > x1 .and. &
                      below(x1, z1, x2, z2, real(is), real(iz))) then
                  nelem = nelem + 1
                  proc(is-1,iz-1) = .true.
              else if (x2 < x1 .and. iz < z2 .and. iz > z1 .and. &
                      rightof(x1, z1, x2, z2, real(is), real(iz))) then
                  nelem = nelem + 1
                  proc(is-1,iz-1) = .true.
              ! third segment
              else if (x3 > x2 .and. is > x2 .and. x2 > x1 .and. &
                      below(x2, z2, x3, z3, real(is), real(iz))) then
                  nelem = nelem + 1
                  proc(is-1,iz-1) = .true.
              else if (x3 > x2 .and. is > x2 .and. x2 < x1 .and. iz > z2 .and. &
                      below(x2, z2, x3, z3, real(is), real(iz))) then
                  nelem = nelem + 1
                  proc(is-1,iz-1) = .true.
              else if (x3 < x2 .and. iz > z2 .and. &
                      rightof(x2, z2, x3, z3, real(is), real(iz))) then
                  nelem = nelem + 1
                  proc(is-1,iz-1) = .true.
              ! lower right corner
              else if (z3 > z2 .and. iz < z2 .and. is > x2 .and. is > x1) then
                  nelem = nelem + 1
                  proc(is-1,iz-1) = .true.
              endif
          endif
      enddo
  enddo
end subroutine nelem_under
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine ascii_print(ndivs, proc, mode)
  integer :: is, iz
  integer, intent(in) :: mode, ndivs
  integer, intent(in) :: proc(0:ndivs-1,0:ndivs-1)

  if (mode == 1) then
      do iz = ndivs - 1, 0, - 1
          do is = 0, ndivs - 1, 1
              write(*,'(i3,$)') proc(is,iz)
          enddo
          print *
      enddo
  else if (mode == 2) then
      do iz = ndivs - 1, 0, - 1
          do is = 0, ndivs - 1, 1
              if (proc(is,iz) == ((proc(is,iz) / 2) * 2)) then
                  write(*,'(A)', advance='no') '0 '
              else
                  write(*,'(A)', advance='no') 'X '
              endif
          enddo
          print *
      enddo
  else if (mode == 3) then
      do iz = ndivs - 1, 0, - 1
          do is = 0, ndivs - 1, 1
              if (proc(is,iz) == ((proc(is,iz) / 2) * 2)) then
                  write(*,'(A)', advance='no') '0'
              else
                  write(*,'(A)', advance='no') '-'
              endif
          enddo
          print *
      enddo
  endif

end subroutine ascii_print
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine ascii_print_markregion(ndivs, proc, px, pz)
  integer :: is, iz
  integer, intent(in) :: ndivs, px, pz
  integer, intent(in) :: proc(0:ndivs-1,0:ndivs-1)

  do iz = ndivs - 1, 0, - 1
      do is = 0, ndivs - 1, 1
          if ((is == px - 2 .or. is == px + 2) &
                  .or. (iz == pz - 2 .or. iz == pz + 2)) then
              write(*,'(A)', advance='no') '#'
          else if (proc(is,iz) == ((proc(is,iz) / 2) * 2)) then
              write(*,'(A)', advance='no') 'x'
          else
              write(*,'(A)', advance='no') '-'
          endif
      enddo
      print *
  enddo

end subroutine ascii_print_markregion
!-----------------------------------------------------------------------------------------

end module parallelization
!=========================================================================================
