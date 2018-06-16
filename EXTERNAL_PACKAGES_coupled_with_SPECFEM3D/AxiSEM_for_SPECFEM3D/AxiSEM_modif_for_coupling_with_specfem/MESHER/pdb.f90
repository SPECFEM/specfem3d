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
module pdb
  use data_mesh
  use data_numbering
  use data_spec
  use data_diag
  use data_pdb
  use data_bkgrdmodel
  use data_grid

  implicit none
  public :: create_pdb
  private
  contains

!-----------------------------------------------------------------------------------------
subroutine create_pdb
  ! Wrapper routine to define everything in the parallel global and solid/fluid!
  ! worlds and dump the database for each processor.
  ! Additionally creating the header mesh_params.h containing all major mesh
  ! sizes for the solver
  ! Whenever DATABASE occurs as a comment below, those parameters will be saved
  ! into the solver database.
  !
  ! Terminology:
  ! global - counting through the entire solid+fluid domain across processors
  ! slobal - counting through the entire solid domain across processors
  ! glocal - counting through a single processor's entire solid+fluid domain
  ! slocal - counting through a single processor's solid domain
  ! accordingly for fluid: flobal, flocal, and colloquially for sflobal, sflocal.

  use data_gllmesh, only: sgll, zgll
  use data_time
  use clocks_mod

  integer   :: nelmax

  if (dump_mesh_info_screen) then
    write(*,*)
    write(*,*)' ||||||||||||||| CREATING THE PARALLEL DATABASE ||||||||||||||||'
    write(*,*)
  endif

  deallocate(iglob)

  write(*,*) '  define glocal numbering....'; call flush(6)
  iclock12 = tick()
  call define_glocal_numbering ! needs sgll,zgll, creates igloc
  iclock12 = tick(id=idold12, since=iclock12)

  ! Solid-fluid distinction
  write(*,*) '  define solflu coordinates....'; call flush(6)
  call define_sflocal_coordinates ! needs sgll, zgll
                                  ! creates procel_solidp     procel_fluidp
                                  !         inv_procel_solidp inv_procel_fluidp
  nelmax = maxval(nel)

  write(*,*) '  define axial elems....'; call flush(6)
  call define_axial_elem ! needs sgll, sgll_solid, sgll_fluid

  write(*,*) '  define solflu numbering....'; call flush(6)
  iclock13 = tick()
  call define_sflocal_numbering   ! needs sgll, zgll
                                  ! creates igloc_solid, igloc_fluid
  iclock13 = tick(id=idold13, since=iclock13)
  deallocate(sgll,zgll)

  write(*,*) '  define search sflobal index....'; call flush(6)
  call define_search_sflobal_index ! needs iglob_solid, iglob_fluid

  write(*,*) '  partition sflobal index....'; call flush(6)
  call partition_sflobal_index

  write(*,*) '  define local bdry elems....'; call flush(6)
  call define_local_bdry_elem


  ! For compliance with solver, we need the control points info for each process
  write(*,*) '  generate processor serendipity....'; call flush(6)
  iclock14 = tick()
  call generate_serendipity_per_proc(sg,zg) ! needs sgp, zgp
  iclock14 = tick(id=idold14, since=iclock14)

  write(*,*) '  define element type....'; call flush(6)
  call define_element_type !

  ! Write out mesh database
  write(*,*) '  write database....'; call flush(6)
  call write_db

  write(*,*) ' create static header mesh_params.h ....';call flush(6)
  call create_static_header

end subroutine create_pdb
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine define_sflocal_coordinates
! NOTE: procel_fluid and procel_solid map the GLOBAL element numbers!!!
! i.e.: procel(1:nel_fluid(iproc),iproc)=procel_fluid(:,iproc)
!       procel(nel_fluid(iproc)+1:nel(iproc),iproc)= procel_solid(:,iproc)
!
! Also re-defining procel_solid here, and defining new mapping
! from solf/flu parallel elem to global parallel elem number
! This will be the output to the solver to become ielsolid and ielfluid!!
!
! procel_solidp: mapping from slocal to glocal elem number
! inv_procel_solidp: mapping from glocal to slocal elem number
!
! DATABASE: procel_solidp,procel_fluidp

  use data_gllmesh
  integer           :: iproc, iel, ielg, jpol, ipol
  integer           :: nelmax_solid, nelmax_fluid
  real(kind=dp)     :: rsol_min, rsol_max, rflu_min, rflu_max

  nelmax_solid = maxval(nel_solid)
  nelmax_fluid = maxval(nel_fluid)

  ! Min/max values
  if (dump_mesh_info_screen) then
     ! initialize to crazy values
     rsol_max = -1e30
     rflu_max = -1e30
     rsol_min = 1e30
     rflu_min = 1e30

     ! Solid
     do iproc = 0, nproc-1
       do iel = 1, nel_solid(iproc)
         ielg = procel_solid(iel,iproc)
          do jpol = 0, npol
            do ipol = 0, npol
              rsol_max = max(rsol_max, sgll(ipol,jpol,ielg)**2 + zgll(ipol,jpol,ielg)**2)
              rsol_min = min(rsol_min, sgll(ipol,jpol,ielg)**2 + zgll(ipol,jpol,ielg)**2)
            enddo
          enddo
       enddo
     enddo

     write(*,*)'HAVE FLUID?',have_fluid

     ! Fluid
     if (have_fluid) then
       do iproc = 0, nproc-1
         do iel = 1, nel_fluid(iproc)
           ielg = procel_fluid(iel,iproc)
            do jpol = 0, npol
              do ipol = 0, npol
                rflu_max = max(rflu_max, sgll(ipol,jpol,ielg)**2 + zgll(ipol,jpol,ielg)**2)
                rflu_min = min(rflu_min, sgll(ipol,jpol,ielg)**2 + zgll(ipol,jpol,ielg)**2)
              enddo
            enddo
         enddo
       enddo
     endif ! have_fluid

     rsol_max = sqrt(rsol_max)
     rsol_min = sqrt(rsol_min)
     rflu_max = sqrt(rflu_max)
     rflu_min = sqrt(rflu_min)

     write(*,*)'Solid-fluid coordinates:'
     if (have_solid) &
        write(*,*)'Min/max radius solid:', rsol_max, rsol_min
     if (have_fluid) &
        write(*,*)'Min/max radius fluid:', rflu_max, rflu_min
  endif

  ! Define processor-specific mapping:
  ! procel_solidp: given slocal el, return glocal el number
  ! inv_procel_solidp, given glocal, return slocal el number

  allocate(procel_solidp(nelmax_solid,0:nproc-1))
  allocate(procel_fluidp(nelmax_fluid,0:nproc-1))
  allocate(inv_procel_solidp(maxval(nel),0:nproc-1))
  allocate(inv_procel_fluidp(maxval(nel),0:nproc-1))

  do iproc=0, nproc-1
     do iel=1, nel_solid(iproc)
        procel_solidp(iel,iproc) = inv_procel(procel_solid(iel,iproc))
        inv_procel_solidp(procel_solidp(iel,iproc),iproc) = iel
     enddo
     do iel=1, nel_fluid(iproc)
        procel_fluidp(iel,iproc) = inv_procel(procel_fluid(iel,iproc))
        inv_procel_fluidp(procel_fluidp(iel,iproc),iproc) = iel
     enddo
  enddo

  ! check procel_solidp and inverses
  do iproc=0, nproc-1
     do iel=1, nel_solid(iproc)
        if (iel /= inv_procel_solidp(procel_solidp(iel,iproc),iproc) ) then
           write(*,*)'PROBLEM in procel_solidp vs. inv_procel_solidp!'
           write(*,*)'iproc,iel,inv(procel_solidp(iel)):',iproc,iel, &
                      inv_procel_solidp(procel_solidp(iel,iproc),iproc)
           stop
        endif

        if (.not. solid(procel( procel_solidp(iel,iproc),iproc) ) ) then
           write(*,*)'PROBLEM with procel or procel_solidp or both!'
           write(*,*)'slocal,glocal,global el:',iel,procel_solidp(iel,iproc), &
                     solid(procel(procel_solidp(iel,iproc),iproc))
           stop
        endif
     enddo
  enddo

  ! check procel_fluidp and inverses
  do iproc=0, nproc-1
     do iel=1,nel_fluid(iproc)
        if (iel /= inv_procel_fluidp(procel_fluidp(iel,iproc),iproc) ) then
           write(*,*)'PROBLEM in procel_fluidp vs. inv_procel_fluidp!'
           write(*,*)'iproc,iel,inv(procel_fluidp(iel)):',iproc,iel, &
                      inv_procel_fluidp(procel_fluidp(iel,iproc),iproc)
           stop
        endif

        if (.not. fluid(procel(procel_fluidp(iel,iproc),iproc)) ) then
           write(*,*)'PROBLEM with procel or procel_fluidp or both!'
           write(*,*)'flocal,glocal,global el:',iel,procel_fluidp(iel,iproc), &
                     fluid(procel(procel_fluidp(iel,iproc),iproc))
           stop
        endif
     enddo
  enddo

  ! check procel_solidp/procel_fluidp and inverses per global search
  do iproc=0, nproc-1
     do iel=1, nel(iproc)
        if (solid(procel(iel,iproc))) then
           if (iel /= procel_solidp(inv_procel_solidp(iel,iproc),iproc) ) then
              write(*,*)'PROBLEM in procel_solidp vs. inv_procel_solidp!'
              write(*,*)'iproc,iel,procel_solidp(inv(iel)):',iproc,iel, &
                   procel_solidp(inv_procel_solidp(iel,iproc),iproc)
              stop
           endif

        else if (fluid(procel(iel,iproc))) then
           if (iel /= procel_fluidp(inv_procel_fluidp(iel,iproc),iproc) ) then
              write(*,*)'PROBLEM in procel_fluidp vs. inv_procel_fluidp!'
              write(*,*)'iproc,iel,procel_fluidp(inv(iel)):',iproc,iel, &
                         procel_fluidp(inv_procel_fluidp(iel,iproc),iproc)
              stop
           endif
        else
           write(*,*) 'If not solid nor fluid, what then??'
           write(*,*) 'iel (glocal),iproc:',iel,iproc
           stop
        endif
     enddo
  enddo

end subroutine define_sflocal_coordinates
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine define_glocal_numbering
! The glocal numbering for each processor's entire domain.
! The glocal number bookkeeping is not passed to the database,
! as we only need solid and fluid global numbers.

  use numbering
  use data_gllmesh, only: sgll, zgll
  use data_time
  use clocks_mod

  !$ use omp_lib

  integer :: nelmax, nelp
  integer :: npointotp, wnglob
  integer :: iproc, ipol, jpol, iel, ipt, ielg
  real(kind=dp), dimension(:), allocatable :: wsgll, wzgll
  integer, dimension(:), allocatable :: wloc, wigloc
  logical, dimension(:), allocatable :: wifseg

  ! valence test
  real(kind=dp), dimension(:), allocatable     :: uglob2
  real(kind=dp), dimension(:,:,:), allocatable :: val
  integer :: idest,i
  integer :: valnum_cent(6),totvalnum_cent
  integer :: valnum_semi(6),totvalnum_semi
  integer :: nthreads = 1

  nelmax = maxval(nel)
  allocate(igloc(nelmax*(npol+1)**2,0:nproc-1))
  allocate(nglobp(0:nproc-1))

  !$ nthreads = max(OMP_get_max_threads() / nproc, 1)
  !$omp parallel do private(iproc, nelp, npointotp, wsgll, wzgll, iel, ielg, ipol, &
  !$omp                     jpol, ipt, wigloc, wifseg, wloc, wnglob, uglob2, val, &
  !$omp                     valnum_cent, valnum_semi, idest, i, totvalnum_semi, &
  !$omp                     totvalnum_cent) &
  !$omp              shared(nglobp, igloc)
  do iproc = 0, nproc-1

     nelp = nel(iproc)
     npointotp = nelp*(npol+1)**2
     allocate(wsgll(npointotp))
     allocate(wzgll(npointotp))


     do iel = 1, nelp
       ielg = procel(iel,iproc)
       do jpol = 0, npol
         do ipol = 0, npol
           ipt = (iel-1)*(npol+1)**2 + jpol*(npol+1) + ipol + 1
           wsgll(ipt) = sgll(ipol,jpol,ielg)
           wzgll(ipt) = zgll(ipol,jpol,ielg)
         enddo
       enddo
     enddo

     allocate(wigloc(npointotp))
     allocate(wifseg(npointotp))
     allocate(wloc(npointotp))

     call get_global(nelp, wsgll, wzgll, wigloc, wloc, wifseg, wnglob, &
                     npointotp, NGLLCUBE, nthreads)

     deallocate(wzgll,wsgll)
     deallocate(wifseg)
     deallocate(wloc)

     do ipt = 1, npointotp
        igloc(ipt,iproc) = wigloc(ipt)
     enddo
     nglobp(iproc) = wnglob

     ! valence: test global numbering
     allocate(uglob2(nglobp(iproc)))
     allocate(val(0:npol,0:npol,nel(iproc)))

     ! valence test, equivalent to how assembly is used in the solver
     ! use script plot_proc_valence.csh to generate GMT valence grids for each proc
     ! glocally and zoomed into r < 0.2 ( denoted as *_central )
     val(:,:,:) = 1.0

     valnum_cent = 0
     valnum_semi = 0

     uglob2(:) = 0.d0
     do iel = 1, nel(iproc)
       do ipol = 0, npol
         do jpol = 0, npol
           ipt = (iel - 1) * (npol + 1)**2 + jpol * (npol + 1) + ipol + 1
           idest = wigloc(ipt)
           uglob2(idest) = uglob2(idest) + val(ipol,jpol,iel)
         enddo
       enddo
     enddo

     do iel = 1, nel(iproc)
       do ipol = 0, npol
         do jpol = 0, npol
           ipt = (iel - 1) * (npol + 1)**2 + jpol * (npol + 1) + ipol + 1
           idest = wigloc(ipt)
           val(ipol,jpol,iel) = uglob2(idest)

           ! Check valence/global number inside central region
           if (eltypeg(procel(iel,iproc)) == 'linear') then
              do i=1,6 !possible valences
                 if (val(ipol,jpol,iel) == i) valnum_cent(i)=valnum_cent(i) + 1
              enddo
           endif

           if (eltypeg(procel(iel,iproc)) == 'semino' .or. &
               eltypeg(procel(iel,iproc)) == 'semiso') then
              do i=1,6 !possible valences
                 if (val(ipol,jpol,iel) == i) valnum_semi(i)=valnum_semi(i) + 1
              enddo
           endif

         enddo
       enddo
     enddo
     deallocate(uglob2, val)
     deallocate(wigloc)

     totvalnum_cent = 0
     totvalnum_semi = 0
     do i=1,6
        totvalnum_cent = totvalnum_cent + valnum_cent(i) / i
        totvalnum_semi = totvalnum_semi + valnum_semi(i) / i
     enddo

     if (dump_mesh_info_screen) then
        write(*,*)
        write(*,*) iproc, 'glocal number      :', nglobp(iproc)
        write(*,*) iproc, 'central region only:', totvalnum_cent
        write(*,*) iproc, 'semicurved els only:', totvalnum_semi
        write(*,*) iproc, 'everywhere else    :', &
                    nglobp(iproc) - totvalnum_cent-totvalnum_semi
        call flush(6)
     endif

  enddo
  !$omp end parallel do

  if (dump_mesh_info_screen) then
    do iproc = 0, nproc-1
       write(*,*) 'proc.,nglob:', iproc, nglobp(iproc)
    enddo
    write(*,*) 'SUM = ',SUM(nglobp(0:nproc-1))
    write(*,*)
  endif

end subroutine define_glocal_numbering
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine define_sflocal_numbering
! The sflocal numbering for each processor's solid & fluid domains.
! This will be passed to the solver as the defining bookkeeping arrays
! for the assembly procedure.
! DATABASE: igloc_solid,igloc_solid

  use numbering
  use data_gllmesh, only: sgll, zgll
  use data_time
  use clocks_mod
  !$ use omp_lib

  integer :: nelmax_solid, nelmax_fluid, nelp_solid, nelp_fluid
  integer :: npointotp_solid, npointotp_fluid, wnglob_solid, wnglob_fluid
  integer :: iproc, ipol, jpol, iel, ipt,  ielg
  real(kind=dp), dimension(:), allocatable :: wsgll_solid, wzgll_solid
  real(kind=dp), dimension(:), allocatable :: wsgll_fluid, wzgll_fluid

  integer, dimension(:), allocatable :: wigloc_solid,wigloc_fluid
  integer, dimension(:), allocatable :: wloc_solid,wloc_fluid
  logical, dimension(:), allocatable :: wifseg_solid,wifseg_fluid

  ! valence test
  integer :: idest,i
  real(kind=dp), dimension(:), allocatable     :: uglob2_solid
  real(kind=dp), dimension(:,:,:), allocatable :: val_solid
  integer :: valnum_cent_solid(6),totvalnum_cent_solid
  integer :: valnum_semi_solid(6),totvalnum_semi_solid

  real(kind=dp), dimension(:), allocatable     :: uglob2_fluid
  real(kind=dp), dimension(:,:,:), allocatable :: val_fluid
  integer :: valnum_fluid(6),totvalnum_fluid
  integer :: nthreads = 1

  nelmax_solid = maxval(nel_solid)
  nelmax_fluid = maxval(nel_fluid)

  allocate(igloc_solid(nelmax_solid*(npol+1)**2,0:nproc-1))
  allocate(igloc_fluid(nelmax_fluid*(npol+1)**2,0:nproc-1))

  allocate(nglobp_solid(0:nproc-1))
  allocate(nglobp_fluid(0:nproc-1))

  nglobp_solid(:) = 0
  nglobp_fluid(:) = 0

  !$ nthreads = max(OMP_get_max_threads() / nproc, 1)
  !$omp parallel do private(iproc, nelp_solid, nelp_fluid, npointotp_solid, npointotp_fluid, &
  !$omp                     wsgll_solid, wsgll_fluid, wzgll_solid, wzgll_fluid, &
  !$omp                     iel, ielg, ipol, jpol, ipt, &
  !$omp                     wigloc_solid, wigloc_fluid, wifseg_solid, wifseg_fluid, &
  !$omp                     wloc_solid, wloc_fluid, wnglob_solid, wnglob_fluid, &
  !$omp                     uglob2_solid, val_solid, valnum_cent_solid, valnum_semi_solid, &
  !$omp                     uglob2_fluid, val_fluid, valnum_fluid, &
  !$omp                     idest, i, totvalnum_semi_solid, totvalnum_cent_solid, &
  !$omp                     totvalnum_fluid) &
  !$omp             shared(nglobp_solid, nglobp_fluid, igloc_solid, igloc_fluid)
  do iproc = 0, nproc-1

     ! Solid
     if (have_solid) then
        nelp_solid = nel_solid(iproc)
        npointotp_solid = nelp_solid*(npol+1)**2
        allocate(wsgll_solid(npointotp_solid))
        allocate(wzgll_solid(npointotp_solid))

        do iel = 1, nelp_solid
           ielg = procel_solid(iel,iproc)
           do jpol = 0, npol
              do ipol = 0, npol
                 ipt = (iel-1)*(npol+1)**2 + jpol*(npol+1) + ipol + 1
                 wsgll_solid(ipt) = sgll(ipol,jpol,ielg)
                 wzgll_solid(ipt) = zgll(ipol,jpol,ielg)
              enddo
           enddo
        enddo

        allocate(wigloc_solid(npointotp_solid))
        allocate(wifseg_solid(npointotp_solid))
        allocate(wloc_solid(npointotp_solid))

        iclock09 = tick()
        call get_global(nelp_solid, wsgll_solid, wzgll_solid, wigloc_solid, &
                        wloc_solid, wifseg_solid, wnglob_solid, npointotp_solid, &
                        NGLLCUBE, nthreads)
        iclock09 = tick(id=idold09, since=iclock09)

        deallocate(wzgll_solid, wsgll_solid)
        deallocate(wloc_solid)
        deallocate(wifseg_solid)

        do ipt = 1, npointotp_solid
           igloc_solid(ipt,iproc) = wigloc_solid(ipt)
        enddo
        nglobp_solid(iproc) = wnglob_solid

        ! SOLID valence: test slocal numbering
        allocate(uglob2_solid(nglobp_solid(iproc)))
        allocate(val_solid(0:npol,0:npol,nel_solid(iproc)))

        ! valence test, equivalent to how assembly is used in the solver
        ! use script plot_proc_valence.csh to generate GMT valence grids for each proc
        ! glocally and zoomed into r < 0.2 ( denoted as *_central )
        val_solid(:,:,:) = 1.0d0

        valnum_cent_solid = 0
        valnum_semi_solid = 0

        uglob2_solid(:) = 0.d0
        do iel = 1, nel_solid(iproc)
           do ipol = 0, npol
              do jpol = 0, npol
                 ipt = (iel-1)*(npol+1)**2 + jpol*(npol+1) + ipol+1
                 idest = wigloc_solid(ipt)
                 uglob2_solid(idest) = uglob2_solid(idest) + val_solid(ipol,jpol,iel)
              enddo
           enddo
        enddo

        do iel = 1, nel_solid(iproc)
          do ipol = 0, npol
             do jpol = 0, npol
                ipt = (iel-1)*(npol+1)**2 + jpol*(npol+1) + ipol+1
                idest = wigloc_solid(ipt)
                val_solid(ipol,jpol,iel) = uglob2_solid(idest)

                ! Check valence/global number inside central region
                if (eltypeg(procel_solid(iel,iproc)) == 'linear') then
                   do i=1,6 !possible valences
                      if (val_solid(ipol,jpol,iel) == i) &
                           valnum_cent_solid(i)=valnum_cent_solid(i)+1
                   enddo
                endif

                if (eltypeg(procel_solid(iel,iproc)) == 'semino' .or. &
                     eltypeg(procel_solid(iel,iproc)) == 'semiso') then
                   do i=1,6 !possible valences
                      if (val_solid(ipol,jpol,iel) == i) &
                           valnum_semi_solid(i)=valnum_semi_solid(i)+1
                   enddo
                endif
             enddo
          enddo
        enddo

        deallocate(uglob2_solid, val_solid)
        deallocate(wigloc_solid)

     endif ! have_solid

     ! Fluid
     if (have_fluid) then
        nelp_fluid = nel_fluid(iproc)
        npointotp_fluid = nelp_fluid*(npol+1)**2
        allocate(wsgll_fluid(npointotp_fluid))
        allocate(wzgll_fluid(npointotp_fluid))

        do iel = 1, nelp_fluid
           ielg = procel_fluid(iel,iproc)
           do jpol = 0, npol
              do ipol = 0, npol
                 ipt = (iel-1)*(npol+1)**2 + jpol*(npol+1) + ipol + 1
                 wsgll_fluid(ipt) = sgll(ipol,jpol,ielg)
                 wzgll_fluid(ipt) = zgll(ipol,jpol,ielg)
              enddo
           enddo
        enddo

        allocate(wigloc_fluid(npointotp_fluid))
        allocate(wifseg_fluid(npointotp_fluid))
        allocate(wloc_fluid(npointotp_fluid))

        iclock09 = tick()
        call get_global(nelp_fluid, wsgll_fluid, wzgll_fluid, wigloc_fluid, &
             wloc_fluid, wifseg_fluid, wnglob_fluid, npointotp_fluid, NGLLCUBE)
        iclock09 = tick(id=idold09, since=iclock09)

        deallocate(wzgll_fluid, wsgll_fluid)
        deallocate(wloc_fluid)
        deallocate(wifseg_fluid)

        do ipt = 1, npointotp_fluid
           igloc_fluid(ipt,iproc) = wigloc_fluid(ipt)
        enddo
        nglobp_fluid(iproc) = wnglob_fluid

        ! FLUID valence: test flocal numbering
        allocate(uglob2_fluid(nglobp_fluid(iproc)))
        allocate(val_fluid(0:npol,0:npol,nel_fluid(iproc)))

        ! valence test, equivalent to how assembly is used in the solver
        ! use script plot_proc_valence.csh to generate GMT valence grids for each proc
        ! glocally and zoomed into r < 0.2 ( denoted as *_central )
        val_fluid(:,:,:) = 1.0d0

        valnum_fluid = 0

        uglob2_fluid(:) = 0.d0
        do iel = 1, nel_fluid(iproc)
          do ipol = 0, npol
            do jpol = 0, npol
              ipt = (iel-1)*(npol+1)**2 + jpol*(npol+1) + ipol+1
              idest = wigloc_fluid(ipt)
              uglob2_fluid(idest) = uglob2_fluid(idest) + val_fluid(ipol,jpol,iel)
            enddo
          enddo
        enddo
        !QT
        do iel = 1, nel_fluid(iproc)
          do ipol = 0, npol
            do jpol = 0, npol
              ipt = (iel-1)*(npol+1)**2 + jpol*(npol+1) + ipol+1
              idest = wigloc_fluid(ipt)
              val_fluid(ipol,jpol,iel) = uglob2_fluid(idest)

              ! Check valence/global number inside central region
              do i=1,6 !possible valences
                 if (val_fluid(ipol,jpol,iel) == i) &
                      valnum_fluid(i)=valnum_fluid(i)+1
              enddo
            enddo
          enddo
        enddo

        deallocate(uglob2_fluid, val_fluid)
        deallocate(wigloc_fluid)

     endif ! have_fluid

     if (have_solid) then
        totvalnum_cent_solid = 0
        totvalnum_semi_solid = 0
        do i=1,6
           totvalnum_cent_solid = totvalnum_cent_solid + valnum_cent_solid(i)/i
           totvalnum_semi_solid = totvalnum_semi_solid + valnum_semi_solid(i)/i
        enddo

        if (dump_mesh_info_screen) then
           write(*,*)
           write(*,*) ' SOLID VALENCE & GLOBAL NUMBERING:'
           write(*,*) iproc,'slocal number      :',nglobp_solid(iproc)
           write(*,*) iproc,'scentral region only:',totvalnum_cent_solid
           write(*,*) iproc,'ssemicurved els only:',totvalnum_semi_solid
           write(*,*) iproc,'severywhere else    :',nglobp_solid(iproc)-&
                totvalnum_cent_solid-totvalnum_semi_solid
        endif
     endif ! have_solid

     if (have_fluid) then
        totvalnum_fluid = 0
        do i=1,6
           totvalnum_fluid = totvalnum_fluid + valnum_fluid(i)/i

        enddo

        if (dump_mesh_info_screen) then
         write(*,*) ' FLUID VALENCE & GLOBAL NUMBERING:'
         write(*,*) iproc,'flocal number      :',nglobp_fluid(iproc)
         write(*,*) iproc,'fluid region only  :',totvalnum_fluid
         write(*,*) iproc,'feverywhere else   :',nglobp_fluid(iproc)-&
                                                 totvalnum_fluid
        endif
     endif ! have_fluid
  enddo ! iproc
  !$omp end parallel do

  if (dump_mesh_info_screen) then
     ! here are used some arrays such as nglobp that are not needed
     ! nglobp will not be defined in subsequent versions of the mesher
     ! this whole printout is bound to be erased
     write(*,*)
     write(*,*)'Amount of unique grid points in each processor:'
     do iproc = 0, nproc-1
      write(*,21) iproc, nglobp_solid(iproc), nglobp_fluid(iproc), nglobp(iproc), &
                  nglobp_solid(iproc) + nglobp_fluid(iproc) - nglobp(iproc)
     enddo

     write(*,*)
     write(*,*) 'SUM # grid points Solid:       ',SUM(nglobp_solid(0:nproc-1))
     write(*,*) 'SUM # grid points Fluid:       ',SUM(nglobp_fluid(0:nproc-1))
     write(*,*) 'SUM # grid points Total:       ',SUM(nglobp(0:nproc-1))
     write(*,*) '    # total grid points:       ',nglobglob
     write(*,*) 'SUM # grid points S/F Boundary:',SUM(nglobp_fluid(0:nproc-1))+ &
                               SUM(nglobp_solid(0:nproc-1))-SUM(nglobp(0:nproc-1))
     write(*,*) 'Predicted S/F Boundary points :', &
                                 2*(npol+1)+(sum(nbelem)/2-nbcnd)*(npol)

  endif

  ! This test is not valid for nproc>1 since the proc boundaries are singular
  ! global points for each proc and hence one would need to subtract the
  ! sum of the proc boundary points here.
  ! Could be done ;)
  !
  !  if (have_fluid .and. 2*(npol+1)+(sum(nbelem)/2-nbcnd)*(npol) /= &
  !          SUM(nglobp_fluid(0:nproc-1))+ &
  !          SUM(nglobp_solid(0:nproc-1))-SUM(nglobp(0:nproc-1)) ) then
  !     write(*,*)'...something wrong with global number of S/F boundary points..'
  !     call flush(6)
  !     stop
  !  endif

  if (dump_mesh_info_screen) write(*,*)

21 format('Proc',i3, ' has',i9, ' solid ',i9,' fluid',i9, &
          ' total and',i6,' S/F boundary pts')
end subroutine define_sflocal_numbering
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine define_global2glocal
! Define bookkeeping array that for a given global number, returns
! the glocal number. Used later on to partition global index into message bins,
! and deallocated there, in partition_global_index
! Both this and partition_global_index might well be redundant after all....

  integer :: ipt, iel, ipol, jpol, iproc, ielg, iptg, igg, igp
  allocate(glob2gloc(nglobglob,0:nproc-1))
  glob2gloc(:,:) = 0

  do iproc = 0, nproc-1
     do iel = 1, nel(iproc)
        ielg = procel(iel,iproc)
        do jpol = 0, npol
           do ipol = 0, npol
              iptg = (ielg-1)*(npol+1)**2 + jpol*(npol+1) + ipol + 1
              ipt  = (iel -1)*(npol+1)**2 + jpol*(npol+1) + ipol + 1
              igg = iglob(iptg)      ! global #
              igp = igloc(ipt,iproc) ! glocal #
              glob2gloc(igg,iproc) = igp
           enddo
        enddo
     enddo
  enddo

  deallocate(igloc)

end subroutine define_global2glocal
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine define_sflobal2sflocal(iproc, solorflu)
! Define bookkeeping array that for a given sflobal number, returns
! the sflocal number. Used later on to partition global index into message
! bins for the solid/fluid subdomains, and should be deallocated there,
! in partition_sflobal_index

  integer, intent(in) :: iproc
  logical, intent(in) :: solorflu
  integer             :: ipt, iel, ipol, jpol, ielg, iptg, igg, igp

  if (dump_mesh_info_screen .and. iproc == 0) then
     write(*,*)'Info on global numbering array sizes/values:'
     write(*,*)'MAX INV_SOL/FLU     :',maxval(inv_ielem_solid), &
                                      maxval(inv_ielem_fluid)
     if (have_solid) write(*,*)'MAX/SIZE IGLOB_SOLID:', &
                        maxval(iglob_solid),size(iglob_solid)
     if (have_fluid) &
         write(*,*)'MAX/SIZ IGLOB_FLUID :',maxval(iglob_fluid), size(iglob_fluid)

     if (have_solid) write(*,*)'MAX/SIZE SOL IGLOC  :', &
                        maxval(igloc_solid),size(igloc_solid)

     if (have_fluid) &
        write(*,*)'MAX/SIZE FLU IGLOC  :',maxval(igloc_fluid),size(igloc_fluid)


     write(*,*)'MAX PROCEL/SOL/FLU  :',maxval(procel),maxval(procel_solid), &
                                      maxval(procel_fluid)
     write(*,*)'SIZ PROCEL/SOL/FLU  :',size(procel),size(procel_solid), &
                                      size(procel_fluid)
     write(*,*)
     call flush(6)
  endif

  ! Solid
  if (solorflu) then
     slob2sloc(:) = 0

     do iel = 1, nel_solid(iproc)
        ielg = procel_solid(iel,iproc)
        do jpol = 0, npol
           do ipol = 0, npol
              iptg = (inv_ielem_solid(ielg)-1)*(npol+1)**2 + jpol*(npol+1) + ipol + 1
              ipt  = (iel -1)*(npol+1)**2 + jpol*(npol+1) + ipol + 1
              igg = iglob_solid(iptg)      ! slobal #
              igp = igloc_solid(ipt,iproc) ! slocal #
              slob2sloc(igg) = igp
           enddo
        enddo
     enddo
    !enddo

  ! Fluid
  else
    if (have_fluid) then
       flob2floc(:) = 0
       do iel = 1, nel_fluid(iproc)
          ielg = procel_fluid(iel,iproc)
          do jpol = 0, npol
             do ipol = 0, npol
                iptg = (inv_ielem_fluid(ielg)-1)*(npol+1)**2 + jpol*(npol+1) + ipol + 1
                ipt  = (iel -1)*(npol+1)**2 + jpol*(npol+1) + ipol + 1
                igg = iglob_fluid(iptg)
                igp = igloc_fluid(ipt,iproc)
                flob2floc(igg) = igp
             enddo
          enddo
       enddo
    endif
  endif

end subroutine define_sflobal2sflocal
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine define_search_sflobal_index
! Define a sflobal index: given GLL coords and proc returns sflobal number.
! Also defined are
! nprocb_solid: Given slobal number returns # of procs it belongs to (valence)
! lprocb_solid: Given slobal number and index of valence returns processor ID
! .... and correspondingly for fluid arrays.
! These are crucial for the partitioning, i.e. the heart of message passing.

  integer :: iproc, ielg, iel, ipol, jpol, ipt, nelmax_solid, nelmax_fluid
  integer :: il, nprocbmax_solid, nprocbmax_fluid
  integer :: nsearch, nneighbors
  integer, allocatable :: nbelong2_solid(:,:), nbelong2_fluid(:,:)

  nelmax_solid = maxval(nel_solid)
  nelmax_fluid = maxval(nel_fluid)

  allocate(nprocb_solid(nglobslob))
  allocate(nprocb_fluid(nglobflob))

  allocate(nbelong_solid(nglobslob))

  if (nradialslices == 1) then
     nneighbors = 2
  else
     nneighbors = 8
  endif
  allocate(nbelong2_solid(1:nneighbors,nglobslob))

  nbelong2_solid(:,:) = -1
  nbelong_solid(:)    =  0
  nprocbmax_solid     =  1

  ! solid
  do iproc = 0, nproc - 1
     do iel = 1, nel_solid(iproc)
        ielg = procel_solid(iel, iproc) ! global element number
        do jpol = 0, npol
           outer: do ipol = 0, npol
               ipt = (inv_ielem_solid(ielg) - 1) * (npol + 1)**2 &
                        + jpol * (npol + 1) + ipol + 1
               nbelong_solid(iglob_solid(ipt)) = nbelong_solid(iglob_solid(ipt)) + 1

               nsearch = 1
               do while (nbelong2_solid(nsearch,iglob_solid(ipt)) /= -1)
                  if (nbelong2_solid(nsearch,iglob_solid(ipt)) == iproc) &
                        cycle outer
                  nsearch = nsearch + 1
                  if (nsearch == nneighbors + 1) then
                     write(*,*) 'ERROR: too many neighbors in solid'
                     stop
                  endif
               enddo

               nbelong2_solid(nsearch,iglob_solid(ipt)) = iproc
               nprocbmax_solid = max(nsearch, nprocbmax_solid)
           enddo outer
        enddo
     enddo
  enddo !iproc

  if (dump_mesh_info_screen) then
      write(*,*) 'Maximum number of neighbors in the solid:', nprocbmax_solid
  endif

  allocate(lprocb_solid(nprocbmax_solid,nglobslob))

  do ipt=1, nglobslob
     do il = 1, nneighbors
        if (nbelong2_solid(il,ipt) /= -1) then
           lprocb_solid(il,ipt) = nbelong2_solid(il,ipt)
           nprocb_solid(ipt) = il
        endif
     enddo
  enddo

  deallocate(nbelong2_solid)


  ! fluid
  if (have_fluid) then
     allocate(nbelong_fluid(nglobflob))
     allocate(nbelong2_fluid(1:nneighbors,nglobflob))

     nbelong2_fluid(:,:) = -1
     nbelong_fluid(:)    = 0
     nprocbmax_fluid     = 1

     do iproc = 0, nproc - 1
        do iel = 1, nel_fluid(iproc)
           ielg = procel_fluid(iel,iproc) ! global element number
           do jpol = 0, npol
              outerfl: do ipol = 0, npol
                 ipt = (inv_ielem_fluid(ielg) - 1) * (npol + 1)**2 &
                        + jpol * (npol + 1) + ipol + 1

                 nbelong_fluid(iglob_fluid(ipt)) = nbelong_fluid(iglob_fluid(ipt)) + 1

                 nsearch = 1
                 do while (nbelong2_fluid(nsearch,iglob_fluid(ipt)) /= -1)
                    if (nbelong2_fluid(nsearch,iglob_fluid(ipt)) == iproc) &
                          cycle outerfl
                    nsearch = nsearch + 1
                    if (nsearch == nneighbors + 1) then
                       write(*,*) 'ERROR: too many neighbors in fluid'
                       stop
                    endif
                 enddo

                 nbelong2_fluid(nsearch,iglob_fluid(ipt)) = iproc
                 nprocbmax_fluid = max(nsearch, nprocbmax_fluid)
              enddo outerfl
           enddo
        enddo
     enddo

     if (dump_mesh_info_screen) then
         write(*,*) 'Maximum number of neighbors in the solid:', nprocbmax_solid
     endif
     allocate(lprocb_fluid(nprocbmax_fluid,nglobflob))

     do ipt=1, nglobflob
        do il = 1, nneighbors
           if (nbelong2_fluid(il,ipt) /= -1) then
              lprocb_fluid(il,ipt) = nbelong2_fluid(il,ipt)
              nprocb_fluid(ipt) = il
           endif
        enddo
     enddo
     deallocate(nbelong2_fluid)

  endif !have_fluid

end subroutine define_search_sflobal_index
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine partition_sflobal_index
! Partitioning the index defined in define_search_sflobal_index, i.e.
! creating bins for the message passing across processors. These inherently
! depend on the sflocal number as the message passing is done during the
! assembly, i.e. sflocal stage in the solver.
!
! SOME TERM EXPLANATIONS (All have _solid appended, left out for brevity):
!
! nglobglob:    global number of points
! nprocb:       Given global numnber, lists the number of procs it belongs to
! lprocb:       given global number, indices 1,..,nprocb list proc it belongs to
! sizebin:      given a proc, lists glocal number
! binp:         given glocal number and processor, lists global number
! sizemsg:      for each proc-proc pair, lists number of global points shared
! global_index_msg: for each proc-proc pair and shared points 1,...,sizemsgmax
!                   lists global number
!
! THE FOLLOWING ENTER THE DATABASE FOR THE SOLVER:
! sizerecvp, sizesendp:         given a proc, lists # other procs to communicate with
! listrecvp, listsendp:         given a proc and all its shared points, lists proc ID
!                               to send/receive that point with
! sizemsgrecvp,sizemsgsendp:    given proc and each communicating proc, list size
!                               of respective messages sent/received
! glocal_index_msg_recvp, glocal_index_msg_sendp:
!                               same as glocal_ind..for glocal
!
! DATABASE: sizerecvp_solid, listrecvp_solid
! DATABASE: sizemsgrecvp_solid
! DATABASE: glocal_index_msg_recvp_solid
! DATABASE: sizesendp_solid, listsendp_solid
! DATABASE: sizemsgsendpg_solid
! DATABASE: glocal_index_msg_sendp_solid
!
! DATABASE: sizerecvp_fluid, listrecvp_fluid
! DATABASE: sizemsgrecvp_fluid
! DATABASE: glocal_index_msg_recvp_fluid
! DATABASE: sizesendp_fluid, listsendp_fluid
! DATABASE: sizemsgsendpg_fluid
! DATABASE: glocal_index_msg_sendp_fluid
  use sorting

  integer :: ipt, iproct, ibp, ibel, ig, ip
  integer :: ipdes, ipsrc, imsg

  integer, dimension(:), allocatable        :: ibin_solid
  integer, dimension(:,:), allocatable      :: myneighbors_solid
  integer, dimension(:,:), allocatable      :: sizemsg_solid
  integer, dimension(:,:), allocatable      :: index_msg_solid
  integer, dimension(:,:,:), allocatable    :: global_index_msg_solid
  integer, dimension(:,:), allocatable      :: binp_solid
  integer, dimension(:), allocatable        :: sizebin_solid
  integer :: sizebinmax_solid, sizemsgmax_solid
  integer :: sizerecvpmax_solid, sizesendpmax_solid

  integer, dimension(:), allocatable        :: ibin_fluid
  integer, dimension(:,:), allocatable      :: myneighbors_fluid
  integer, dimension(:,:), allocatable      :: sizemsg_fluid
  integer, dimension(:,:), allocatable      :: index_msg_fluid
  integer, dimension(:,:,:), allocatable    :: global_index_msg_fluid
  integer, dimension(:,:), allocatable      :: binp_fluid
  integer, dimension(:), allocatable        :: sizebin_fluid

  integer :: sizebinmax_fluid, sizemsgmax_fluid
  integer :: sizerecvpmax_fluid, sizesendpmax_fluid
  integer :: nneighbors, inbr
  real(kind=dp), dimension(:), allocatable  :: sort_buf


  if (dump_mesh_info_screen) then
     write(*,*)
     write(*,*)'****************************************************************'
     write(*,*)'******************* CREATING MESSAGING ARRAYS ******************'
     write(*,*)'****************************************************************'
  endif

  if (nradialslices == 1) then
     nneighbors = 2
  else
     nneighbors = 8
  endif

  if (dump_mesh_info_screen) write(*,*) 'nneighbors', nneighbors

  ! SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS SOLID SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS

  if (dump_mesh_info_screen) write(*,*) 'Creating solid bins '
  allocate(sizebin_solid(0:nproc-1))

  sizebin_solid(0:nproc-1) = 0
  do ipt = 1, nglobslob
     do ibel = 1, nprocb_solid(ipt)
        iproct = lprocb_solid(ibel,ipt)
        sizebin_solid(iproct) = sizebin_solid(iproct) + 1
     enddo
  enddo

  sizebinmax_solid = maxval(sizebin_solid(:))

  if (dump_mesh_info_screen) then
     do iproct = 0, nproc - 1
        write(*,*) 'proc:', iproct, 'size solid bin:', sizebin_solid(iproct)
     enddo
     write(*,*) 'sum sizebin,slobal:', sum(sizebin_solid),nglobslob
     write(*,*) ' sizebinmax_solid = ' , sizebinmax_solid
  endif

  allocate(binp_solid(sizebinmax_solid,0:nproc-1))
  binp_solid = 0

  allocate(ibin_solid(0:nproc-1))
  ibin_solid(0:nproc-1) = 0

  do ipt = 1, nglobslob
     do ibel = 1, nprocb_solid(ipt)
        iproct = lprocb_solid(ibel,ipt)
        ibin_solid(iproct) = ibin_solid(iproct) + 1
        ibp = ibin_solid(iproct)
        binp_solid(ibp,iproct) = ipt
     enddo
  enddo
  deallocate(ibin_solid)

  allocate(sizemsg_solid(0:nproc-1,nneighbors))
  sizemsg_solid = 0
  allocate(myneighbors_solid(0:nproc-1,nneighbors))
  myneighbors_solid = -1

  do iproct = 0, nproc - 1
     do ipt = 1, sizebin_solid(iproct)
        ig = binp_solid(ipt,iproct)
        if (nprocb_solid(ig) > 1) then
           do ibel = 1, nprocb_solid(ig)
              ipdes = lprocb_solid(ibel,ig)
              if (ipdes /= iproct) then
                 ! find first empty neighbor location
                 do inbr=1, nneighbors
                    if (myneighbors_solid(iproct,inbr) == ipdes .or. &
                        myneighbors_solid(iproct,inbr) == -1) exit
                 enddo
                 if (inbr > 8) then
                    write(*,*) 'ERORR: having more then 8 neighbors (+myself)'
                    write(*,*) '       check mesh decomposition in the solid)'
                    write(*,*) '       iproct = ', iproct
                    stop
                 endif
                 if (myneighbors_solid(iproct,inbr) == -1) myneighbors_solid(iproct,inbr) = ipdes
                 sizemsg_solid(iproct,inbr) = sizemsg_solid(iproct,inbr) + 1
              endif
           enddo
        endif
     enddo
  enddo

  ! sort neighbors and messaging array according to processor number
  allocate(sort_buf(nneighbors))
  do iproct = 0, nproc - 1
     if (any(myneighbors_solid(iproct,:) == -1)) then
        inbr = minval(minloc(myneighbors_solid(iproct,:))) - 1
     else
        inbr = nneighbors
     endif
     sort_buf(1:inbr) = dble(myneighbors_solid(iproct,1:inbr))
     call mergesort_3(sort_buf(1:inbr), il=myneighbors_solid(iproct,1:inbr), &
                      il2=sizemsg_solid(iproct,1:inbr), p=1)
     if (dump_mesh_info_screen) &
         write(*,'(100(i4))') myneighbors_solid(iproct,1:inbr)
  enddo
  deallocate(sort_buf)

  if (dump_mesh_info_screen) then
     write(*,*) 'Size of solid messages for each proc-proc pair:'
     do iproct = 0,  nproc-1
        ! MvD: 100 beeing the maximum number of procs??? yes, but is debuggin only
        write(*,'(100(i4))') sizemsg_solid(iproct,:)
     enddo
     write(*,*) 'Total solid messages size:',SUM(SUM(sizemsg_solid,DIM=1))
  endif

  sizemsgmax_solid = maxval(maxval(sizemsg_solid,DIM=1))
  if (dump_mesh_info_screen) write(*,*) 'size msg max solid is ' , sizemsgmax_solid

  allocate(index_msg_solid(0:nproc-1,nneighbors))
  index_msg_solid = 0

  allocate(global_index_msg_solid(sizemsgmax_solid,0:nproc-1,nneighbors))
  global_index_msg_solid = -1

  do iproct = 0, nproc-1
     do ipt = 1, sizebin_solid(iproct)
        ig = binp_solid(ipt,iproct)
        if (nprocb_solid(ig) > 1) then
           do ibel = 1, nprocb_solid(ig)
              ipdes = lprocb_solid(ibel,ig)
              if (ipdes /= iproct) then

                 do inbr=1, nneighbors
                    if (ipdes == myneighbors_solid(iproct,inbr)) exit
                 enddo

                 index_msg_solid(iproct,inbr) = index_msg_solid(iproct,inbr) + 1
                 imsg = index_msg_solid(iproct,inbr)
                 global_index_msg_solid(imsg,iproct,inbr) = ig
              endif
           enddo
        endif
     enddo
  enddo

  deallocate(index_msg_solid)
  deallocate(sizebin_solid)
  deallocate(binp_solid)

  ! How many messages am I sending / receiving ?
  allocate(sizerecvp_solid(0:nproc-1))
  allocate(sizesendp_solid(0:nproc-1))
  sizerecvp_solid(0:nproc-1) = 0
  sizesendp_solid(0:nproc-1) = 0

  do iproct = 0, nproc-1
     do inbr=1, nneighbors
        ipsrc = myneighbors_solid(iproct,inbr)
        if (ipsrc == -1) exit
     enddo
     sizerecvp_solid(iproct) = inbr - 1
  enddo

  ! somewhat redundant, but makes sure symmetry of the communication
  do iproct = 0, nproc -1
     do ipdes = 0, nproc-1
        do inbr=1, nneighbors
           ipsrc = myneighbors_solid(ipdes,inbr)
           if (ipsrc == iproct) sizesendp_solid(iproct) = sizesendp_solid(iproct) + 1
        enddo
     enddo
  enddo

  if (dump_mesh_info_screen) then
     do iproct = 0, nproc -1
       write(*,'("Proc", i3, " receives solid stuff from", i3, " procs and sends to", i3, " procs")') &
            iproct, sizerecvp_solid(iproct), sizesendp_solid(iproct)
     enddo
  endif

  sizerecvpmax_solid = maxval(sizerecvp_solid(:))
  sizesendpmax_solid = maxval(sizesendp_solid(:))

  if (dump_mesh_info_screen) then
     write(*,*)'max size recv solid:', sizerecvpmax_solid
     write(*,*)'max size send solid:', sizesendpmax_solid
  endif

  ! To which processors ?
  allocate(listrecvp_solid(sizerecvpmax_solid,0:nproc-1))
  allocate(listsendp_solid(sizesendpmax_solid,0:nproc-1))
  listrecvp_solid = -1
  listsendp_solid = -1

  do iproct = 0, nproc-1
     ip = 0
     do inbr=1, nneighbors
        ipsrc = myneighbors_solid(iproct,inbr)
        if (ipsrc == -1) exit
        ip = ip + 1
        listrecvp_solid(ip,iproct) = ipsrc
     enddo
  enddo

  ! somewhat redundant, but makes sure symmetry of the communication
  do iproct = 0, nproc -1
     ip = 0
     do ipdes = 0, nproc-1
        do inbr=1, nneighbors
           ipsrc = myneighbors_solid(ipdes,inbr)
           if (ipsrc == iproct) then
              ip = ip + 1
              listsendp_solid(ip,iproct) = ipdes
           endif
        enddo
     enddo
  enddo

  if (dump_mesh_info_screen .and. nproc > 1) then
     do iproct = 0, nproc-1
        write(*,'("Proc", i3, " will receive ", i2, " solid messages from procs ", 20(i3,1x))') &
              iproct, sizerecvp_solid(iproct), &
              (listrecvp_solid(ip,iproct),ip=1,sizerecvp_solid(iproct))
        write(*,'("Proc", i3, " will send    ", i2, " solid messages to   procs ", 20(i3,1x))') &
              iproct, sizesendp_solid(iproct), &
              (listsendp_solid(ip,iproct),ip=1,sizesendp_solid(iproct))
     enddo
  endif


  ! What size ?
  allocate(sizemsgrecvp_solid(sizerecvpmax_solid,0:nproc-1))
  allocate(sizemsgsendp_solid(sizesendpmax_solid,0:nproc-1))

  sizemsgrecvp_solid = 0
  sizemsgsendp_solid = 0

  do iproct = 0, nproc-1
     do inbr=1, nneighbors
        ipsrc = myneighbors_solid(iproct,inbr)
        if (ipsrc == -1) exit
        sizemsgrecvp_solid(inbr,iproct) = sizemsg_solid(iproct,inbr)
     enddo
  enddo

  ! somewhat redundant, but makes sure symmetry of the communication
  do iproct = 0, nproc -1
     do ipdes = 0, nproc-1
        do inbr=1, nneighbors
           ipsrc = myneighbors_solid(ipdes,inbr)
           if (ipsrc == iproct) then
              sizemsgsendp_solid(inbr,ipdes) = sizemsg_solid(ipdes,inbr)
           endif
        enddo
     enddo
  enddo

  ! OUTPUT message size
  if (dump_mesh_info_screen) then
     write(*,*)
     do iproct = 0, nproc-1
        if (sizerecvp_solid(iproct) > 0) then
           do ip = 1, sizerecvp_solid(iproct)
              write(*,'("Proc",i3," receiving solid message from",i3," sized",i6)') &
                    iproct, listrecvp_solid(ip,iproct), sizemsgrecvp_solid(ip,iproct)
           enddo
        endif

        if (sizesendp_solid(iproct) > 0) then
           do ip = 1, sizesendp_solid(iproct)
              write(*,'("Proc",i3," sending   solid message to  ",i3," sized",i6)') &
                    iproct, listsendp_solid(ip,iproct), sizemsgsendp_solid(ip,iproct)
           enddo
        endif
     enddo
  endif

  ! NOW CREATE GLOCAL INDEX FOR MESSAGES
  ! Which glocal indices?
  allocate(glocal_index_msg_recvp_solid(sizemsgmax_solid,sizerecvpmax_solid,0:nproc-1))
  allocate(glocal_index_msg_sendp_solid(sizemsgmax_solid,sizesendpmax_solid,0:nproc-1))

  glocal_index_msg_recvp_solid = 0
  glocal_index_msg_sendp_solid = 0

  allocate(slob2sloc(nglobslob))

  do iproct = 0, nproc-1

     call define_sflobal2sflocal(iproct, .true.)

     if (sizerecvp_solid(iproct) > 0) then
        do ip=1, sizerecvp_solid(iproct)
           ipsrc = listrecvp_solid(ip,iproct)
           do ipt = 1, sizemsgrecvp_solid(ip,iproct)
              do inbr=1, nneighbors
                 if (iproct == myneighbors_solid(ipsrc,inbr)) exit
              enddo
              if (inbr > nneighbors) exit
              ig = global_index_msg_solid(ipt,ipsrc,inbr)
              glocal_index_msg_recvp_solid(ipt,ip,iproct) = slob2sloc(ig)
           enddo
        enddo
     endif

     if (sizesendp_solid(iproct) > 0) then
        do ip = 1, sizesendp_solid(iproct)
           ipdes = listsendp_solid(ip,iproct)
           do ipt = 1, sizemsgsendp_solid(ip,iproct)
              do inbr=1, nneighbors
                 if (ipdes == myneighbors_solid(iproct,inbr)) exit
              enddo
              if (inbr > nneighbors) exit
              ig = global_index_msg_solid(ipt,iproct,inbr)
              glocal_index_msg_sendp_solid(ipt,ip,iproct) = slob2sloc(ig)
           enddo
        enddo
     endif
     call flush(6)
  enddo

  deallocate(slob2sloc)
  deallocate(myneighbors_solid)
  deallocate(global_index_msg_solid)

  if (any(glocal_index_msg_sendp_solid /= glocal_index_msg_recvp_solid)) then
      write(*,*) 'ERROR: Index Array for send and recv should be identical in'
      write(*,*) '       the new communication scheme, but are not!'
      stop
  endif

  if (any(sizesendp_solid /= sizerecvp_solid)) then
      write(*,*) 'ERROR: Messages for send and recv should have same size in'
      write(*,*) '       the new communication scheme, but are not!'
      stop
  endif

  write(*,*)'End of solid messaging'
  write(*,*)
  call flush(6)


  !FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF FLUID FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF

  if (have_fluid) then
     if (dump_mesh_info_screen) write(*,*) ' Creating fluid bins '
     allocate(sizebin_fluid(0:nproc-1))
     sizebin_fluid(0:nproc-1) = 0
     do ipt = 1, nglobflob
        do ibel = 1, nprocb_fluid(ipt)
           iproct = lprocb_fluid(ibel,ipt)
           sizebin_fluid(iproct) = sizebin_fluid(iproct) + 1
        enddo
     enddo

     if (dump_mesh_info_screen) then
        do iproct = 0, nproc - 1
           write(*,*) 'proc:', iproct, 'size fluid bin:', sizebin_fluid(iproct)
        enddo
     endif

     if (dump_mesh_info_screen) write(*,*) sum(sizebin_fluid),nglobflob
     sizebinmax_fluid = maxval(sizebin_fluid(:))

     if (dump_mesh_info_screen) &
     write(*,*) ' sizebinmax_fluid = ' , sizebinmax_fluid

     allocate(binp_fluid(sizebinmax_fluid,0:nproc-1))
     binp_fluid(1:sizebinmax_fluid,0:nproc-1) = 0

     allocate(ibin_fluid(0:nproc-1))
     ibin_fluid(0:nproc-1) = 0

     do ipt = 1, nglobflob
        do ibel = 1, nprocb_fluid(ipt)
           iproct =  lprocb_fluid(ibel,ipt)
           ibin_fluid(iproct) = ibin_fluid(iproct) + 1
           ibp = ibin_fluid(iproct)
           binp_fluid(ibp,iproct) = ipt
        enddo
     enddo
     deallocate(ibin_fluid)

     allocate(sizemsg_fluid(0:nproc-1,nneighbors))
     sizemsg_fluid = 0
     allocate(myneighbors_fluid(0:nproc-1,nneighbors))
     myneighbors_fluid = -1

     do iproct = 0, nproc - 1
        do ipt = 1, sizebin_fluid(iproct)
           ig = binp_fluid(ipt,iproct)
           if (nprocb_fluid(ig) > 1) then
              do ibel = 1, nprocb_fluid(ig)
                 ipdes = lprocb_fluid(ibel,ig)
                 if (ipdes /= iproct) then
                    ! find first empty neighbor location
                    do inbr=1, nneighbors
                       if (myneighbors_fluid(iproct,inbr) == ipdes .or. &
                           myneighbors_fluid(iproct,inbr) == -1) exit
                    enddo
                    if (inbr > 8) then
                       write(*,*) 'ERORR: having more then 8 neighbors (+myself)'
                       write(*,*) '       check mesh decomposition in the fluid)'
                       write(*,*) '       iproct = ', iproct
                       stop
                    endif
                    if (myneighbors_fluid(iproct,inbr) == -1) myneighbors_fluid(iproct,inbr) = ipdes
                    sizemsg_fluid(iproct,inbr) = sizemsg_fluid(iproct,inbr) + 1
                 endif
              enddo
           endif
        enddo
     enddo

     ! sort neighbors and messaging array according to processor number
     allocate(sort_buf(nneighbors))
     do iproct = 0, nproc - 1
        if (any(myneighbors_fluid(iproct,:) == -1)) then
           inbr = minval(minloc(myneighbors_fluid(iproct,:))) - 1
        else
           inbr = nneighbors
        endif
        sort_buf(1:inbr) = dble(myneighbors_fluid(iproct,1:inbr))
        call mergesort_3(sort_buf(1:inbr), il=myneighbors_fluid(iproct,1:inbr), &
                         il2=sizemsg_fluid(iproct,1:inbr), p=1)
        if (dump_mesh_info_screen) &
           write(*,'(100(i4))') myneighbors_fluid(iproct,1:inbr)
     enddo
     deallocate(sort_buf)

     if (dump_mesh_info_screen) then
        write(*,*) 'Size of fluid messages for each proc-proc pair:'
        do iproct = 0,  nproc-1
           ! MvD: 100 beeing the maximum number of procs??? yes, but is debuggin only
           write(*,'(100(i4))') sizemsg_fluid(iproct,:)
        enddo
        write(*,*) 'Total fluid messages size:',SUM(SUM(sizemsg_fluid,DIM=1))
     endif

     sizemsgmax_fluid = maxval(maxval(sizemsg_fluid,DIM=1))

     if (dump_mesh_info_screen) write(*,*) 'size msg max fluid is ' , sizemsgmax_fluid

     allocate(index_msg_fluid(0:nproc-1,nneighbors))
     index_msg_fluid = 0
     allocate(global_index_msg_fluid(sizemsgmax_fluid,0:nproc-1,nneighbors))
     global_index_msg_fluid = -1

     do iproct = 0, nproc-1
        do ipt = 1, sizebin_fluid(iproct)
           ig = binp_fluid(ipt,iproct)
           if (nprocb_fluid(ig) > 1) then
              do ibel = 1, nprocb_fluid(ig)
                 ipdes = lprocb_fluid(ibel,ig)
                 if (ipdes /= iproct) then

                    do inbr=1, nneighbors
                       if (ipdes == myneighbors_fluid(iproct,inbr)) exit
                    enddo

                    index_msg_fluid(iproct,inbr) = index_msg_fluid(iproct,inbr) + 1
                    imsg = index_msg_fluid(iproct,inbr)
                    global_index_msg_fluid(imsg,iproct,inbr) = ig
                 endif
              enddo
           endif
        enddo
     enddo

     deallocate(index_msg_fluid)
     deallocate(sizebin_fluid)
     deallocate(binp_fluid)

     ! How many messages am I sending / receiving ?
     allocate(sizerecvp_fluid(0:nproc-1))
     allocate(sizesendp_fluid(0:nproc-1))
     sizerecvp_fluid = 0
     sizesendp_fluid = 0

     do iproct = 0, nproc-1
        do inbr=1, nneighbors
           ipsrc = myneighbors_fluid(iproct,inbr)
           if (ipsrc == -1) exit
        enddo
        sizerecvp_fluid(iproct) = inbr - 1
     enddo

     ! somewhat redundant, but makes sure symmetry of the communication
     do iproct = 0, nproc -1
        do ipdes = 0, nproc-1
           do inbr=1, nneighbors
              ipsrc = myneighbors_fluid(ipdes,inbr)
              if (ipsrc == iproct) sizesendp_fluid(iproct) = sizesendp_fluid(iproct) + 1
           enddo
        enddo
     enddo

     if (dump_mesh_info_screen) then
        do iproct = 0, nproc -1
           write(*,'("Proc", i3, " receives fluid stuff from",i3, &
                    &" procs and sends to", i3, " procs")') &
                 iproct, sizerecvp_fluid(iproct), sizesendp_fluid(iproct)
        enddo
     endif


     sizerecvpmax_fluid = maxval(sizerecvp_fluid(:))
     sizesendpmax_fluid = maxval(sizesendp_fluid(:))

     if (dump_mesh_info_screen) then
        write(*,*) 'max size recv fluid:', sizerecvpmax_fluid
        write(*,*) 'max size send fluid:', sizesendpmax_fluid
     endif

     ! To which processor ?
     allocate(listrecvp_fluid(sizerecvpmax_fluid,0:nproc-1))
     allocate(listsendp_fluid(sizesendpmax_fluid,0:nproc-1))
     listrecvp_fluid = -1
     listsendp_fluid = -1

     do iproct = 0, nproc-1
        ip = 0
        do inbr=1, nneighbors
           ipsrc = myneighbors_fluid(iproct,inbr)
           if (ipsrc == -1) exit
           ip = ip + 1
           listrecvp_fluid(ip,iproct) = ipsrc
        enddo
     enddo

     ! somewhat redundant, but makes sure symmetry of the communication
     do iproct = 0, nproc -1
        ip = 0
        do ipdes = 0, nproc-1
           do inbr=1, nneighbors
              ipsrc = myneighbors_fluid(ipdes,inbr)
              if (ipsrc == iproct) then
                 ip = ip + 1
                 listsendp_fluid(ip,iproct) = ipdes
              endif
           enddo
        enddo
     enddo

     if (dump_mesh_info_screen .and. nproc > 1) then
        do iproct = 0, nproc-1
           write(*,'("Proc", i2, " will receive ", i2, &
                    &" fluid messages from procs ", 20(i3,1x))') &
                 iproct, sizerecvp_fluid(iproct), &
                 (listrecvp_fluid(ip,iproct), ip=1,sizerecvp_fluid(iproct))
           write(*,'("Proc", i2, " will send    ", i2, &
                    &" fluid messages to   procs ", 20(i3,1x))') &
                 iproct, sizesendp_fluid(iproct), &
                 (listsendp_fluid(ip,iproct), ip=1,sizesendp_fluid(iproct))
        enddo
     endif


     ! What size ?
     allocate(sizemsgrecvp_fluid(sizerecvpmax_fluid,0:nproc-1))
     allocate(sizemsgsendp_fluid(sizesendpmax_fluid,0:nproc-1))
     sizemsgrecvp_fluid = 0
     sizemsgsendp_fluid = 0

     do iproct = 0, nproc-1
        do inbr=1, nneighbors
           ipsrc = myneighbors_fluid(iproct,inbr)
           if (ipsrc == -1) exit
           sizemsgrecvp_fluid(inbr,iproct) = sizemsg_fluid(iproct,inbr)
        enddo
     enddo

     ! somewhat redundant, but makes sure symmetry of the communication
     do iproct = 0, nproc -1
        do ipdes = 0, nproc-1
           do inbr=1, nneighbors
              ipsrc = myneighbors_fluid(ipdes,inbr)
              if (ipsrc == iproct) then
                 sizemsgsendp_fluid(inbr,ipdes) = sizemsg_fluid(ipdes,inbr)
              endif
           enddo
        enddo
     enddo

     ! OUTPUT message size
     if (dump_mesh_info_screen .and. nproc > 1) then
        write(*,*)
        do iproct = 0, nproc -1
           if (sizerecvp_fluid(iproct) > 0) then
              do ip = 1, sizerecvp_fluid(iproct)
                 write(*,'("Proc",i3," receiving fluid message from",i3," sized",i6)') &
                       iproct, listrecvp_fluid(ip,iproct), sizemsgrecvp_fluid(ip,iproct)
              enddo
           endif

           if (sizesendp_fluid(iproct) > 0) then
              do ip = 1, sizesendp_fluid(iproct)
                 write(*,'("Proc",i3," sending   fluid message to  ",i3," sized",i6)') &
                       iproct, listsendp_fluid(ip,iproct), sizemsgsendp_fluid(ip,iproct)
              enddo
           endif
        enddo
     endif

     ! NOW CREATE GLOCAL INDEX FOR MESSAGES
     ! Which glocal indices?
     allocate(glocal_index_msg_recvp_fluid(sizemsgmax_fluid,sizerecvpmax_fluid,0:nproc-1))
     allocate(glocal_index_msg_sendp_fluid(sizemsgmax_fluid,sizesendpmax_fluid,0:nproc-1))

     glocal_index_msg_recvp_fluid = 0
     glocal_index_msg_sendp_fluid = 0

     allocate(flob2floc(nglobflob))

     do iproct = 0, nproc-1

        call define_sflobal2sflocal(iproct, .false.)

        if (sizerecvp_fluid(iproct) > 0) then
           do ip=1, sizerecvp_fluid(iproct)
              ipsrc = listrecvp_fluid(ip,iproct)
              do ipt = 1, sizemsgrecvp_fluid(ip,iproct)
                 do inbr=1, nneighbors
                    if (iproct == myneighbors_fluid(ipsrc,inbr)) exit
                 enddo
                 if (inbr > nneighbors) exit
                 ig = global_index_msg_fluid(ipt,ipsrc,inbr)
                 glocal_index_msg_recvp_fluid(ipt,ip,iproct) = flob2floc(ig)
              enddo
           enddo
        endif

        if (sizesendp_fluid(iproct) > 0) then
           do ip = 1, sizesendp_fluid(iproct)
              ipdes = listsendp_fluid(ip,iproct)
              do ipt = 1, sizemsgsendp_fluid(ip,iproct)
                 do inbr=1, nneighbors
                    if (ipdes == myneighbors_fluid(iproct,inbr)) exit
                 enddo
                 if (inbr > nneighbors) exit
                 ig = global_index_msg_fluid(ipt,iproct,inbr)
                 glocal_index_msg_sendp_fluid(ipt,ip,iproct) = flob2floc(ig)
              enddo
           enddo
        endif
        call flush(6)
     enddo

     deallocate(myneighbors_fluid)
     deallocate(global_index_msg_fluid)

     if (any(glocal_index_msg_sendp_fluid /= glocal_index_msg_recvp_fluid)) then
         write(*,*) 'ERROR: Index Array for send and recv should be identical in'
         write(*,*) '       the new communication scheme, but are not!'
         stop
     endif

     if (any(sizesendp_fluid /= sizerecvp_fluid)) then
         write(*,*) 'ERROR: Messages for send and recv should have same size in'
         write(*,*) '       the new communication scheme, but are not!'
         stop
     endif

  endif ! have_fluid

  if (allocated(flob2floc)) deallocate(flob2floc)

  if (dump_mesh_info_screen .and. nproc > 1) then
     write(*,*)
     write(*,*) '----------------------------------------------------------------'
     write(*,*) 'Sum over all solid message sizes  :', SUM(SUM(sizemsg_solid,DIM=1))
     write(*,*) 'Total global points in solid,ratio:', nglobslob, &
                              real(SUM(SUM(sizemsg_solid,DIM=1))) / real(nglobslob)
     write(*,*) '----------------------------------------------------------------'

     if (have_fluid) then
        write(*,*) '----------------------------------------------------------------'
        write(*,*) 'Sum over all fluid message sizes  :', SUM(SUM(sizemsg_fluid,DIM=1))
        write(*,*) 'Total global points in fluid,ratio:', nglobflob, &
                                real(SUM(SUM(sizemsg_fluid,DIM=1))) / real(nglobflob)
        write(*,*) '----------------------------------------------------------------'
        write(*,*) '----------------------------------------------------------------'
        write(*,*) 'Sum over all s/f message sizes:', SUM(SUM(sizemsg_fluid,DIM=1)) + &
                                                      SUM(SUM(sizemsg_solid,DIM=1))
        write(*,*) 'Total global points, ratio    :', nglobglob, &
                      real(SUM(SUM(sizemsg_fluid,DIM=1)) + &
                           SUM(SUM(sizemsg_solid,DIM=1))) / real(nglobglob)
        write(*,*) '----------------------------------------------------------------'
     endif !have_fluid

     write(*,*)
     write(*,*) '****************************************************************'
     write(*,*) '******************* END OF MESSAGING ARRAYS ********************'
     write(*,*) '****************************************************************'
  endif

end subroutine partition_sflobal_index
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine define_local_bdry_elem

  use global_parameters
  use data_gllmesh, only: sgll, zgll

  integer                              :: iel, iproc, j
  integer, dimension(0:nproc-1)        :: solid_count, fluid_count
  integer, allocatable, dimension(:,:) :: tmpsolid, tmpfluid, jpolsol, jpolflu
  real(kind=dp)                        :: rbound
  integer                              :: myproc, herproc, myielglob, herielglob

  allocate(nbdry_el(0:nproc-1))
  allocate(have_bdry_elemp(0:nproc-1))
  allocate(belemp(sum(nbelem),0:nproc-1))
  allocate(tmpsolid(sum(nbelem)/2,0:nproc-1))
  allocate(tmpfluid(sum(nbelem)/2,0:nproc-1))
  allocate(jpolsol(sum(nbelem)/2,0:nproc-1))
  allocate(jpolflu(sum(nbelem)/2,0:nproc-1))

  tmpsolid = -1
  tmpfluid = -1

  solid_count(:) = 0
  fluid_count(:) = 0
  nbdry_el(:) = 0

  do j=1, nbcnd

     if (mod(j,2) /= 0) then ! upper boundary of fluid region
        rbound = discont(idom_fluid(j)) / router
     else  ! lower boundary of fluid region
        rbound = discont(idom_fluid(j-1)+1) / router
     endif

     do iel=1, nbelem(j)
        nbdry_el(el2proc(belem(iel,j))) = nbdry_el(el2proc(belem(iel,j))) + 1

        ! belemp: given element counter on boundary for a proc, return glocal element number
        belemp(nbdry_el(el2proc(belem(iel,j))),el2proc(belem(iel,j))) = &
             inv_procel(belem(iel,j))

        if (solid(belem(iel,j))) then

           myielglob = belem(iel,j)
           herielglob = belem(my_neighbor(iel,j),j)

           myproc =  el2proc(myielglob)
           herproc = el2proc(herielglob)

           if (herproc /= myproc) then
              write(*,*)
              write(*,*) 'PROBLEM: changing processors across solid/fluid boundary!'
              write(*,*) 'This case is not implemented at this point as it requires'
              write(*,*) 'message passing when adding the boundary term on both sides'
              write(*,*) 'and therefore seems highly ineffective. Sorry...'
              write(*,*) 'solid domain proc:', myproc
              write(*,*) 'fluid domain proc:', herproc
              write(*,*) 'boundary,global element num:',j,iel,myielglob
              write(*,*) 'One possible reaon: NRADIAL_SLICES > 8 in inparam_mesh might lead to this problem'
              stop
           endif

           solid_count(myproc) = solid_count(myproc) + 1

           ! Determine jpol at the boundary: 0,npol depending on above/below, North/South
           if (rcom(myielglob) > rbound) then !solid above solid-fluid boundary
              if (zcom(myielglob) >= 0.d0) then ! North
                 jpolsol(solid_count(myproc),myproc) = 0
                 jpolflu(solid_count(myproc),myproc) = npol
              else ! South
                 jpolsol(solid_count(myproc),myproc) = npol
                 jpolflu(solid_count(myproc),myproc) = 0
              endif
           else ! solid below solid-fluid boundary
              if (zcom(myielglob) >= 0.d0) then ! North
                 jpolsol(solid_count(myproc),myproc) = npol
                 jpolflu(solid_count(myproc),myproc) = 0
              else ! South
                 jpolsol(solid_count(myproc),myproc) = 0
                 jpolflu(solid_count(myproc),myproc) = npol
              endif
              ! crude way to accomodate the case of having the buffer layer below the ICB.
              ! ...i.e., the jpol indices are switched for 45 < theta < 135 deg
              if (eltypeg(myielglob) /= 'curved') then
                 if ( scom(myielglob) > abs(zcom(myielglob)) ) then
                    jpolsol(solid_count(myproc),myproc)=&
                         abs(jpolsol(solid_count(myproc),myproc)-npol)
                 endif
              endif

           endif

           tmpsolid(solid_count(myproc),myproc) = &
              inv_procel_solidp( inv_procel(myielglob), myproc)

           tmpfluid(solid_count(myproc),herproc) = &
              inv_procel_fluidp( inv_procel(herielglob), herproc)

        else
           fluid_count(myproc)=fluid_count(myproc)+1
        endif
     enddo
  enddo

  if (dump_mesh_info_screen) write(*,*) 'ended solid-fluid boundary loop'

  if (sum(solid_count) /= sum(fluid_count) ) then
     write(*,*) 'Something wrong with the # solid/fluid bdry elements:'
     write(*,*) 'Solid count:',sum(solid_count)
     write(*,*) 'Fluid count:',sum(fluid_count)
     stop
  endif

  nbdry_el = solid_count

  allocate(bdry_solid_elp(maxval(nbdry_el),0:nproc-1))
  allocate(bdry_fluid_elp(maxval(nbdry_el),0:nproc-1))

  bdry_solid_elp = -1
  bdry_fluid_elp = -1

  do iproc=0, nproc-1
     if (minval(tmpsolid(1:nbdry_el(iproc),iproc)) < 1) then
        write(*,*) 'Problem with bdry_solid count!'
        write(*,*) minval(tmpsolid(1:nbdry_el(iproc),iproc))
        write(*,*) 'Unassigned elements...'
        stop
     endif

     if (minval(tmpfluid(1:nbdry_el(iproc),iproc)) < 1) then
        write(*,*) 'Problem with bdry_fluid count!'
        write(*,*) minval(tmpfluid(1:nbdry_el(iproc),iproc))
        write(*,*) 'Unassigned elements...'
        stop
     endif

     if (nbdry_el(iproc) > 0) then
        bdry_solid_elp(1:nbdry_el(iproc),iproc) = tmpsolid(1:nbdry_el(iproc),iproc)
        bdry_fluid_elp(1:nbdry_el(iproc),iproc) = tmpfluid(1:nbdry_el(iproc),iproc)
     endif

     if (dump_mesh_info_screen) then
        write(*,'("Proc", i3, " has ", i5, " S/F boundary elements (on one side)")') &
              iproc, nbdry_el(iproc)
        if (maxval(nbdry_el) >= 1) &
            write(*,*) bdry_solid_elp(1,iproc), bdry_fluid_elp(1,iproc)
     endif
  enddo



  if (dump_mesh_info_screen) then
     write(*,*) 'max bdry_solid_elp,max loc el:', &
                maxval(bdry_solid_elp),maxval(nel_solid)
     if (have_fluid) &
        write(*,*) 'max bdry_fluid_elp,max loc el:', &
                    maxval(bdry_fluid_elp),maxval(nel_fluid)
  endif

  if (have_fluid) then
     allocate(bdry_jpol_solidp(maxval(nbdry_el),0:nproc-1))
     bdry_jpol_solidp(:,:) = jpolsol(1:maxval(nbdry_el),:)

     allocate(bdry_jpol_fluidp(maxval(nbdry_el),0:nproc-1))
     bdry_jpol_fluidp(:,:) = jpolflu(1:maxval(nbdry_el),:)
  endif

  ! boolean array to check whether proc has solid-fluid boundary
  have_bdry_elemp(:) = .false.
  do iproc=0, nproc-1
     if (nbdry_el(iproc) > 0 ) have_bdry_elemp(iproc) = .true.
  enddo


  if (dump_mesh_info_screen .and. have_fluid) then
     write(*,*) 'BOUNDARY TERMS:'
     do iproc=0, nproc-1
        if (have_bdry_elemp(iproc)) then
           write(*,*) iproc, 'minmax bdry sol:', &
                      minval(bdry_solid_elp(1:nbdry_el(iproc),iproc)), &
                      maxval(bdry_solid_elp(1:nbdry_el(iproc),iproc))
           write(*,*) iproc, 'minmax bdry flu:', &
                      minval(bdry_fluid_elp(1:nbdry_el(iproc),iproc)), &
                      maxval(bdry_fluid_elp(1:nbdry_el(iproc),iproc))
        endif
     enddo
  endif
end subroutine define_local_bdry_elem
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine define_axial_elem

  use data_gllmesh

  integer :: iel, iproc, ielg
  integer :: dummyax_elp(maxval(nel),0:nproc-1)
  integer :: dummyax_el_solidp(maxval(nel_solid),0:nproc-1)
  integer :: dummyax_el_fluidp(maxval(nel_fluid),0:nproc-1)

  allocate(naxelp(0:nproc-1))
  allocate(naxel_solidp(0:nproc-1))
  allocate(naxel_fluidp(0:nproc-1))
  allocate(axis(maxval(nel),0:nproc-1))
  allocate(axis_solid(maxval(nel_solid),0:nproc-1))
  allocate(axis_fluid(maxval(nel_fluid),0:nproc-1))

  axis = 0
  axis_solid = 0
  axis_fluid = 0
  dummyax_elp = 0
  dummyax_el_solidp = 0
  dummyax_el_fluidp = 0
  naxelp = 0
  naxel_fluidp = 0
  naxel_solidp = 0

  if (dump_mesh_info_screen) then
    write(*,*)
    write(*,*) 'Defining axial elements:'
  endif

  do iproc=0, nproc-1

     ! glocal domain
     do iel=1, nel(iproc)
        ! jpol=npol is random choice
        ielg = procel(iel,iproc)

        if ( sgll(0,npol,ielg) < min_distance_nondim ) then
           axis(iel,iproc) = 1
           naxelp(iproc) = naxelp(iproc) + 1
           dummyax_elp(naxelp(iproc),iproc) = iel
        endif
     enddo ! iel

     !slocal domain
     do iel=1, nel_solid(iproc)
        ielg = procel_solid(iel,iproc)
        ! jpol=npol is random choice
        if ( sgll(0,npol,ielg) < min_distance_nondim ) then
           axis_solid(iel,iproc) = 1
           naxel_solidp(iproc) = naxel_solidp(iproc) + 1
           dummyax_el_solidp(naxel_solidp(iproc),iproc) = iel
        endif
     enddo ! iel

     !flocal domain
     do iel=1, nel_fluid(iproc)
        ielg = procel_fluid(iel,iproc)
        ! jpol=npol is random choice
        if ( sgll(0,npol,ielg) < min_distance_nondim ) then
           axis_fluid(iel,iproc) = 1
           naxel_fluidp(iproc) = naxel_fluidp(iproc) + 1
           dummyax_el_fluidp(naxel_fluidp(iproc),iproc) = iel
        endif
     enddo ! iel

  enddo !iproc

  write(*,*) 'NAXEL:', naxelp(0), naxel_fluidp(0), min_distance_nondim

  allocate(ax_elp(maxval(naxelp),0:nproc-1))
  allocate(ax_el_solidp(maxval(naxel_solidp),0:nproc-1))
  allocate(ax_el_fluidp(maxval(naxel_fluidp),0:nproc-1))

  do iproc=0, nproc-1
     if (dump_mesh_info_screen) &
        write(*,'(" Proc", i3, " has ", i4," total, ",i4, " solid, and ",i4, " fluid axial elements")') &
             iproc, naxelp(iproc), naxel_solidp(iproc), naxel_fluidp(iproc)

     ax_elp(1:naxelp(iproc),iproc) = dummyax_elp(1:naxelp(iproc),iproc)
     ax_el_solidp(1:naxel_solidp(iproc),iproc) = &
                 dummyax_el_solidp(1:naxel_solidp(iproc),iproc)
     ax_el_fluidp(1:naxel_fluidp(iproc),iproc) = &
                 dummyax_el_fluidp(1:naxel_fluidp(iproc),iproc)

     if (naxelp(iproc) /= naxel_solidp(iproc)+naxel_fluidp(iproc)) then
        write(*,*) 'Processor', iproc, ': PROBLEM in counting axial elements:'
        write(*,*) 'Solid axial elements:', naxel_solidp(iproc)
        write(*,*) 'Fluid axial elements:', naxel_fluidp(iproc)
        write(*,*) 'Total axial elements:', naxelp(iproc)
        stop
     endif

  enddo

end subroutine define_axial_elem
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine generate_serendipity_per_proc(sg, zg)

  use numbering
  use data_time
  use clocks_mod

  real(kind=dp), dimension(4,neltot), intent(in) :: sg, zg

  integer :: iproc, nelp, iel, ielg
  integer :: npointotp, nelpmax, nglobmeshpmax
  integer :: ncp
  integer :: wnglob, iptcp, inode

  real(kind=dp), dimension(:,:,:), allocatable   :: sgp, zgp
  real(kind=dp), dimension(:,:), allocatable     :: sgpw, zgpw
  integer, dimension(:), allocatable             :: wiglob, wloc
  logical, dimension(:), allocatable             :: wifseg

  if (dump_mesh_info_screen) then
    write(*,*)
    write(*,*) ' SERENDIPITY DB '
    write(*,*)
  endif

  ncp = 8 ! exclusively using the serendipity quadrilateral element topology

  ! MvD: as far as I can see, this is somewhat of an overkill: e.g. the
  !      spheroidal elements (which is the bulk of the mantel) do not need all
  !      the 8 points, but only the 4 corners. Compare e.g.
  !      SOLVER/analytic_spheroid_mapping.f90

  ! FIRST DETERMINE NUMBER OF CONTROL POINTS FOR A GIVEN PROCESSOR
  allocate(nglobmeshp(0:nproc-1))
  nelpmax = maxval(nel(:))

  ! elemental physical coordinates of serendipity control nodes, parallel
  allocate(sgp(8,nelpmax,0:nproc-1))
  allocate(zgp(8,nelpmax,0:nproc-1))
  sgp(8,nelpmax,0:nproc-1) = 0.d0
  zgp(8,nelpmax,0:nproc-1) = 0.d0

  ! Global number of control nodes, parallel
  allocate(iglobcp(nelpmax*8,0:nproc-1))
  iglobcp(nelpmax*8,0:nproc-1) = 0

  do iproc = 0, nproc-1
     nelp = nel(iproc)
     npointotp = nelp*8
     do iel = 1, nelp
        ielg = procel(iel,iproc)
        sgp(1,iel,iproc) = sg(1,ielg)
        zgp(1,iel,iproc) = zg(1,ielg)

        sgp(3,iel,iproc) = sg(2,ielg)
        zgp(3,iel,iproc) = zg(2,ielg)

        sgp(5,iel,iproc) = sg(3,ielg)
        zgp(5,iel,iproc) = zg(3,ielg)

        sgp(7,iel,iproc) = sg(4,ielg)
        zgp(7,iel,iproc) = zg(4,ielg)

        sgp(2,iel,iproc) = .5d0*(sgp(1,iel,iproc) + sgp(3,iel,iproc))
        sgp(4,iel,iproc) = .5d0*(sgp(3,iel,iproc) + sgp(5,iel,iproc))
        sgp(6,iel,iproc) = .5d0*(sgp(5,iel,iproc) + sgp(7,iel,iproc))
        sgp(8,iel,iproc) = .5d0*(sgp(7,iel,iproc) + sgp(1,iel,iproc))
        zgp(2,iel,iproc) = .5d0*(zgp(1,iel,iproc) + zgp(3,iel,iproc))
        zgp(4,iel,iproc) = .5d0*(zgp(3,iel,iproc) + zgp(5,iel,iproc))
        zgp(6,iel,iproc) = .5d0*(zgp(5,iel,iproc) + zgp(7,iel,iproc))
        zgp(8,iel,iproc) = .5d0*(zgp(7,iel,iproc) + zgp(1,iel,iproc))
     enddo

     allocate(sgpw(8,nelp))
     allocate(zgpw(8,nelp))

     do iel=1, nelp
        sgpw(:,iel) = sgp(:,iel,iproc)
        zgpw(:,iel) = zgp(:,iel,iproc)
     enddo

     allocate(wiglob(npointotp))
     allocate(wloc(npointotp))
     allocate(wifseg(npointotp))

     wiglob(1:npointotp) = 0
     wloc(1:npointotp) = 0

     iclock09 = tick()
     call get_global(nelp, sgpw, zgpw, wiglob, wloc, wifseg, wnglob, npointotp, ncp)
     iclock09 = tick(id=idold09, since=iclock09)

     nglobmeshp(iproc) = wnglob
     if (dump_mesh_info_screen) &
        write(*,*) ' iproc = ', iproc, ' nglobmesh = ', nglobmeshp(iproc)
     do iel = 1, nelp
        do inode = 1, 8
           iptcp = (iel - 1) * 8 + inode
           iglobcp(iptcp,iproc) = wiglob(iptcp)
        enddo
     enddo

     deallocate(wifseg, wloc, wiglob)
     deallocate(zgpw, sgpw)
  enddo

  ! COORDINATES OF THESE CONTROL POINTS (GLOCALLY NUMBERED)
  ! TOPOLOGY OF THE ELEMENTS
  nglobmeshpmax = maxval(nglobmeshp)
  allocate(scpp(nglobmeshpmax,0:nproc-1))
  allocate(zcpp(nglobmeshpmax,0:nproc-1))
  scpp(nglobmeshpmax,0:nproc-1) = 0
  zcpp(nglobmeshpmax,0:nproc-1) = 0
  allocate(lnodescp(nelpmax,8,0:nproc-1))

  do iproc = 0, nproc - 1
     nelp = nel(iproc)
     do iel = 1, nelp
        do inode = 1, 8
           iptcp = (iel-1)*8 + inode

           ! glocally numbered physical coordinates of serendipity control nodes, parallel
           scpp(iglobcp(iptcp,iproc),iproc) = sgp(inode,iel,iproc)
           zcpp(iglobcp(iptcp,iproc),iproc) = zgp(inode,iel,iproc)

           ! global number of control nodes stored by element, parallel
           lnodescp(iel,inode,iproc) = iglobcp(iptcp,iproc)
        enddo
     enddo
  enddo

end subroutine generate_serendipity_per_proc
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine define_element_type

  integer :: nelp, nelpmax
  integer :: iproc, iel, ielg

  if (dump_mesh_info_screen) write(*,*) ' DEFINING ELEMENT TYPE ARRAY FOR EACH PROCESSOR '
  nelpmax = maxval(nel)

  allocate(eltypep(nelpmax,0:nproc-1))
  allocate(coarsingp(nelpmax,0:nproc-1))

  do iproc = 0, nproc-1
     nelp = nel(iproc)
     do iel = 1, nelp
        ielg = procel(iel,iproc)
        eltypep(iel,iproc) = eltypeg(ielg)
        coarsingp(iel,iproc) = coarsing(ielg)
     enddo
  enddo

  ! Fluid
  nelpmax = maxval(nel_fluid)
  allocate(eltypep_fluid(nelpmax,0:nproc-1))

  do iproc = 0, nproc-1
     do iel = 1, nel_fluid(iproc)
        ielg = procel_fluid(iel,iproc)
        eltypep_fluid(iel,iproc) = eltypeg(ielg)
     enddo
  enddo

  ! Solid
  nelpmax = maxval(nel_solid)
  allocate(eltypep_solid(nelpmax,0:nproc-1))

  do iproc = 0, nproc-1
     do iel = 1, nel_solid(iproc)
        ielg = procel_solid(iel,iproc)
        eltypep_solid(iel,iproc) = eltypeg(ielg)
     enddo
  enddo

end subroutine define_element_type
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine write_db
! Writes out a database file to be read by the solver for each processor.

  use data_gllmesh
  use data_spec
  use background_models, only: override_ext_q

  integer           :: iproc, iptp, npointotp, ipsrc, imsg, iel, inode, ielg, idom
  character(len=4)  :: appiproc
  character(len=80) :: dbname
  integer           :: lfdbname

  do iproc=0, nproc-1

     call define_io_appendix(appiproc,iproc)
     dbname='meshdb.dat'//appiproc
     lfdbname=index(dbname,' ')-1
     if (dump_mesh_info_screen) write(*,*) 'WRITING OUT DATABASE TO ',dbname(1:lfdbname)
     open(10,file=dbname(1:lfdbname), FORM="UNFORMATTED", STATUS="REPLACE")

     npointotp = nel(iproc)*(npol+1)**2

     !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     ! Write basic mesh parameters. What was once in mesh_params.h
     write(10) nproc                        ! nproc_mesh
     write(10) npol                         ! npol in SOLVER
     write(10) nel(iproc)                   ! nelem in SOLVER
     write(10) nel(iproc)*(npol+1)**2       ! npoint in SOLVER
     write(10) nel_solid(iproc)             ! nel_solid
     write(10) nel_fluid(iproc)             ! nel_fluid
     write(10) nel_solid(iproc)*(npol+1)**2 ! npoint_solid
     write(10) nel_fluid(iproc)*(npol+1)**2 ! npoint_fluid
     write(10) nglobp_solid(iproc)          ! nglob_solid
     write(10) nglobp_fluid(iproc)          ! nglob_fluid
     write(10) nbdry_el(iproc)              ! nel_bdry
     write(10) ndisc                        ! ndisc
     write(10) lfbkgrdmodel                 ! lfbkgrdmodel

     ! write spectral stuff
     write(10) xi_k
     write(10) eta
     write(10) dxi
     write(10) wt
     write(10) wt_axial_k
     write(10) G0
     write(10) G1
     write(10) G1T
     write(10) G2
     write(10) G2T

     ! Coordinates of control points
     if (dump_mesh_info_screen) &
        write(*,*)'PARALLEL DATABASE: writing coordinates/control points...',iproc

     ! global number of control nodes (slightly differs for each processor!)
     write(10) nglobmeshp(iproc)

     write(10) router * scpp(:,iproc)
     write(10) router * zcpp(:,iproc)

     ! Topology of control points
     do iel = 1, nel(iproc)
        ielg = procel(iel,iproc)
        if (zcom(ielg) >= 0.) then ! NORTH
           write(10) (lnodescp(iel,inode,iproc),inode=1,8)
        else ! SOUTH
           write(10) (lnodescp(iel,inode,iproc),inode=7,1,-1),lnodescp(iel,8,iproc)
        endif
     enddo

     !! Number of global distinct points (slightly differs for each processor!)
     write(10) nglobp(iproc)

     ! Element types
     write(10) eltypep(:,iproc)
     write(10) coarsingp(:,iproc)

     !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     ! Solid-Fluid distinction
     if (dump_mesh_info_screen) &
        write(*,*)'PARALLEL DATABASE: writing solid/fluid domain info...',iproc

     ! mapping from sol/flu to global element numbers
     write(10) procel_solidp(:,iproc)
     write(10) procel_fluidp(:,iproc)

     ! slocal numbering
     npointotp = nel_solid(iproc)*(npol+1)**2
     write(10) igloc_solid(1:npointotp,iproc)

     ! flocal numbering
     npointotp = nel_fluid(iproc)*(npol+1)**2
     write(10) igloc_fluid(1:npointotp,iproc)

     !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     ! Solid-Fluid boundary
     if (dump_mesh_info_screen) &
        write(*,*)'PARALLEL DATABASE: writing solid/fluid boundary info...',iproc
     write(10) have_bdry_elemp(iproc)

     if (have_bdry_elemp(iproc)) then
        write(10) bdry_solid_elp(1:nbdry_el(iproc),iproc)
        write(10) bdry_fluid_elp(1:nbdry_el(iproc),iproc)
        write(10) bdry_jpol_solidp(1:nbdry_el(iproc),iproc)
        write(10) bdry_jpol_fluidp(1:nbdry_el(iproc),iproc)
     endif

     !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     ! General numerical input/output parameters
     if (dump_mesh_info_screen) &
        write(*,*)'PARALLEL DATABASE: writing numerical parameters...',iproc
     write(10)pts_wavelngth,period,courant,dt

     !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     ! Background model
     if (dump_mesh_info_screen) &
        write(*,*)'PARALLEL DATABASE: writing background model info...',iproc
     write(10) bkgrdmodel(1:lfbkgrdmodel)
     write(10) override_ext_q
     write(10) router,have_fluid
     do idom=1,ndisc
        write(10) discont(idom),solid_domain(idom),idom_fluid(idom)
     enddo
     write(10)rmin,minh_ic,maxh_ic,maxh_icb

     !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     ! Min/max grid spacing
     write(10)hmin_glob,hmax_glob
     write(10)min_distance_dim,min_distance_nondim

     !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     ! critical ratios h/v min/max and locations
     write(10)char_time_max,char_time_max_globel
     write(10)char_time_max_rad,char_time_max_theta
     write(10)char_time_min,char_time_min_globel
     write(10)char_time_min_rad,char_time_min_theta

     !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     ! Axial element arrays
     write(10)naxelp(iproc),naxel_solidp(iproc),naxel_fluidp(iproc)

     write(10) ax_elp(1:naxelp(iproc),iproc)
     write(10) ax_el_solidp(1:naxel_solidp(iproc),iproc)
     write(10) ax_el_fluidp(1:naxel_fluidp(iproc),iproc)

     !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     if (dump_mesh_info_screen) write(*,*)'PARALLEL DATABASE: writing communication info...',iproc

     ! SSSSSSSSSSSS SOLID MESSAGING SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
     ! number of processors to communicate with
     write(10) sizerecvp_solid(iproc)

     if (sizerecvp_solid(iproc) > 0 ) then
        ! list of processors to communicate with
        write(10) listrecvp_solid(1:sizerecvp_solid(iproc),iproc)

        ! size of messages coming from each of these processors
        write(10) sizemsgrecvp_solid(1:sizerecvp_solid(iproc),iproc)

        ! glocal_index corresponding to each of these messages
        do imsg = 1, sizerecvp_solid(iproc)
           ipsrc = listrecvp_solid(imsg,iproc) ! To be deleted
           do iptp = 1, sizemsgrecvp_solid(imsg,iproc)
              write(10) glocal_index_msg_recvp_solid(iptp,imsg,iproc)
           enddo
        enddo
     endif

     ! FFFFFFFFFFFFF FLUID MESSAGING FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
     ! Only written if there is a fluid!
     if (have_fluid) then
        if (dump_mesh_info_screen) write(*,*)'writing fluid messaging'
        ! number of processors to communicate with
        write(10) sizerecvp_fluid(iproc)

        if (sizerecvp_fluid(iproc) > 0 ) then
           ! list of processors to communicate with
           write(10) listrecvp_fluid(1:sizerecvp_fluid(iproc),iproc)

           ! size of messages for each of these processors
           write(10) sizemsgrecvp_fluid(1:sizerecvp_fluid(iproc),iproc)

           ! glocal_index corresponding to each of these messages
           do imsg = 1, sizerecvp_fluid(iproc)
              ipsrc = listrecvp_fluid(imsg,iproc) ! To be deleted
              do iptp = 1, sizemsgrecvp_fluid(imsg,iproc)
                 write(10) glocal_index_msg_recvp_fluid(iptp,imsg,iproc)
              enddo
           enddo
        endif

     endif ! have_fluid

     close(10)
     write(*,*)'....Wrote database for processor',iproc; call flush(6)

  enddo

end subroutine write_db
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine create_static_header

    integer           :: iproc
    character(len=8)  :: mydate
    character(len=10) :: mytime
    character(len=80) :: dbname

    call date_and_time(mydate,mytime)

    dbname = 'mesh_params.h'
    open(97,file=trim(dbname), STATUS="REPLACE")
    iproc = 0
    write(97,10) nproc
    write(97,11) mydate(5:6),mydate(7:8),mydate(1:4),mytime(1:2),mytime(3:4)
    write(97,*)''
    write(97,29)
    write(97,12)'Background model     :',bkgrdmodel(1:lfbkgrdmodel)
    write(97,14)'Dominant period [s]  :',period
    write(97,14)'Elements/wavelength  :',pts_wavelngth
    write(97,14)'Courant number       :',courant
    write(97,15)'Coarsening levels    :',nc_init
    write(97,30)
    write(97,*)''
    write(97,9)'npol',npol,'polynomial order'
    write(97,9)'nelem',nel(iproc),'proc. els'
    write(97,9)'npoint',nel(iproc)*(npol+1)**2,'proc. all pts'
    write(97,9)'nel_solid',nel_solid(iproc),'proc. solid els'
    write(97,9)'nel_fluid',nel_fluid(iproc),'proc. fluid els'
    write(97,9)'npoint_solid',nel_solid(iproc)*(npol+1)**2,'proc. solid pts'
    write(97,9)'npoint_fluid',nel_fluid(iproc)*(npol+1)**2,'proc. fluid pts'
    write(97,9)'nglob_fluid',nglobp_fluid(iproc),'proc. flocal pts'
    write(97,9)'nel_bdry',nbdry_el(iproc),'proc. solid-fluid bndry els'
    write(97,9)'ndisc',ndisc,'# disconts in bkgrd model'
    write(97,9)'nproc_mesh',nproc,'number of processors'
    write(97,9)'lfbkgrdmodel',lfbkgrdmodel,'length of bkgrdmodel name'

    write(97,*)''
    write(97,31)
    write(97,14)'Time step [s]        :',dt
    write(97,16)'Min(h/vp),dt/courant :',minhvp,dt/courant*real(npol)
    write(97,16)'max(h/vs),T0/wvlngth :',maxhvs,period/pts_wavelngth
    write(97,14)'Inner core r_min [km]:',rmin/1.d3
    write(97,14)'Max(h) r/ns(icb) [km]:',maxhnsicb/1.d3
    write(97,14)'Max(h) precalc.  [km]:',maxh_icb/1.d3
    write(97,30)
    write(97,*)''
    close(97)
    write(*,*)'wrote parameters for static solver into ', trim(dbname)

9 format(' integer, parameter :: ',A12,' =',i10,'  ! ',A27)
10 format('! Proc ',i3,': Header for mesh information to run static solver')
11 format('! created by the mesher on ', &
            A2,'/',A2,'/',A4,', at ',A2,'h ',A2,'min')
29 format('!:::::::::::::::::::: Input parameters :::::::::::::::::::::::::::')
12 format('!  ',A23,A20)
14 format('!  ',A23,1f10.4)
16 format('!  ',A23,2(f10.4))
15 format('!  ',A23,I10)
30 format('!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::')

31 format('!:::::::::::::::::::: Output parameters ::::::::::::::::::::::::::')


end subroutine create_static_header
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine define_io_appendix(app, iproc)
! Defines the 4 digit character string appended to any
! data or io file related to process myid.

  implicit none
  integer, intent(in)           :: iproc
  character(len=4), intent(out) :: app

  write(app,"(I4.4)") iproc

end subroutine define_io_appendix
!-----------------------------------------------------------------------------------------

end module pdb
!=========================================================================================
