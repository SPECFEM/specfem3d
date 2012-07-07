!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  2 . 1
!               ---------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Princeton University, USA and CNRS / INRIA / University of Pau
! (c) Princeton University / California Institute of Technology and CNRS / INRIA / University of Pau
!                             July 2012
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


! define sets of colors that contain disconnected elements for the CUDA solver.
! also split the elements into two subsets: inner and outer elements, in order
! to be able to compute the outer elements first in the solver and then
! start non-blocking MPI calls and overlap them with the calculation of the inner elements
! (which works fine because there are always far more inner elements than outer elements)

! note:    these are modified routines to use element domain flags given in ispec_is_d, thus
!             coloring only acoustic or elastic (or..) elements in one run, then repeat run for other domains.
!             also, the permutation re-starts at 1 for outer and for inner elements,
!             making it usable for the phase_ispec_inner_** arrays for acoustic and elastic elements.

  subroutine get_perm_color_faster(is_on_a_slice_edge,ispec_is_d, &
                                  ibool,perm,color, &
                                  nspec,nglob, &
                                  nb_colors_outer_elements,nb_colors_inner_elements, &
                                  nspec_outer,nspec_inner,nspec_domain, &
                                  first_elem_number_in_this_color, &
                                  myrank)

  implicit none

  include "constants.h"

  integer, intent(in) :: nspec, nglob
  logical, dimension(nspec), intent(in) :: is_on_a_slice_edge
  logical, dimension(nspec), intent(in) :: ispec_is_d

  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec), intent(in) :: ibool
  integer, dimension(nspec),intent(inout) :: perm

  integer, dimension(nspec),intent(inout) :: color
  integer, dimension(MAX_NUMBER_OF_COLORS+1),intent(inout) :: first_elem_number_in_this_color
  integer, intent(out) :: nb_colors_outer_elements,nb_colors_inner_elements

  integer, intent(out) :: nspec_outer,nspec_inner,nspec_domain
  integer, intent(in) :: myrank

  ! local variables
  integer :: nb_colors

  ! coloring algorithm w/ Droux
  call get_color_faster(ibool, is_on_a_slice_edge, ispec_is_d, &
                        myrank, nspec, nglob, &
                        color, nb_colors_outer_elements, nb_colors_inner_elements, &
                        nspec_outer,nspec_inner,nspec_domain)

  !debug output
  if(myrank == 0) then
    write(IMAIN,*) '     colors:'
    write(IMAIN,*) '     number of colors for inner elements = ',nb_colors_inner_elements
    write(IMAIN,*) '     number of colors for outer elements = ',nb_colors_outer_elements
    write(IMAIN,*) '     total number of colors (sum of both) = ', nb_colors_inner_elements + nb_colors_outer_elements
    write(IMAIN,*) '     elements:'
    write(IMAIN,*) '     number of elements for outer elements  = ',nspec_outer
    write(IMAIN,*) '     number of elements for inner elements  = ',nspec_inner
    write(IMAIN,*) '     total number of elements for domain elements  = ',nspec_domain
  endif

  ! total number of colors used
  nb_colors = nb_colors_inner_elements+nb_colors_outer_elements
  first_elem_number_in_this_color(:) = 0

  ! gets element permutation depending on colors
  call get_final_perm(color,perm,first_elem_number_in_this_color(1:nb_colors), &
                     nspec,nb_colors,nb_colors_outer_elements, &
                     ispec_is_d,nspec_domain)


  end subroutine get_perm_color_faster

!
!-------------------------------------------------------------------------------------------------
!

  subroutine get_color_faster(ibool, is_on_a_slice_edge, ispec_is_d, &
                             myrank, nspec, nglob, &
                             color, nb_colors_outer_elements, nb_colors_inner_elements, &
                             nspec_outer,nspec_inner,nspec_domain)

  implicit none

  include "constants.h"

  integer nspec,nglob
  logical, dimension(nspec) :: is_on_a_slice_edge,ispec_is_d

  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool
  integer, dimension(nspec) :: color
  integer :: nb_colors_outer_elements,nb_colors_inner_elements,myrank

  integer :: nspec_outer,nspec_inner,nspec_domain

  ! local variables
  integer :: ispec
  logical, dimension(:), allocatable :: mask_ibool
  integer :: icolor, nb_already_colored
  integer :: iglob1,iglob2,iglob3,iglob4,iglob5,iglob6,iglob7,iglob8
  integer :: ier
  logical :: conflict_found_need_new_color
  ! Droux
  logical :: try_Droux_coloring
  logical :: fail_safe
  ! valence
  integer :: maxval_count_ibool_outer,maxval_count_ibool_inner

  ! display absolute minimum possible number of colors, i.e., maximum valence (for information only)
  ! beware: this wastes memory (needs an additional array called "count_ibool")
  logical, parameter :: DISPLAY_MIN_POSSIBLE_COLORS = .false.

  ! user output
  if( myrank == 0 ) then
    if( USE_DROUX_OPTIMIZATION ) then
      write(IMAIN,*) '     fast coloring mesh algorithm w/ Droux optimization'
    else if( BALANCE_COLORS_SIMPLE_ALGO ) then
      write(IMAIN,*) '     fast coloring mesh algorithm w/ color balancing'
    else
      write(IMAIN,*) '     fast coloring mesh algorithm'
    endif
  endif

  ! counts number of elements for inner, outer and total domain
  nspec_outer = 0
  nspec_inner = 0
  nspec_domain = 0
  do ispec=1,nspec
    ! domain elements
    if(ispec_is_d(ispec)) then
      ! outer/inner elements
      if(is_on_a_slice_edge(ispec)) then
        nspec_outer=nspec_outer+1
      else
        nspec_inner=nspec_inner+1
      endif
      nspec_domain=nspec_domain+1
    endif
  enddo

  ! debug
  !if(myrank == 0) then
  !  print *
  !  print *,'----------------------------------'
  !  print *,'coloring the mesh'
  !  print *,'----------------------------------'
  !  print *
  !endif

  ! Droux optimization
  try_Droux_coloring = USE_DROUX_OPTIMIZATION

  if(BALANCE_COLORS_SIMPLE_ALGO .and. USE_DROUX_OPTIMIZATION ) then
    if( myrank == 0 ) then
      print *,'noticed a problem with mesh coloring options: '
      print *,'  cannot set both USE_DROUX_OPTIMAL_ALGO and BALANCE_COLORS_SIMPLE_ALGO'
      print *,'  -> this run will use only BALANCE_COLORS_SIMPLE_ALGO'
      print *,'please check parameter settings in constants.h...'
    endif
    try_Droux_coloring = .false.
  endif

  ! gives a lower bound for the number of colors needed
  if(DISPLAY_MIN_POSSIBLE_COLORS .or. try_Droux_coloring) then
    ! gets maximum values of valence for inner and outer element points
    call count_mesh_valence(ibool,is_on_a_slice_edge,ispec_is_d, &
                           myrank, nspec, nglob, &
                           maxval_count_ibool_outer,maxval_count_ibool_inner)
  endif

  ! allocates mask
  allocate(mask_ibool(nglob),stat=ier)
  if( ier /= 0 ) stop 'error allocating mask_ibool array'

  ! entry point for fail-safe mechanism when Droux 1993 fails
  999 continue

  ! first set color of all elements to 0,
  ! to use it as a flag to detect elements not yet colored
  color(:) = 0
  icolor = 0
  nb_already_colored = 0

  ! colors outer elements
  do while( nb_already_colored < nspec_outer )

    333 continue
    icolor = icolor + 1

    ! debug: user output
    !if(myrank == 0) then
    !  print *,'  analyzing color ',icolor,' - outer elements'
    !endif

    ! resets flags
    mask_ibool(:) = .false.
    conflict_found_need_new_color = .false.

    ! finds un-colored elements
    do ispec = 1,nspec
      ! domain elements only
      if( ispec_is_d(ispec) ) then
        ! outer elements
        if( is_on_a_slice_edge(ispec) ) then
          if(color(ispec) == 0) then
            ! the eight corners of the current element
            iglob1=ibool(1,1,1,ispec)
            iglob2=ibool(NGLLX,1,1,ispec)
            iglob3=ibool(NGLLX,NGLLY,1,ispec)
            iglob4=ibool(1,NGLLY,1,ispec)
            iglob5=ibool(1,1,NGLLZ,ispec)
            iglob6=ibool(NGLLX,1,NGLLZ,ispec)
            iglob7=ibool(NGLLX,NGLLY,NGLLZ,ispec)
            iglob8=ibool(1,NGLLY,NGLLZ,ispec)

            if(mask_ibool(iglob1) .or. mask_ibool(iglob2) .or. mask_ibool(iglob3) .or. mask_ibool(iglob4) .or. &
               mask_ibool(iglob5) .or. mask_ibool(iglob6) .or. mask_ibool(iglob7) .or. mask_ibool(iglob8)) then
              ! if element of this color has a common point with another element of that same color
              ! then we need to create a new color, i.e., increment the color of the current element
              conflict_found_need_new_color = .true.
            else
              color(ispec) = icolor
              nb_already_colored = nb_already_colored + 1
              mask_ibool(iglob1) = .true.
              mask_ibool(iglob2) = .true.
              mask_ibool(iglob3) = .true.
              mask_ibool(iglob4) = .true.
              mask_ibool(iglob5) = .true.
              mask_ibool(iglob6) = .true.
              mask_ibool(iglob7) = .true.
              mask_ibool(iglob8) = .true.
            endif
          endif
        endif
      endif
    enddo

    ! debug: user output
    !if(myrank == 0) then
    !  print *,'  done ',(100.0*nb_already_colored)/nspec_domain,'% of ',nspec_domain,'elements'
    !endif

    if(conflict_found_need_new_color) then
      if( icolor >= MAX_NUMBER_OF_COLORS ) stop 'error MAX_NUMBER_OF_COLORS too small'
      goto 333
    endif
  enddo

  nb_colors_outer_elements = icolor

  ! colors inner elements
  do while(nb_already_colored < nspec_domain)

    334 continue
    icolor = icolor + 1

    ! debug: user output
    !if(myrank == 0) then
    !  print *,'  analyzing color ',icolor,' - inner elements'
    !endif

    ! resets flags
    mask_ibool(:) = .false.
    conflict_found_need_new_color = .false.

    do ispec = 1,nspec
      ! domain elements only
      if(ispec_is_d(ispec)) then
        ! inner elements
        if (.not. is_on_a_slice_edge(ispec)) then
          if(color(ispec) == 0) then
            ! the eight corners of the current element
            iglob1=ibool(1,1,1,ispec)
            iglob2=ibool(NGLLX,1,1,ispec)
            iglob3=ibool(NGLLX,NGLLY,1,ispec)
            iglob4=ibool(1,NGLLY,1,ispec)
            iglob5=ibool(1,1,NGLLZ,ispec)
            iglob6=ibool(NGLLX,1,NGLLZ,ispec)
            iglob7=ibool(NGLLX,NGLLY,NGLLZ,ispec)
            iglob8=ibool(1,NGLLY,NGLLZ,ispec)

            if(mask_ibool(iglob1) .or. mask_ibool(iglob2) .or. mask_ibool(iglob3) .or. mask_ibool(iglob4) .or. &
               mask_ibool(iglob5) .or. mask_ibool(iglob6) .or. mask_ibool(iglob7) .or. mask_ibool(iglob8)) then
              ! if element of this color has a common point with another element of that same color
              ! then we need to create a new color, i.e., increment the color of the current element
              conflict_found_need_new_color = .true.
            else
              color(ispec) = icolor
              nb_already_colored = nb_already_colored + 1
              mask_ibool(iglob1) = .true.
              mask_ibool(iglob2) = .true.
              mask_ibool(iglob3) = .true.
              mask_ibool(iglob4) = .true.
              mask_ibool(iglob5) = .true.
              mask_ibool(iglob6) = .true.
              mask_ibool(iglob7) = .true.
              mask_ibool(iglob8) = .true.
            endif
          endif
        endif
      endif
    enddo

    ! debug user output
    !if(myrank == 0) then
    !  print *,'  done ',(100.0*nb_already_colored)/nspec_domain,'% of ',nspec_domain,'elements'
    !endif

    if(conflict_found_need_new_color) then
      if( icolor >= MAX_NUMBER_OF_COLORS ) stop 'error MAX_NUMBER_OF_COLORS too small'
      goto 334
    endif
  enddo

  nb_colors_inner_elements = icolor - nb_colors_outer_elements

  ! Droux optimization:
  ! added this to create more balanced colors according to JJ Droux (1993)
  ! note: this might not find an optimial solution.
  !          we will probably have to try a few times with increasing colors
  if( try_Droux_coloring ) then
    ! initializes fail-safe mechanism
    fail_safe = .false.

    ! tries to find a balanced coloring
    call balance_colors_Droux(ibool,is_on_a_slice_edge,ispec_is_d, &
                              myrank, nspec, nglob, &
                              color,nb_colors_outer_elements,nb_colors_inner_elements, &
                              nspec_outer,nspec_inner,maxval_count_ibool_inner, &
                              mask_ibool,fail_safe)

    ! in case it fails go back to simple coloring algorithm
    if( fail_safe ) then
      try_Droux_coloring = .false.
      if(myrank == 0) write(IMAIN,*) '     giving up on Droux 1993 algorithm, calling fail-safe mechanism'
      goto 999
    endif
  endif ! of if(try_Droux_coloring)

  ! balances colors using a simple algorithm (if Droux was not used)
  if( BALANCE_COLORS_SIMPLE_ALGO ) then
    call balance_colors_simple(ibool,is_on_a_slice_edge,ispec_is_d, &
                              myrank, nspec, nglob, &
                              color,nb_colors_outer_elements,nb_colors_inner_elements, &
                              nspec_outer,nspec_inner,mask_ibool)
  endif

  ! checks that all the color sets are independent
  do icolor = 1,maxval(color)
    mask_ibool(:) = .false.
    do ispec = 1,nspec
      ! domain elements only
      if(ispec_is_d(ispec)) then
        if(color(ispec) == icolor ) then
          ! the eight corners of the current element
          iglob1=ibool(1,1,1,ispec)
          iglob2=ibool(NGLLX,1,1,ispec)
          iglob3=ibool(NGLLX,NGLLY,1,ispec)
          iglob4=ibool(1,NGLLY,1,ispec)
          iglob5=ibool(1,1,NGLLZ,ispec)
          iglob6=ibool(NGLLX,1,NGLLZ,ispec)
          iglob7=ibool(NGLLX,NGLLY,NGLLZ,ispec)
          iglob8=ibool(1,NGLLY,NGLLZ,ispec)

          if(mask_ibool(iglob1) .or. mask_ibool(iglob2) .or. mask_ibool(iglob3) .or. mask_ibool(iglob4) .or. &
             mask_ibool(iglob5) .or. mask_ibool(iglob6) .or. mask_ibool(iglob7) .or. mask_ibool(iglob8)) then
            ! if element of this color has a common point with another element of that same color
            ! then there is a problem, the color set is not correct
            print*,'error check color:',icolor
            stop 'error detected: found a common point inside a color set'
          else
            mask_ibool(iglob1) = .true.
            mask_ibool(iglob2) = .true.
            mask_ibool(iglob3) = .true.
            mask_ibool(iglob4) = .true.
            mask_ibool(iglob5) = .true.
            mask_ibool(iglob6) = .true.
            mask_ibool(iglob7) = .true.
            mask_ibool(iglob8) = .true.
          endif
        endif
      endif
    enddo

    !debug output
    !if(myrank == 0) print *,'  color ',icolor,' has disjoint elements only and is therefore OK'
    !if(myrank == 0) print *,'  color ',icolor,' contains ',count(color == icolor),' elements'
  enddo
  ! debug output
  !if(myrank == 0) then
  !  print*, '     the ',maxval(color),' color sets are OK'
  !endif

  deallocate(mask_ibool)

  end subroutine get_color_faster

!
!-------------------------------------------------------------------------------------------------
!

  subroutine count_mesh_valence(ibool,is_on_a_slice_edge,ispec_is_d, &
                                myrank, nspec, nglob, &
                                maxval_count_ibool_outer,maxval_count_ibool_inner)

  implicit none

  include "constants.h"

  integer :: nspec,nglob

  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool

  logical, dimension(nspec) :: is_on_a_slice_edge,ispec_is_d

  integer :: myrank
  integer :: maxval_count_ibool_outer,maxval_count_ibool_inner

  ! local parameters
  integer, dimension(:), allocatable :: count_ibool
  integer :: ispec
  integer :: iglob1,iglob2,iglob3,iglob4,iglob5,iglob6,iglob7,iglob8
  integer :: ier

  ! allocates count array
  allocate(count_ibool(nglob),stat=ier)
  if( ier /= 0 ) stop 'error allocating count_ibool array'

  ! valence numbers of the mesh
  maxval_count_ibool_outer = 0
  maxval_count_ibool_inner = 0

  ! valence for outer elements
  count_ibool(:) = 0
  do ispec = 1,nspec
    ! domain elements only
    if(ispec_is_d(ispec)) then
      ! outer elements
      if (is_on_a_slice_edge(ispec)) then
        ! the eight corners of the current element
        iglob1=ibool(1,1,1,ispec)
        iglob2=ibool(NGLLX,1,1,ispec)
        iglob3=ibool(NGLLX,NGLLY,1,ispec)
        iglob4=ibool(1,NGLLY,1,ispec)
        iglob5=ibool(1,1,NGLLZ,ispec)
        iglob6=ibool(NGLLX,1,NGLLZ,ispec)
        iglob7=ibool(NGLLX,NGLLY,NGLLZ,ispec)
        iglob8=ibool(1,NGLLY,NGLLZ,ispec)

        count_ibool(iglob1) = count_ibool(iglob1) + 1
        count_ibool(iglob2) = count_ibool(iglob2) + 1
        count_ibool(iglob3) = count_ibool(iglob3) + 1
        count_ibool(iglob4) = count_ibool(iglob4) + 1
        count_ibool(iglob5) = count_ibool(iglob5) + 1
        count_ibool(iglob6) = count_ibool(iglob6) + 1
        count_ibool(iglob7) = count_ibool(iglob7) + 1
        count_ibool(iglob8) = count_ibool(iglob8) + 1
      endif
    endif
  enddo
  maxval_count_ibool_outer = maxval(count_ibool)

  ! valence for inner elements
  count_ibool(:) = 0
  do ispec = 1,nspec
    ! domain elements only
    if(ispec_is_d(ispec)) then
      ! inner elements
      if (.not. is_on_a_slice_edge(ispec)) then
        ! the eight corners of the current element
        iglob1=ibool(1,1,1,ispec)
        iglob2=ibool(NGLLX,1,1,ispec)
        iglob3=ibool(NGLLX,NGLLY,1,ispec)
        iglob4=ibool(1,NGLLY,1,ispec)
        iglob5=ibool(1,1,NGLLZ,ispec)
        iglob6=ibool(NGLLX,1,NGLLZ,ispec)
        iglob7=ibool(NGLLX,NGLLY,NGLLZ,ispec)
        iglob8=ibool(1,NGLLY,NGLLZ,ispec)

        count_ibool(iglob1) = count_ibool(iglob1) + 1
        count_ibool(iglob2) = count_ibool(iglob2) + 1
        count_ibool(iglob3) = count_ibool(iglob3) + 1
        count_ibool(iglob4) = count_ibool(iglob4) + 1
        count_ibool(iglob5) = count_ibool(iglob5) + 1
        count_ibool(iglob6) = count_ibool(iglob6) + 1
        count_ibool(iglob7) = count_ibool(iglob7) + 1
        count_ibool(iglob8) = count_ibool(iglob8) + 1
      endif
    endif
  enddo
  maxval_count_ibool_inner = maxval(count_ibool)

  ! debug outupt
  if( myrank == 0 ) then
    write(IMAIN,*) '     maximum valence (i.e. minimum possible nb of colors) for outer = ',maxval_count_ibool_outer
    write(IMAIN,*) '     maximum valence (i.e. minimum possible nb of colors) for inner = ',maxval_count_ibool_inner
  endif

  deallocate(count_ibool)

  end subroutine count_mesh_valence

!
!-------------------------------------------------------------------------------------------------
!

  subroutine balance_colors_Droux(ibool,is_on_a_slice_edge,ispec_is_d, &
                                  myrank, nspec, nglob, &
                                  color, nb_colors_outer_elements, nb_colors_inner_elements, &
                                  nspec_outer,nspec_inner,maxval_count_ibool_inner, &
                                  mask_ibool,fail_safe)

  implicit none

  include "constants.h"

  integer :: nspec,nglob

  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool

  logical, dimension(nspec) :: is_on_a_slice_edge,ispec_is_d
  logical, dimension(nglob) :: mask_ibool

  integer, dimension(nspec) :: color

  integer :: myrank
  integer :: nb_colors_outer_elements,nb_colors_inner_elements

  integer :: nspec_outer,nspec_inner
  integer :: maxval_count_ibool_inner

  logical :: fail_safe

  ! local parameters
  logical, dimension(:), allocatable :: icolor_conflict_found
  integer, dimension(:), allocatable :: nb_elems_in_this_color
  integer :: ispec,ispec2,icolor,ncolors,icolormin,icolormax,icolor_chosen,nb_elems_in_color_chosen
  integer :: nb_tries_of_Droux_1993,last_ispec_studied
  integer :: ier

  ! debug outupt
  if( myrank == 0 ) then
    write(IMAIN,*) '     balancing colors: Droux algorithm'
    write(IMAIN,*) '       initial number of outer element colors = ',nb_colors_outer_elements
    write(IMAIN,*) '       initial number of inner element colors = ',nb_colors_inner_elements
    write(IMAIN,*) '       initial number of total colors = ',nb_colors_outer_elements + nb_colors_inner_elements
  endif

  ! initial guess of number of colors needed
  if( maxval_count_ibool_inner > 0 .and. maxval_count_ibool_inner < nb_colors_inner_elements ) then
    ! uses maximum valence to estimate number of colors for Droux
    nb_colors_inner_elements = maxval_count_ibool_inner
  endif

  !! DK DK do it for inner elements only for now
  ! Droux optimization run
  nb_tries_of_Droux_1993 = 1

  ! entry point to re-try Droux
  765 continue

  ! initial guess of number of colors needed
  ncolors = nb_colors_outer_elements + nb_colors_inner_elements

  ! debug output
  if( myrank == 0 ) then
    write(IMAIN,*) '     Droux optimization: try = ',nb_tries_of_Droux_1993,'colors = ',ncolors
  endif

  icolormin = nb_colors_outer_elements + 1
  icolormax = ncolors

  ! allocates temporary arrays
  allocate(nb_elems_in_this_color(ncolors), &
          icolor_conflict_found(ncolors),stat=ier)
  if( ier /= 0 ) stop 'error allocating nb_elems_in_this_color arrays'

  nb_elems_in_this_color(:) = 0
  mask_ibool(:) = .false.
  last_ispec_studied = -1

  do ispec = 1,nspec
    ! domain elements only
    if(ispec_is_d(ispec)) then

      ! only inner elements
      if (is_on_a_slice_edge(ispec)) cycle

      ! unmark the eight corners of the previously marked element
      if(last_ispec_studied > 0) then
        mask_ibool(ibool(1,1,1,last_ispec_studied)) = .false.
        mask_ibool(ibool(NGLLX,1,1,last_ispec_studied)) = .false.
        mask_ibool(ibool(NGLLX,NGLLY,1,last_ispec_studied)) = .false.
        mask_ibool(ibool(1,NGLLY,1,last_ispec_studied)) = .false.
        mask_ibool(ibool(1,1,NGLLZ,last_ispec_studied)) = .false.
        mask_ibool(ibool(NGLLX,1,NGLLZ,last_ispec_studied)) = .false.
        mask_ibool(ibool(NGLLX,NGLLY,NGLLZ,last_ispec_studied)) = .false.
        mask_ibool(ibool(1,NGLLY,NGLLZ,last_ispec_studied)) = .false.
      endif
      icolor_conflict_found(icolormin:icolormax) = .false.

      ! mark the eight corners of the current element
      mask_ibool(ibool(1,1,1,ispec)) = .true.
      mask_ibool(ibool(NGLLX,1,1,ispec)) = .true.
      mask_ibool(ibool(NGLLX,NGLLY,1,ispec)) = .true.
      mask_ibool(ibool(1,NGLLY,1,ispec)) = .true.
      mask_ibool(ibool(1,1,NGLLZ,ispec)) = .true.
      mask_ibool(ibool(NGLLX,1,NGLLZ,ispec)) = .true.
      mask_ibool(ibool(NGLLX,NGLLY,NGLLZ,ispec)) = .true.
      mask_ibool(ibool(1,NGLLY,NGLLZ,ispec)) = .true.
      last_ispec_studied = ispec

      if(ispec > 1) then
        do ispec2 = 1,ispec - 1
          ! domain elements only
          if(ispec_is_d(ispec2)) then

            ! only inner elements
            if (is_on_a_slice_edge(ispec2)) cycle

            ! if conflict already found previously with this color, no need to test again
            if (icolor_conflict_found(color(ispec2))) cycle

            ! test the eight corners of the current element for a common point with element under study
            if (mask_ibool(ibool(1,1,1,ispec2)) .or. &
                mask_ibool(ibool(NGLLX,1,1,ispec2)) .or. &
                mask_ibool(ibool(NGLLX,NGLLY,1,ispec2)) .or. &
                mask_ibool(ibool(1,NGLLY,1,ispec2)) .or. &
                mask_ibool(ibool(1,1,NGLLZ,ispec2)) .or. &
                mask_ibool(ibool(NGLLX,1,NGLLZ,ispec2)) .or. &
                mask_ibool(ibool(NGLLX,NGLLY,NGLLZ,ispec2)) .or. &
                mask_ibool(ibool(1,NGLLY,NGLLZ,ispec2))) &
              icolor_conflict_found(color(ispec2)) = .true.

          endif ! domain elements
        enddo
      endif

      ! check if the Droux 1993 algorithm found a solution
      if (all(icolor_conflict_found(icolormin:icolormax))) then
        ! user output
        !if(myrank == 0) write(IMAIN,*) '     Droux 1993 algorithm did not find any solution for ncolors = ',ncolors

        ! try with one more color
        if(nb_tries_of_Droux_1993 < MAX_NB_TRIES_OF_DROUX_1993) then
          nb_colors_inner_elements = nb_colors_inner_elements + 1
          deallocate(nb_elems_in_this_color)
          deallocate(icolor_conflict_found)
          nb_tries_of_Droux_1993 = nb_tries_of_Droux_1993 + 1
          goto 765
        else
          ! fail-safe mechanism: if Droux 1993 still fails after all the tries with one more color,
          ! then go back to my original simple and fast coloring algorithm
          fail_safe = .true.
          return
        endif
      endif

      ! loop on all the colors to determine the color with the smallest number
      ! of elements and for which there is no conflict
      nb_elems_in_color_chosen = 2147000000 ! start with extremely large unrealistic value
      icolor_chosen = 0
      do icolor = icolormin,icolormax
        if (.not. icolor_conflict_found(icolor) .and. nb_elems_in_this_color(icolor) < nb_elems_in_color_chosen) then
          icolor_chosen = icolor
          nb_elems_in_color_chosen = nb_elems_in_this_color(icolor)
        endif
      enddo

      ! store the color finally chosen
      color(ispec) = icolor_chosen
      nb_elems_in_this_color(icolor_chosen) = nb_elems_in_this_color(icolor_chosen) + 1

    endif ! domain elements
  enddo

  ! debug output
  if(myrank == 0) then
    write(IMAIN,*) '     created a total of ',maxval(color),' colors in this domain' ! 'for all the domain elements of the mesh'
    if( nb_colors_outer_elements > 0 ) &
      write(IMAIN,*) '     typical nb of elements per color for outer elements should be ', &
        nspec_outer / nb_colors_outer_elements
    if( nb_colors_inner_elements > 0 ) &
      write(IMAIN,*) '     typical nb of elements per color for inner elements should be ', &
        nspec_inner / nb_colors_inner_elements
  endif


  end subroutine balance_colors_Droux

!
!-------------------------------------------------------------------------------------------------
!

  subroutine balance_colors_simple(ibool,is_on_a_slice_edge,ispec_is_d, &
                                  myrank, nspec, nglob, &
                                  color, nb_colors_outer_elements, nb_colors_inner_elements, &
                                  nspec_outer,nspec_inner,mask_ibool)

  implicit none

  include "constants.h"

  integer :: nspec,nglob

  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool

  logical, dimension(nspec) :: is_on_a_slice_edge,ispec_is_d
  logical, dimension(nglob) :: mask_ibool

  integer, dimension(nspec) :: color

  integer :: myrank
  integer :: nb_colors_outer_elements,nb_colors_inner_elements

  integer :: nspec_outer,nspec_inner

  ! local parameters
  logical, dimension(:), allocatable :: icolor_conflict_found
  integer, dimension(:), allocatable :: nb_elems_in_this_color
  integer :: ispec,ispec2,icolor,ncolors,icolormin,icolormax,icolor_chosen,nb_elems_in_color_chosen
  integer :: last_ispec_studied
  integer :: target_nb_elems_per_color,icolor_target
  integer :: ier

  ! debug outupt
  if( myrank == 0 ) then
    write(IMAIN,*) '     balancing colors: simple algorithm'
    write(IMAIN,*) '       number of outer element colors = ',nb_colors_outer_elements
    write(IMAIN,*) '       number of inner element colors = ',nb_colors_inner_elements
    write(IMAIN,*) '       number of total colors = ',nb_colors_outer_elements + nb_colors_inner_elements
  endif

  ! balances colors in postprocess if Droux (1993) is not used
  ncolors = nb_colors_outer_elements + nb_colors_inner_elements

  ! allocates temporary arrays
  allocate(nb_elems_in_this_color(ncolors), &
          icolor_conflict_found(ncolors),stat=ier)
  if( ier /= 0 ) stop 'error allocating nb_elems_in_this_color arrays'

  !! DK DK do it for outer elements
  icolormin = 1
  icolormax = nb_colors_outer_elements

  ! ideal value if all colors are perfectly balanced
  if( nb_colors_outer_elements > 0 ) then
    target_nb_elems_per_color = nspec_outer / nb_colors_outer_elements + 1
  else
    target_nb_elems_per_color = 1
  endif

  ! print *,'nspec_outer,target_nb_elems_per_color = ',nspec_outer,target_nb_elems_per_color

  ! count the initial number of elements in each color
  nb_elems_in_this_color(:) = 0
  do icolor = icolormin,icolormax
    nb_elems_in_this_color(icolor) = count(color == icolor)
  enddo

  ! do not balance the last one, because it will be balanced automatically by the others
  do icolor = icolormin,icolormax-1

    ! if color is already balanced, do nothing
    ! (this works because in the initial set the number of elements per color decreases when the color number increases)
    if(nb_elems_in_this_color(icolor) <= target_nb_elems_per_color) cycle

    mask_ibool(:) = .false.
    last_ispec_studied = -1

    do ispec = 1,nspec
      ! domain elements only
      if(ispec_is_d(ispec)) then

        ! only outer elements
        if (.not. is_on_a_slice_edge(ispec)) cycle

        ! only elements of this color
        if (color(ispec) /= icolor) cycle

        ! if color is now balanced because we have moved enough elements then stop searching
        if(nb_elems_in_this_color(icolor) <= target_nb_elems_per_color) exit

        ! unmark the eight corners of the previously marked element
        if(last_ispec_studied > 0) then
          mask_ibool(ibool(1,1,1,last_ispec_studied)) = .false.
          mask_ibool(ibool(NGLLX,1,1,last_ispec_studied)) = .false.
          mask_ibool(ibool(NGLLX,NGLLY,1,last_ispec_studied)) = .false.
          mask_ibool(ibool(1,NGLLY,1,last_ispec_studied)) = .false.
          mask_ibool(ibool(1,1,NGLLZ,last_ispec_studied)) = .false.
          mask_ibool(ibool(NGLLX,1,NGLLZ,last_ispec_studied)) = .false.
          mask_ibool(ibool(NGLLX,NGLLY,NGLLZ,last_ispec_studied)) = .false.
          mask_ibool(ibool(1,NGLLY,NGLLZ,last_ispec_studied)) = .false.
        endif
        icolor_conflict_found(icolormin:icolormax) = .false.
        icolor_conflict_found(icolor) = .true. ! cannot move element to the color it already has

        ! mark the eight corners of the current element
        mask_ibool(ibool(1,1,1,ispec)) = .true.
        mask_ibool(ibool(NGLLX,1,1,ispec)) = .true.
        mask_ibool(ibool(NGLLX,NGLLY,1,ispec)) = .true.
        mask_ibool(ibool(1,NGLLY,1,ispec)) = .true.
        mask_ibool(ibool(1,1,NGLLZ,ispec)) = .true.
        mask_ibool(ibool(NGLLX,1,NGLLZ,ispec)) = .true.
        mask_ibool(ibool(NGLLX,NGLLY,NGLLZ,ispec)) = .true.
        mask_ibool(ibool(1,NGLLY,NGLLZ,ispec)) = .true.
        last_ispec_studied = ispec

        ! test if we can move this element to another color
        do ispec2 = 1,nspec
          ! domain elements only
          if(ispec_is_d(ispec2)) then

            ! do not test that element itself
            if (ispec2 == ispec) cycle

            ! only outer elements
            if (.not. is_on_a_slice_edge(ispec2)) cycle

            ! if conflict already found previously with this color, no need to test again
            if (icolor_conflict_found(color(ispec2))) cycle

            ! test the eight corners of the current element for a common point with element under study
            if (mask_ibool(ibool(1,1,1,ispec2)) .or. &
              mask_ibool(ibool(NGLLX,1,1,ispec2)) .or. &
              mask_ibool(ibool(NGLLX,NGLLY,1,ispec2)) .or. &
              mask_ibool(ibool(1,NGLLY,1,ispec2)) .or. &
              mask_ibool(ibool(1,1,NGLLZ,ispec2)) .or. &
              mask_ibool(ibool(NGLLX,1,NGLLZ,ispec2)) .or. &
              mask_ibool(ibool(NGLLX,NGLLY,NGLLZ,ispec2)) .or. &
              mask_ibool(ibool(1,NGLLY,NGLLZ,ispec2))) &
                icolor_conflict_found(color(ispec2)) = .true.
          endif ! domain elements
        enddo

        ! if color is already above target size for a balanced set, do not move to it
        do icolor_target = icolormin,icolormax
          if(nb_elems_in_this_color(icolor_target) >= target_nb_elems_per_color) &
            icolor_conflict_found(icolor_target) = .true.
        enddo

        ! if cannot find any other color to move this element to
        if (all(icolor_conflict_found(icolormin:icolormax))) cycle

        ! loop on all the colors to determine the color with the smallest number of elements
        ! and for which there is no conflict
        nb_elems_in_color_chosen = 2147000000 ! start with extremely large unrealistic value
        icolor_chosen = 0
        do icolor_target = icolormin,icolormax
          if (.not. icolor_conflict_found(icolor_target) .and. &
             nb_elems_in_this_color(icolor_target) < nb_elems_in_color_chosen) then
            icolor_chosen = icolor_target
            nb_elems_in_color_chosen = nb_elems_in_this_color(icolor_target)
          endif
        enddo

        ! move the element to that new color
        ! remove element from its current color
        nb_elems_in_this_color(color(ispec)) = nb_elems_in_this_color(color(ispec)) - 1
        color(ispec) = icolor_chosen
        ! and add it to the new color
        nb_elems_in_this_color(icolor_chosen) = nb_elems_in_this_color(icolor_chosen) + 1

      endif ! domain elements
    enddo

  enddo ! icolor

!! DK DK do it for inner elements

  icolormin = nb_colors_outer_elements + 1
  icolormax = ncolors

  ! ideal value if all colors are perfectly balanced
  if( nb_colors_inner_elements > 0 ) then
    target_nb_elems_per_color = nspec_inner / nb_colors_inner_elements + 1
  else
    target_nb_elems_per_color = 1
  endif
  ! print *,'nspec_inner,target_nb_elems_per_color = ',nspec_inner,target_nb_elems_per_color

  ! count the initial number of elements in each color
  nb_elems_in_this_color(:) = 0
  do icolor = icolormin,icolormax
    nb_elems_in_this_color(icolor) = count(color == icolor)
  enddo

  ! do not balance the last one, because it will be balanced automatically by the others
  do icolor = icolormin,icolormax-1

    ! if color is already balanced, do nothing
    ! (this works because in the initial set the number of elements per color decreases when the color number increases)
    if(nb_elems_in_this_color(icolor) <= target_nb_elems_per_color) cycle

    mask_ibool(:) = .false.
    last_ispec_studied = -1

    do ispec = 1,nspec
      ! domain elements only
      if(ispec_is_d(ispec)) then

        ! only inner elements
        if (is_on_a_slice_edge(ispec)) cycle

        ! only elements of this color
        if (color(ispec) /= icolor) cycle

        ! if color is now balanced because we have moved enough elements then stop searching
        if(nb_elems_in_this_color(icolor) <= target_nb_elems_per_color) exit

        ! unmark the eight corners of the previously marked element
        if(last_ispec_studied > 0) then
          mask_ibool(ibool(1,1,1,last_ispec_studied)) = .false.
          mask_ibool(ibool(NGLLX,1,1,last_ispec_studied)) = .false.
          mask_ibool(ibool(NGLLX,NGLLY,1,last_ispec_studied)) = .false.
          mask_ibool(ibool(1,NGLLY,1,last_ispec_studied)) = .false.
          mask_ibool(ibool(1,1,NGLLZ,last_ispec_studied)) = .false.
          mask_ibool(ibool(NGLLX,1,NGLLZ,last_ispec_studied)) = .false.
          mask_ibool(ibool(NGLLX,NGLLY,NGLLZ,last_ispec_studied)) = .false.
          mask_ibool(ibool(1,NGLLY,NGLLZ,last_ispec_studied)) = .false.
        endif
        icolor_conflict_found(icolormin:icolormax) = .false.
        icolor_conflict_found(icolor) = .true. ! cannot move element to the color it already has

        ! mark the eight corners of the current element
        mask_ibool(ibool(1,1,1,ispec)) = .true.
        mask_ibool(ibool(NGLLX,1,1,ispec)) = .true.
        mask_ibool(ibool(NGLLX,NGLLY,1,ispec)) = .true.
        mask_ibool(ibool(1,NGLLY,1,ispec)) = .true.
        mask_ibool(ibool(1,1,NGLLZ,ispec)) = .true.
        mask_ibool(ibool(NGLLX,1,NGLLZ,ispec)) = .true.
        mask_ibool(ibool(NGLLX,NGLLY,NGLLZ,ispec)) = .true.
        mask_ibool(ibool(1,NGLLY,NGLLZ,ispec)) = .true.
        last_ispec_studied = ispec

        ! test if we can move this element to another color
        do ispec2 = 1,nspec
          ! domain elements only
          if(ispec_is_d(ispec2)) then

            ! do not test that element itself
            if (ispec2 == ispec) cycle

            ! only inner elements
            if (is_on_a_slice_edge(ispec2)) cycle

            ! if conflict already found previously with this color, no need to test again
            if (icolor_conflict_found(color(ispec2))) cycle

            ! test the eight corners of the current element for a common point with element under study
            if (mask_ibool(ibool(1,1,1,ispec2)) .or. &
              mask_ibool(ibool(NGLLX,1,1,ispec2)) .or. &
              mask_ibool(ibool(NGLLX,NGLLY,1,ispec2)) .or. &
              mask_ibool(ibool(1,NGLLY,1,ispec2)) .or. &
              mask_ibool(ibool(1,1,NGLLZ,ispec2)) .or. &
              mask_ibool(ibool(NGLLX,1,NGLLZ,ispec2)) .or. &
              mask_ibool(ibool(NGLLX,NGLLY,NGLLZ,ispec2)) .or. &
              mask_ibool(ibool(1,NGLLY,NGLLZ,ispec2))) &
                icolor_conflict_found(color(ispec2)) = .true.

          endif ! domain elements
        enddo

        ! if color is already above target size for a balanced set, do not move to it
        do icolor_target = icolormin,icolormax
          if(nb_elems_in_this_color(icolor_target) >= target_nb_elems_per_color) &
            icolor_conflict_found(icolor_target) = .true.
        enddo

        ! if cannot find any other color to move this element to
        if (all(icolor_conflict_found(icolormin:icolormax))) cycle

        ! loops on all the colors to determine the color with the smallest number of elements
        ! and for which there is no conflict
        nb_elems_in_color_chosen = 2147000000 ! start with extremely large unrealistic value
        icolor_chosen = 0
        do icolor_target = icolormin,icolormax
          if (.not. icolor_conflict_found(icolor_target) .and. &
            nb_elems_in_this_color(icolor_target) < nb_elems_in_color_chosen) then
            icolor_chosen = icolor_target
            nb_elems_in_color_chosen = nb_elems_in_this_color(icolor_target)
          endif
        enddo

        ! moves the element to that new color
        ! remove element from its current color
        nb_elems_in_this_color(color(ispec)) = nb_elems_in_this_color(color(ispec)) - 1
        color(ispec) = icolor_chosen
        ! and add it to the new color
        nb_elems_in_this_color(icolor_chosen) = nb_elems_in_this_color(icolor_chosen) + 1

      endif ! domain elements
    enddo

  enddo ! icolor

  end subroutine balance_colors_simple

!
!-------------------------------------------------------------------------------------------------
!

  subroutine get_final_perm(color,perm,first_elem_number_in_this_color, &
                            nspec,nb_colors,nb_colors_outer_elements, &
                            ispec_is_d,nspec_domain)

  integer, intent(in) :: nspec,nb_colors

  integer,dimension(nspec), intent(in) :: color
  integer,dimension(nspec), intent(inout) :: perm

  integer, intent(inout) :: first_elem_number_in_this_color(nb_colors)

  logical,dimension(nspec),intent(in) :: ispec_is_d

  integer,intent(in) :: nb_colors_outer_elements,nspec_domain

  ! local parameters
  integer :: ispec,icolor,icounter,counter_outer

  ! note: permutations are only valid within each domain
  !          also, the counters start at 1 for each inner/outer element range

  ! outer elements first ( note: inner / outer order sensitive)
  icounter = 1
  do icolor = 1, nb_colors_outer_elements
    first_elem_number_in_this_color(icolor) = icounter
    do ispec = 1, nspec
      ! elements in this domain only
      if( ispec_is_d(ispec) ) then
        if(color(ispec) == icolor) then
          perm(ispec) = icounter
          icounter = icounter + 1
        endif
      endif
    enddo
  enddo
  counter_outer = icounter - 1

  ! inner elements second
  icounter = 1
  do icolor = nb_colors_outer_elements+1, nb_colors
    first_elem_number_in_this_color(icolor) = icounter + counter_outer
    do ispec = 1, nspec
      ! elements in this domain only
      if( ispec_is_d(ispec) ) then
        ! outer elements
        if(color(ispec) == icolor) then
          perm(ispec) = icounter
          icounter = icounter + 1
        endif
      endif
    enddo
  enddo

  ! checks
  if( counter_outer + icounter -1 /= nspec_domain ) then
    print*,'error: perm: ',nspec_domain,counter_outer,icounter,counter_outer+icounter-1
    stop 'error get_final_perm: counter incomplete'
  endif

  end subroutine get_final_perm


!
!-------------------------------------------------------------------------------------------------
!
! PERMUTATIONS
!
!-------------------------------------------------------------------------------------------------
!

! implement permutation of elements for arrays of real (CUSTOM_REAL) type

  subroutine permute_elements_real(array_to_permute,temp_array,perm,nspec)

  implicit none

  include "constants.h"

  integer, intent(in) :: nspec
  integer, intent(in), dimension(nspec) :: perm

  real(kind=CUSTOM_REAL), intent(inout), dimension(NGLLX,NGLLY,NGLLZ,nspec) :: &
    array_to_permute,temp_array

  integer old_ispec,new_ispec

  ! copy the original array
  temp_array(:,:,:,:) = array_to_permute(:,:,:,:)

  do old_ispec = 1,nspec
    new_ispec = perm(old_ispec)
    array_to_permute(:,:,:,new_ispec) = temp_array(:,:,:,old_ispec)
  enddo

  end subroutine permute_elements_real

!
!-------------------------------------------------------------------------------------------------
!

! implement permutation of elements for arrays of integer type

  subroutine permute_elements_integer(array_to_permute,temp_array,perm,nspec)

  implicit none

  include "constants.h"

  integer, intent(in) :: nspec
  integer, intent(in), dimension(nspec) :: perm

  integer, intent(inout), dimension(NGLLX,NGLLY,NGLLZ,nspec) :: &
    array_to_permute,temp_array

  integer old_ispec,new_ispec

  ! copy the original array
  temp_array(:,:,:,:) = array_to_permute(:,:,:,:)

  do old_ispec = 1,nspec
    new_ispec = perm(old_ispec)
    array_to_permute(:,:,:,new_ispec) = temp_array(:,:,:,old_ispec)
  enddo

  end subroutine permute_elements_integer

!
!-------------------------------------------------------------------------------------------------
!

! implement permutation of elements for arrays of double precision type

  subroutine permute_elements_dble(array_to_permute,temp_array,perm,nspec)

  implicit none

  include "constants.h"

  integer, intent(in) :: nspec
  integer, intent(in), dimension(nspec) :: perm

  double precision, intent(inout), dimension(NGLLX,NGLLY,NGLLZ,nspec) :: &
    array_to_permute,temp_array

  integer old_ispec,new_ispec

  ! copy the original array
  temp_array(:,:,:,:) = array_to_permute(:,:,:,:)

  do old_ispec = 1,nspec
    new_ispec = perm(old_ispec)
    array_to_permute(:,:,:,new_ispec) = temp_array(:,:,:,old_ispec)
  enddo

  end subroutine permute_elements_dble

!
!-------------------------------------------------------------------------------------------------
!

! implement permutation of elements for arrays of double precision type

  subroutine permute_elements_logical1D(array_to_permute,temp_array,perm,nspec)

  implicit none

  include "constants.h"

  integer, intent(in) :: nspec
  integer, intent(in), dimension(nspec) :: perm

  logical, intent(inout), dimension(nspec) :: array_to_permute,temp_array

  integer old_ispec,new_ispec

  ! copy the original array
  temp_array(:) = array_to_permute(:)

  do old_ispec = 1,nspec
    new_ispec = perm(old_ispec)
    array_to_permute(new_ispec) = temp_array(old_ispec)
  enddo

  end subroutine permute_elements_logical1D
