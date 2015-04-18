
! define sets of colors that contain disconnected elements for the CUDA solver.
! also split the elements into two subsets: inner and outer elements, in order
! to be able to compute the outer elements first in the solver and then
! start non-blocking MPI calls and overlap them with the calculation of the inner elements
! (which works fine because there are always far more inner elements than outer elements)

  subroutine get_perm_color(is_on_a_slice_edge,ibool,perm,nspec,nglob, &
     nb_colors_outer_elements,nb_colors_inner_elements,nspec_outer,first_elem_number_in_this_color,myrank)

  implicit none

  include "constants.h"

  logical, dimension(nspec) :: is_on_a_slice_edge

  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool
  integer, dimension(nspec) :: perm
  integer, dimension(nspec) :: color
  integer, dimension(MAX_NUMBER_OF_COLORS) :: first_elem_number_in_this_color
  integer :: nb_colors_outer_elements,nb_colors_inner_elements,nspec_outer,myrank

! local variables
  integer nspec,nglob_GLL_full

! a neighbor of a hexahedral node is a hexahedron that shares a face with it -> max degree of a node = 6
  integer, parameter :: MAX_NUMBER_OF_NEIGHBORS = 100

! global corner numbers that need to be created
  integer, dimension(nglob) :: global_corner_number

  integer mn(nspec*NGNOD_HEXAHEDRA),mp(nspec+1)
  integer, dimension(:), allocatable :: ne,np,adj
  integer xadj(nspec+1)

  logical maskel(nspec)

  integer i,istart,istop,number_of_neighbors

  integer nglob_eight_corners_only,nglob

! only count the total size of the array that will be created, or actually create it
  logical count_only
  integer total_size_ne,total_size_adj

!
!-----------------------------------------------------------------------
!

! total number of points in the mesh
    nglob_GLL_full = nglob

!---- call Charbel Farhat's routines
    if(myrank == 0) &
      write(IMAIN,*) 'calling form_elt_connectivity_foelco to perform mesh coloring and inner/outer element splitting'
    call form_elt_connectivity_foelco(mn,mp,nspec,global_corner_number,nglob_GLL_full,ibool,nglob_eight_corners_only)
    do i=1,nspec
      istart = mp(i)
      istop = mp(i+1) - 1
    enddo

! count only, to determine the size needed for the array
    allocate(np(nglob_eight_corners_only+1))
    count_only = .true.
    total_size_ne = 1
    if(myrank == 0) write(IMAIN,*) 'calling form_node_connectivity_fonoco to determine the size of the table'
    allocate(ne(total_size_ne))
    call form_node_connectivity_fonoco(mn,mp,ne,np,nglob_eight_corners_only,nspec,count_only,total_size_ne)
    deallocate(ne)

! allocate the array with the right size
    allocate(ne(total_size_ne))

! now actually generate the array
    count_only = .false.
    if(myrank == 0) write(IMAIN,*) 'calling form_node_connectivity_fonoco to actually create the table'
    call form_node_connectivity_fonoco(mn,mp,ne,np,nglob_eight_corners_only,nspec,count_only,total_size_ne)
    do i=1,nglob_eight_corners_only
      istart = np(i)
      istop = np(i+1) - 1
    enddo

! count only, to determine the size needed for the array
    count_only = .true.
    total_size_adj = 1
    if(myrank == 0) write(IMAIN,*) 'calling create_adjacency_table_adjncy to determine the size of the table'
    allocate(adj(total_size_adj))
    call create_adjacency_table_adjncy(mn,mp,ne,np,adj,xadj,maskel,nspec,nglob_eight_corners_only,&
    count_only,total_size_ne,total_size_adj,.false.)
    deallocate(adj)

! allocate the array with the right size
    allocate(adj(total_size_adj))

! now actually generate the array
    count_only = .false.
    if(myrank == 0) write(IMAIN,*) 'calling create_adjacency_table_adjncy again to actually create the table'
    call create_adjacency_table_adjncy(mn,mp,ne,np,adj,xadj,maskel,nspec,nglob_eight_corners_only,&
    count_only,total_size_ne,total_size_adj,.false.)

    do i=1,nspec
      istart = xadj(i)
      istop = xadj(i+1) - 1
      number_of_neighbors = istop-istart+1
      if(number_of_neighbors < 1 .or. number_of_neighbors > MAX_NUMBER_OF_NEIGHBORS) stop 'incorrect number of neighbors'
    enddo

    deallocate(ne,np)

    call get_color(adj,xadj,color,nspec,total_size_adj,is_on_a_slice_edge, &
       nb_colors_outer_elements,nb_colors_inner_elements,nspec_outer)

    if(myrank == 0) then
      write(IMAIN,*) 'number of colors of the graph for inner elements = ',nb_colors_inner_elements
      write(IMAIN,*) 'number of colors of the graph for outer elements = ',nb_colors_outer_elements
      write(IMAIN,*) 'total number of colors of the graph (sum of both) = ', &
       nb_colors_inner_elements + nb_colors_outer_elements
      write(IMAIN,*) 'number of elements of the graph for outer elements = ',nspec_outer
    endif

    deallocate(adj)

    if(myrank == 0) write(IMAIN,*) 'generating the final colors'
    first_elem_number_in_this_color(:) = -1
    call get_final_perm(color,perm,first_elem_number_in_this_color,nspec,nb_colors_inner_elements+nb_colors_outer_elements)

    if(myrank == 0) write(IMAIN,*) 'done with mesh coloring and inner/outer element splitting'

  end subroutine get_perm_color

!------------------------------------------------------------------

subroutine get_final_perm(color,perm,first_elem_number_in_this_color,nspec,nb_color)

  integer, intent(in) :: nspec,nb_color
  integer, intent(in) :: color(nspec)
  integer, intent(inout) :: perm(nspec)
  integer, intent(inout) :: first_elem_number_in_this_color(nb_color)
  integer :: ielem,icolor,counter

  counter = 1
  do icolor = 1, nb_color
    first_elem_number_in_this_color(icolor) = counter
    do ielem = 1, nspec
      if(color(ielem) == icolor) then
        perm(ielem) = counter
        counter = counter + 1
      endif
    enddo
  enddo

end subroutine get_final_perm

!------------------------------------------------------------------

subroutine get_color(adj,xadj,color,nspec,total_size_adj,is_on_a_slice_edge, &
     nb_colors_outer_elements,nb_colors_inner_elements,nspec_outer)

  integer, intent(in) :: nspec,total_size_adj
  integer, intent(in) :: adj(total_size_adj),xadj(nspec+1)
  integer :: color(nspec)
  integer :: this_color,nb_already_colored,ispec,ixadj,ok
  logical, dimension(nspec) :: is_on_a_slice_edge
  integer :: nb_colors_outer_elements,nb_colors_inner_elements,nspec_outer
  logical :: is_outer_element(nspec)

  nspec_outer = 0

  is_outer_element(:) = .false.

  do ispec=1,nspec
    if (is_on_a_slice_edge(ispec)) then
      is_outer_element(ispec) = .true.
      nspec_outer=nspec_outer+1
    endif
  enddo

! outer elements
  color(:) = 0
  this_color = 0
  nb_already_colored = 0
  do while(nb_already_colored<nspec_outer)
    this_color = this_color + 1
    do ispec = 1, nspec
      if (is_outer_element(ispec)) then
        if (color(ispec) == 0) then
          ok = 1
          do ixadj = xadj(ispec), (xadj(ispec+1)-1)
            if (is_outer_element(adj(ixadj)) .and. color(adj(ixadj)) == this_color) ok = 0
          enddo
          if (ok /= 0) then
            color(ispec) = this_color
            nb_already_colored = nb_already_colored + 1
          endif
        endif
      endif
    enddo
  enddo
  nb_colors_outer_elements = this_color

! inner elements
  do while(nb_already_colored<nspec)
    this_color = this_color + 1
    do ispec = 1, nspec
      if (.not. is_outer_element(ispec)) then
        if (color(ispec) == 0) then
          ok = 1
          do ixadj = xadj(ispec), (xadj(ispec+1)-1)
            if (.not. is_outer_element(adj(ixadj)) .and. color(adj(ixadj)) == this_color) ok = 0
          enddo
          if (ok /= 0) then
            color(ispec) = this_color
            nb_already_colored = nb_already_colored + 1
          endif
        endif
      endif
    enddo
  enddo

  nb_colors_inner_elements = this_color - nb_colors_outer_elements

end subroutine get_color

!------------------------------------------------------------------

!=======================================================================
!
!  Charbel Farhat's FEM topology routines
!
!  Dimitri Komatitsch, February 1996 - Code based on Farhat's original version from 1987
!
!  modified and adapted by Dimitri Komatitsch, May 2006
!
!=======================================================================

  subroutine form_elt_connectivity_foelco(mn,mp,nspec,global_corner_number,&
nglob_GLL_full,ibool,nglob_eight_corners_only)

!-----------------------------------------------------------------------
!
!   Forms the MN and MP arrays
!
!     Input :
!     -------
!           ibool    Array needed to build the element connectivity table
!           nspec    Number of elements in the domain
!           NGNOD_HEXAHEDRA    number of nodes per hexahedron (brick with 8 corners)
!
!     Output :
!     --------
!           MN, MP   This is the element connectivity array pair.
!                    Array MN contains the list of the element
!                    connectivity, that is, the nodes contained in each
!                    element. They are stored in a stacked fashion.
!
!                    Pointer array MP stores the location of each
!                    element list. Its length is equal to the number
!                    of elements plus one.
!
!-----------------------------------------------------------------------

  implicit none

  include "constants.h"

  integer nspec,nglob_GLL_full

! arrays with mesh parameters per slice
  integer, intent(in), dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool

! global corner numbers that need to be created
  integer, intent(out), dimension(nglob_GLL_full) :: global_corner_number
  integer, intent(out) :: mn(nspec*NGNOD_HEXAHEDRA),mp(nspec+1)
  integer, intent(out) :: nglob_eight_corners_only

  integer ninter,nsum,ispec,node,k,inumcorner,ix,iy,iz

  ninter = 1
  nsum = 1
  mp(1) = 1

!---- define topology of the elements in the mesh
!---- we need to define adjacent numbers from the sub-mesh consisting of the corners only
  nglob_eight_corners_only = 0
  global_corner_number(:) = -1

  do ispec=1,nspec

    inumcorner = 0
    do iz = 1,NGLLZ,NGLLZ-1
      do iy = 1,NGLLY,NGLLY-1
        do ix = 1,NGLLX,NGLLX-1

          inumcorner = inumcorner + 1
          if(inumcorner > NGNOD_HEXAHEDRA) stop 'corner number too large'

! check if this point was already assigned a number previously, otherwise create one and store it
          if(global_corner_number(ibool(ix,iy,iz,ispec)) == -1) then
            nglob_eight_corners_only = nglob_eight_corners_only + 1
            global_corner_number(ibool(ix,iy,iz,ispec)) = nglob_eight_corners_only
          endif

          node = global_corner_number(ibool(ix,iy,iz,ispec))
            do k=nsum,ninter-1
              if(node == mn(k)) goto 200
            enddo

            mn(ninter) = node
            ninter = ninter + 1
  200 continue

        enddo
      enddo
    enddo

      nsum = ninter
      mp(ispec + 1) = nsum

  enddo

  end subroutine form_elt_connectivity_foelco

!
!----------------------------------------------------
!

  subroutine form_node_connectivity_fonoco(mn,mp,ne,np,nglob_eight_corners_only,&
nspec,count_only,total_size_ne)

!-----------------------------------------------------------------------
!
!   Forms the NE and NP arrays
!
!     Input :
!     -------
!           MN, MP, nspec
!           nglob_eight_corners_only    Number of nodes in the domain
!
!     Output :
!     --------
!           NE, NP   This is the node-connected element array pair.
!                    Integer array NE contains a list of the
!                    elements connected to each node, stored in stacked fashion.
!
!                    Array NP is the pointer array for the
!                    location of a node's element list in the NE array.
!                    Its length is equal to the number of points plus one.
!
!-----------------------------------------------------------------------

  implicit none

  include "constants.h"

! only count the total size of the array that will be created, or actually create it
  logical count_only
  integer total_size_ne

  integer nglob_eight_corners_only,nspec

  integer, intent(in) ::  mn(nspec*NGNOD_HEXAHEDRA),mp(nspec+1)

  integer, intent(out) ::  ne(total_size_ne),np(nglob_eight_corners_only+1)

  integer nsum,inode,ispec,j

  nsum = 1
  np(1) = 1

  do inode=1,nglob_eight_corners_only
      do 200 ispec=1,nspec

            do j=mp(ispec),mp(ispec + 1) - 1
                  if (mn(j) == inode) then
                        if(count_only) then
                          total_size_ne = nsum
                        else
                          ne(nsum) = ispec
                        endif
                        nsum = nsum + 1
                        goto 200
                  endif
            enddo
  200 continue

      np(inode + 1) = nsum

  enddo

  end subroutine form_node_connectivity_fonoco

!
!----------------------------------------------------
!

  subroutine create_adjacency_table_adjncy(mn,mp,ne,np,adj,xadj,maskel,nspec,nglob_eight_corners_only,&
                                           count_only,total_size_ne,total_size_adj,face)

!-----------------------------------------------------------------------
!
!   Establishes the element adjacency information of the mesh
!   Two elements are considered adjacent if they share a face.
!
!     Input :
!     -------
!           MN, MP, NE, NP, nspec
!           MASKEL    logical mask (length = nspec)
!
!     Output :
!     --------
!           ADJ, XADJ This is the element adjacency array pair. Array
!                     ADJ contains the list of the elements adjacent to
!                     element i. They are stored in a stacked fashion.
!                     Pointer array XADJ stores the location of each element list.
!
!-----------------------------------------------------------------------

  implicit none

  include "constants.h"

! only count the total size of the array that will be created, or actually create it
  logical count_only,face
  integer total_size_ne,total_size_adj

  integer nglob_eight_corners_only

  integer, intent(in) :: mn(nspec*NGNOD_HEXAHEDRA),mp(nspec+1),ne(total_size_ne),np(nglob_eight_corners_only+1)

  integer, intent(out) :: adj(total_size_adj),xadj(nspec+1)

  logical maskel(nspec)
  integer countel(nspec)

  integer nspec,iad,ispec,istart,istop,ino,node,jstart,jstop,nelem,jel

  xadj(1) = 1
  iad = 1

  do ispec=1,nspec

! reset mask
  maskel(:) = .false.

! mask current element
  maskel(ispec) = .true.
  if (face) countel(:) = 0

  istart = mp(ispec)
  istop = mp(ispec+1) - 1
    do ino=istart,istop
      node = mn(ino)
      jstart = np(node)
      jstop = np(node + 1) - 1
        do 120 jel=jstart,jstop
            nelem = ne(jel)
            if(maskel(nelem)) goto 120
            if (face) then
              ! if 2 elements share at least 3 corners, therefore they share a face
              countel(nelem) = countel(nelem) + 1
              if (countel(nelem)>=3) then
                if(count_only) then
                  total_size_adj = iad
                else
                  adj(iad) = nelem
                endif
                maskel(nelem) = .true.
                iad = iad + 1
              endif
            else
              if(count_only) then
                total_size_adj = iad
              else
                adj(iad) = nelem
              endif
              maskel(nelem) = .true.
              iad = iad + 1
            endif
  120   continue
    enddo

    xadj(ispec+1) = iad

  enddo

  end subroutine create_adjacency_table_adjncy


  subroutine permute_elements_real(array_to_permute,temp_array,perm,nspec)

  implicit none

  include "constants.h"

  integer, intent(in) :: nspec
  integer, intent(in), dimension(nspec) :: perm

  real(kind=CUSTOM_REAL), intent(inout), dimension(NGLLX,NGLLY,NGLLZ,nspec) :: array_to_permute,temp_array

  integer old_ispec,new_ispec

! copy the original array
  temp_array(:,:,:,:) = array_to_permute(:,:,:,:)

  do old_ispec = 1,nspec
    new_ispec = perm(old_ispec)
    array_to_permute(:,:,:,new_ispec) = temp_array(:,:,:,old_ispec)
  enddo

  end subroutine permute_elements_real

!
!-----------------------------------------------------------------------
!

! implement permutation of elements for arrays of integer type
  subroutine permute_elements_integer(array_to_permute,temp_array,perm,nspec)

  implicit none

  include "constants.h"

  integer, intent(in) :: nspec
  integer, intent(in), dimension(nspec) :: perm

  integer, intent(inout), dimension(NGLLX,NGLLY,NGLLZ,nspec) :: array_to_permute,temp_array

  integer old_ispec,new_ispec

! copy the original array
  temp_array(:,:,:,:) = array_to_permute(:,:,:,:)

  do old_ispec = 1,nspec
    new_ispec = perm(old_ispec)
    array_to_permute(:,:,:,new_ispec) = temp_array(:,:,:,old_ispec)
  enddo

  end subroutine permute_elements_integer

!
!-----------------------------------------------------------------------
!

! implement permutation of elements for arrays of double precision type
  subroutine permute_elements_dble(array_to_permute,temp_array,perm,nspec)

  implicit none

  include "constants.h"

  integer, intent(in) :: nspec
  integer, intent(in), dimension(nspec) :: perm

  double precision, intent(inout), dimension(NGLLX,NGLLY,NGLLZ,nspec) :: array_to_permute,temp_array

  integer old_ispec,new_ispec

! copy the original array
  temp_array(:,:,:,:) = array_to_permute(:,:,:,:)

  do old_ispec = 1,nspec
    new_ispec = perm(old_ispec)
    array_to_permute(:,:,:,new_ispec) = temp_array(:,:,:,old_ispec)
  enddo

  end subroutine permute_elements_dble

