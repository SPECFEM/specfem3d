!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  1 . 4
!               ---------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!         (c) California Institute of Technology September 2006
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

module PML_par

  use constants,only: CUSTOM_REAL

  !--------------------------------------------------------------------
  ! USER PARAMETERS
  
  ! damping profile coefficients: 
  !   R: theoretical reflection coefficient after discretization
  real(kind=CUSTOM_REAL),parameter:: PML_damp_R = 1.e-3 

  ! number of element layers for PML region
  ! default is 2 element layers
  integer :: PML_LAYERS = 2

  ! additional absorbing, Sommerfeld (^Stacey) condition at the boundaries
  logical,parameter:: PML_USE_SOMMERFELD = .false.
  
  !--------------------------------------------------------------------

  real(kind=CUSTOM_REAL):: PML_width
  real(kind=CUSTOM_REAL):: PML_width_min,PML_width_max
  
  ! PML element type flag
  integer,dimension(:),allocatable :: ispec_is_PML_inum

  ! PML global points
  integer,dimension(:),allocatable :: iglob_is_PML

  ! PML spectral elements
  integer,dimension(:),allocatable :: PML_ispec
  integer :: num_PML_ispec
  
  ! PML normal for each PML spectral element
  real(kind=CUSTOM_REAL),dimension(:,:),allocatable :: PML_normal
  ! PML damping coefficients d & dprime
  real(kind=CUSTOM_REAL),dimension(:,:,:,:),allocatable :: PML_damping_d,PML_damping_dprime

  !real(kind=CUSTOM_REAL),dimension(:),allocatable :: PML_damping_d_global
  
  ! PML interface
  integer,dimension(:),allocatable :: iglob_is_PML_interface
  
  ! mask ibool needed for time marching
  logical,dimension(:),allocatable :: PML_mask_ibool  
  
  ! PML damping flag
  logical:: PML = .false.

end module PML_par

!--------

module PML_par_acoustic

  ! potentials split into 4 terms plus temporary potential:
  ! chi = chi1 + chi2 + chi3 + chi4
  ! temporary: chi2_t = (\partial_t + d ) chi2

  use constants,only: CUSTOM_REAL
  
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
                        chi1,chi2,chi2_t,chi3,chi4
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
                        chi1_dot,chi2_t_dot,chi3_dot,chi4_dot
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
                        chi1_dot_dot,chi2_t_dot_dot,chi3_dot_dot,chi4_dot_dot

end module PML_par_acoustic

!
!-------------------------------------------------------------------------------------------------
!

subroutine PML_damping_profile_l(d,x,vp,delta)

! calculates damping coefficient value d  for a given 
!   x: distance x and 
!   vp: p-velocity alpha
!   delta: PML width
!
! returns: d damping coefficients
  use PML_par,only: CUSTOM_REAL,PML_damp_R
  implicit none
  real(kind=CUSTOM_REAL),intent(out):: d
  real(kind=CUSTOM_REAL),intent(in):: x,vp,delta

  ! damping profile coefficients: 
  !   d : damping function of (x)
  !   vp:  P-velocity
  !   delta: width of PML layer 
  !   R: theoretical reflection coefficient after discretization
  
  ! damping profile function: d = f(x)
  ! Komatitsch & Tromp, 2003: eq. 24 page 150
  d = 3.0*vp/(2.0*delta)*log(1.0/PML_damp_R)*x*x/(delta*delta) 
  
end subroutine PML_damping_profile_l

!
!-------------------------------------------------------------------------------------------------
!

subroutine PML_initialize()

  use specfem_par,only: NGLOB_AB,NSPEC_AB,myrank, &
                        ibool,xstore,ystore,zstore,&
                        model_speed_max,hdur
  use PML_par
  use PML_par_acoustic
  use constants,only: FIX_UNDERFLOW_PROBLEM,VERYSMALLVAL,IMAIN,&
                      NGLLX,NGLLY,NGLLZ,TINYVAL
  use specfem_par_acoustic,only: ACOUSTIC_SIMULATION
  implicit none

  ! local parameters
  real(kind=CUSTOM_REAL):: d,dprime,d_glob,dprime_glob
  real(kind=CUSTOM_REAL) :: dominant_wavelength,hdur_max
  integer :: count,ilayer,sign

  ! sets flag
  PML = .true.

  ! user output
  if( myrank == 0 ) then
    write(IMAIN,*)
    write(IMAIN,*) 'incorporating PML  '
    write(IMAIN,*)
  endif  
  
  ! PML element type array: 1 = face, 2 = edge, 3 = corner
  allocate(ispec_is_PML_inum(NSPEC_AB))  
  num_PML_ispec = 0
  
  ! PML interface points between PML and "regular" region
  allocate(iglob_is_PML_interface(NGLOB_AB))
  iglob_is_PML_interface(:) = 0
  
  ! PML global points
  allocate(iglob_is_PML(NGLOB_AB))
  iglob_is_PML(:) = 0

  ! PML ibool mask
  allocate(PML_mask_ibool(NGLOB_AB))
  PML_mask_ibool(:) = .false.
  
  ! determines dominant wavelength based on maximum model speed 
  ! and source half time duration
  hdur_max = maxval(hdur(:))
  if( hdur_max > 0.0 ) then
    dominant_wavelength = model_speed_max * 2.0 * hdur_max
  else
    dominant_wavelength = 0._CUSTOM_REAL
  endif

  ! for multiple PML element layers
  ilayer = 0
  do while( ilayer < PML_LAYERS  )
    ilayer = ilayer + 1

    if( ilayer == 1 ) then
      ! sets ispec occurrences for first element layer in PML region based on absorbing boundary elements
      call PML_set_firstlayer()
    else
      ! adds an additional element layer based on adjacent elements on PML interface points
      call PML_add_layer()    
    endif
    
    ! update global interface points of PML region to "regular" domain
    call PML_determine_interfacePoints()  
    
    ! optional? update PML width according to dominant wavelength
    !call PML_get_width()
    ! checks with wavelength criteria
    !if( dominant_wavelength > 0.0 ) then    
    !  if( PML_width > dominant_wavelength/2.0 ) then
    !    PML_LAYERS = ilayer
    !    exit
    !  else
    !    PML_LAYERS = ilayer + 1
    !  endif
    !endif
  enddo
  
  ! checks PML normals at edges and corners, 
  ! tries to gather elements at edges & corners
  do ilayer=1,PML_LAYERS-1
    call PML_update_normals(ilayer)
  enddo

  ! updates statistics global PML width
  call PML_get_width()

  ! pre-calculates damping profiles on PML points
  ! damping coefficients
  call PML_set_local_dampingcoeff()

  ! pre-calculates derivatives of damping coefficients
  call PML_determine_dprime()  

  ! wavefield array initialization
  allocate(chi1(NGLLX,NGLLY,NGLLZ,num_PML_ispec),&
          chi2(NGLLX,NGLLY,NGLLZ,num_PML_ispec),&
          chi2_t(NGLLX,NGLLY,NGLLZ,num_PML_ispec),&
          chi3(NGLLX,NGLLY,NGLLZ,num_PML_ispec),&
          chi4(NGLLX,NGLLY,NGLLZ,num_PML_ispec))
  allocate(chi1_dot(NGLLX,NGLLY,NGLLZ,num_PML_ispec),&
          chi2_t_dot(NGLLX,NGLLY,NGLLZ,num_PML_ispec),&
          chi3_dot(NGLLX,NGLLY,NGLLZ,num_PML_ispec),&
          chi4_dot(NGLLX,NGLLY,NGLLZ,num_PML_ispec))
  allocate(chi1_dot_dot(NGLLX,NGLLY,NGLLZ,num_PML_ispec),&
          chi2_t_dot_dot(NGLLX,NGLLY,NGLLZ,num_PML_ispec),&
          chi3_dot_dot(NGLLX,NGLLY,NGLLZ,num_PML_ispec),&
          chi4_dot_dot(NGLLX,NGLLY,NGLLZ,num_PML_ispec))

  ! potentials
  chi1 = 0._CUSTOM_REAL
  chi2 = 0._CUSTOM_REAL
  chi2_t = 0._CUSTOM_REAL
  chi3 = 0._CUSTOM_REAL
  chi4 = 0._CUSTOM_REAL

  ! "velocity" potential
  chi1_dot = 0._CUSTOM_REAL
  chi2_t_dot = 0._CUSTOM_REAL
  chi3_dot = 0._CUSTOM_REAL
  chi4_dot = 0._CUSTOM_REAL

  ! "acceleration"/pressure potential
  chi1_dot_dot = 0._CUSTOM_REAL
  chi2_t_dot_dot = 0._CUSTOM_REAL
  chi3_dot_dot = 0._CUSTOM_REAL
  chi4_dot_dot = 0._CUSTOM_REAL    
  if(FIX_UNDERFLOW_PROBLEM) then 
    chi1_dot_dot = VERYSMALLVAL
    chi2_t_dot_dot = VERYSMALLVAL
    chi3_dot_dot = VERYSMALLVAL
    chi4_dot_dot = VERYSMALLVAL    
  endif

  ! statistics user output    
  d = maxval(abs(PML_damping_d(:,:,:,:)))
  if( d > TINYVAL ) then
    sign = maxval(PML_damping_d(:,:,:,:)) / maxval(abs(PML_damping_d(:,:,:,:)))
  else
    sign = 1.0
  endif
  dprime = maxval(abs(PML_damping_dprime(:,:,:,:)))
  call max_all_cr(d,d_glob)
  call max_all_cr(dprime,dprime_glob)
  call sum_all_i(num_PML_ispec,count)  
  if( myrank == 0 ) then
    write(IMAIN,*)
    write(IMAIN,*) 'PML region: '
    write(IMAIN,*) '    total spectral elements:',count
    write(IMAIN,*) '    number of layers : ',PML_LAYERS
    write(IMAIN,*) '    dominant wavelength max: ',dominant_wavelength
    write(IMAIN,*) '    width min / max:',PML_width_min,PML_width_max
    write(IMAIN,*) '    reflection coefficient:',PML_damp_R
    write(IMAIN,*) '    maximum d : ',sign*d_glob
    write(IMAIN,*) '    maximum dprime : ',sign*dprime_glob
    write(IMAIN,*)
  endif
  
  ! VTK file output
  call PML_output_VTKs()
    
end subroutine PML_initialize

!
!-------------------------------------------------------------------------------------------------
!

subroutine PML_set_firstlayer()

! sets ispec occurrences for first element layer in PML region based on absorbing boundary elements

  use PML_par
  use specfem_par,only: NSPEC_AB,NGLOB_AB, &
                        abs_boundary_ispec,abs_boundary_normal,num_abs_boundary_faces,&
                        abs_boundary_ijk,ibool,myrank
  use constants,only: NDIM,TINYVAL,NGNOD,NGLLX,NGLLY,NGLLZ,NGLLSQUARE
  implicit none
  ! local parameters
  real(kind=CUSTOM_REAL),dimension(:,:),allocatable:: temp_ispec_pml_normal
  integer,dimension(:),allocatable:: temp_is_pml_elem  
  integer:: iface,count,new_elemts,ispec,icorner,igll,iglobf
  integer:: i,j,k,iglobcount,iglobcorners(NGNOD)
  integer,dimension(3,NGNOD),parameter :: ielem_corner_ijk = &
       reshape((/ 1,1,1, 1,NGLLY,1, 1,NGLLY,NGLLZ, 1,1,NGLLZ, &
              NGLLX,1,1, NGLLX,NGLLY,1, NGLLX,NGLLY,NGLLZ, NGLLX,1,NGLLZ /),(/3,NGNOD/))
  ! temporary arrays
  allocate(temp_is_pml_elem(NSPEC_AB))
  allocate(temp_ispec_pml_normal(NDIM,NSPEC_AB))

  temp_is_pml_elem(:) = 0
  temp_ispec_pml_normal(:,:) = 0._CUSTOM_REAL

  count = 0
  do iface=1,num_abs_boundary_faces
    ! gets spectral elements with boundary face
    ispec = abs_boundary_ispec(iface)
           
    ! counts new PML elements
    if( temp_is_pml_elem(ispec) == 0 ) count = count + 1
    
    ! counts number of occurrences
    !  1 : element with 1 face to regular one,
    !  2 : element with 2 faces (elements at edges)
    !  3 : element with 3 faces (elements at corners)
    temp_is_pml_elem(ispec) = temp_is_pml_elem(ispec) + 1    
    
    ! adds contribution to element normal
    temp_ispec_pml_normal(:,ispec) = temp_ispec_pml_normal(:,ispec) + abs_boundary_normal(:,1,iface)
  enddo
  new_elemts = count

  ! doubling layers might have elements with only an edge on the absorbing surface
  ! poses problems if not accounted for
  count = 0
  do ispec = 1,NSPEC_AB
    ! only elements not recognized so far
    if( temp_is_pml_elem(ispec) > 0 ) cycle
    
    ! stores global indices of element corners
    do icorner=1,NGNOD
      i = ielem_corner_ijk(1,icorner)
      j = ielem_corner_ijk(2,icorner)
      k = ielem_corner_ijk(3,icorner)      
      iglobcorners(icorner) = ibool(i,j,k,ispec)      
    enddo
    
    ! checks if element has an edge (two corner points) on a absorbing boundary    
    ! (refers mainly to elements in doubling layers)
    do iface=1,num_abs_boundary_faces
      ! checks if already encountered this element
      if( abs_boundary_ispec(iface) == ispec ) exit
          
      ! loops over face points
      iglobcount = 0
      do igll=1,NGLLSQUARE
        i = abs_boundary_ijk(1,igll,iface)
        j = abs_boundary_ijk(2,igll,iface)
        k = abs_boundary_ijk(3,igll,iface)    
        iglobf = ibool(i,j,k,abs_boundary_ispec(iface))
        ! checks with corners
        do icorner=1,NGNOD
          if( iglobcorners(icorner) == iglobf ) iglobcount = iglobcount + 1
        enddo
      enddo
      
      ! adds as pml element
      if( iglobcount >= 2 ) then
        ! counter        
        if( temp_is_pml_elem(ispec) == 0 ) count = count + 1
        temp_is_pml_elem(ispec) = temp_is_pml_elem(ispec) + 1
        ! updates normal
        temp_ispec_pml_normal(:,ispec) = temp_ispec_pml_normal(:,ispec) &
                              + abs_boundary_normal(:,1,iface)
        exit
      endif
    enddo ! iface
    
  enddo
  new_elemts = new_elemts + count

  ! stores PML element indices and resulting normal
  call PML_set_elements(temp_is_pml_elem,temp_ispec_pml_normal,new_elemts)
  
  deallocate( temp_is_pml_elem)
  deallocate( temp_ispec_pml_normal)  
  
end subroutine PML_set_firstlayer

!
!-------------------------------------------------------------------------------------------------
!

subroutine PML_set_elements(temp_is_pml_elem,temp_ispec_pml_normal,new_elemts)

! adds new elements to PML region

  use PML_par
  use specfem_par,only: NSPEC_AB,myrank
  use constants,only: NDIM,TINYVAL
  implicit none
  
  integer:: temp_is_pml_elem(NSPEC_AB)
  real(kind=CUSTOM_REAL):: temp_ispec_pml_normal(NDIM,NSPEC_AB)
  integer:: new_elemts
  
  ! local parameters
  real(kind=CUSTOM_REAL) :: length
  integer :: ispec,ispecPML

  ! sets new element type flags
  ispec_is_PML_inum(:) = temp_is_pml_elem(:)

  ! sets new number of elements
  num_PML_ispec = new_elemts    

  ! re-allocates arrays
  if( allocated(PML_normal) ) deallocate(PML_normal)
  if( allocated(PML_ispec) ) deallocate(PML_ispec)
  allocate(PML_ispec(num_PML_ispec))
  allocate(PML_normal(NDIM,num_PML_ispec))
  
  ! stores PML elements flags and normals
  ispecPML = 0
  do ispec=1,NSPEC_AB
    if( ispec_is_PML_inum(ispec) > 0 ) then
      ! stores indices
      ispecPML = ispecPML+1
      PML_ispec(ispecPML) = ispec   
          
      ! gets resulting element normal
      PML_normal(:,ispecPML) = temp_ispec_pml_normal(:,ispec)

      ! normalizes normal
      length = sqrt( PML_normal(1,ispecPML)**2 &
                   + PML_normal(2,ispecPML)**2 &
                   + PML_normal(3,ispecPML)**2 )
      if( length < TINYVAL ) then
        print*,'error set elements: normal length:',length
        print*,'elem:',ispec,ispecPML
        print*,'num_pml_ispec:',num_PML_ispec
        call exit_mpi(myrank,'error PML normal length')
      else
        ! normalizes normal
        PML_normal(:,ispecPML) = PML_normal(:,ispecPML)/length
      endif      
    endif
  enddo
  if( ispecPML /= num_PML_ispec) call exit_mpi(myrank,'PML add layer count error')

end subroutine PML_set_elements

!
!-------------------------------------------------------------------------------------------------
!

subroutine PML_determine_interfacePoints()

! finds global interface points of PML region to "regular" domain

  use specfem_par,only: ibool,myrank,NGLOB_AB,NSPEC_AB, &
                        num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                        nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                        my_neighbours_ext_mesh,NPROC
  use PML_par
  use PML_par_acoustic
  use constants,only: NGLLX,NGLLY,NGLLZ
  use specfem_par_acoustic,only: ispec_is_acoustic,ACOUSTIC_SIMULATION
  implicit none

  ! local parameters
  integer,dimension(:),allocatable:: temp_regulardomain
  integer:: i,j,k,ispec,iglob

  ! PML interface points array
  iglob_is_PML_interface(:) = 0
  
  ! temporary arrays
  allocate(temp_regulardomain(NGLOB_AB))    
  temp_regulardomain(:) = 0
  
  ! global PML points
  iglob_is_PML(:) = 0
  
  ! sets flags on PML and regular domain points
  do ispec=1,NSPEC_AB
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX
          iglob = ibool(i,j,k,ispec)
          ! sets flag for PML/regular domain
          if( ispec_is_PML_inum(ispec) > 0 ) then
            ! global points
            iglob_is_PML(iglob) = iglob_is_PML(iglob) + 1                        
          else
            ! not a PML point
            temp_regulardomain(iglob) = temp_regulardomain(iglob) + 1
          endif
        enddo
      enddo
    enddo
  enddo
  
  ! assemble on MPI interfaces
  call assemble_MPI_scalar_i_ext_mesh(NPROC,NGLOB_AB,iglob_is_PML, &
                        num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                        nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh,&
                        my_neighbours_ext_mesh)  
  call assemble_MPI_scalar_i_ext_mesh(NPROC,NGLOB_AB,temp_regulardomain, &
                        num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                        nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh,&
                        my_neighbours_ext_mesh)  
                        
  ! stores interface points
  do ispec=1,NSPEC_AB
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX
          iglob = ibool(i,j,k,ispec)
          ! checks if it belongs to both, PML and regular domain
          if( temp_regulardomain(iglob) > 0 .and. iglob_is_PML(iglob) > 0 ) then
            ! increases flag on global point
            iglob_is_PML_interface(iglob) = iglob_is_PML_interface(iglob) + 1            
          endif
        enddo
      enddo
    enddo
  enddo

  deallocate(temp_regulardomain)

end subroutine PML_determine_interfacePoints

!
!-------------------------------------------------------------------------------------------------
!


subroutine PML_get_width()

! calculates PML width for statistics

  use specfem_par,only: abs_boundary_ispec,abs_boundary_normal,abs_boundary_ijk,&
                        num_abs_boundary_faces,&
                        ibool,xstore,ystore,zstore,myrank, &
                        NGLOB_AB
  use PML_par
  use constants,only: NGLLSQUARE,TINYVAL,HUGEVAL
  implicit none

  ! local parameters
  real(kind=CUSTOM_REAL) :: length,width
  integer:: i,j,k,ispec,iglob,iface,igll,iglobf

  ! determines global PML width
  ! loops over domain surface
  width = HUGEVAL
  do iface=1,num_abs_boundary_faces
  
    ispec = abs_boundary_ispec(iface)

    ! avoids taking corner or edge elements for width
    if( ispec_is_PML_inum(ispec) > 1 ) cycle
    
    ! determines smallest distance to interface points
    do iglob=1,NGLOB_AB
      if( iglob_is_PML_interface(iglob) > 0 ) then                    
        ! loops over face points
        do igll=1,NGLLSQUARE
          i = abs_boundary_ijk(1,igll,iface)
          j = abs_boundary_ijk(2,igll,iface)
          k = abs_boundary_ijk(3,igll,iface)
    
          ! takes distance between two points
          iglobf = ibool(i,j,k,ispec)
          length =  sqrt((xstore(iglobf)-xstore(iglob))**2 &
                       + (ystore(iglobf)-ystore(iglob))**2 &
                       + (zstore(iglobf)-zstore(iglob))**2 )
          
          ! checks length
          if( length < TINYVAL ) then
            print*,'PML:',myrank,'length:',length
            print*,'  ijk:',i,j,k,ispec,'face:',iface,'iglob:',iglobf
            print*,'  ijk xyz:',xstore(iglobf),ystore(iglobf),zstore(iglobf)
            print*,'  iglob interface',iglob
            print*,'  iglob xyz:',xstore(iglob),ystore(iglob),zstore(iglob)
            call exit_mpi(myrank,'PML length zero error')
          endif
                
          ! updates minimum width      
          if( length < width ) width = length
          
        enddo        
      endif      
    enddo
  enddo
  
  ! determines maximum width on all MPI processes
  ! all process gets overall maximum
  call max_all_all_cr(width,PML_width_max)
  call min_all_all_cr(width,PML_width_min)
  
  ! sets PML width
  if( PML_width_min > TINYVAL ) then
    PML_width = PML_width_min
  else
    PML_width = PML_width_max
  endif
    
  ! checks
  if( PML_width < TINYVAL ) call exit_mpi(myrank,'PML width error: width too small')

end subroutine PML_get_width

!
!-------------------------------------------------------------------------------------------------
!

subroutine PML_set_local_dampingcoeff()

! calculates damping profiles on PML points

  use specfem_par,only: ibool,xstore,ystore,zstore,myrank, &
                        kappastore,mustore,NGLOB_AB,&
                        abs_boundary_ispec,abs_boundary_ijk,num_abs_boundary_faces                        
  use specfem_par_acoustic,only: ACOUSTIC_SIMULATION,rhostore
  use specfem_par_elastic,only: ELASTIC_SIMULATION,rho_vp
  use PML_par
  use constants,only: NGLLX,NGLLY,NGLLZ,HUGEVAL,FOUR_THIRDS,NGLLSQUARE,TINYVAL
  implicit none

  ! local parameters
  real(kind=CUSTOM_REAL) :: length
  real(kind=CUSTOM_REAL) :: dist,vp
  real(kind=CUSTOM_REAL) :: d
  real(kind=CUSTOM_REAL) :: width

  integer:: i,j,k,ispec,iglob,ispecPML,iglobf
  integer:: ispecB,igll,iface
  
  ! stores damping coefficient
  allocate( PML_damping_d(NGLLX,NGLLY,NGLLZ,num_PML_ispec))    
  PML_damping_d(:,:,:,:) = 0._CUSTOM_REAL    
  
  ! loops over all PML elements             
  do ispecPML=1,num_PML_ispec
  
    ispec = PML_ispec(ispecPML)

    ! determines smallest distance to interface points
    ! and determines smallest distance to absorbing boundary points 
    ! (note: MPI partitioning not considered here yet; might be a problem)
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX
          ! global index
          iglobf = ibool(i,j,k,ispec)

          ! ensures that PML interface points have zero damping coefficients
          if( iglob_is_PML_interface(iglobf) > 0 ) then
            PML_damping_d(i,j,k,ispecPML) = 0._CUSTOM_REAL
            cycle
          endif
          
          ! distance to PML interface points
          dist = HUGEVAL
          do iglob=1,NGLOB_AB
            if( iglob_is_PML_interface(iglob) > 0 ) then                    
              ! distance to interface
              length =  (xstore(iglobf)-xstore(iglob))**2 &
                      + (ystore(iglobf)-ystore(iglob))**2 &
                      + (zstore(iglobf)-zstore(iglob))**2              
              if( length < dist ) dist = length 
            endif                    
          enddo !iglob
          !dist = distances(i,j,k) 
          if( dist >= HUGEVAL ) then
            dist = PML_width_max
          else
            dist = sqrt( dist ) 
          endif          

          ! distance to boundary points
          width = HUGEVAL
          do iface=1,num_abs_boundary_faces
            ispecB = abs_boundary_ispec(iface)      
            do igll=1,NGLLSQUARE
              iglob = ibool(abs_boundary_ijk(1,igll,iface),&
                             abs_boundary_ijk(2,igll,iface),&
                             abs_boundary_ijk(3,igll,iface),ispecB)
              ! distance to boundary
              length =  (xstore(iglobf)-xstore(iglob))**2 &
                      + (ystore(iglobf)-ystore(iglob))**2 &
                      + (zstore(iglobf)-zstore(iglob))**2 
              if( length < width ) width = length
            enddo
          enddo ! iface
          ! apparent width of PML for this point
          if( width >= HUGEVAL ) then
            width = PML_width_max
          else
            width = sqrt( width ) + dist
          endif          
          
          ! checks width 
          if( width < TINYVAL ) then
            print*,'error: pml width ',width
            print*,'ijk:',ispec,i,j,k
            print*,'xyz:',xstore(ibool(i,j,k,ispec)),ystore(ibool(i,j,k,ispec)),zstore(ibool(i,j,k,ispec))
            print*,'dist:',dist
            print*,'pml min/max:',PML_width_max,PML_width_min
            call exit_mpi(myrank,'PML error getting width')
          endif          
              
          ! P-velocity
          if( ACOUSTIC_SIMULATION ) then
            vp = sqrt( kappastore(i,j,k,ispec)/rhostore(i,j,k,ispec) )
          else if( ELASTIC_SIMULATION ) then
            vp = (FOUR_THIRDS * mustore(i,j,k,ispec) + kappastore(i,j,k,ispec)) &
                        / rho_vp(i,j,k,ispec)
          else
            call exit_mpi(myrank,'PML error getting vp')
          endif          
          
          ! gets damping coefficient
          call PML_damping_profile_l(d,dist,vp,width)
              
          ! stores d & dprime for this element's GLL points              
          PML_damping_d(i,j,k,ispecPML) = d          
          
        enddo
      enddo
    enddo
  enddo !ispecPML

end subroutine PML_set_local_dampingcoeff

!
!-------------------------------------------------------------------------------------------------
!

subroutine PML_determine_dprime()

! calculates derivatives dprime of damping coefficients on GLL points

  use PML_par
  use PML_par_acoustic
  use constants,only: NGLLX,NGLLY,NGLLZ
  use specfem_par,only: xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz,&
                        hprime_xx,hprime_yy,hprime_zz
  implicit none
  ! local parameters
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLZ,NGLLZ) :: dprime_elem
  real(kind=CUSTOM_REAL) :: xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl
  real(kind=CUSTOM_REAL) :: nx,ny,nz
  real(kind=CUSTOM_REAL) :: d_dx,d_dy,d_dz,tempd_dx,tempd_dy,tempd_dz
  integer :: ispec,i,j,k,l,ispecPML 

  ! dprime derivatives
  allocate( PML_damping_dprime(NGLLX,NGLLY,NGLLZ,num_PML_ispec))  
  PML_damping_dprime(:,:,:,:) = 0._CUSTOM_REAL  

  ! loops over all PML elements           
  do ispecPML=1,num_PML_ispec
  
    ispec = PML_ispec(ispecPML)
    
    ! PML normal 
    nx = PML_normal(1,ispecPML)
    ny = PML_normal(2,ispecPML)
    nz = PML_normal(3,ispecPML)

    ! calculates terms:
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX
          ! derivative along x, y, z
          ! first double loop over GLL points to compute and store gradients
          ! we can merge the loops because NGLLX == NGLLY == NGLLZ
          tempd_dx = 0._CUSTOM_REAL
          tempd_dy = 0._CUSTOM_REAL
          tempd_dz = 0._CUSTOM_REAL          
          do l = 1,NGLLX
            tempd_dx = tempd_dx + PML_damping_d(l,j,k,ispecPML)*hprime_xx(i,l)
            tempd_dy = tempd_dy + PML_damping_d(i,l,k,ispecPML)*hprime_yy(j,l)
            tempd_dz = tempd_dz + PML_damping_d(i,j,l,ispecPML)*hprime_zz(k,l)
          enddo 

          ! get derivatives of potential with respect to x, y and z
          xixl = xix(i,j,k,ispec)
          xiyl = xiy(i,j,k,ispec)
          xizl = xiz(i,j,k,ispec)
          etaxl = etax(i,j,k,ispec)
          etayl = etay(i,j,k,ispec)
          etazl = etaz(i,j,k,ispec)
          gammaxl = gammax(i,j,k,ispec)
          gammayl = gammay(i,j,k,ispec)
          gammazl = gammaz(i,j,k,ispec)
          
          ! derivatives dprime
          d_dx = xixl*tempd_dx + etaxl*tempd_dy + gammaxl*tempd_dz
          d_dy = xiyl*tempd_dx + etayl*tempd_dy + gammayl*tempd_dz
          d_dz = xizl*tempd_dx + etazl*tempd_dy + gammazl*tempd_dz
          dprime_elem(i,j,k) = d_dx*nx + d_dy*ny + d_dz*nz

        enddo
      enddo
    enddo

    ! stores dprime coefficients
    PML_damping_dprime(:,:,:,ispecPML) = dprime_elem(:,:,:)

  enddo

end subroutine PML_determine_dprime

!
!-------------------------------------------------------------------------------------------------
!

subroutine PML_add_layer()

! adds an element layer to the PML region

  use PML_par
  use specfem_par,only: NSPEC_AB,NGLOB_AB, &
                        abs_boundary_ispec,abs_boundary_normal,num_abs_boundary_faces,&
                        ibool,myrank,&
                        num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                        nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                        my_neighbours_ext_mesh,NPROC                        
  use constants,only: NDIM,TINYVAL,NGLLX,NGLLY,NGLLZ,NGNOD2D
  implicit none

  ! local parameters
  real(kind=CUSTOM_REAL),dimension(:,:),allocatable:: iglob_pml_normal
  real(kind=CUSTOM_REAL),dimension(:,:),allocatable:: ispec_pml_normal
  integer,dimension(:),allocatable:: is_pml_elem
  integer:: i,j,k,iglob,count,ispecPML,ispec,new_elemts
  integer :: iface,icorner,ipmlcorners
  
  integer,dimension(3,4),parameter :: iface1_corner_ijk = &
       reshape((/ 1,1,1, 1,NGLLY,1, 1,NGLLY,NGLLZ, 1,1,NGLLZ /),(/3,4/)) ! xmin  
  integer,dimension(3,4),parameter :: iface2_corner_ijk = &
       reshape((/ NGLLX,1,1, NGLLX,NGLLY,1, NGLLX,NGLLY,NGLLZ, NGLLX,1,NGLLZ  /),(/3,4/)) ! xmax  
  integer,dimension(3,4),parameter :: iface3_corner_ijk = &
       reshape((/ 1,1,1, 1,1,NGLLZ, NGLLX,1,NGLLZ, NGLLX,1,1  /),(/3,4/)) ! ymin  
  integer,dimension(3,4),parameter :: iface4_corner_ijk = &
       reshape((/ 1,NGLLY,1, NGLLX,NGLLY,1, NGLLX,NGLLY,NGLLZ, 1,NGLLY,NGLLZ /),(/3,4/)) ! ymax  
  integer,dimension(3,4),parameter :: iface5_corner_ijk = &
       reshape((/ 1,1,1, 1,NGLLY,1, NGLLX,NGLLY,1, NGLLX,1,1 /),(/3,4/)) ! bottom    
  integer,dimension(3,4),parameter :: iface6_corner_ijk = &
       reshape((/ 1,1,NGLLZ, NGLLX,1,NGLLZ, NGLLX,NGLLY,NGLLZ, 1,NGLLY,NGLLZ  /),(/3,4/)) ! top  
  integer,dimension(3,4,6),parameter :: iface_all_corner_ijk = &
       reshape((/ iface1_corner_ijk,iface2_corner_ijk, &
                  iface3_corner_ijk,iface4_corner_ijk, &
                  iface5_corner_ijk,iface6_corner_ijk /),(/3,4,6/)) ! all faces
  ! midpoint indices for each face (xmin,xmax,ymin,ymax,zmin,zmax)               
  integer,dimension(3,6),parameter :: iface_all_midpointijk = &
             reshape( (/ 1,2,2, NGLLX,2,2, 2,1,2, 2,NGLLY,2, 2,2,1, 2,2,NGLLZ  /),(/3,6/))
  logical :: is_done
  
  ! temporary arrays
  allocate(is_pml_elem(NSPEC_AB))
  allocate(iglob_pml_normal(NDIM,NGLOB_AB))
  allocate(ispec_pml_normal(NDIM,NSPEC_AB))
  
  iglob_pml_normal(:,:) = 0._CUSTOM_REAL
  ispec_pml_normal(:,:) = 0._CUSTOM_REAL

  ! sets pml normals on PML interface, global points  
  do ispecPML=1,num_PML_ispec

    ispec = PML_ispec(ispecPML)
    ! checks
    if( ispec_is_PML_inum(ispec) < 1 ) call exit_mpi(myrank,'PML error add ispec layer')
    
    ! starts from first layer elements 
    ! stores normal information on temporary global points
    if( ispec_is_PML_inum(ispec) >= 1 ) then          
      ! stores PML normal on interface points
      do k=1,NGLLZ
        do j=1,NGLLY
          do i=1,NGLLX
            iglob = ibool(i,j,k,ispec)            
            if( iglob_is_PML_interface(iglob) > 0 ) then     
              iglob_pml_normal(:,iglob) = iglob_pml_normal(:,iglob) + PML_normal(:,ispecPML)            
            endif  
          enddo
        enddo
      enddo
    endif
    
  enddo

  ! assembles with other MPI processes
  call assemble_MPI_vector_ext_mesh(NPROC,NGLOB_AB,iglob_pml_normal, &
                        num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                        nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                        my_neighbours_ext_mesh)


  ! adds new elements sharing PML interface 
  count = 0
  is_pml_elem(:) = 0
  do ispec=1,NSPEC_AB
  
    ! checks if we already have this element set as pml element in first layer
    is_done = .false.
    do ispecPML=1,num_PML_ispec
      if( PML_ispec(ispecPML) == ispec ) then
        ! adds as pml element
        if(is_pml_elem(ispec) == 0) count = count + 1        
        ! copies normal
        ispec_pml_normal(:,ispec) = PML_normal(:,ispecPML)
        ! copies element type flag
        is_pml_elem(ispec) = ispec_is_PML_inum(ispec)

        is_done = .true.
        exit
      endif
    enddo  
    if( is_done ) cycle
    
    ! loops over element faces
    do iface=1,6
      ipmlcorners = 0
      do icorner=1,NGNOD2D
        i = iface_all_corner_ijk(1,icorner,iface)
        j = iface_all_corner_ijk(2,icorner,iface)
        k = iface_all_corner_ijk(3,icorner,iface)
        iglob = ibool(i,j,k,ispec)
        if( iglob_is_PML_interface(iglob) > 0 ) ipmlcorners = ipmlcorners + 1
      enddo
    
      ! face is pml interface
      if( ipmlcorners == NGNOD2D ) then              
        ! counts new pml elements
        if(is_pml_elem(ispec) == 0) count = count + 1
        
        ! increments flag
        is_pml_elem(ispec) = is_pml_elem(ispec) + 1            
        
        ! sets normal    
        ! reference midpoint on face
        i = iface_all_midpointijk(1,iface)
        j = iface_all_midpointijk(2,iface)
        k = iface_all_midpointijk(3,iface)      
        iglob = ibool(i,j,k,ispec)
        if( iglob_is_PML_interface(iglob) < 1 ) call exit_mpi(myrank,'PML error midpoint interface')  
        
        ! checks new normal
        if( sqrt(iglob_pml_normal(1,iglob)**2+iglob_pml_normal(2,iglob)**2 &
                +iglob_pml_normal(3,iglob)**2) < TINYVAL ) then
          print*,'error add layer: normal length zero: iglob',iglob
          print*,'face ',iface,ipmlcorners
          print*,'ijk ispec',i,j,k,ispec
          call exit_mpi(myrank,'PML add layer has new normal length error')
        endif
        
        ! adds contribution to normal 
        ispec_pml_normal(:,ispec) = ispec_pml_normal(:,ispec) + iglob_pml_normal(:,iglob)        
      endif
      
    enddo ! iface    
  enddo ! ispec
  new_elemts = count
  
  ! adds new pml elements to PML region
  call PML_set_elements(is_pml_elem,ispec_pml_normal,new_elemts)

end subroutine PML_add_layer

!
!-------------------------------------------------------------------------------------------------
!

subroutine PML_update_normals(ilayer)

! updates normal's directions for elements in PML region

  use PML_par
  use specfem_par,only: NSPEC_AB,NGLOB_AB, &
                        ibool,myrank,&
                        num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                        nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                        my_neighbours_ext_mesh,NPROC                        
  use constants,only: NGNOD2D,NGLLX,NGLLY,NGLLZ
  implicit none
  integer :: ilayer

  ! local parameters
  integer::  iglob,ispecPML,ispec
  integer :: iface,icorner
  integer :: ipmlcorners,ipmledges,ipmlsngl
  integer :: ipmlcorners_tot,ipmledges_tot,ipmlsngl_tot
  
  integer,dimension(3,4),parameter :: iface1_corner_ijk = &
       reshape((/ 1,1,1, 1,NGLLY,1, 1,NGLLY,NGLLZ, 1,1,NGLLZ /),(/3,4/)) ! xmin  
  integer,dimension(3,4),parameter :: iface2_corner_ijk = &
       reshape((/ NGLLX,1,1, NGLLX,NGLLY,1, NGLLX,NGLLY,NGLLZ, NGLLX,1,NGLLZ  /),(/3,4/)) ! xmax  
  integer,dimension(3,4),parameter :: iface3_corner_ijk = &
       reshape((/ 1,1,1, 1,1,NGLLZ, NGLLX,1,NGLLZ, NGLLX,1,1  /),(/3,4/)) ! ymin  
  integer,dimension(3,4),parameter :: iface4_corner_ijk = &
       reshape((/ 1,NGLLY,1, NGLLX,NGLLY,1, NGLLX,NGLLY,NGLLZ, 1,NGLLY,NGLLZ /),(/3,4/)) ! ymax  
  integer,dimension(3,4),parameter :: iface5_corner_ijk = &
       reshape((/ 1,1,1, 1,NGLLY,1, NGLLX,NGLLY,1, NGLLX,1,1 /),(/3,4/)) ! bottom  
  integer,dimension(3,4),parameter :: iface6_corner_ijk = &
       reshape((/ 1,1,NGLLZ, NGLLX,1,NGLLZ, NGLLX,NGLLY,NGLLZ, 1,NGLLY,NGLLZ  /),(/3,4/)) ! top    
  integer,dimension(3,4,6),parameter :: iface_all_corner_ijk = &
       reshape((/ iface1_corner_ijk,iface2_corner_ijk, &
                  iface3_corner_ijk,iface4_corner_ijk, &
                  iface5_corner_ijk,iface6_corner_ijk /),(/3,4,6/)) ! all faces
  integer:: ispecngb,iadj,ipmlinterface
  integer :: ispecPMLngb_corner,ispecPMLngb_edge,ispecPMLngb_sngl
  integer,dimension(:),allocatable :: iglob_nadj,ispec_is_PML_inum_org
  integer,dimension(:,:,:),allocatable :: iglob_adj


  ! checks normals for elements adjacent to edge/corner elements
  ! assigns element information to each global point
  ! (note: mpi partitioning/interface between elements not considered yet)    
  allocate(iglob_nadj(NGLOB_AB),iglob_adj(2,32,NGLOB_AB))
  iglob_nadj(:) = 0
  iglob_adj(:,:,:) = 0
  do ispecPML=1,num_PML_ispec
    ispec = PML_ispec(ispecPML)    
    ! sets element corners
    do iface=1,2
      do icorner=1,NGNOD2D
        iglob = ibool(iface_all_corner_ijk(1,icorner,iface),&
                      iface_all_corner_ijk(2,icorner,iface),&
                      iface_all_corner_ijk(3,icorner,iface),ispec)
        ! number of occurrences
        iglob_nadj(iglob) = iglob_nadj(iglob) + 1
        ! first parameter is assigned element id ispec
        iglob_adj(1,iglob_nadj(iglob),iglob) = ispec
        ! second parameter is corresponding pml element id ispecPML
        iglob_adj(2,iglob_nadj(iglob),iglob) = ispecPML
      enddo
    enddo
  enddo
  if( maxval(iglob_nadj(:)) > 32 ) then
    print*,'info neighbors:',myrank
    print*,'max number of adjacents:',maxval(iglob_nadj(:)),maxloc(iglob_nadj(:))
    call exit_mpi(myrank,'error iglob number of adj')
  endif
    
  ! finds neighbors based on common nodes  and changes type and normal accordingly
  ! for edges and corners
  allocate(ispec_is_PML_inum_org(NSPEC_AB))
  ispec_is_PML_inum_org(:) = ispec_is_PML_inum(:)
  do ispecPML=1,num_PML_ispec
    ispec = PML_ispec(ispecPML)
    
    ! only non-corner elements
    if( ispec_is_PML_inum_org(ispec) <= 2 ) then
      ipmlsngl_tot = 0
      ipmlcorners_tot = 0
      ipmledges_tot = 0
      ipmlinterface = 0
      ! loops over element corners
      do iface=1,2
        ! checks corner neighbors
        do icorner=1,NGNOD2D
          iglob = ibool(iface_all_corner_ijk(1,icorner,iface),&
                       iface_all_corner_ijk(2,icorner,iface),&
                       iface_all_corner_ijk(3,icorner,iface),ispec)
          ! adjacent elements
          ipmlsngl = 0
          ipmlcorners = 0
          ipmledges = 0          
          do iadj=1,iglob_nadj(iglob)
            ispecngb = iglob_adj(1,iadj,iglob)
            if( ispecngb /= ispec ) then
              ! counts single normal neighbors
              if( ispec_is_PML_inum_org(ispecngb) == 1 ) then
                ipmlsngl = ipmlsngl + 1
                ispecPMLngb_sngl = iglob_adj(2,iadj,iglob)
              endif
              ! counts corner neighbors
              if( ispec_is_PML_inum_org(ispecngb) == 3 ) then
                ipmlcorners = ipmlcorners + 1
                ispecPMLngb_corner = iglob_adj(2,iadj,iglob)
              endif
              ! counts edge neighbors
              if( ispec_is_PML_inum_org(ispecngb) == 2 ) then
                ipmledges = ipmledges + 1
                ispecPMLngb_edge = iglob_adj(2,iadj,iglob)
              endif            
            endif
          enddo  
          if( ipmlsngl > 0 ) ipmlsngl_tot = ipmlsngl_tot + 1        
          if( ipmlcorners > 0 ) ipmlcorners_tot = ipmlcorners_tot + 1
          if( ipmledges > 0 ) ipmledges_tot = ipmledges_tot + 1
          
          ! interface points
          if( iglob_is_PML_interface(iglob) > 0 ) ipmlinterface = ipmlinterface + 1
          
        enddo !icorner        
      enddo
    
      ! elements inside PML
      if( ipmlinterface < 4 ) then
      
        ! shares two faces with edge elements, so it becomes an edge element too
        if( ispec_is_PML_inum_org(ispec) == 1 ) then
          if( ipmledges_tot >= 6 ) then
            ispec_is_PML_inum(ispec) = 2
            PML_normal(:,ispecPML) = PML_normal(:,ispecPMLngb_edge)
          endif
          if( ipmlcorners_tot >= 5 ) then
            ispec_is_PML_inum(ispec) = 3
            PML_normal(:,ispecPML) = PML_normal(:,ispecPMLngb_corner)
          endif        
        else if( ispec_is_PML_inum_org(ispec) == 2 ) then
        
        ! shares at least a face and a face edge with a corner element, 
        ! so it becomes a corner element too
          if( ipmlcorners_tot >= 5 ) then
            ispec_is_PML_inum(ispec) = 3
            PML_normal(:,ispecPML) = PML_normal(:,ispecPMLngb_corner)
          endif
        endif            
      endif
      ! avoid elements between two edges and next to corner to become edge elements
      if( ispec_is_PML_inum(ispec) == 2 .and. ilayer > 1 ) then
        if( ipmlsngl_tot == 8 .and. ipmlcorners_tot == 2 ) then 
          ispec_is_PML_inum(ispec) = 1
          PML_normal(:,ispecPML) = PML_normal(:,ispecPMLngb_sngl)
        endif
      endif
      
    endif  
  enddo
  deallocate(iglob_adj,iglob_nadj)
  deallocate(ispec_is_PML_inum_org)
  
end subroutine PML_update_normals

!
!-------------------------------------------------------------------------------------------------
!

subroutine PML_output_VTKs()

! outputs informations about PML elements 

  use PML_par
  use specfem_par,only: NGLOB_AB,NSPEC_AB,myrank, &
                        prname,ibool,xstore,ystore,zstore
  use constants,only: NGLLX,NGLLY,NGLLZ,IMAIN                        
  implicit none
  ! local parameters
  real(kind=CUSTOM_REAL),dimension(:,:),allocatable:: ispec_normal
  real(kind=CUSTOM_REAL),dimension(:,:,:,:),allocatable:: temp_gllvalues
  integer,dimension(:),allocatable :: temp_iglob
  integer :: count,iglob,ispecPML,ispec
  character(len=256) :: vtkfilename

  ! element type flags
  if( .false. ) then
    vtkfilename = prname(1:len_trim(prname))//'PML_ispec_inum'
    call write_VTK_data_elem_i(NSPEC_AB,NGLOB_AB,xstore,ystore,zstore,ibool,&
                          ispec_is_PML_inum,vtkfilename)
  endif  
  
  ! interface points
  if( .false. ) then
    ! puts global points in a temporary array for plotting
    count = 0
    do iglob=1,NGLOB_AB
      if( iglob_is_PML_interface(iglob) > 0 ) then
        count = count+1
      endif
    enddo      
    allocate(temp_iglob(count))
    count = 0
    do iglob=1,NGLOB_AB
      if( iglob_is_PML_interface(iglob) > 0 ) then
        count = count+1
        temp_iglob(count) = iglob
      endif
    enddo
    vtkfilename = prname(1:len_trim(prname))//'PML_interface_points'
    call write_VTK_data_points(NGLOB_AB,xstore,ystore,zstore, &
                          temp_iglob,count,vtkfilename)
    deallocate(temp_iglob)
  endif

  ! pml normals
  if( .false. ) then
    allocate(ispec_normal(3,NSPEC_AB) )
    ispec_normal(:,:) = 0._CUSTOM_REAL
    do ispecPML=1,num_PML_ispec
      ispec = PML_ispec(ispecPML)
      ispec_normal(:,ispec) = PML_normal(:,ispecPML)
    enddo
    vtkfilename = prname(1:len_trim(prname))//'PML_normals'
    call write_VTK_data_elem_vectors(NSPEC_AB,NGLOB_AB,xstore,ystore,zstore,ibool, &
                          ispec_normal,vtkfilename)  
    deallocate(ispec_normal)                        
  endif  

  ! pml damping coefficients
  if( .false. ) then
    allocate(temp_gllvalues(NGLLX,NGLLY,NGLLZ,NSPEC_AB))
    temp_gllvalues = 0._CUSTOM_REAL
    do ispecPML=1,num_PML_ispec
      ispec = PML_ispec(ispecPML)
      temp_gllvalues(:,:,:,ispec) = PML_damping_d(:,:,:,ispecPML)
    enddo
    vtkfilename = prname(1:len_trim(prname))//'PML_damping_d'
    call write_VTK_data_gll_cr(NSPEC_AB,NGLOB_AB, &
              xstore,ystore,zstore,ibool, &
              temp_gllvalues,vtkfilename)
    deallocate(temp_gllvalues)    
  endif ! VTK file output

  if(myrank == 0) write(IMAIN,*)

end subroutine PML_output_VTKs