!
!----
!

subroutine get_element_face_id(ispec,xcoord,ycoord,zcoord,&
                              ibool,nspec,nglob, &
                              xstore_dummy,ystore_dummy,zstore_dummy, &
                              iface_id )

! returns iface_id of face in reference element, determined by corner locations xcoord/ycoord/zcoord;

  implicit none
  
  include "constants.h"
                     
  integer :: ispec,nspec,nglob,iface_id
  
! face corner locations
  real(kind=CUSTOM_REAL),dimension(NGNOD2D) :: xcoord,ycoord,zcoord

! index array
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool
  
! global point locations          
  real(kind=CUSTOM_REAL) :: xstore_dummy(nglob),ystore_dummy(nglob),zstore_dummy(nglob)
  
! local parameters  
  real(kind=CUSTOM_REAL),dimension(NGNOD2D) :: xcoord_face,ycoord_face,zcoord_face
  real(kind=CUSTOM_REAL) :: midpoint_faces(NDIM,6),midpoint(NDIM),midpoint_distances(6)
  
! corners indices of reference cube faces
  ! shapes of arrays below
  integer,dimension(2),parameter :: face_shape = (/3,4/)
  integer,dimension(3),parameter :: all_faces_shape = (/3,4,6/)

  ! xmin
  integer,dimension(3,4),parameter :: iface1_corner_ijk = &
       reshape((/ 1,1,1, 1,NGLLY,1, 1,NGLLY,NGLLZ, 1,1,NGLLZ /),face_shape)
  ! xmax
  integer,dimension(3,4),parameter :: iface2_corner_ijk = &
       reshape((/ NGLLX,1,1, NGLLX,NGLLY,1, NGLLX,NGLLY,NGLLZ, NGLLX,1,NGLLZ  /),face_shape)
  ! ymin
  integer,dimension(3,4),parameter :: iface3_corner_ijk = &
       reshape((/ 1,1,1, 1,1,NGLLZ, NGLLX,1,NGLLZ, NGLLX,1,1  /),face_shape)
  ! ymax
  integer,dimension(3,4),parameter :: iface4_corner_ijk = &
       reshape((/ 1,NGLLY,1, NGLLX,NGLLY,1, NGLLX,NGLLY,NGLLZ, 1,NGLLY,NGLLZ /),face_shape)
  ! bottom
  integer,dimension(3,4),parameter :: iface5_corner_ijk = &
       reshape((/ 1,1,1, 1,NGLLY,1, NGLLX,NGLLY,1, NGLLX,1,1 /),face_shape)
  ! top  
  integer,dimension(3,4),parameter :: iface6_corner_ijk = &
       reshape((/ 1,1,NGLLZ, NGLLX,1,NGLLZ, NGLLX,NGLLY,NGLLZ, 1,NGLLY,NGLLZ  /),face_shape)
  ! all faces
  integer,dimension(3,4,6),parameter :: iface_all_corner_ijk = &
       reshape((/ iface1_corner_ijk,iface2_corner_ijk, &
                  iface3_corner_ijk,iface4_corner_ijk, &
                  iface5_corner_ijk,iface6_corner_ijk /),all_faces_shape)
                 
! face orientation
  !real(kind=CUSTOM_REAL) :: face_n(3),face_ntmp(3),tmp
  integer  :: ifa,icorner,i,j,k,iglob,iloc(1)

! initializes
  iface_id = -1
  
! gets face midpoint by its corners 
  midpoint(:) = 0.0
  do icorner=1,NGNOD2D
    midpoint(1) = midpoint(1) + xcoord(icorner)
    midpoint(2) = midpoint(2) + ycoord(icorner)
    midpoint(3) = midpoint(3) + zcoord(icorner)      
  enddo
  midpoint(:) = midpoint(:) / 4.0

  ! checks: this holds only for planar face
  !if( midpoint(1) /= (xcoord(1)+xcoord(3))/2.0 .or. midpoint(1) /= (xcoord(2)+xcoord(4))/2.0  ) then
  !  print*,'error midpoint x:',midpoint(1),(xcoord(1)+xcoord(3))/2.0,(xcoord(2)+xcoord(4))/2.0
  !endif
  !if( midpoint(2) /= (ycoord(1)+ycoord(3))/2.0 .or. midpoint(2) /= (ycoord(2)+ycoord(4))/2.0  ) then
  !  print*,'error midpoint y:',midpoint(1),(ycoord(1)+ycoord(3))/2.0,(ycoord(2)+ycoord(4))/2.0
  !endif
  !if( midpoint(3) /= (zcoord(1)+zcoord(3))/2.0 .or. midpoint(3) /= (zcoord(2)+zcoord(4))/2.0  ) then
  !  print*,'error midpoint z:',midpoint(1),(zcoord(1)+zcoord(3))/2.0,(zcoord(2)+zcoord(4))/2.0
  !endif
     
! determines element face by minimum distance of midpoints
  midpoint_faces(:,:) = 0.0
  do ifa=1,6
    ! face corners
    do icorner = 1,NGNOD2D
      i = iface_all_corner_ijk(1,icorner,ifa)
      j = iface_all_corner_ijk(2,icorner,ifa)
      k = iface_all_corner_ijk(3,icorner,ifa)
      !print*,'corner:',i,j,k,ispec
      
      ! coordinates
      iglob = ibool(i,j,k,ispec)
      xcoord_face(icorner) = xstore_dummy(iglob)
      ycoord_face(icorner) = ystore_dummy(iglob)
      zcoord_face(icorner) = zstore_dummy(iglob)
      
      ! face midpoint coordinates
      midpoint_faces(1,ifa) =  midpoint_faces(1,ifa) + xcoord_face(icorner)
      midpoint_faces(2,ifa) =  midpoint_faces(2,ifa) + ycoord_face(icorner)
      midpoint_faces(3,ifa) =  midpoint_faces(3,ifa) + zcoord_face(icorner)
    enddo
    midpoint_faces(:,ifa) = midpoint_faces(:,ifa) / 4.0
    
    ! distance
    midpoint_distances(ifa) = (midpoint(1)-midpoint_faces(1,ifa))**2 &
                            + (midpoint(2)-midpoint_faces(2,ifa))**2 &
                            + (midpoint(3)-midpoint_faces(3,ifa))**2 
  enddo 

! gets closest point, which determines face
  iloc = minloc(midpoint_distances)

  ! checks that found midpoint is close enough  
  !print*,'face:', midpoint_distances(iloc(1))
  if( midpoint_distances(iloc(1)) > 1.e-5 * &
          (   (xcoord(1)-xcoord(2))**2 &
            + (ycoord(1)-ycoord(2))**2 &
            + (zcoord(1)-zcoord(2))**2 ) ) then
    print*,'error element face midpoint distance:',midpoint_distances(iloc(1)),(xcoord(1)-xcoord(2))**2
    ! corner locations 
    do icorner=1,NGNOD2D      
      i = iface_all_corner_ijk(1,icorner,iloc(1))
      j = iface_all_corner_ijk(2,icorner,iloc(1))
      k = iface_all_corner_ijk(3,icorner,iloc(1))
      iglob = ibool(i,j,k,ispec)    
      print*,'error corner:',icorner,'xyz:',xstore_dummy(iglob),ystore_dummy(iglob),zstore_dummy(iglob)
    enddo
    ! stop
    stop 'error element face midpoint'
  else
    iface_id = iloc(1)

    !print*,'face:',iface_id
    !do icorner=1,NGNOD2D      
    !  i = iface_all_corner_ijk(1,icorner,iloc(1))
    !  j = iface_all_corner_ijk(2,icorner,iloc(1))
    !  k = iface_all_corner_ijk(3,icorner,iloc(1))
    !  iglob = ibool(i,j,k,ispec)    
    !  print*,'corner:',icorner,'xyz:',sngl(xstore_dummy(iglob)), &
    !            sngl(ystore_dummy(iglob)),sngl(zstore_dummy(iglob))
    !enddo

  endif

end subroutine get_element_face_id

!
!----
!

subroutine get_element_face_gll_indices(iface,ijk_face,NGLLA,NGLLB )

! returns local indices in ijk_face for specified face

  implicit none
  
  include "constants.h"
                     
  integer :: iface !,nspec,nglob
  
! gll point indices i,j,k for face, format corresponds to ijk_face(1,*) = i, ijk_face(2,*) = j, ijk_face(3,*) = k
  integer :: NGLLA,NGLLB
  integer,dimension(3,NGLLA,NGLLB) :: ijk_face
  
!  integer  :: icorner,i,j,k,iglob,iloc(1)
  integer :: i,j,k
  integer :: ngll,i_gll,j_gll,k_gll
 
! sets i,j,k indices of GLL points on boundary face
  ngll = 0
  select case( iface )
  
  ! reference xmin face
  case(1)
    if( NGLLA /= NGLLY .or. NGLLB /= NGLLZ ) stop 'error absorbing face 1 indexing'
    i_gll = 1
    do k=1,NGLLZ
      do j=1,NGLLY
        ngll = ngll + 1
        ijk_face(1,j,k) = i_gll
        ijk_face(2,j,k) = j
        ijk_face(3,j,k) = k          
      enddo
    enddo
    
  ! reference xmax face
  case(2)
    if( NGLLA /= NGLLY .or. NGLLB /= NGLLZ ) stop 'error absorbing face 2 indexing'
    i_gll = NGLLX
    do k=1,NGLLZ
      do j=1,NGLLY
        ngll = ngll + 1
        ijk_face(1,j,k) = i_gll
        ijk_face(2,j,k) = j
        ijk_face(3,j,k) = k          
      enddo
    enddo

  ! reference ymin face
  case(3)
    if( NGLLA /= NGLLX .or. NGLLB /= NGLLZ ) stop 'error absorbing face 3 indexing'
    j_gll = 1
    do k=1,NGLLZ
      do i=1,NGLLX
        ngll = ngll + 1
        ijk_face(1,i,k) = i
        ijk_face(2,i,k) = j_gll
        ijk_face(3,i,k) = k          
      enddo
    enddo
    
  ! reference ymax face
  case(4)
    if( NGLLA /= NGLLX .or. NGLLB /= NGLLZ ) stop 'error absorbing face 4 indexing'  
    j_gll = NGLLY
    do k=1,NGLLZ
      do i=1,NGLLX
        ngll = ngll + 1
        ijk_face(1,i,k) = i
        ijk_face(2,i,k) = j_gll
        ijk_face(3,i,k) = k          
      enddo
    enddo
    
  ! reference bottom face
  case(5)
    if( NGLLA /= NGLLX .or. NGLLB /= NGLLY ) stop 'error absorbing face 5 indexing'  
    k_gll = 1
    do j=1,NGLLY
      do i=1,NGLLX
        ngll = ngll + 1
        ijk_face(1,i,j) = i
        ijk_face(2,i,j) = j
        ijk_face(3,i,j) = k_gll 
      enddo
    enddo
    
  ! reference bottom face
  case(6)
    if( NGLLA /= NGLLX .or. NGLLB /= NGLLY ) stop 'error absorbing face 6 indexing'  
    k_gll = NGLLZ
    do j=1,NGLLY
      do i=1,NGLLX
        ngll = ngll + 1
        ijk_face(1,i,j) = i
        ijk_face(2,i,j) = j
        ijk_face(3,i,j) = k_gll
      enddo
    enddo    
    
  case default
    stop 'error element face not found'
    
  end select

  ! checks number of gll points set on face
  if( ngll /= NGLLA*NGLLB ) then
    print*,'error element face ngll:',ngll,NGLLA,NGLLB
    stop 'error element face ngll'
  endif
!
!! corner locations 
!  do icorner=1,NGNOD2D      
!    i = iface_all_corner_ijk(1,icorner,iface)
!    j = iface_all_corner_ijk(2,icorner,iface)
!    k = iface_all_corner_ijk(3,icorner,iface)
!    iglob = ibool(i,j,k,ispec)    
!    xcoord_iboun(icorner) = xstore_dummy(iglob)
!    ycoord_iboun(icorner) = ystore_dummy(iglob) 
!    zcoord_iboun(icorner) = zstore_dummy(iglob)       
!    ! looks at values
!    !print*,'corner:',icorner,'xyz:',sngl(xcoord_iboun(icorner)),sngl(ycoord_iboun(icorner)),sngl(zcoord_iboun(icorner))      
!  enddo
!
!! determines initial orientation given by three corners of the face 
!  ! (CUBIT orders corners such that normal points outwards of element)
!  ! cross-product of vectors from corner 1 to corner 2 and from corner 1 to corner 3
!  face_n(1) =   (ycoord(2)-ycoord(1))*(zcoord(3)-zcoord(1)) - (zcoord(2)-zcoord(1))*(ycoord(3)-ycoord(1))
!  face_n(2) = - (xcoord(2)-xcoord(1))*(zcoord(3)-zcoord(1)) + (zcoord(2)-zcoord(1))*(xcoord(3)-xcoord(1))
!  face_n(3) =   (xcoord(2)-xcoord(1))*(ycoord(3)-ycoord(1)) - (ycoord(2)-ycoord(1))*(xcoord(3)-xcoord(1))
!  face_n(:) = face_n(:)/(sqrt( face_n(1)**2 + face_n(2)**2 + face_n(3)**2) )
!
!! checks that this normal direction is outwards of element: 
!  ! takes additional corner out of face plane and determines scalarproduct to normal
!  select case( iface )
!  case(1) ! opposite to xmin face
!    iglob = ibool(NGLLX,1,1,ispec)      
!  case(2) ! opposite to xmax face
!    iglob = ibool(1,1,1,ispec)      
!  case(3) ! opposite to ymin face
!    iglob = ibool(1,NGLLY,1,ispec)      
!  case(4) ! opposite to ymax face
!    iglob = ibool(1,1,1,ispec)        
!  case(5) ! opposite to bottom
!    iglob = ibool(1,1,NGLLZ,ispec)      
!  case(6) ! opposite to top
!    iglob = ibool(1,1,1,ispec)      
!  end select
!  ! vector from corner 1 to this opposite one
!  xcoord(4) = xstore_dummy(iglob) - xcoord(1)
!  ycoord(4) = ystore_dummy(iglob) - ycoord(1)
!  zcoord(4) = zstore_dummy(iglob) - zcoord(1)
!  
!  ! scalarproduct
!  tmp = xcoord(4)*face_n(1) + ycoord(4)*face_n(2) + zcoord(4)*face_n(3)
!  
!  ! makes sure normal points outwards, that is points away from this additional corner and scalarproduct is negative
!  if( tmp > 0.0 ) then
!    face_n(:) = - face_n(:)
!  endif  
!  !print*,'face ',iface,'scalarproduct:',tmp
!  
!! determines orientation of gll corner locations and sets it such that normal points outwards
!  ! cross-product
!  face_ntmp(1) =   (ycoord_iboun(2)-ycoord_iboun(1))*(zcoord_iboun(3)-zcoord_iboun(1)) &
!                     - (zcoord_iboun(2)-zcoord_iboun(1))*(ycoord_iboun(3)-ycoord_iboun(1))
!  face_ntmp(2) = - (xcoord_iboun(2)-xcoord_iboun(1))*(zcoord_iboun(3)-zcoord_iboun(1)) &
!                      + (zcoord_iboun(2)-zcoord_iboun(1))*(xcoord_iboun(3)-xcoord_iboun(1))
!  face_ntmp(3) =   (xcoord_iboun(2)-xcoord_iboun(1))*(ycoord_iboun(3)-ycoord_iboun(1))&
!                       - (ycoord_iboun(2)-ycoord_iboun(1))*(xcoord_iboun(3)-xcoord_iboun(1))
!  face_ntmp(:) = face_ntmp(:)/(sqrt( face_ntmp(1)**2 + face_ntmp(2)**2 + face_ntmp(3)**2) )
!  if( abs( (face_n(1)-face_ntmp(1))**2+(face_n(2)-face_ntmp(2))**2+(face_n(3)-face_ntmp(3))**2) > 0.1 ) then
!    !print*,'error orientation face 1:',ispec,face_n(:)
!    !swap corners 2 and 4 ( switches clockwise / anti-clockwise )
!    tmp = xcoord_iboun(2)
!    xcoord_iboun(2) = xcoord_iboun(4)
!    xcoord_iboun(4) = tmp
!    tmp = ycoord_iboun(2)
!    ycoord_iboun(2) = ycoord_iboun(4)
!    ycoord_iboun(4) = tmp
!    tmp = zcoord_iboun(2)
!    zcoord_iboun(2) = zcoord_iboun(4)
!    zcoord_iboun(4) = tmp      
!  endif

end subroutine get_element_face_gll_indices                  

!
!----
!

subroutine get_element_face_normal(ispec,iface,xcoord,ycoord,zcoord, &
                                ibool,nspec,nglob, &
                                xstore_dummy,ystore_dummy,zstore_dummy, &
                                normal)

! only changes direction of normal to point outwards of element

  implicit none
  
  include "constants.h"
                     
  integer :: ispec,iface,nspec,nglob
  
! face corner locations
  real(kind=CUSTOM_REAL),dimension(NGNOD2D) :: xcoord,ycoord,zcoord

! index array
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool
  
! global point locations          
  real(kind=CUSTOM_REAL),dimension(nglob) :: xstore_dummy,ystore_dummy,zstore_dummy
  
! face normal  
  real(kind=CUSTOM_REAL),dimension(NDIM) :: normal
  
! local parameters  
  real(kind=CUSTOM_REAL) :: face_n(3),tmp,v_tmp(3)
  integer :: iglob
 
! determines initial orientation given by three corners on the face 
  ! cross-product of vectors from corner 1 to corner 2 and from corner 1 to corner 3
  face_n(1) =   (ycoord(2)-ycoord(1))*(zcoord(3)-zcoord(1)) - (zcoord(2)-zcoord(1))*(ycoord(3)-ycoord(1))
  face_n(2) = - (xcoord(2)-xcoord(1))*(zcoord(3)-zcoord(1)) + (zcoord(2)-zcoord(1))*(xcoord(3)-xcoord(1))
  face_n(3) =   (xcoord(2)-xcoord(1))*(ycoord(3)-ycoord(1)) - (ycoord(2)-ycoord(1))*(xcoord(3)-xcoord(1))
  tmp = sqrt( face_n(1)*face_n(1) + face_n(2)*face_n(2) + face_n(3)*face_n(3) ) 
  if( abs(tmp) < TINYVAL ) then
    print*,'error get face normal: length',tmp
    print*,'normal:',face_n(:)
    call exit_mpi(0,'error get element face normal')
  endif
  face_n(:) = face_n(:)/tmp

! checks that this normal direction is outwards of element: 
  ! takes additional corner out of face plane and determines scalarproduct to normal
  select case( iface )
  case(1) ! opposite to xmin face
    iglob = ibool(NGLLX,1,1,ispec)      
  case(2) ! opposite to xmax face
    iglob = ibool(1,1,1,ispec)      
  case(3) ! opposite to ymin face
    iglob = ibool(1,NGLLY,1,ispec)      
  case(4) ! opposite to ymax face
    iglob = ibool(1,1,1,ispec)        
  case(5) ! opposite to bottom
    iglob = ibool(1,1,NGLLZ,ispec)      
  case(6) ! opposite to top
    iglob = ibool(1,1,1,ispec)      
  end select
  ! vector from corner 1 to this opposite one
  v_tmp(1) = xstore_dummy(iglob) - xcoord(1)
  v_tmp(2) = ystore_dummy(iglob) - ycoord(1)
  v_tmp(3) = zstore_dummy(iglob) - zcoord(1)
  
  ! scalarproduct
  tmp = v_tmp(1)*face_n(1) + v_tmp(2)*face_n(2) + v_tmp(3)*face_n(3)
  
  ! makes sure normal points outwards, that is points away from this additional corner and scalarproduct is negative
  if( tmp > 0.0 ) then
    face_n(:) = - face_n(:)
  endif  
 
! in case given normal has zero length, sets it to computed face normal
  if( ( normal(1)**2 + normal(2)**2 + normal(3)**2 ) < TINYVAL ) then
    normal(:) = face_n(:)
    return
  endif
   
! otherwise determines orientation of normal and flips direction such that normal points outwards
  tmp = face_n(1)*normal(1) + face_n(2)*normal(2) + face_n(3)*normal(3)
  if( tmp < 0.0 ) then
    !print*,'element face normal: orientation ',ispec,iface,tmp
    !print*,'face normal: ',face_n(:)
    !print*,'     normal: ',normal(:)
    !swap 
    normal(:) = - normal(:)      
  endif
  !print*,'face ',iface,'scalarproduct:',tmp

end subroutine get_element_face_normal      

!
!----
!

subroutine get_element_face_normal_idirect(ispec,iface,xcoord,ycoord,zcoord, &
                                ibool,nspec,nglob, &
                                xstore_dummy,ystore_dummy,zstore_dummy, &
                                normal,idirect)

! returns direction of normal: 
!   idirect = 1 to point outwards of/away from element
!   idirect = 2 to point into element

  implicit none
  
  include "constants.h"
                     
  integer :: ispec,iface,nspec,nglob
  
! face corner locations
  real(kind=CUSTOM_REAL),dimension(NGNOD2D) :: xcoord,ycoord,zcoord

! index array
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool
  
! global point locations          
  real(kind=CUSTOM_REAL) :: xstore_dummy(nglob),ystore_dummy(nglob),zstore_dummy(nglob)
  
! face normal  
  real(kind=CUSTOM_REAL),dimension(NDIM) :: normal
  
! direction type
  integer, intent(out) :: idirect
  
! local parameters  
  real(kind=CUSTOM_REAL) :: face_n(3),tmp,v_tmp(3)
  integer :: iglob
 
! initializes 
  idirect = 0
 
! determines initial orientation given by three corners on the face 
  ! cross-product of vectors from corner 1 to corner 2 and from corner 1 to corner 3
  face_n(1) =   (ycoord(2)-ycoord(1))*(zcoord(3)-zcoord(1)) - (zcoord(2)-zcoord(1))*(ycoord(3)-ycoord(1))
  face_n(2) = - (xcoord(2)-xcoord(1))*(zcoord(3)-zcoord(1)) + (zcoord(2)-zcoord(1))*(xcoord(3)-xcoord(1))
  face_n(3) =   (xcoord(2)-xcoord(1))*(ycoord(3)-ycoord(1)) - (ycoord(2)-ycoord(1))*(xcoord(3)-xcoord(1))
  tmp = sqrt( face_n(1)**2 + face_n(2)**2 + face_n(3)**2 ) 
  if( abs(tmp) < TINYVAL ) then
    print*,'error get face normal: length',tmp
    print*,'normal:',face_n(:)
    call exit_mpi(0,'error get element face normal')
  endif
  face_n(:) = face_n(:)/tmp

! checks that this normal direction is outwards of element: 
  ! takes additional corner out of face plane and determines scalarproduct to normal
  select case( iface )
  case(1) ! opposite to xmin face
    iglob = ibool(NGLLX,1,1,ispec)      
  case(2) ! opposite to xmax face
    iglob = ibool(1,1,1,ispec)      
  case(3) ! opposite to ymin face
    iglob = ibool(1,NGLLY,1,ispec)      
  case(4) ! opposite to ymax face
    iglob = ibool(1,1,1,ispec)        
  case(5) ! opposite to bottom
    iglob = ibool(1,1,NGLLZ,ispec)      
  case(6) ! opposite to top
    iglob = ibool(1,1,1,ispec)      
  end select
  ! vector from corner 1 to this opposite one
  v_tmp(1) = xstore_dummy(iglob) - xcoord(1)
  v_tmp(2) = ystore_dummy(iglob) - ycoord(1)
  v_tmp(3) = zstore_dummy(iglob) - zcoord(1)
  
  ! scalarproduct
  tmp = v_tmp(1)*face_n(1) + v_tmp(2)*face_n(2) + v_tmp(3)*face_n(3)
  
  ! makes sure normal points outwards, that is points away from this additional corner and scalarproduct is negative
  if( tmp > 0.0 ) then
    face_n(:) = - face_n(:)
  endif  
 
! in case given normal has zero length, exit
  if( ( normal(1)**2 + normal(2)**2 + normal(3)**2 ) < TINYVAL ) then    
    print*,'problem: given normal is zero'
    return
  endif
   
! otherwise determines orientation of normal 
  tmp = face_n(1)*normal(1) + face_n(2)*normal(2) + face_n(3)*normal(3)
  if( tmp < 0.0 ) then
    ! points into element
    idirect = 2
  else
    ! points away from element/ outwards
    idirect = 1
  endif

end subroutine get_element_face_normal_idirect
   
