  module  post_processing 

    contains 
      
      subroutine connect_points()
        use global_parameters
        use interp_mod

        real(kind=CUSTOM_REAL) xi,eta,scur,zcur
        real(kind=CUSTOM_REAL) smin,smax,zmin,zmax
        real(kind=CUSTOM_REAL) nodes_crd(NGNOD,2),z_axi(4)
        real(kind=CUSTOM_REAL), parameter :: eps=1.e-3
        integer  IGRIDs(8),IGRIDz(8),i_axi(4)
        integer irec,iel,inode
        integer nmax

        ! conversion GLL point to 8-node control elements
        IGRIDs(1)=0
        IGRIDs(2)=2
        IGRIDs(3)=4
        IGRIDs(4)=4
        IGRIDs(5)=4
        IGRIDs(6)=2
        IGRIDs(7)=0
        IGRIDs(8)=0
        
        IGRIDz(1)=0
        IGRIDz(2)=0
        IGRIDz(3)=0
        IGRIDz(4)=2
        IGRIDz(5)=4
        IGRIDz(6)=4
        IGRIDz(7)=4
        IGRIDz(8)=2
        
        allocate(xi_rec(nbrec),eta_rec(nbrec))
        allocate(rec2elm(nbrec))
        allocate(ele_candidate(4,nbrec))
        allocate(ele_seen(nbrec))
        ele_seen=0
        ele_candidate=-1
        rec2elm=-1
        ! CONNECTION POINT <-> MESH-----------------------------------------------------------
        do irec=1,nbrec
           scur=reciever_cyl(1,irec)
           zcur=reciever_cyl(3,irec)
           do iel = 1, NEL
              
              ! element
              smin=1d40
              smax=-1d40
              zmin=smin
              zmax=smax
              do inode=1,NGNOD
                 nodes_crd(inode,1)=scoor(IGRIDs(inode),IGRIDz(inode),iel)
                 nodes_crd(inode,2)=zcoor(IGRIDs(inode),IGRIDz(inode),iel)
                 !write(*,*) nodes_crd(inode,1), nodes_crd(inode,2)
                 smin=min(smin, nodes_crd(inode,1))
                 smax=max(smax, nodes_crd(inode,1))
                 zmin=min(zmin, nodes_crd(inode,2))
                 zmax=max(zmax, nodes_crd(inode,2))
              end do
              if ( scur > smin-eps .and. scur < smax + eps .and. zcur > zmin-eps .and. zcur < zmax + eps) then
                 call find_xix_eta(nodes_crd,xi,eta,scur,zcur)
                 if (xi > -1.05 .and. xi < 1.05 .and. eta > -1.05 .and. eta < 1.05) then
                    ele_seen(irec) =  ele_seen(irec) + 1
                    ele_candidate(ele_seen(irec),irec) = iel 
                    !!rec2elm(irec)=iel
                    !!xi_rec(irec)=xi
                    !!eta_rec(irec)=eta
                    !exit
                 end if
              end if
           end do
        end do

        !---------------- choice of candidate points --------------
        do irec=1,nbrec

           scur=reciever_cyl(1,irec)
           zcur=reciever_cyl(3,irec)

           z_axi=1.e30
           i_axi=-1
           nmax=ele_seen(irec)
           if (nmax == 0) then 
              write(*,*) 'pb with ', irec,reciever_cyl(1,irec),reciever_cyl(3,irec)
             
           end if
           do iel=1,ele_seen(irec)
              z_axi(iel) = depth_ele(ele_candidate(iel,irec))
              i_axi(iel) = ele_candidate(iel,irec)
           end do

           call sort_real_array(z_axi,i_axi)

           if (up(irec)) then ! then max value
              rec2elm(irec)=i_axi(nmax)
           else ! min value
              rec2elm(irec)=i_axi(1)
           end if

           iel = rec2elm(irec)
           do inode=1,NGNOD
              nodes_crd(inode,1)=scoor(IGRIDs(inode),IGRIDz(inode),iel)
              nodes_crd(inode,2)=zcoor(IGRIDs(inode),IGRIDz(inode),iel)
              !write(*,*) nodes_crd(inode,1), nodes_crd(inode,2)
              smin=min(smin, nodes_crd(inode,1))
              smax=max(smax, nodes_crd(inode,1))
              zmin=min(zmin, nodes_crd(inode,2))
              zmax=max(zmax, nodes_crd(inode,2))
           end do

           call find_xix_eta(nodes_crd,xi,eta,scur,zcur)
           xi_rec(irec)=xi
           eta_rec(irec)=eta

        end do
        call check_rec2elm
        ! END CONNECTION POINT <-> MESH ------------------------------------------------------

      end subroutine connect_points

      subroutine check_rec2elm()
        use global_parameters
        integer irec,forgot_point
        forgot_point=0
        do irec=1,nbrec
           if(rec2elm(irec)==-1) then
              forgot_point= forgot_point+1
           end if
        end do
        if (forgot_point > 0) write(*,*) 'forgot ', forgot_point,' points'
      end subroutine check_rec2elm

      !---------------------------------------------------------------------------------------

      !----------------------------------------------------------------------------------------
      subroutine sort_real_array(r,i)
        use global_parameters
        !implicit none
        real(kind=CUSTOM_REAL) r(4),tmp
        integer i(4)
        integer ii,jj,itmp
        
        do ii = 3, 1,-1
           do jj = 1 , ii
              if (r(jj) > r(jj+1)) then
                 tmp = r(jj)
                 r(jj)=r(jj+1)
                 r(jj+1) = tmp
                 
                 itmp = i(jj)
                 i(jj)=i(jj+1)
                 i(jj+1)=itmp
                 
              end if
           end do
        end do
        
      end subroutine sort_real_array
      
      subroutine interpol_field(ifield)
        
        use global_parameters
        use interp_mod
        integer irec,iel,ifield
        real(kind=CUSTOM_REAL) xi,eta,interp_value
        real(kind=CUSTOM_REAL) field(NGLLX,NGLLY)

        data_rec(:,ifield) = 0.
        !data_rec_temp(:)=0.
        do irec=irecmin , irecmax

           xi=xi_rec(irec)
           eta=eta_rec(irec)
           iel=rec2elm(irec)
           field(:,:)=data_read(:,:,iel)

           call interpole_field(xi,eta,field,interp_value) 

           data_rec(irec,ifield)=interp_value
           !data_rec_temp(irec)=interp_value
           !write(*,*) field
           !write(*,*) interp_value
           !write(*,*)
           !phi=reciever_cyl(3,irec)
           !if (ifield == 3 .and. irec == 10) then
           !   write(*,*) iel,xi,eta
           !   write(*,*) data_read(:,:,iel)
           !   write(*,*) interp_value
           !   write(*,*)  
           !end if
           ! prefactors depending on source type
           !call compute_prefactor(f1,f2,phi,src_type(isim,1),src_type(isim,2))
           
           ! expand 2D wave field in 3D axisymetric field
           

        end do
        !data_reduce(:)=0.
        !call mpi_reduce(data_rec_temp,data_reduce,)
        !data_rec(:,ifield)=data_deduce(:)

      end subroutine interpol_field

      subroutine interpol_stress(ifield)
        
        use global_parameters
        use interp_mod
        integer irec,iel,ifield
        real(kind=CUSTOM_REAL) xi,eta,interp_value
        real(kind=CUSTOM_REAL) field(NGLLX,NGLLY)
        
        stress_rec(:,ifield)=0.
        do irec = irecmin, irecmax

           xi=xi_rec(irec)
           eta=eta_rec(irec)
           iel=rec2elm(irec)
           field(:,:)=data_read(:,:,iel)

           call interpole_field(xi,eta,field,interp_value) 

           stress_rec(irec,ifield)=interp_value

           !write(*,*) interp_value

           !phi=reciever_cyl(3,irec)
           
           ! prefactors depending on source type
           !call compute_prefactor(f1,f2,phi,src_type(isim,1),src_type(isim,2))
           
           ! expand 2D wave field in 3D axisymetric field
           

        end do
        !* call mpi_reduce()

      end subroutine interpol_stress

      subroutine compute_prefactor(src_type,Mcomp)
        use global_parameters , only : SINGLE_REAL,f1,f2,phi,nbrec
        !real(kind=SINGLE_REAL) f1,f2,phi
        character(len=10) src_type,Mcomp
        integer irec


        select case (trim(src_type))
        case('monopole')
           f1=1.
           f2=0.

        case('dipole')

           select case (trim(Mcomp))
           case('mtr')
              do irec=1,nbrec
                 f1(irec)=cos(phi(irec))
                 f2(irec)=-sin(phi(irec))
              end do
           case ('thetaforce')
              do irec=1,nbrec
                 f1(irec)=cos(phi(irec))
                 f2(irec)=-sin(phi(irec))
              end do
           case('mpr')
              do irec=1,nbrec
                 f1(irec)=sin(phi(irec))
                 f2(irec)=cos(phi(irec))
              end do
           case('phiforce')
              do irec=1,nbrec
                 f1(irec)=sin(phi(irec))
                 f2(irec)=cos(phi(irec))
              end do
           end select

        case('quadpole')
           select case (trim(Mcomp))
           case ('mtt_m_mpp')
              do irec=1,nbrec
                 f1(irec)=cos(2.*phi(irec))
                 f2(irec)=-sin(2.*phi(irec))
              end do
           case('mtp')
              do irec=1,nbrec
                 f1(irec)=sin(2*phi(irec))
                 f2(irec)=cos(2*phi(irec))
              end do
           end select

        end select

      end subroutine compute_prefactor


      subroutine compute_3D_cyl()
        use global_parameters

        integer irec

        do irec=irecmin,irecmax
           data_rec(irec,1)=f1(irec)*data_rec(irec,1)
           data_rec(irec,2)=f2(irec)*data_rec(irec,2)
           data_rec(irec,3)=f1(irec)*data_rec(irec,3)
           !if (irec==10) write(*,*) data_rec(irec,3)
        end do

      end subroutine compute_3D_cyl

      subroutine compute_stress_3D_cyl()
        use global_parameters

        integer irec

        do irec=irecmin,irecmax
           stress_rec(irec,1)=f1(irec)*stress_rec(irec,1)
           stress_rec(irec,2)=f1(irec)*stress_rec(irec,2)
           stress_rec(irec,3)=f1(irec)*stress_rec(irec,3)
           stress_rec(irec,4)=f2(irec)*stress_rec(irec,4)
           stress_rec(irec,5)=f1(irec)*stress_rec(irec,5)
           stress_rec(irec,6)=f2(irec)*stress_rec(irec,6)
           !write(*,*) stress_rec(irec,:)
        end do

      end subroutine compute_stress_3D_cyl


      subroutine rotate2cartesian_with_source_in_pole_stress()
        use global_parameters

        integer irec,i,j,k
        real(kind=SINGLE_REAL) tmp(6,6),tmp1(3,3),B(3,3),st(3,3) 

        ! compute B*st*Bt
        do irec=irecmin,irecmax
          
           ! rotation matrix 
           B(1,1)=  cos(phi(irec))
           B(1,2)= - sin(phi(irec))
           B(1,3)= 0. 

           B(2,1)=  sin(phi(irec))
           B(2,2)=  cos(phi(irec))
           B(2,3)=  0.

           B(3,1)= 0. 
           B(3,2)= 0. 
           B(3,3)= 1. 
           !write(*,*) '1 ', stress_rec(irec,:)
           !write(*,*) ' matrix ', B
           ! stress in cylindical coordinates
           st(1,1)=stress_rec(irec,1)
           st(1,2)=stress_rec(irec,4)
           st(1,3)=stress_rec(irec,5)

           st(2,1)=stress_rec(irec,4)
           st(2,2)=stress_rec(irec,2)
           st(2,3)=stress_rec(irec,6)

           st(3,1)=stress_rec(irec,5)
           st(3,2)=stress_rec(irec,6)
           st(3,3)=stress_rec(irec,3)


          
           ! st*Bt
           tmp=0.
           do j=1,3
              do i=1,3
                 do k=1,3
                    tmp(i,j)=tmp(i,j)+st(i,k)*B(j,k)
                 end do
              end do
           end do

           ! B*st*Bt
           tmp1=0.
           do j=1,3
              do i=1,3
                 do k=1,3
                    tmp1(i,j)=tmp1(i,j)+B(i,k)*tmp(k,j)
                 end do
              end do
           end do
           
           ! stress in cartesian coordinates 
           stress_rec(irec,1)=tmp1(1,1)
           stress_rec(irec,2)=tmp1(2,2)
           stress_rec(irec,3)=tmp1(3,3)
           stress_rec(irec,4)=tmp1(1,2) 
           stress_rec(irec,5)=tmp1(1,3)
           stress_rec(irec,6)=tmp1(2,3)
           !write(*,*) '2 ', stress_rec(irec,:)
        end do
        
      end subroutine rotate2cartesian_with_source_in_pole_stress

!!$      subroutine rotate_back_source()
!!$        use global_parameters
!!$
!!$        integer irec
!!$        real(kind=SINGLE_REAL) tmp1,tmp2,tmp3
!!$
!!$        do irec=1,nbrec
!!$
!!$           tmp1 = rot_mat(1,1)*data_rec(irec,1) + rot_mat(1,2)*data_rec(irec,2) + rot_mat(1,3)*data_rec(irec,3)
!!$           tmp2 = rot_mat(2,1)*data_rec(irec,1) + rot_mat(2,2)*data_rec(irec,2) + rot_mat(2,3)*data_rec(irec,3)
!!$           tmp3 = rot_mat(3,1)*data_rec(irec,1) + rot_mat(3,2)*data_rec(irec,2) + rot_mat(3,3)*data_rec(irec,3)
!!$
!!$           data_rec(irec,1) = tmp1
!!$           data_rec(irec,2) = tmp2
!!$           data_rec(irec,3) = tmp3
!!$
!!$        end do
!!$
!!$
!!$      end subroutine rotate_back_source

     subroutine rotate_back_source()
        use global_parameters

        integer irec,i,j,k
        real(kind=SINGLE_REAL) tmp(3),veloc(3)

        do irec=irecmin,irecmax
           
           ! veloc in cylindical coordinates
           veloc(1)=data_rec(irec,1)
           veloc(2)=data_rec(irec,2)
           veloc(3)=data_rec(irec,3)

          
           !
           ! R*veloc
           tmp=0.  
           do i=1,3
              do k=1,3
                 tmp(i)=tmp(i)+veloc(k)*rot_mat(i,k)
              end do
           end do
          
           ! valocity in cartesian
           data_rec(irec,1)=tmp(1)
           data_rec(irec,2)=tmp(2)
           data_rec(irec,3)=tmp(3)
          

        end do


      end subroutine rotate_back_source

      

     subroutine rotate_back_to_local_cart()
        use global_parameters

        integer irec,i,j,k
        real(kind=SINGLE_REAL) tmp(3),veloc(3)

        do irec=irecmin,irecmax
           
           ! veloc in global coordinates
           veloc(1)=data_rec(irec,1)
           veloc(2)=data_rec(irec,2)
           veloc(3)=data_rec(irec,3)

          
           !
           ! Rt*veloc
           tmp=0.  
           do i=1,3
              do k=1,3
                 tmp(i)=tmp(i)+veloc(k)*trans_rot_mat_mesh(i,k)  
              end do
           end do
          
           ! valocity in cartesian
           data_rec(irec,1)=tmp(1)
           data_rec(irec,2)=tmp(2)
           data_rec(irec,3)=tmp(3)
          

        end do


      end subroutine rotate_back_to_local_cart


      subroutine rotate_back_source_stress()
        use global_parameters

        integer irec,i,j,k
        real(kind=SINGLE_REAL) tmp(3,3),tmp1(3,3),st(3,3)

        do irec=irecmin,irecmax
           
            ! stress in cylindical coordinates
           st(1,1)=stress_rec(irec,1)
           st(1,2)=stress_rec(irec,4)
           st(1,3)=stress_rec(irec,5)

           st(2,1)=stress_rec(irec,4)
           st(2,2)=stress_rec(irec,2)
           st(2,3)=stress_rec(irec,6)

           st(3,1)=stress_rec(irec,5)
           st(3,2)=stress_rec(irec,6)
           st(3,3)=stress_rec(irec,3)

           !
           ! 1 -> tmp =st*Rt
           tmp=0.
           do j=1,3
              do i=1,3
                 do k=1,3
                    tmp(i,j)=tmp(i,j)+st(i,k)*trans_rot_mat(k,j)
                 end do
              end do
           end do

           ! R*(st*Rt) =R*tmp
           tmp1=0.
           do j=1,3
              do i=1,3
                 do k=1,3
                    tmp1(i,j)=tmp1(i,j)+tmp(k,j)*rot_mat(i,k)
                 end do
              end do
           end do
       
           ! stress in cartesian
           stress_rec(irec,1)=tmp1(1,1)
           stress_rec(irec,2)=tmp1(2,2)
           stress_rec(irec,3)=tmp1(3,3)
           stress_rec(irec,4)=tmp1(1,2) 
           stress_rec(irec,5)=tmp1(1,3)
           stress_rec(irec,6)=tmp1(2,3)

        end do


      end subroutine rotate_back_source_stress

      subroutine rotate_back_to_local_cart_stress()
        use global_parameters

        integer irec,i,j,k
        real(kind=SINGLE_REAL) tmp(3,3),tmp1(3,3),st(3,3)

        do irec=irecmin,irecmax
           
            ! stress in cylindical coordinates
           st(1,1)=stress_rec(irec,1)
           st(1,2)=stress_rec(irec,4)
           st(1,3)=stress_rec(irec,5)

           st(2,1)=stress_rec(irec,4)
           st(2,2)=stress_rec(irec,2)
           st(2,3)=stress_rec(irec,6)

           st(3,1)=stress_rec(irec,5)
           st(3,2)=stress_rec(irec,6)
           st(3,3)=stress_rec(irec,3)

           !
           ! st*R
           tmp=0.
           do j=1,3
              do i=1,3
                 do k=1,3
                    tmp(i,j)=tmp(i,j)+st(i,k)*rot_mat_mesh(k,j)
                 end do
              end do
           end do

           ! Rt*st*R
           tmp1=0.
           do j=1,3
              do i=1,3
                 do k=1,3
                    tmp1(i,j)=tmp1(i,j)+trans_rot_mat_mesh(i,k)*tmp(k,j) 
                 end do
              end do
           end do
       
           ! stress in cartesian
           stress_rec(irec,1)=tmp1(1,1)
           stress_rec(irec,2)=tmp1(2,2)
           stress_rec(irec,3)=tmp1(3,3)
           stress_rec(irec,4)=tmp1(1,2) 
           stress_rec(irec,5)=tmp1(1,3)
           stress_rec(irec,6)=tmp1(2,3)

        end do


      end subroutine rotate_back_to_local_cart_stress


      subroutine rotate2cartesian_with_source_in_pole()
        use global_parameters
        integer irec,i,k
        real(kind=SINGLE_REAL) tmp(3),B(3,3),veloc(3)

        do irec=irecmin,irecmax
           
           ! rotation matrix
           B(1,1)=  cos(phi(irec))
           B(1,2)= - sin(phi(irec))
           B(1,3)= 0. 

           B(2,1)=  sin(phi(irec))
           B(2,2)=  cos(phi(irec))
           B(2,3)=  0.

           B(3,1)= 0. 
           B(3,2)= 0. 
           B(3,3)= 1. 

           ! veloc in cylindical coordinates
           veloc(1)=data_rec(irec,1)
           veloc(2)=data_rec(irec,2)
           veloc(3)=data_rec(irec,3)

          
           !
           ! B*veloc
           tmp=0.  
           do i=1,3
              do k=1,3
                 tmp(i)=tmp(i)+veloc(k)*B(i,k)
              end do
           end do
          
           ! valocity in cartesian
           data_rec(irec,1)=tmp(1)
           data_rec(irec,2)=tmp(2)
           data_rec(irec,3)=tmp(3)

        end do
      end subroutine rotate2cartesian_with_source_in_pole


    end module post_processing
    
