
! Below are my projection routines. These routines are working, I used them for the submitted GJI paper on inversion.
! I did some changes: to avoid projections artifacts
! I am minimizing the number of projections SEM mesh <-> tomo grid. So now I process like this  (-> means projection) :

! 1/ model SEM -> model TOMO  : this model is used for BFGS and regularisation term.
! 2/ gradient SEM -> gradient TOMO : used by BFGS

! 3/ Computing BFGS perturbation in tomo gird and then project it back in SEM mesh.

! 4/ Update model in SEM mesh.

! Previously I updated the model in TOMO grid and then project back the model in SEM mesh.
! I think it is better to update the model directly in SEM mesh because I found that the model projection SEM->TOMO->SEM->TOMO...
! generate some artifacts in model. (The result is not identical but it should be).
! Projecting only the perturbation form TOMO->SEM reduce the magnitudes of artifacts
! because perturbations are generaly small. Projecting the model TOMO->SEM generate bigger artifacts.

! So, you can see in the attached file the way I am doing the projections now:

! interpole_one_field + project_one_field : the two routines are used to project SEM->TOMO
! project_sem_interp : the routine used for projection TOMO->SEM

! I am working for integrating thoses routines in the current devel version of specfem in the directory ./inverse_problem
! created by Dimitri for this purpose, but I need some time (probably a few months).
! Since you want to perform inversion as soon as possible,
! probably the temporary solution is that you update your projection subroutines to use mine.

! Vadim Monteiller (March 2015).

module project_tomo_grid_mod


  use specfem_par
  implicit none
  include 'mpif.h'
  include 'precision.h'
  integer :: ier
  integer :: nx,ny,nz,nnx,nny,nnz,ngrid_local
  double precision, parameter :: R_EARTH=6371000.d0
  ! real array for data
  real(kind=CUSTOM_REAL),dimension(:,:,:,:),allocatable :: dat
  real(kind=CUSTOM_REAL),dimension(:,:,:),allocatable ::  vol_data,data_glob
  real(kind=CUSTOM_REAL),dimension(:,:,:),allocatable ::  data_tomo1,data_tomo2,data_tomo3,data_tomo4,data_tomo
  real(kind=CUSTOM_REAL),dimension(:,:,:),allocatable :: rho_grid,kappa_grid,mu_grid,alpha_grid,beta_grid
   real(kind=CUSTOM_REAL),dimension(:,:,:),allocatable ::  rho_grid_prior,alpha_grid_prior,beta_grid_prior
  integer, allocatable :: indice(:,:,:),indice_glob(:,:,:)
  integer, allocatable :: indice_grille_i(:,:,:,:),indice_grille_j(:,:,:,:),indice_grille_k(:,:,:,:)
  integer, allocatable :: indice_grille_i1(:,:,:,:),indice_grille_j1(:,:,:,:),indice_grille_k1(:,:,:,:)
  integer, allocatable :: indice_spec(:),ni(:)
  double precision, allocatable :: xstore_tomo_grid(:,:,:),xstore_tomo_grid_interp_tri(:,:,:)
  double precision, allocatable :: ystore_tomo_grid(:,:,:),ystore_tomo_grid_interp_tri(:,:,:),work_array(:,:,:)
  double precision, allocatable :: zstore_tomo_grid(:,:,:),zstore_tomo_grid_interp_tri(:,:,:),profondeur_tomo_grid(:,:,:)
  double precision, allocatable :: x_grid_gll(:,:,:,:), y_grid_gll(:,:,:,:),z_grid_gll(:,:,:,:)
  double precision, allocatable :: wk_reduce(:,:,:,:)
  real(kind=CUSTOM_REAL), allocatable :: valeur_integration1(:,:,:,:),valeur_integration2(:,:,:,:),valeur_integration(:,:,:,:)
  real(kind=CUSTOM_REAL), allocatable :: valeur_integration3(:,:,:,:),valeur_integration4(:,:,:,:)
  real(kind=CUSTOM_REAL), allocatable :: volume_integration(:,:,:,:),valeur_int_glob(:,:,:,:)
  real(kind=CUSTOM_REAL), allocatable :: sum_grad_alpha_kl(:,:,:,:),sum_grad_beta_kl(:,:,:,:)
  integer, allocatable :: indice_integration(:,:,:,:),indice_int_glob(:,:,:,:),indice_grid(:,:,:)
  double precision xmin,xmax,ymin,ymax,zmin,zmax
  ! 3D shape functions and their derivatives
  double precision shape3D(NGNOD,NGLLX,NGLLY,NGLLZ)
  double precision dershape3D(NDIM,NGNOD,NGLLX,NGLLY,NGLLZ)
  double precision valeur_fonction(NGLLX,NGLLY,NGLLZ)
  double precision xgll(NGLLX,NGLLY,NGLLZ),ygll(NGLLX,NGLLY,NGLLZ),zgll(NGLLX,NGLLY,NGLLZ)
  double precision,dimension(NGNOD) :: xstore_local, ystore_local, zstore_local
  double precision,dimension(NGNOD) :: xstore_local_int, ystore_local_int, zstore_local_int

  double precision  xmin_local,xmax_local,ymin_local,ymax_local,zmin_local,zmax_local
  double precision r_bot,r_cur,depth_box,pmin,pmax,px,py,pz,HLAGRANGE,valeur
  double precision  xi,ksi,eta,gamma,ratio_eta,ratio_xi,ksimin,ksimax,etamin,etamax
  double precision volume,volume_cellule
  double precision start, finish
  real(kind=CUSTOM_REAL) hx,hy,hz

  !! debug
  integer i666
  character(len=10) debug_file

  contains

    subroutine setup_tomo_grid()



      real(kind=CUSTOM_REAL)  xmins,xmaxs,ymins,ymaxs,zmins,zmaxs


      ! verifier que ce numero n'est pas deja utilise
      open(10,file='Par_for_projection_in_grid_tomo.par')
      read(10,*) depth_box
      read(10,*) nx,ny,nz
      close(10)

      ! integration grid
      allocate(indice(nx,ny,nz))
      allocate(indice_glob(nx,ny,nz))
      allocate(vol_data(nx,ny,nz))

      ! integration grid GLL points
      allocate(xstore_tomo_grid(nx,ny,nz))
      allocate(ystore_tomo_grid(nx,ny,nz))
      allocate(zstore_tomo_grid(nx,ny,nz))
      allocate(xstore_tomo_grid_interp_tri(nx+1,ny+1,nz+1))
      allocate(ystore_tomo_grid_interp_tri(nx+1,ny+1,nz+1))
      allocate(zstore_tomo_grid_interp_tri(nx+1,ny+1,nz+1))
      allocate(profondeur_tomo_grid(nx,ny,nz))
      allocate(work_array(nx,ny,nz))
      allocate(data_glob(nx+1,ny+1,nz+1))


      !! ces tableaux doivent etres distribues pour gagner de la memoire.
      !! il faut compter le nombre de points que chaque tranche mpi a dans la grille
      !allocate(valeur_integration(NGLLX,NGLLY,NGLLZ,nx*ny*nz))
      !allocate(valeur_integration1(NGLLX,NGLLY,NGLLZ,nx*ny*nz))
      !allocate(valeur_integration2(NGLLX,NGLLY,NGLLZ,nx*ny*nz))
      !allocate(valeur_integration3(NGLLX,NGLLY,NGLLZ,nx*ny*nz))
      !allocate(valeur_integration4(NGLLX,NGLLY,NGLLZ,nx*ny*nz))
      !allocate(volume_integration(NGLLX,NGLLY,NGLLZ,nx*ny*nz))
      !allocate(indice_integration(NGLLX,NGLLY,NGLLZ,nx*ny*nz))
      !allocate(valeur_int_glob(NGLLX,NGLLY,NGLLZ,nx*ny*nz))
      !allocate(indice_int_glob(NGLLX,NGLLY,NGLLZ,nx*ny*nz))


      ! model
      allocate(rho_grid(nx,ny,nz),kappa_grid(nx,ny,nz),mu_grid(nx,ny,nz))
      allocate(alpha_grid(nx,ny,nz),beta_grid(nx,ny,nz))
      allocate(rho_grid_prior(nx,ny,nz),alpha_grid_prior(nx,ny,nz),beta_grid_prior(nx,ny,nz))

      ! data projected onto tomographic grid
      allocate(data_tomo(nx,ny,nz))
      !allocate(data_tomo1(nx,ny,nz))
      !allocate(data_tomo2(nx,ny,nz))
      !allocate(data_tomo3(nx,ny,nz))
      !allocate(data_tomo4(nx,ny,nz))

      allocate(indice_spec(nx*ny*nz))
      indice_spec(:)=0
      allocate(indice_grid(3,1000,nspec_ab))
      indice_grid(:,:,:)=0
      allocate(ni(nspec_ab))
      ni(:)=0

      !allocate(x_grid_gll(NGLLX,NGLLY,NGLLZ,nx*ny*nz))
      !allocate(y_grid_gll(NGLLX,NGLLY,NGLLZ,nx*ny*nz))
      !allocate(z_grid_gll(NGLLX,NGLLY,NGLLZ,nx*ny*nz))
      !allocate(wk_reduce(NGLLX,NGLLY,NGLLZ,nx*ny*nz))

      allocate(dat(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
      if( ier /= 0 ) stop 'error allocating dat array'
      allocate(indice_grille_i(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
      if( ier /= 0 ) stop 'error allocating indice array'
      allocate(indice_grille_j(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
      if( ier /= 0 ) stop 'error allocating indice array'
      allocate(indice_grille_k(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
      if( ier /= 0 ) stop 'error allocating indice array'
      allocate(indice_grille_i1(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
      if( ier /= 0 ) stop 'error allocating indice array'
      allocate(indice_grille_j1(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
      if( ier /= 0 ) stop 'error allocating indice array'
      allocate(indice_grille_k1(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
      if( ier /= 0 ) stop 'error allocating indice array'

      allocate(sum_grad_alpha_kl(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
      if( ier /= 0 ) stop 'error allocating indice array'
      allocate(sum_grad_beta_kl(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
      if( ier /= 0 ) stop 'error allocating indice array'

      r_bot=  R_EARTH - depth_box

      ! on ajoute 1 pour le nombre de "piquets"
      nnx = nx + 1
      nny = ny + 1
      nnz = nz + 1

      call mpi_allreduce(minval(xstore(:)), xmins, 1, CUSTOM_MPI_TYPE, MPI_MIN, MPI_COMM_WORLD,ier)
      call mpi_allreduce(maxval(xstore(:)), xmaxs, 1, CUSTOM_MPI_TYPE, MPI_MAX, MPI_COMM_WORLD,ier)
      call mpi_allreduce(minval(ystore(:)), ymins, 1, CUSTOM_MPI_TYPE, MPI_MIN, MPI_COMM_WORLD,ier)
      call mpi_allreduce(maxval(ystore(:)), ymaxs, 1, CUSTOM_MPI_TYPE, MPI_MAX, MPI_COMM_WORLD,ier)
      call mpi_allreduce(minval(zstore(:)), zmins, 1, CUSTOM_MPI_TYPE, MPI_MIN, MPI_COMM_WORLD,ier)
      call mpi_allreduce(maxval(zstore(:)), zmaxs, 1, CUSTOM_MPI_TYPE, MPI_MAX, MPI_COMM_WORLD,ier)

      ! casting pour se mettre dans le "bon" type
      xmin=xmins
      xmax=xmaxs
      ymin=ymins
      ymax=ymaxs
      zmin=zmins
      zmax=zmaxs

!!$      ! weight GLL quadrature
!!$      call zwgljd(xigll,wxgll,NGLLX,GAUSSALPHA,GAUSSBETA)
!!$      call zwgljd(yigll,wygll,NGLLY,GAUSSALPHA,GAUSSBETA)
!!$      call zwgljd(zigll,wzgll,NGLLZ,GAUSSALPHA,GAUSSBETA)
      ! fonction shape
      call get_shape3D(0,shape3D,dershape3D,xigll,yigll,zigll)

      ! create cartesain tomo grid for inversion
      if (.not. COUPLING_WITH_DSM) then
         call create_box_grid_projection()
         if (myrank==0) write(*,*) ' box grid'
      else  ! chunk domain
         call create_chunk_grid_projection()
         if (myrank==0) write(*,*) ' chunk grid'
      endif
      call check_projection_grid()
    end subroutine setup_tomo_grid

!###################################################################################################

    subroutine  project_tomo_grid(dat_to_project,ires)
      implicit none
      integer ispec,ires,ix,iy,iz
      integer i,j,k,igll,jgll,kgll
      integer imin,imax,jmin,jmax,kmin,kmax
      integer index_cellule,ierr,np
      double precision  ANGULAR_WIDTH_ETA_RAD,ANGULAR_WIDTH_XI_RAD,deg2rad
      double precision, parameter :: R_EARTH=6371000.d0
      double precision profondeur_courante
      double precision rx,ry,rz
      double precision x,y,z,x_centred,y_centred,z_centred
      double precision x0,x1,y0,y1,p0,p1,p,volume,volume_cellule
      real(kind=CUSTOM_REAL)   dat_to_project(NGLLX,NGLLY,NGLLZ,NSPEC_AB)
      integer :: iii,jjj,kkk


      ! il faut definir dat (ie grad ou modele)
      write(debug_file,'(a3,i4.4)') 'dbp',myrank
      open(666,file=trim(debug_file))
      deg2rad = 3.141592653589793d0/180.d0

      ! initialization
      !! DEBUG
      !dat_to_project(:,:,:,:) = 10._CUSTOM_REAL
      !!
      data_tomo(:,:,:)=0._CUSTOM_REAL
      vol_data(:,:,:)=0._CUSTOM_REAL
      indice(:,:,:)=0
      indice_glob(:,:,:)=0
      indice_int_glob(:,:,:,:)=0
      valeur_int_glob(:,:,:,:)=0._CUSTOM_REAL
      valeur_integration(:,:,:,:) = 0._CUSTOM_REAL
      volume_integration(:,:,:,:) = 1._CUSTOM_REAL
      indice_integration(:,:,:,:) = 0
      indice_grille_i(:,:,:,:)=0
      indice_grille_j(:,:,:,:)=0
      indice_grille_k(:,:,:,:)=0

      wk_reduce(:,:,:,:)=0.d0
      xstore_tomo_grid(:,:,:)=0.d0
      ystore_tomo_grid(:,:,:)=0.d0
      zstore_tomo_grid(:,:,:)=0.d0
      profondeur_tomo_grid(:,:,:)=0.d0

!======== on se met dans la sphere cubique et on definit le pas des grilles

      ANGULAR_WIDTH_XI_RAD =  2.d0* dasin(dabs(xmin*1.0000)/R_EARTH)
      ANGULAR_WIDTH_ETA_RAD = 2.d0* dasin(dabs(ymin*1.0000)/R_EARTH)

      hx = ANGULAR_WIDTH_XI_RAD / (nx)
      hy = ANGULAR_WIDTH_ETA_RAD / (ny)
      hz = (r_earth-r_bot) / (nz)

      if (myrank==0) then
         !print *, 'Total number of points: ', np
         !print *, 'domain :'
         !print *, 'Xmin, Xmax :', xmin,xmax
         !print *, 'Ymin, Ymax :', ymin,ymax
         !print *, 'Zmin, Zmax :', zmin,zmax
         !print *, ' '
         !print *, 'angular width xi : ',ANGULAR_WIDTH_XI_RAD/deg2rad
         !print *, 'angular width eta : ',ANGULAR_WIDTH_ETA_RAD/deg2rad
         !print *, 'hx , hy, hz : ' , hx/deg2rad, hy/deg2rad, hz/1000.

         ! informations relatives a la grille
         open(10,file='grille.par')
         write(10,*)  nx,ny,nz
         write(10,*)  ANGULAR_WIDTH_XI_RAD, ANGULAR_WIDTH_ETA_RAD
         write(10,*)  r_earth,r_bot,hz
         write(10,*)
         close(10)
      endif

      if (myrank==0)   call cpu_time(start)
  !================================== 1/ interpolation grille tomographique  =========================================
      !write(666,*) 'nb el spectraux :' ,NSPEC_AB
      do ispec = 1, NSPEC_AB   !===== boucle element ispec

       xmin_local =  HUGEVAL
       xmax_local = -HUGEVAL
       ymin_local =  HUGEVAL
       ymax_local = -HUGEVAL
       zmin_local =  HUGEVAL
       zmax_local = -HUGEVAL

      !======================= Sommets de l'element ============
!1
       xstore_local(1)=xstore(ibool(1,1,1,ispec))
       ystore_local(1)=ystore(ibool(1,1,1,ispec))
       zstore_local(1)=zstore(ibool(1,1,1,ispec))
!2
       xstore_local(2)=xstore(ibool(NGLLX,1,1,ispec))
       ystore_local(2)=ystore(ibool(NGLLX,1,1,ispec))
       zstore_local(2)=zstore(ibool(NGLLX,1,1,ispec))
!3
       xstore_local(3)=xstore(ibool(NGLLX,NGLLY,1,ispec))
       ystore_local(3)=ystore(ibool(NGLLX,NGLLY,1,ispec))
       zstore_local(3)=zstore(ibool(NGLLX,NGLLY,1,ispec))
!4
       xstore_local(4)=xstore(ibool(1,NGLLY,1,ispec))
       ystore_local(4)=ystore(ibool(1,NGLLY,1,ispec))
       zstore_local(4)=zstore(ibool(1,NGLLY,1,ispec))
!5
       xstore_local(5)=xstore(ibool(1,1,NGLLZ,ispec))
       ystore_local(5)=ystore(ibool(1,1,NGLLZ,ispec))
       zstore_local(5)=zstore(ibool(1,1,NGLLZ,ispec))
!6
       xstore_local(6)=xstore(ibool(NGLLX,1,NGLLZ,ispec))
       ystore_local(6)=ystore(ibool(NGLLX,1,NGLLZ,ispec))
       zstore_local(6)=zstore(ibool(NGLLX,1,NGLLZ,ispec))
!7
       xstore_local(7)=xstore(ibool(NGLLX,NGLLY,NGLLZ,ispec))
       ystore_local(7)=ystore(ibool(NGLLX,NGLLY,NGLLZ,ispec))
       zstore_local(7)=zstore(ibool(NGLLX,NGLLY,NGLLZ,ispec))
!8
       xstore_local(8)=xstore(ibool(1,NGLLY,NGLLZ,ispec))
       ystore_local(8)=ystore(ibool(1,NGLLY,NGLLZ,ispec))
       zstore_local(8)=zstore(ibool(1,NGLLY,NGLLZ,ispec))
       !write(*,*) zstore_local
       !=============== on cherche le pave circonscrit a l'element en se basant sur les point GLL
       do iz=1,NGLLZ
         do iy=1,NGLLY
           do ix=1,NGLLX
             xmin_local=min(xmin_local,xstore(ibool(ix,iy,iz,ispec)))
             xmax_local=max(xmin_local,xstore(ibool(ix,iy,iz,ispec)))
             ymin_local=min(ymin_local,ystore(ibool(ix,iy,iz,ispec)))
             ymax_local=max(ymin_local,ystore(ibool(ix,iy,iz,ispec)))
             zmin_local=min(zmin_local,zstore(ibool(ix,iy,iz,ispec)))
             zmax_local=max(zmax_local,zstore(ibool(ix,iy,iz,ispec)))
           enddo
         enddo
       enddo

       ! =========== calcul du rayon min et max correspondant au pave
       call rayon_min_max(R_EARTH-zmax,xmin_local,xmax_local,ymin_local,ymax_local,zmin_local,zmax_local,pmin,pmax)

       ! indice sur grille fine d'integration (sphere cubique) !! ajout de +/-2 pour deborder sur l'element adjacent
       kmin = 1+ floor((pmin-r_bot)/(r_earth-r_bot)*(nnz-1) )  !- 2
       kmax = 1+ floor((pmax-r_bot)/(r_earth-r_bot)*(nnz-1) )  !+ 2

       call ksi_eta_min_max(R_EARTH-zmax,xmin_local,xmax_local,ymin_local,ymax_local,zmin_local,zmax_local,&
                             ksimin,ksimax,etamin,etamax)

       ! indice sur grille fine d'integration (sphere cubique)
       imin = floor(1. + (NNX - 1) * (ksimin + 0.5d0 *  ANGULAR_WIDTH_XI_RAD)  / ANGULAR_WIDTH_XI_RAD) !- 2
       jmin = floor(1. + (NNY - 1) * (etamin + 0.5d0 * ANGULAR_WIDTH_ETA_RAD) / ANGULAR_WIDTH_ETA_RAD) !- 2

       ! indice sur grille fine d'integration (sphere cubique)
       imax = floor(1. + (NNX - 1) * (ksimax + 0.5d0 *  ANGULAR_WIDTH_XI_RAD)  / ANGULAR_WIDTH_XI_RAD) !+ 2
       jmax = floor(1. + (NNY - 1) * (etamax + 0.5d0 * ANGULAR_WIDTH_ETA_RAD) / ANGULAR_WIDTH_ETA_RAD) !+ 2
       !write(666,*) ksimax ,  ANGULAR_WIDTH_XI_RAD , ((ksimax + 0.5d0 * ANGULAR_WIDTH_XI_RAD) / ANGULAR_WIDTH_XI_RAD)

       imin=max(imin,1)
       imax=min(imax,nx)
       jmin=max(jmin,1)
       jmax=min(jmax,ny)
       kmin=max(kmin,1)
       kmax=min(kmax,nz)

       !================calcul de la valeur du noyau au centre de chaque cellule de la grille tomo
       !write(666,*) 'imin ', imin,ksimin
       !write(666,*) 'imax ', imax,ksimax,nx-1
       !write(666,*) xmin_local,xmax_local,ymin_local,ymax_local,zmin_local,zmax_local
       !write(666,*) 'jmin ', jmin
       !write(666,*) 'jmax ', jmax
       !write(666,*) 'kmin ', kmin
       !write(666,*) 'kmax ', kmax

       do k=kmin,kmax

         z_centred = r_bot + (k-0.5d0)*hz  ! valeur au centre de la cellule (sphere cubique)
         profondeur_courante = R_EARTH  - z_centred
         !write(*,*) R_EARTH, z_centred, profondeur_courante
         do j=jmin,jmax

          ratio_eta = (dble(j)-0.5d0) / dble(NY)
          y_centred = 2.d0*ratio_eta-1
          y_centred = tan((ANGULAR_WIDTH_ETA_RAD/2.d0) * y_centred)


          do i=imin,imax

            index_cellule = i + (j-1)*(nnx-1) + (k-1)*(nnx-1)*(nny-1)

            ratio_xi = (dble(i)-0.5d0) / dble(NX)
            x_centred = 2.d0*ratio_xi-1
            x_centred = tan((ANGULAR_WIDTH_XI_RAD/2.d0) * x_centred)
!!$            write(*,*) i,j,k
!!$            write(*,*) ratio_xi,dble(NX)
!!$            write(*,*) ANGULAR_WIDTH_XI_RAD/2.d0
!!$            write(*,*) i, 2.d0*ratio_xi-1,x
!!$            write(*,*) i,j,k
!!$            write(*,*) ratio_eta,dble(NY)
!!$            write(*,*) ANGULAR_WIDTH_ETA_RAD/2.d0
!!$            write(*,*) j, 2.d0*ratio_eta-1,y


            ! on repasse en cartesien pour travailler dans le repere SEM
            px = x_centred*z_centred
            py = y_centred*z_centred
            pz =  -(R_EARTH - z_centred/dsqrt(1.d0 + y_centred**2 + x_centred**2) - zmax)
!!$            write(*,*) px/1000.,py/1000.,pz/1000.
!!$            read(*,*) jjj
            ! on cherche le xix, eta et gamma du point central qui est px,py,pz
            call Find_xix_eta_gamma(xstore_local,ystore_local,zstore_local,xi,eta,gamma,px,py,pz)

            ! on regarde si on est bien dans l'element ispec pour faire la correspondance SEM <=> tomo
            if (dabs(xi)<=1.05d0.and.dabs(eta)<=1.05d0.and.dabs(gamma)<=1.05d0) then

!!$               do kgll = 1,NGLLZ
!!$                 do jgll = 1,NGLLY
!!$                   do igll = 1,NGLLX
!!$
!!$                      indice_grille_i(igll,jgll,kgll,ispec) = i
!!$                      indice_grille_j(igll,jgll,kgll,ispec) = j
!!$                      indice_grille_k(igll,jgll,kgll,ispec) = k
!!$
!!$                   enddo
!!$                 enddo
!!$               enddo
               !write(666,*) 'found ',i,j,k
               !write(666,*)
               !if (j==ny) write(*,*) 'py ',py/1000.
               !if (i==nx) write(*,*) 'px ',px/1000.
               indice(i,j,k)=1 ! on a bien visite la cellule de la grille
               xstore_tomo_grid(i,j,k)=px
               ystore_tomo_grid(i,j,k)=py
               zstore_tomo_grid(i,j,k)=pz
               profondeur_tomo_grid(i,j,k)=profondeur_courante
!!$               IF (dabs(gamma)>1.) THEN
!!$                  write(*,*)
!!$                  write(*,'(3f12.3,4x,4i8)') px/1000.,py/1000.,pz/1000.,i,j,k,ispec
!!$                  do jjj=1,8
!!$                     write(*,'(3f12.3)') xstore_local(jjj)/1000.,ystore_local(jjj)/1000.,zstore_local(jjj)/1000.
!!$                  enddo
!!$               endif
!!$               read(*,*) jjj
            else
               !write(666,*) xi,eta,gamma
               !write(666,*) px,py,pz
               !write(666,*) 'not found ',i,j,k

               !do i666=1,8
               ! write(666,*)xstore_local(i666),ystore_local(i666),zstore_local(i666)
               !enddo
               !write(666,*)
            endif

            ! == calcul des 8 sommets de la cellule d'integration =====================
            ! 1
            rx=1.d0;ry=1.d0;rz=1.d0

            pz = r_bot + (dble(k)-rz)*hz

            ratio_eta = (dble(j)-ry) / dble(NY)
            y = 2.d0*ratio_eta-1;y = tan((ANGULAR_WIDTH_ETA_RAD/2.d0) * y)

            ratio_xi = (dble(i)-rx) / dble(NX)
            x = 2.d0*ratio_xi-1;x = tan((ANGULAR_WIDTH_XI_RAD/2.d0) * x)

            px = x*pz; py = y*pz
            pz =  -(R_EARTH - pz/dsqrt(1.d0 + y*y + x*x) - zmax)



            xstore_local_int(1) = px
            ystore_local_int(1) = py
            zstore_local_int(1) = pz

            ! 2
            rx=0.d0;ry=1.d0;rz=1.d0

            pz = r_bot + (dble(k)-rz)*hz

            ratio_eta = (dble(j)-ry) / dble(NY)
            y = 2.d0*ratio_eta-1;y = tan((ANGULAR_WIDTH_ETA_RAD/2.d0) * y)

            ratio_xi = (dble(i)-rx) / dble(NX)
            x = 2.d0*ratio_xi-1;x = tan((ANGULAR_WIDTH_XI_RAD/2.d0) * x)

            px = x*pz; py = y*pz
            pz =  -(R_EARTH - pz/dsqrt(1.d0 + y*y + x*x) - zmax)



            xstore_local_int(2) = px
            ystore_local_int(2) = py
            zstore_local_int(2) = pz

            ! 3
            rx=0.d0;ry=0.d0;rz=1.d0

            pz = r_bot + (dble(k)-rz)*hz

            ratio_eta = (dble(j)-ry) / dble(NY)
            y = 2.d0*ratio_eta-1;y = tan((ANGULAR_WIDTH_ETA_RAD/2.d0) * y)

            ratio_xi = (dble(i)-rx) / dble(NX)
            x = 2.d0*ratio_xi-1;x = tan((ANGULAR_WIDTH_XI_RAD/2.d0) * x)

            px = x*pz; py = y*pz
            pz =  -(R_EARTH - pz/dsqrt(1.d0 + y*y + x*x) - zmax)



            xstore_local_int(3) = px
            ystore_local_int(3) = py
            zstore_local_int(3) = pz


            ! 4
            rx=1.d0;ry=0.d0;rz=1.d0

            pz = r_bot + (dble(k)-rz)*hz

            ratio_eta = (dble(j)-ry) / dble(NY)
            y = 2.d0*ratio_eta-1;y = tan((ANGULAR_WIDTH_ETA_RAD/2.d0) * y)

            ratio_xi = (dble(i)-rx) / dble(NX)
            x = 2.d0*ratio_xi-1;x = tan((ANGULAR_WIDTH_XI_RAD/2.d0) * x)

            px = x*pz; py = y*pz
            pz =  -(R_EARTH - pz/dsqrt(1.d0 + y*y + x*x) - zmax)



            xstore_local_int(4) = px
            ystore_local_int(4) = py
            zstore_local_int(4) = pz

           ! 5
            rx=1.d0;ry=1.d0;rz=0.d0

            pz = r_bot + (dble(k)-rz)*hz

            ratio_eta = (dble(j)-ry) / dble(NY)
            y = 2.d0*ratio_eta-1;y = tan((ANGULAR_WIDTH_ETA_RAD/2.d0) * y)

            ratio_xi = (dble(i)-rx) / dble(NX)
            x = 2.d0*ratio_xi-1;x = tan((ANGULAR_WIDTH_XI_RAD/2.d0) * x)

            px = x*pz; py = y*pz
            pz =  -(R_EARTH - pz/dsqrt(1.d0 + y*y + x*x) - zmax)



            xstore_local_int(5) = px
            ystore_local_int(5) = py
            zstore_local_int(5) = pz



            ! 6
            rx=0.d0;ry=1.d0;rz=0.d0

            pz = r_bot + (dble(k)-rz)*hz

            ratio_eta = (dble(j)-ry) / dble(NY)
            y = 2.d0*ratio_eta-1;y = tan((ANGULAR_WIDTH_ETA_RAD/2.d0) * y)

            ratio_xi = (dble(i)-rx) / dble(NX)
            x = 2.d0*ratio_xi-1;x = tan((ANGULAR_WIDTH_XI_RAD/2.d0) * x)

            px = x*pz; py = y*pz
            pz =  -(R_EARTH - pz/dsqrt(1.d0 + y*y + x*x) - zmax)



            xstore_local_int(6) = px
            ystore_local_int(6) = py
            zstore_local_int(6) = pz

            ! 7
            rx=0.d0;ry=0.d0;rz=0.d0

            pz = r_bot + (dble(k)-rz)*hz

            ratio_eta = (dble(j)-ry) / dble(NY)
            y = 2.d0*ratio_eta-1;y = tan((ANGULAR_WIDTH_ETA_RAD/2.d0) * y)

            ratio_xi = (dble(i)-rx) / dble(NX)
            x = 2.d0*ratio_xi-1;x = tan((ANGULAR_WIDTH_XI_RAD/2.d0) * x)

            px = x*pz; py = y*pz
            pz =  -(R_EARTH - pz/dsqrt(1.d0 + y*y + x*x) - zmax)



            xstore_local_int(7) = px
            ystore_local_int(7) = py
            zstore_local_int(7) = pz


            ! 8
            rx=1.d0;ry=0.d0;rz=0.d0

            pz = r_bot + (dble(k)-rz)*hz

            ratio_eta = (dble(j)-ry) / dble(NY)
            y = 2.d0*ratio_eta-1;y = tan((ANGULAR_WIDTH_ETA_RAD/2.d0) * y)

            ratio_xi = (dble(i)-rx) / dble(NX)
            x = 2.d0*ratio_xi-1;x = tan((ANGULAR_WIDTH_XI_RAD/2.d0) * x)

            px = x*pz; py = y*pz
            pz =  -(R_EARTH - pz/dsqrt(1.d0 + y*y + x*x) - zmax)



            xstore_local_int(8) = px
            ystore_local_int(8) = py
            zstore_local_int(8) = pz

            ! ========== calcul des points GLL
            call calcule_point_gll(xgll,ygll,zgll,xstore_local_int,ystore_local_int,zstore_local_int,shape3D)

            ! ========== interpolation de la fonction dat aux points GLL
            do kgll = 1,NGLLZ
               do jgll = 1,NGLLY
                  do igll = 1,NGLLX



                    px = xgll(igll,jgll,kgll)
                    py = ygll(igll,jgll,kgll)
                    pz = zgll(igll,jgll,kgll)

                    call Find_xix_eta_gamma(xstore_local,ystore_local,zstore_local,xi,eta,gamma,px,py,pz)

                    if (dabs(xi)<=1.05d0.and.dabs(eta)<=1.05d0.and.dabs(gamma)<=1.05d0) then

                       call interpolation_valeur(valeur,px,py,pz,dat_to_project,xstore_local,ystore_local,zstore_local,&
                                                 ispec,NSPEC_AB,xigll,yigll,zigll)
                       valeur_integration(igll,jgll,kgll,index_cellule)=valeur
                       indice_integration(igll,jgll,kgll,index_cellule)=1
                       !if (xi > 1.) then
!!$                       if (i==18 .and. j==21 .and. k==13) then
!!$                          write(*,*) igll,jgll,kgll,index_cellule
!!$                          write(*,*) xi,eta,gamma
!!$                          write(*,*) px/1000.,py/1000.,pz/1000.
!!$                          do i666=1,8
!!$                             write(*,*)xstore_local(i666)/1000.,ystore_local(i666)/1000.,zstore_local(i666)/1000.
!!$                          enddo
!!$                          write(*,*)
!!$                          write(*,*) '----- > : ', valeur
!!$                          do kkk=1,5
!!$                             do jjj=1,5
!!$                                do iii=1,5
!!$                                   write(*,*) iii,jjj,kkk,ispec,dat_to_project(iii,jjj,kkk,ispec)
!!$                                enddo
!!$                             enddo
!!$                          enddo
!!$                          write(*,*) '------------------------------------------'

!!$                          read(*,*) jjj
!!$                       endif
                    else
                      !write(666,*) xi,eta,gamma
                      !do i666=1,8
                      !  write(666,*)xstore_local(i666),ystore_local(i666),zstore_local(i666)
                      !enddo
                      !write(666,*)


                    endif

                  enddo
               enddo
            enddo
            ! ========================= on cherche les points GLL qui sont dans la cellule de la grille tomo
             call calcule_point_gll(xgll,ygll,zgll,xstore_local,ystore_local,zstore_local,shape3D)
             do kgll = 1,NGLLZ
                do jgll = 1,NGLLY
                   do igll = 1,NGLLX



                    px = xgll(igll,jgll,kgll)
                    py = ygll(igll,jgll,kgll)
                    pz = zgll(igll,jgll,kgll)
                    call Find_xix_eta_gamma(xstore_local_int,ystore_local_int,zstore_local_int,xi,eta,gamma,px,py,pz)

                    if (dabs(xi)<=1.05d0.and.dabs(eta)<=1.05d0.and.dabs(gamma)<=1.05d0) then

                       indice_grille_i(igll,jgll,kgll,ispec) = i
                       indice_grille_i1(igll,jgll,kgll,ispec)= i

                       indice_grille_j(igll,jgll,kgll,ispec) = j
                       indice_grille_j1(igll,jgll,kgll,ispec)= j

                       indice_grille_k(igll,jgll,kgll,ispec) = k
                       indice_grille_k1(igll,jgll,kgll,ispec)= k

                       if (xi > 0. ) then
                          if (i<nx) indice_grille_i1(igll,jgll,kgll,ispec) = i+1
                       else
                          if (i>1)  indice_grille_i1(igll,jgll,kgll,ispec) = i-1
                       endif

                       if (eta > 0. ) then
                          if (j<ny) indice_grille_j1(igll,jgll,kgll,ispec) = j+1
                       else
                          if (j>1)  indice_grille_j1(igll,jgll,kgll,ispec) = j-1
                       endif

                       if (gamma > 0. ) then
                          if (k<nz) indice_grille_k1(igll,jgll,kgll,ispec) = k+1
                       else
                          if (k>1)  indice_grille_k1(igll,jgll,kgll,ispec) = k-1
                       endif

                    endif

                 enddo
              enddo
           enddo
            ! ===========================================================


          enddo
         enddo
       enddo

    enddo !==== boucle ispec
    if (myrank==0) then
       call cpu_time(finish)
       write(*,*) 'time // proj 0: ', finish-start
       call cpu_time(start)
    endif
    !! === verif des correspondances (gll,ispec) -> tomo
    do ispec=1,nspec_ab
       do kgll = 1,NGLLZ
          do jgll = 1,NGLLY
             do igll = 1,NGLLX

                 !! patch ad hoc pour effet de bords
                 if (indice_grille_i(igll,jgll,kgll,ispec)==0 .and. igll==1) then
                    indice_grille_i(igll,jgll,kgll,ispec)=indice_grille_i(igll+1,jgll,kgll,ispec)
                    indice_grille_j(igll,jgll,kgll,ispec)=indice_grille_j(igll+1,jgll,kgll,ispec)
                    indice_grille_k(igll,jgll,kgll,ispec)=indice_grille_k(igll+1,jgll,kgll,ispec)
                 endif

                  if (indice_grille_i(igll,jgll,kgll,ispec)==0 .and. igll==5) then
                    indice_grille_i(igll,jgll,kgll,ispec)=indice_grille_i(igll-1,jgll,kgll,ispec)
                    indice_grille_j(igll,jgll,kgll,ispec)=indice_grille_j(igll-1,jgll,kgll,ispec)
                    indice_grille_k(igll,jgll,kgll,ispec)=indice_grille_k(igll-1,jgll,kgll,ispec)
                 endif


                 if (indice_grille_i(igll,jgll,kgll,ispec)==0 .and. jgll==1) then
                    indice_grille_i(igll,jgll,kgll,ispec)=indice_grille_i(igll,jgll+1,kgll,ispec)
                    indice_grille_j(igll,jgll,kgll,ispec)=indice_grille_j(igll,jgll+1,kgll,ispec)
                    indice_grille_k(igll,jgll,kgll,ispec)=indice_grille_k(igll,jgll+1,kgll,ispec)
                 endif


                 if (indice_grille_i(igll,jgll,kgll,ispec)==0 .and. jgll==5) then
                    indice_grille_i(igll,jgll,kgll,ispec)=indice_grille_i(igll,jgll-1,kgll,ispec)
                    indice_grille_j(igll,jgll,kgll,ispec)=indice_grille_j(igll,jgll-1,kgll,ispec)
                    indice_grille_k(igll,jgll,kgll,ispec)=indice_grille_k(igll,jgll-1,kgll,ispec)
                 endif

                 if (indice_grille_i(igll,jgll,kgll,ispec)==0 .and. kgll==1) then
                    indice_grille_i(igll,jgll,kgll,ispec)=indice_grille_i(igll,jgll,kgll+1,ispec)
                    indice_grille_j(igll,jgll,kgll,ispec)=indice_grille_j(igll,jgll,kgll+1,ispec)
                    indice_grille_k(igll,jgll,kgll,ispec)=indice_grille_k(igll,jgll,kgll+1,ispec)
                 endif


                 if (indice_grille_i(igll,jgll,kgll,ispec)==0 .and. kgll==5) then
                    indice_grille_i(igll,jgll,kgll,ispec)=indice_grille_i(igll,jgll,kgll-1,ispec)
                    indice_grille_j(igll,jgll,kgll,ispec)=indice_grille_j(igll,jgll,kgll-1,ispec)
                    indice_grille_k(igll,jgll,kgll,ispec)=indice_grille_k(igll,jgll,kgll-1,ispec)
                 endif


                 if (indice_grille_i(igll,jgll,kgll,ispec)==0 .and. jgll==1 .and. kgll==1) then
                    indice_grille_i(igll,jgll,kgll,ispec)=indice_grille_i(igll,jgll+1,kgll+1,ispec)
                    indice_grille_j(igll,jgll,kgll,ispec)=indice_grille_j(igll,jgll+1,kgll+1,ispec)
                    indice_grille_k(igll,jgll,kgll,ispec)=indice_grille_k(igll,jgll+1,kgll+1,ispec)
                 endif


                 if (indice_grille_i(igll,jgll,kgll,ispec)==0 .and. jgll==5 .and. kgll==5) then
                    indice_grille_i(igll,jgll,kgll,ispec)=indice_grille_i(igll,jgll-1,kgll-1,ispec)
                    indice_grille_j(igll,jgll,kgll,ispec)=indice_grille_j(igll,jgll-1,kgll-1,ispec)
                    indice_grille_k(igll,jgll,kgll,ispec)=indice_grille_k(igll,jgll-1,kgll-1,ispec)
                 endif

                 if (indice_grille_i(igll,jgll,kgll,ispec)==0 .and. igll==1 .and. jgll==1) then
                    indice_grille_i(igll,jgll,kgll,ispec)=indice_grille_i(igll+1,jgll+1,kgll,ispec)
                    indice_grille_j(igll,jgll,kgll,ispec)=indice_grille_j(igll+1,jgll+1,kgll,ispec)
                    indice_grille_k(igll,jgll,kgll,ispec)=indice_grille_k(igll+1,jgll+1,kgll,ispec)
                 endif


                 if (indice_grille_i(igll,jgll,kgll,ispec)==0 .and. igll==5 .and. jgll==5) then
                    indice_grille_i(igll,jgll,kgll,ispec)=indice_grille_i(igll-1,jgll-1,kgll,ispec)
                    indice_grille_j(igll,jgll,kgll,ispec)=indice_grille_j(igll-1,jgll-1,kgll,ispec)
                    indice_grille_k(igll,jgll,kgll,ispec)=indice_grille_k(igll-1,jgll-1,kgll,ispec)
                 endif


                 write(666,'(3i2,i6,3i5)') igll,jgll,kgll,ispec,&
                      indice_grille_i(igll,jgll,kgll,ispec),&
                      indice_grille_j(igll,jgll,kgll,ispec),&
                      indice_grille_k(igll,jgll,kgll,ispec)

             enddo
          enddo
       enddo
    enddo
    if (myrank==0) then
       call cpu_time(finish)
       write(*,*) 'time // proj : ', finish-start
    endif
    call mpi_reduce(valeur_integration,valeur_int_glob,NGLLX*NGLLY*NGLLZ*(nnx-1)*(nny-1)*(nnz-1),CUSTOM_MPI_TYPE, MPI_SUM,0,MPI_COMM_WORLD,ier)
    call mpi_reduce(indice_integration,indice_int_glob,NGLLX*NGLLY*NGLLZ*(nnx-1)*(nny-1)*(nnz-1),MPI_INTEGER, MPI_SUM,0,MPI_COMM_WORLD,ier)
    call mpi_reduce(indice,indice_glob,(nnx-1)*(nny-1)*(nnz-1),MPI_INTEGER ,MPI_SUM,0, MPI_COMM_WORLD, ierr)

    call  mpi_reduce(xstore_tomo_grid,work_array,(nnx-1)*(nny-1)*(nnz-1),MPI_DOUBLE_PRECISION ,MPI_SUM,0, MPI_COMM_WORLD, ierr)
    xstore_tomo_grid(:,:,:)=work_array(:,:,:)

    call  mpi_reduce(ystore_tomo_grid,work_array,(nnx-1)*(nny-1)*(nnz-1),MPI_DOUBLE_PRECISION ,MPI_SUM,0, MPI_COMM_WORLD, ierr)
    ystore_tomo_grid(:,:,:)=work_array(:,:,:)

    call  mpi_reduce(zstore_tomo_grid,work_array,(nnx-1)*(nny-1)*(nnz-1),MPI_DOUBLE_PRECISION ,MPI_SUM,0, MPI_COMM_WORLD, ierr)
    zstore_tomo_grid(:,:,:)=work_array(:,:,:)

    call  mpi_reduce(profondeur_tomo_grid,work_array,(nnx-1)*(nny-1)*(nnz-1),MPI_DOUBLE_PRECISION ,MPI_SUM,0, MPI_COMM_WORLD, ierr)
    profondeur_tomo_grid(:,:,:)=work_array(:,:,:)

    indice(:,:,:) = indice_glob(:,:,:)
    indice_integration(:,:,:,:) = indice_int_glob(:,:,:,:)
    valeur_integration(:,:,:,:) = valeur_int_glob(:,:,:,:)

!================== 2/ projection sur grille tomo =====================================

!================== 2/ projection sur grille tomo =====================================

    if (myrank==0) then  !======================== myrank=0
       call  cpu_time(start)
       index_cellule=0
       do k=1,nnz-1
          do j=1,nny-1
             do i=1,nnx-1
                index_cellule=index_cellule+1

                if (indice(i,j,k)/=0)  then
                   xstore_tomo_grid(i,j,k) =  xstore_tomo_grid(i,j,k) / real(indice(i,j,k),8)
                   ystore_tomo_grid(i,j,k) =  ystore_tomo_grid(i,j,k) / real(indice(i,j,k),8)
                   zstore_tomo_grid(i,j,k) =  zstore_tomo_grid(i,j,k) / real(indice(i,j,k),8)
                   !write(*,*) profondeur_tomo_grid(i,j,k) ,profondeur_tomo_grid(i,j,k) / real(indice(i,j,k),8)
                   profondeur_tomo_grid(i,j,k) = profondeur_tomo_grid(i,j,k) / real(indice(i,j,k),8)
                else
                   write(*,*) 'point oublie :'  ,i,j,k
                endif


                !! on reagarde si les points on ete visites et combien de fois?
                do kgll=1,NGLLZ
                   do jgll=1,NGLLY
                      do igll=1,NGLLX

                         if (indice_integration(igll,jgll,kgll,index_cellule) == 0) then  !! le point n'a pas ete visite
                            !write(666,*) 'point oublie :',igll,jgll,kgll,i,j,k,index_cellule
                            !write(667,*) i,j,k
                            !! chercher le point le plus proche qui est a ete visite
                            !  premier patch : devrait marcher pour mon cas cartesien
                            if (kgll==1 .and. indice_integration(igll,jgll,kgll+1,index_cellule) > 0) then
                               valeur_integration(igll,jgll,kgll,index_cellule) = valeur_integration(igll,jgll,kgll+1,index_cellule) /&
                                    real(indice_integration(igll,jgll,kgll+1,index_cellule),CUSTOM_REAL)
                            endif
                         else  ! le point a ete visite "indice_integration" fois
!!$                            if (j==21 .and. k==13) then
!!$                               write(*,*) i,valeur_integration(igll,jgll,kgll,index_cellule)
!!$                            endif
                            valeur_integration(igll,jgll,kgll,index_cellule) =  valeur_integration(igll,jgll,kgll,index_cellule)/&
                                 real(indice_integration(igll,jgll,kgll,index_cellule),CUSTOM_REAL)
                            ! on divise par indice pour avoir la valeur moyenne car on a fait un MPI_SUM
                         endif

                      enddo
                   enddo
                enddo
!!$                if (j==21 .and. k==13) then
!!$                   write(*,*)
!!$                   read(*,*) jjj
!!$                endif
             enddo
          enddo
       enddo


  !===== boucle grille tomo : verifier qu'on fait exactement la meme boucle que precedemment
  np=0
  do k=1,nnz-1  ! grille integration (sphere cubique)


    p0= r_bot + (k-1)*hz
    p1= r_bot + k*hz

    do j=2,nny ! grille integration (sphere cubique)


       ratio_eta = (dble(j-2)) / dble(NY)
       y0 = 2.d0*ratio_eta-1
       y0 = tan((ANGULAR_WIDTH_ETA_RAD/2.d0) * y0)

       ratio_eta = (dble(j-1)) / dble(NY)
       y1 = 2.d0*ratio_eta-1
       y1 = tan((ANGULAR_WIDTH_ETA_RAD/2.d0) * y1)


       do i=2,nnx ! grille integration (sphere cubique)

          index_cellule=i-1+(j-2)*(nnx-1)+(k-1)*(nnx-1)*(nny-1)

          ratio_xi = (dble(i-2)) / dble(NX)
          x0 = 2.d0*ratio_xi-1
          x0 = tan((ANGULAR_WIDTH_XI_RAD/2.d0) * x0)

          ratio_xi = (dble(i-1)) / dble(NX)
          x1 = 2.d0*ratio_xi-1
          x1 = tan((ANGULAR_WIDTH_XI_RAD/2.d0) * x1)

          !===== calcul des 8 sommets de la cellule d'integration dans le repere cartesien SEM

          ! 1
          x=x0;y=y0;p=p0
          xstore_local(1) = x*p
          ystore_local(1) = y*p
          zstore_local(1) =  -(R_EARTH - p/dsqrt(1.d0 + y*y + x*x) - zmax)

          ! 2
          x=x1;y=y0;p=p0
          xstore_local(2) = x*p
          ystore_local(2) = y*p
          zstore_local(2) =  -(R_EARTH - p/dsqrt(1.d0 + y*y + x*x) - zmax)

          ! 3
          x=x1;y=y1;p=p0
          xstore_local(3) = x*p
          ystore_local(3) = y*p
          zstore_local(3) =  -(R_EARTH - p/dsqrt(1.d0 + y*y + x*x) - zmax)

          ! 4
          x=x0;y=y1;p=p0
          xstore_local(4) = x*p
          ystore_local(4) = y*p
          zstore_local(4) =  -(R_EARTH - p/dsqrt(1.d0 + y*y + x*x) - zmax)

          ! 5
          x=x0;y=y0;p=p1
          xstore_local(5) = x*p
          ystore_local(5) = y*p
          zstore_local(5) =  -(R_EARTH - p/dsqrt(1.d0 + y*y + x*x) - zmax)

          ! 6
          x=x1;y=y0;p=p1
          xstore_local(6) = x*p
          ystore_local(6) = y*p
          zstore_local(6) =  -(R_EARTH - p/dsqrt(1.d0 + y*y + x*x) - zmax)

          ! 7
          x=x1;y=y1;p=p1
          xstore_local(7) = x*p
          ystore_local(7) = y*p
          zstore_local(7) =  -(R_EARTH - p/dsqrt(1.d0 + y*y + x*x) - zmax)


          ! 8
          x=x0;y=y1;p=p1
          xstore_local(8) = x*p
          ystore_local(8) = y*p
          zstore_local(8) =  -(R_EARTH - p/dsqrt(1.d0 + y*y + x*x) - zmax)
          ! === calcul de l'integrale de volume
          call compute_volume_int(volume,xstore_local,ystore_local,zstore_local,wxgll,wygll,wzgll,&
               shape3D,dershape3D,valeur_integration,index_cellule,(nnx-1)*(nny-1)*(nnz-1))
!!$          write(*,*)
!!$          do i666=1,8
!!$             write(*,*) xstore_local(i666),ystore_local(i666),zstore_local(i666)
!!$          enddo
!!$          write(*,*)
!!$          write(*,*) xstore_tomo_grid(i-1,j-1,k),ystore_tomo_grid(i-1,j-1,k), zstore_tomo_grid(i-1,j-1,k)
!!$          write(*,*)
!!$          read(*,*) jjj
          if (ires==1)  then
             ! === si on veut faire juste la projection du modele, on ne doit pas normaliser par
             ! === le volume de la cellule, on prend plutot la valeur moyenne.
             call compute_volume_int(volume_cellule,xstore_local,ystore_local,zstore_local,wxgll,wygll,wzgll,&
                  shape3D,dershape3D,volume_integration,index_cellule,(nnx-1)*(nny-1)*(nnz-1))
             volume=volume/volume_cellule
          endif

          if (indice(i-1,j-1,k)==0) then
           ! on a oublie le point i,j,k
             np=np+1
           !write(30,*) i,j,k
          else
             !zscaling(k) = zscaling(k) + volume*vol_data(i-1,j-1,k)
             data_tomo(i-1,j-1,k) = volume !data_tomo(i_tomo,j_tomo,k_tomo) + volume*vol_data(i-1,j-1,k)
             !zscaling_tomo(k_tomo) = zscaling_tomo(k_tomo) + volume*vol_data(i-1,j-1,k)
          endif

         enddo
      enddo

    enddo
    call cpu_time(finish)
    write(*,*) 'time seq. ',finish-start
!! pour debbug
    open(667,file='profondeur.bin',access='direct',recl=8*nx*ny*nz)
    write(667,rec=1) profondeur_tomo_grid
    close(667)

    open(667,file='valeur_int.bin',access='direct',recl=CUSTOM_REAL*NGLLX*NGLLY*NGLLZ*nx*ny*nz)
    write(667,rec=1) valeur_integration
    close(667)

    open(667,file='xsem.bin',access='direct',recl=CUSTOM_REAL*NGLOB_AB)
    write(667,rec=1) xstore
    close(667)

    open(667,file='ysem.bin',access='direct',recl=CUSTOM_REAL*NGLOB_AB)
    write(667,rec=1) ystore
    close(667)

    open(667,file='zsem.bin',access='direct',recl=CUSTOM_REAL*NGLOB_AB)
    write(667,rec=1) zstore
    close(667)

    open(667,file='xtomo.bin',access='direct',recl=8*nx*ny*nz)
    write(667,rec=1) xstore_tomo_grid
    close(667)

    open(667,file='ytomo.bin',access='direct',recl=8*nx*ny*nz)
    write(667,rec=1) ystore_tomo_grid
    close(667)

    open(667,file='ztomo.bin',access='direct',recl=8*nx*ny*nz)
    write(667,rec=1) zstore_tomo_grid
    close(667)

    open(667,file='ibool.bin',access='direct',recl=4*NGLLX*NGLLY*NGLLZ*NSPEC_AB)
    write(667,rec=1) ibool
    close(667)

    open(667,file='indicei.bin',access='direct',recl=4*NGLLX*NGLLY*NGLLZ*NSPEC_AB)
    write(667,rec=1)  indice_grille_i
    close(667)

    open(667,file='indicej.bin',access='direct',recl=4*NGLLX*NGLLY*NGLLZ*NSPEC_AB)
    write(667,rec=1)  indice_grille_j
    close(667)

    open(667,file='indicek.bin',access='direct',recl=4*NGLLX*NGLLY*NGLLZ*NSPEC_AB)
    write(667,rec=1)  indice_grille_k
    close(667)

   endif
   call mpi_bcast(data_tomo,nx*ny*nz,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ier)
   call mpi_bcast(profondeur_tomo_grid,nx*ny*nz,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
   close(666)

end subroutine project_tomo_grid


subroutine interpolation_valeur(valeur,x,y,z,dat,xstore_local,ystore_local,zstore_local,&
                                ispec,NSPEC_AB,xigll,yigll,zigll)
   implicit none
   include "constants.h"
   integer NSPEC_AB,ispec
   real(kind=CUSTOM_REAL) dat(NGLLX,NGLLY,NGLLZ,NSPEC_AB)
   double precision xstore_local(NGNOD),ystore_local(NGNOD),zstore_local(NGNOD)
   double precision x,y,z,valeur
   double precision xi,eta,gamma,hlagrange
   double precision :: hxir(NGLLX),hetar(NGLLY),hgammar(NGLLZ)
   double precision :: hpxir(NGLLX),hpetar(NGLLY),hpgammar(NGLLZ)
   double precision, dimension(NGLLX) :: xigll!,wxgll
   double precision, dimension(NGLLY) :: yigll!,wygll
   double precision, dimension(NGLLZ) :: zigll!,wzgll
   integer igll,jgll,kgll

   call Find_xix_eta_gamma(xstore_local,ystore_local,zstore_local,xi,eta,gamma,x,y,z)


   call lagrange_any(xi,NGLLX,xigll,hxir,hpxir)
   call lagrange_any(eta,NGLLY,yigll,hetar,hpetar)
   call lagrange_any(gamma,NGLLZ,zigll,hgammar,hpgammar)

   valeur=0.d0

   do kgll = 1,NGLLZ
     do jgll = 1,NGLLY
       do igll = 1,NGLLX

          ! ploynome de lagrange
          hlagrange = hxir(igll)*hetar(jgll)*hgammar(kgll)

          ! interpolation de la valeur
          valeur = valeur + dble(dat(igll,jgll,kgll,ispec))*hlagrange

       enddo
     enddo
   enddo


end subroutine interpolation_valeur
!=============================================================
subroutine calcule_point_gll(xstore,ystore,zstore,xelm,yelm,zelm,shape3D)
  implicit none

  include "constants.h"
  double precision xmesh,ymesh,zmesh
  double precision shape3D(NGNOD,NGLLX,NGLLY,NGLLZ)
  double precision xelm(NGNOD),yelm(NGNOD),zelm(NGNOD)
  double precision xstore(NGLLX,NGLLY,NGLLZ)
  double precision ystore(NGLLX,NGLLY,NGLLZ)
  double precision zstore(NGLLX,NGLLY,NGLLZ)
  integer i,j,k,ia

  do k=1,NGLLZ
   do j=1,NGLLY
     do i=1,NGLLX

       xmesh = ZERO
       ymesh = ZERO
       zmesh = ZERO
       do ia=1,NGNOD
         xmesh = xmesh + shape3D(ia,i,j,k)*xelm(ia)
         ymesh = ymesh + shape3D(ia,i,j,k)*yelm(ia)
         zmesh = zmesh + shape3D(ia,i,j,k)*zelm(ia)
       enddo

       xstore(i,j,k) = xmesh
       ystore(i,j,k) = ymesh
       zstore(i,j,k) = zmesh

     enddo
   enddo
  enddo


end subroutine calcule_point_gll
!=============================================================
subroutine compute_volume_int(volume,xelm,yelm,zelm,wxgll,wygll,wzgll,shape3D,dershape3D,valeur_integration,index_cellule,NBCELL)
  implicit none

  include "constants.h"
  integer NBCELL,index_cellule

! Gauss-Lobatto-Legendre points of integration
  double precision wxgll(NGLLX)
  double precision wygll(NGLLY)
  double precision wzgll(NGLLZ)

! 3D shape functions and their derivatives
  double precision shape3D(NGNOD,NGLLX,NGLLY,NGLLZ)
  double precision dershape3D(NDIM,NGNOD,NGLLX,NGLLY,NGLLZ)
  real(kind=CUSTOM_REAL) valeur_integration(NGLLX,NGLLY,NGLLZ,NBCELL)
!
  double precision xelm(NGNOD), yelm(NGNOD), zelm(NGNOD)
!  double precision xstore_local(NGLLX,NGLLY,NGLLZ)  !!! faire passer ces tableaux en parametre et surtout changer leur nom.
!  double precision ystore_local(NGLLX,NGLLY,NGLLZ)
!  double precision zstore_local(NGLLX,NGLLY,NGLLZ)

  integer i,j,k,ia
  double precision xxi,xeta,xgamma,yxi,yeta,ygamma,zxi,zeta,zgamma
  double precision xmesh,ymesh,zmesh
  double precision xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz
  double precision jacobian,weight,volume
  integer jjj
  !  CALCUL DU VOLUME ------
 volume=0.d0

 do k=1,NGLLZ
    do j=1,NGLLY
      do i=1,NGLLX

      weight = wxgll(i)*wygll(j)*wzgll(k)

      xxi = ZERO
      xeta = ZERO
      xgamma = ZERO
      yxi = ZERO
      yeta = ZERO
      ygamma = ZERO
      zxi = ZERO
      zeta = ZERO
      zgamma = ZERO
      xmesh = ZERO
      ymesh = ZERO
      zmesh = ZERO

      do ia=1,NGNOD
        xxi = xxi + dershape3D(1,ia,i,j,k)*xelm(ia)
        xeta = xeta + dershape3D(2,ia,i,j,k)*xelm(ia)
        xgamma = xgamma + dershape3D(3,ia,i,j,k)*xelm(ia)
        yxi = yxi + dershape3D(1,ia,i,j,k)*yelm(ia)
        yeta = yeta + dershape3D(2,ia,i,j,k)*yelm(ia)
        ygamma = ygamma + dershape3D(3,ia,i,j,k)*yelm(ia)
        zxi = zxi + dershape3D(1,ia,i,j,k)*zelm(ia)
        zeta = zeta + dershape3D(2,ia,i,j,k)*zelm(ia)
        zgamma = zgamma + dershape3D(3,ia,i,j,k)*zelm(ia)
        xmesh = xmesh + shape3D(ia,i,j,k)*xelm(ia)
        ymesh = ymesh + shape3D(ia,i,j,k)*yelm(ia)
        zmesh = zmesh + shape3D(ia,i,j,k)*zelm(ia)
      enddo
      ! GLL points
!      xstore_local(i,j,k)=xmesh
!      ystore_local(i,j,k)=ymesh
!      zstore_local(i,j,k)=zmesh
      !
      jacobian = xxi*(yeta*zgamma-ygamma*zeta) - &
             xeta*(yxi*zgamma-ygamma*zxi) + &
             xgamma*(yxi*zeta-yeta*zxi)

! can ignore negative jacobian in mesher if needed when debugging code
      if(jacobian <= ZERO) then
         write(*,*) xelm
         write(*,*) yelm
         write(*,*) zelm
         call exit_MPI(0,'Warning 3D Jacobian undefined in compute_volume_int')
      endif
!     invert the relation (Fletcher p. 50 vol. 2)
      xix = (yeta*zgamma-ygamma*zeta) / jacobian
      xiy = (xgamma*zeta-xeta*zgamma) / jacobian
      xiz = (xeta*ygamma-xgamma*yeta) / jacobian
      etax = (ygamma*zxi-yxi*zgamma) / jacobian
      etay = (xxi*zgamma-xgamma*zxi) / jacobian
      etaz = (xgamma*yxi-xxi*ygamma) / jacobian
      gammax = (yxi*zeta-yeta*zxi) / jacobian
      gammay = (xeta*zxi-xxi*zeta) / jacobian
      gammaz = (xxi*yeta-xeta*yxi) / jacobian

!     compute and store the jacobian for the solver
      jacobian = 1. / (xix*(etay*gammaz-etaz*gammay) &
                      -xiy*(etax*gammaz-etaz*gammax) &
                      +xiz*(etax*gammay-etay*gammax))

      volume = volume + weight * jacobian * valeur_integration(i,j,k,index_cellule)
      !write(*,*) valeur_integration(i,j,k,index_cellule)
      !read(*,*) jjj
      enddo
    enddo
  enddo

end subroutine compute_volume_int

!=============================================================

  subroutine radius_comp(rbot,x0,y0,z0,p)
    implicit none
    double precision rbot,x,y,z,x0,y0,z0,p
    x=x0
    y=y0
    z=z0
    p = sqrt(x*x + y*y + (z+rbot)*(z+rbot))
  end subroutine radius_comp

  subroutine rayon_min_max(rbot,xmin,xmax,ymin,ymax,zmin,zmax,pmin,pmax)

   implicit none
   include 'constants.h'
   double precision  rbot,xmin,xmax,ymin,ymax,zmin,zmax,pmin,pmax
   double precision  x,y,z,p

   pmax=TINYVAL
   pmin=HUGEVAL

! 1
   x=xmin
   y=ymin
   z=zmin
   p = sqrt(x*x + y*y + (z+rbot)*(z+rbot))
   pmin=min(p,pmin)
   pmax=max(p,pmax)

! 2
   x=xmax
   y=ymin
   z=zmin
   p = sqrt(x*x + y*y + (z+rbot)*(z+rbot))
   pmin=min(p,pmin)
   pmax=max(p,pmax)

! 3
   x=xmax
   y=ymax
   z=zmin
   p = sqrt(x*x + y*y + (z+rbot)*(z+rbot))
   pmin=min(p,pmin)
   pmax=max(p,pmax)

! 4
   x=xmin
   y=ymax
   z=zmin
   p = sqrt(x*x + y*y + (z+rbot)*(z+rbot))
   pmin=min(p,pmin)
   pmax=max(p,pmax)

! 5
   x=xmin
   y=ymin
   z=zmax
   p = sqrt(x*x + y*y + (z+rbot)*(z+rbot))
   pmin=min(p,pmin)
   pmax=max(p,pmax)

! 6
   x=xmax
   y=ymin
   z=zmax
   p = sqrt(x*x + y*y + (z+rbot)*(z+rbot))
   pmin=min(p,pmin)
   pmax=max(p,pmax)

! 7
   x=xmax
   y=ymax
   z=zmax
   p = sqrt(x*x + y*y + (z+rbot)*(z+rbot))
   pmin=min(p,pmin)
   pmax=max(p,pmax)

! 8
   x=xmin
   y=ymax
   z=zmax
   p = sqrt(x*x + y*y + (z+rbot)*(z+rbot))
   pmin=min(p,pmin)
   pmax=max(p,pmax)



 end subroutine rayon_min_max





subroutine Find_xix_eta_gamma(xstore,ystore,zstore,xi,eta,gamma,x_target,y_target,z_target)
  implicit none
  double precision xstore(8),ystore(8),zstore(8),xi,eta,gamma,x,y,z
  double precision x_target,y_target,z_target
  double precision dist,distmin
  double precision xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz
  double precision dx,dy,dz,final_distance,dxi,deta,dgamma
  integer i,k0,iter

  ! initial geuss : sommet le plus proche
  distmin = 1.d20
  do i = 1, 8
     dist = (x_target - xstore(i))**2 + (y_target - ystore(i))**2 + (z_target - zstore(i))**2
     if (dist < distmin) then
        distmin = dist
        k0 = i
     endif
  enddo

  k0=7
  x=xstore(k0)
  y=ystore(k0)
  z=zstore(k0)
  xi=1.   ! xigll(ix_initial_guess(irec))  !!!!!!!!!!!!! pb ici
  eta=1.  ! yigll(iy_initial_guess(irec))
  gamma=1.! zigll(iz_initial_guess(irec))
  !
  if (myrank==0) then
     !write(*,*)
     !write(*,'(a13,2x,3f20.5)') 'first guess :',x,y,z
     !write(*,'(a7,2x,3f20.5)') 'target ',x_target,y_target,z_target
  endif
  do iter=1,6

     ! recompute jacobian for the new point
     call recompute_jacobian(xstore,ystore,zstore,xi,eta,gamma,x,y,z, &
          xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz)
     ! compute distance to target location
     dx = - (x - x_target)
     dy = - (y - y_target)
     dz = - (z - z_target)
     !if (myrank==0) write(*,'(6f20.5)') x,y,z,dx,dy,dz
     ! compute increments
     !
     dxi  = xix*dx + xiy*dy + xiz*dz
     deta = etax*dx + etay*dy + etaz*dz
     dgamma = gammax*dx + gammay*dy + gammaz*dz
     ! update values
     xi = xi + dxi
     eta = eta + deta
     gamma = gamma + dgamma


  enddo
  ! compute final coordinates of point found
  call recompute_jacobian(xstore,ystore,zstore,xi,eta,gamma,x,y,z, &
       xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz)

  ! compute final distance between asked and found (converted to km)
  final_distance = dsqrt((x_target-x)**2 + &
       (y_target-y)**2 + (z_target-z)**2)
  !if (dabs(xi)<1.d0.and.dabs(eta)<1.d0.and.dabs(gamma)<1.d0) then
    !write(*,*) 'dist :',final_distance,xi,eta,gamma
    !write(*,*) x/1000,y/1000,z/1000
    !write(*,*) xi,eta,gamma
    !write(*,*)
  !endif
end subroutine Find_xix_eta_gamma





 subroutine recompute_jacobian(xelm,yelm,zelm,xi,eta,gamma,x,y,z, &
                   xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz)

  implicit none

  !include "constants.h"
  integer, parameter :: NGNOD=8,NDIM=3
  double precision x,y,z,xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz
  double precision xi,eta,gamma,jacobian

! coordinates of the control points
  double precision xelm(NGNOD),yelm(NGNOD),zelm(NGNOD)

! 3D shape functions and their derivatives at receiver
  double precision shape3D(NGNOD)
  double precision dershape3D(NDIM,NGNOD)

  double precision xxi,yxi,zxi
  double precision xeta,yeta,zeta
  double precision xgamma,ygamma,zgamma
  double precision ra1,ra2,rb1,rb2,rc1,rc2

  integer ia

! for 8-node element
  double precision, parameter :: ONE_EIGHTH = 0.125d0
  double precision, parameter :: one=1.d0, zero=0.d0
! recompute jacobian for any (xi,eta,gamma) point, not necessarily a GLL point

! check that the parameter file is correct
  if(NGNOD /= 8) stop 'elements should have 8 control nodes'

! ***
! *** create the 3D shape functions and the Jacobian for an 8-node element
! ***

!--- case of an 8-node 3D element (Dhatt-Touzot p. 115)

  ra1 = one + xi
  ra2 = one - xi

  rb1 = one + eta
  rb2 = one - eta

  rc1 = one + gamma
  rc2 = one - gamma

  shape3D(1) = ONE_EIGHTH*ra2*rb2*rc2
  shape3D(2) = ONE_EIGHTH*ra1*rb2*rc2
  shape3D(3) = ONE_EIGHTH*ra1*rb1*rc2
  shape3D(4) = ONE_EIGHTH*ra2*rb1*rc2
  shape3D(5) = ONE_EIGHTH*ra2*rb2*rc1
  shape3D(6) = ONE_EIGHTH*ra1*rb2*rc1
  shape3D(7) = ONE_EIGHTH*ra1*rb1*rc1
  shape3D(8) = ONE_EIGHTH*ra2*rb1*rc1

  dershape3D(1,1) = - ONE_EIGHTH*rb2*rc2
  dershape3D(1,2) = ONE_EIGHTH*rb2*rc2
  dershape3D(1,3) = ONE_EIGHTH*rb1*rc2
  dershape3D(1,4) = - ONE_EIGHTH*rb1*rc2
  dershape3D(1,5) = - ONE_EIGHTH*rb2*rc1
  dershape3D(1,6) = ONE_EIGHTH*rb2*rc1
  dershape3D(1,7) = ONE_EIGHTH*rb1*rc1
  dershape3D(1,8) = - ONE_EIGHTH*rb1*rc1

  dershape3D(2,1) = - ONE_EIGHTH*ra2*rc2
  dershape3D(2,2) = - ONE_EIGHTH*ra1*rc2
  dershape3D(2,3) = ONE_EIGHTH*ra1*rc2
  dershape3D(2,4) = ONE_EIGHTH*ra2*rc2
  dershape3D(2,5) = - ONE_EIGHTH*ra2*rc1
  dershape3D(2,6) = - ONE_EIGHTH*ra1*rc1
  dershape3D(2,7) = ONE_EIGHTH*ra1*rc1
  dershape3D(2,8) = ONE_EIGHTH*ra2*rc1

  dershape3D(3,1) = - ONE_EIGHTH*ra2*rb2
  dershape3D(3,2) = - ONE_EIGHTH*ra1*rb2
  dershape3D(3,3) = - ONE_EIGHTH*ra1*rb1
  dershape3D(3,4) = - ONE_EIGHTH*ra2*rb1
  dershape3D(3,5) = ONE_EIGHTH*ra2*rb2
  dershape3D(3,6) = ONE_EIGHTH*ra1*rb2
  dershape3D(3,7) = ONE_EIGHTH*ra1*rb1
  dershape3D(3,8) = ONE_EIGHTH*ra2*rb1

! compute coordinates and jacobian matrix
  x=ZERO
  y=ZERO
  z=ZERO
  xxi=ZERO
  xeta=ZERO
  xgamma=ZERO
  yxi=ZERO
  yeta=ZERO
  ygamma=ZERO
  zxi=ZERO
  zeta=ZERO
  zgamma=ZERO

  do ia=1,NGNOD
    x=x+shape3D(ia)*xelm(ia)
    y=y+shape3D(ia)*yelm(ia)
    z=z+shape3D(ia)*zelm(ia)

    xxi=xxi+dershape3D(1,ia)*xelm(ia)
    xeta=xeta+dershape3D(2,ia)*xelm(ia)
    xgamma=xgamma+dershape3D(3,ia)*xelm(ia)
    yxi=yxi+dershape3D(1,ia)*yelm(ia)
    yeta=yeta+dershape3D(2,ia)*yelm(ia)
    ygamma=ygamma+dershape3D(3,ia)*yelm(ia)
    zxi=zxi+dershape3D(1,ia)*zelm(ia)
    zeta=zeta+dershape3D(2,ia)*zelm(ia)
    zgamma=zgamma+dershape3D(3,ia)*zelm(ia)
  enddo

  jacobian = xxi*(yeta*zgamma-ygamma*zeta) - xeta*(yxi*zgamma-ygamma*zxi) + xgamma*(yxi*zeta-yeta*zxi)

  if(jacobian <= ZERO) then
    !if (myrank==0) then
     !write(*,*) myrank
     !write(*,'(3f20.5)') x,y,z
     !write(*,'(3f20.5)') xi,eta,gamma
     !write(*,'(8f20.5)') xelm
     !write(*,'(8f20.5)') yelm
     !write(*,'(8f20.5)') zelm
    !endif
    stop '3D Jacobian undefined in recompute jacobian'
    !endif
  endif
! invert the relation (Fletcher p. 50 vol. 2)
  xix=(yeta*zgamma-ygamma*zeta)/jacobian
  xiy=(xgamma*zeta-xeta*zgamma)/jacobian
  xiz=(xeta*ygamma-xgamma*yeta)/jacobian
  etax=(ygamma*zxi-yxi*zgamma)/jacobian
  etay=(xxi*zgamma-xgamma*zxi)/jacobian
  etaz=(xgamma*yxi-xxi*ygamma)/jacobian
  gammax=(yxi*zeta-yeta*zxi)/jacobian
  gammay=(xeta*zxi-xxi*zeta)/jacobian
  gammaz=(xxi*yeta-xeta*yxi)/jacobian

  end subroutine recompute_jacobian

  subroutine inv_mapping_sph_cub(x0,y0,z0,ksi,eta,rbot)
    implicit none
    double precision x0,y0,z0,x,y,z,rbot
    double precision ksi,eta
    x=x0
    y=y0
    z=rbot+z0
    ksi=datan2(x,z)
    eta=datan2(y,z)
  end subroutine inv_mapping_sph_cub



  subroutine ksi_eta_min_max(rbot,xmin,xmax,ymin,ymax,zmin,zmax,ksimin,ksimax,etamin,etamax)
   implicit none
   include 'constants.h'
   double precision rbot,xmin,xmax,ymin,ymax,zmin,zmax,ksimin,ksimax,etamin,etamax
   double precision p,x,y,z

   ksimax=-HUGEVAL
   ksimin=HUGEVAL
   etamin=HUGEVAL
   etamax=-HUGEVAL

! 1
   x=xmin
   y=ymin
   z=zmin+rbot
   ksimin=min(datan2(x,z),ksimin)
   ksimax=max(datan2(x,z),ksimax)
   etamin=min(datan2(y,z),etamin)
   etamax=max(datan2(y,z),etamax)


! 2
   x=xmax
   y=ymin
   z=zmin+rbot
   ksimin=min(datan2(x,z),ksimin)
   ksimax=max(datan2(x,z),ksimax)
   etamin=min(datan2(y,z),etamin)
   etamax=max(datan2(y,z),etamax)


! 3
   x=xmax
   y=ymax
   z=zmin+rbot
   ksimin=min(datan2(x,z),ksimin)
   ksimax=max(datan2(x,z),ksimax)
   etamin=min(datan2(y,z),etamin)
   etamax=max(datan2(y,z),etamax)

! 4
   x=xmin
   y=ymax
   z=zmin+rbot
   ksimin=min(datan2(x,z),ksimin)
   ksimax=max(datan2(x,z),ksimax)
   etamin=min(datan2(y,z),etamin)
   etamax=max(datan2(y,z),etamax)


! 5
   x=xmin
   y=ymin
   z=zmax+rbot
   ksimin=min(datan2(x,z),ksimin)
   ksimax=max(datan2(x,z),ksimax)
   etamin=min(datan2(y,z),etamin)
   etamax=max(datan2(y,z),etamax)


! 6
   x=xmax
   y=ymin
   z=zmax+rbot
   ksimin=min(datan2(x,z),ksimin)
   ksimax=max(datan2(x,z),ksimax)
   etamin=min(datan2(y,z),etamin)
   etamax=max(datan2(y,z),etamax)


! 7
   x=xmax
   y=ymax
   z=zmax+rbot
   ksimin=min(datan2(x,z),ksimin)
   ksimax=max(datan2(x,z),ksimax)
   etamin=min(datan2(y,z),etamin)
   etamax=max(datan2(y,z),etamax)


! 8
   x=xmin
   y=ymax
   z=zmax+rbot
   ksimin=min(datan2(x,z),ksimin)
   ksimax=max(datan2(x,z),ksimax)
   etamin=min(datan2(y,z),etamin)
   etamax=max(datan2(y,z),etamax)




 end subroutine ksi_eta_min_max


 subroutine project_sem_interp(data_in_tomo_grid,data_in_sem_grid)
   implicit none
   real(kind=CUSTOM_REAL) data_in_tomo_grid(nx,ny,nz),data_in_sem_grid(NGLLX,NGLLY,NGLLZ,NSPEC_AB)
   integer ispec, i,j,k,i1,j1,k1,igll,jgll,kgll

   real(kind=CUSTOM_REAL) f1,f2,f3,f4,f5,f6,f7,f8
   real(kind=CUSTOM_REAL) x1,x2,x3,x4,x5,x6,x7,x8
   real(kind=CUSTOM_REAL) y1,y2,y3,y4,y5,y6,y7,y8
   real(kind=CUSTOM_REAL) z1,z2,z3,z4,z5,z6,z7,z8
   real(kind=CUSTOM_REAL) x,y,z,f

   data_glob(1:nx,1:ny,1:nz)=data_in_tomo_grid(:,:,:)

   data_glob(nx+1,1:ny,1:nz)=data_in_tomo_grid(nx,:,:)
   data_glob(1:nx,ny+1,1:nz)=data_in_tomo_grid(:,ny,:)
   data_glob(1:nx,1:ny,nz+1)=data_in_tomo_grid(:,:,nz)

   data_glob(nx+1,ny+1,1:nz)=data_in_tomo_grid(nx,ny,:)
   data_glob(nx+1,1:ny,nz+1)=data_in_tomo_grid(nx,:,nz)
   data_glob(1:nx,ny+1,nz+1)=data_in_tomo_grid(:,ny,nz)

   data_glob(nx+1,ny+1,nz+1)=data_in_tomo_grid(nx,ny,nz)


    do ispec = 1, NSPEC_AB
       do kgll = 1,NGLLZ
         do jgll = 1,NGLLY
           do igll = 1,NGLLX

                  i =  indice_grille_i(igll,jgll,kgll,ispec)
                  j =  indice_grille_j(igll,jgll,kgll,ispec)
                  k =  indice_grille_k(igll,jgll,kgll,ispec)

                  i1 =  indice_grille_i1(igll,jgll,kgll,ispec)
                  j1 =  indice_grille_j1(igll,jgll,kgll,ispec)
                  k1 =  indice_grille_k1(igll,jgll,kgll,ispec)

                  x= xstore(ibool(igll,jgll,kgll,ispec))
                  y= ystore(ibool(igll,jgll,kgll,ispec))
                  z= zstore(ibool(igll,jgll,kgll,ispec))

                  x1=xstore_tomo_grid_interp_tri(i,j,k) !xstore_tomo_grid(i,j,k)-0.5*hx
                  y1=ystore_tomo_grid_interp_tri(i,j,k) !ystore_tomo_grid(i,j,k)-0.5*hy
                  z1=zstore_tomo_grid_interp_tri(i,j,k) !zstore_tomo_grid(i,j,k)-0.5*hz

                  x2=xstore_tomo_grid_interp_tri(i1,j,k) !xstore_tomo_grid(i1,j,k)-0.5*hx
                  y2=ystore_tomo_grid_interp_tri(i1,j,k) !ystore_tomo_grid(i1,j,k)-0.5*hy
                  z2=zstore_tomo_grid_interp_tri(i1,j,k) !zstore_tomo_grid(i1,j,k)-0.5*hz

!!$                  if (x>=50000.) then
!!$                     write(*,*) i,i1
!!$                     write(*,*) x1,x,x2
!!$                     write(*,*) '---------------'
!!$                  endif
!!$                  if (x-x1 < 0. ) then
!!$                     write(*,*) x1,x,x2
!!$                     write(*,*) i,j,k
!!$                     write(*,*) i1,j1,k1
!!$                     write(*,*) '---------------'
!!$                  endif
                  x3=xstore_tomo_grid_interp_tri(i1,j1,k) !xstore_tomo_grid(i1,j1,k)-0.5*hx
                  y3=ystore_tomo_grid_interp_tri(i1,j1,k) !ystore_tomo_grid(i1,j1,k)-0.5*hy
                  z3=zstore_tomo_grid_interp_tri(i1,j1,k) !zstore_tomo_grid(i1,j1,k)-0.5*hz

                  x4=xstore_tomo_grid_interp_tri(i,j1,k) !xstore_tomo_grid(i,j1,k)-0.5*hx
                  y4=ystore_tomo_grid_interp_tri(i,j1,k) !ystore_tomo_grid(i,j1,k)-0.5*hy
                  z4=zstore_tomo_grid_interp_tri(i,j1,k) !zstore_tomo_grid(i,j1,k)-0.5*hz

                  x5=xstore_tomo_grid_interp_tri(i,j,k1) !xstore_tomo_grid(i,j,k1)-0.5*hx
                  y5=ystore_tomo_grid_interp_tri(i,j,k1) !ystore_tomo_grid(i,j,k1)-0.5*hy
                  z5=zstore_tomo_grid_interp_tri(i,j,k1) !zstore_tomo_grid(i,j,k1)-0.5*hz

                  x6=xstore_tomo_grid_interp_tri(i1,j,k1) !xstore_tomo_grid(i1,j,k1)-0.5*hx
                  y6=ystore_tomo_grid_interp_tri(i1,j,k1) !ystore_tomo_grid(i1,j,k1)-0.5*hy
                  z6=zstore_tomo_grid_interp_tri(i1,j,k1) !zstore_tomo_grid(i1,j,k1)-0.5*hz

                  x7=xstore_tomo_grid_interp_tri(i1,j1,k1) !xstore_tomo_grid(i1,j1,k1)-0.5*hx
                  y7=ystore_tomo_grid_interp_tri(i1,j1,k1) !ystore_tomo_grid(i1,j1,k1)-0.5*hy
                  z7=zstore_tomo_grid_interp_tri(i1,j1,k1) !zstore_tomo_grid(i1,j1,k1)-0.5*hz

                  x8=xstore_tomo_grid_interp_tri(i,j1,k1) !xstore_tomo_grid(i,j1,k1)-0.5*hx
                  y8=ystore_tomo_grid_interp_tri(i,j1,k1) !ystore_tomo_grid(i,j1,k1)-0.5*hy
                  z8=zstore_tomo_grid_interp_tri(i,j1,k1) !zstore_tomo_grid(i,j1,k1)-0.5*hz

                  f1=data_glob(i,j,k)    !data_in_tomo_grid(i,j,k)
                  f2=data_glob(i1,j,k)   !data_in_tomo_grid(i1,j,k)
                  f3=data_glob(i1,j1,k)  !data_in_tomo_grid(i1,j1,k)
                  f4=data_glob(i,j1,k)   !data_in_tomo_grid(i,j1,k)
                  f5=data_glob(i,j,k1)   !data_in_tomo_grid(i,j,k1)
                  f6=data_glob(i1,j,k1)  !data_in_tomo_grid(i1,j,k1)
                  f7=data_glob(i1,j1,k1) !data_in_tomo_grid(i1,j1,k1)
                  f8=data_glob(i,j1,k1)  !data_in_tomo_grid(i,j1,k1)
!!$                  if (myrank==1) then
!!$                     write(*,*) i,j,k
!!$                     write(*,*) i1,j1,k1
!!$                     write(*,*) x1,x2
!!$                     write(*,*) y1,y2
!!$                     write(*,*) z1,z2
!!$                     write(*,*) x,y,z
!!$                     write(*,*) '-------------'
!!$                  endif
                  ! interpolation tri-lineaire:
                  call interpoll3d_tri(&
                       x1,x2,x3,x4,x5,x6,x7,x8,&
                       y1,y2,y3,y4,y5,y6,y7,y8,&
                       z1,z2,z3,z4,z5,z6,z7,z8,&
                       f1,f2,f2,f4,f5,f6,f7,f8,&
                       x,y,z,f)

                  data_in_sem_grid(igll,jgll,kgll,ispec) = f

           enddo
        enddo
      enddo
    enddo


  end subroutine project_sem_interp

  subroutine project_sem(data_in_tomo_grid,data_in_sem_grid)
    implicit none
    real(kind=CUSTOM_REAL) data_in_tomo_grid(nx,ny,nz),data_in_sem_grid(NGLLX,NGLLY,NGLLZ,NSPEC_AB)
    integer ispec, i,j,k,i1,j1,k1,igll,jgll,kgll

    real(kind=CUSTOM_REAL) f1,f2,f3,f4,f5,f6,f7,f8
    real(kind=CUSTOM_REAL) x1,x2,x3,x4,x5,x6,x7,x8
    real(kind=CUSTOM_REAL) y1,y2,y3,y4,y5,y6,y7,y8
    real(kind=CUSTOM_REAL) z1,z2,z3,z4,z5,z6,z7,z8
    real(kind=CUSTOM_REAL) x,y,z,f

    do ispec = 1, NSPEC_AB
       do kgll = 1,NGLLZ
          do jgll = 1,NGLLY
             do igll = 1,NGLLX

                  i =  indice_grille_i(igll,jgll,kgll,ispec)
                  j =  indice_grille_j(igll,jgll,kgll,ispec)
                  k =  indice_grille_k(igll,jgll,kgll,ispec)

                  data_in_sem_grid(igll,jgll,kgll,ispec)=data_in_tomo_grid(i,j,k)
               enddo
            enddo
         enddo
      enddo
    end subroutine project_sem

 subroutine interpoll3d_tri(&
      x1,x2,x3,x4,x5,x6,x7,x8,&
      y1,y2,y3,y4,y5,y6,y7,y8,&
      z1,z2,z3,z4,z5,z6,z7,z8,&
      f1,f2,f3,f4,f5,f6,f7,f8,&
      x,y,z,f)
   !! todo utiliser les 8 sommets pour faire une interpolation de type gll
   implicit none
   real(kind=CUSTOM_REAL) f1,f2,f3,f4,f5,f6,f7,f8
   real(kind=CUSTOM_REAL) x1,x2,x3,x4,x5,x6,x7,x8
   real(kind=CUSTOM_REAL) y1,y2,y3,y4,y5,y6,y7,y8
   real(kind=CUSTOM_REAL) z1,z2,z3,z4,z5,z6,z7,z8
   real(kind=CUSTOM_REAL) x,y,z,f
   real(kind=CUSTOM_REAL) t,u,w

   if (x1==x2) then
      t=0.
   else
      t=(x-x1)/(x2-x1)
   endif

   if (y1==y3) then
      u=0.
   else
      u=(y-y1)/(y3-y1)
   endif

   if (z1==z5) then
      w=0.
   else
      w=(z-z1)/(z5-z1)
   endif

   !! check for debbuging
1000 format(a21,4f30.15)
   !if (t < -0.01) write(*,1000) 'warning interpol t<0',t,x1,x,x2
   !if (u < -0.01) write(*,1000) 'warning interpol u<0',u,y1,y,y3
   !if (w < -0.01) write(*,1000) 'warning interpol w<0',w,z1,z,z5
   !if (t >  1.01) write(*,1000) 'warning interpol t>1',t,x1,x,x2
   !if (u >  1.01) write(*,1000) 'warning interpol u>1',u,y1,y,y3
   !if (w >  1.01) write(*,1000) 'warning interpol w>1',w,z1,z,z5

   f= (1.d0-t) * (1.d0-u) * (1.d0-w)*  f1+&
         t     * (1.d0-u) * (1.d0-w)*  f2+&
         t     *    u     * (1.d0-w)*  f3+&
      (1.d0-t) *    u     * (1.d0-w)*  f4+&
      (1.d0-t) * (1.d0-u) *   w     *  f5+&
         t     * (1.d0-u) *   w     *  f6+&
         t     *   u      *   w     *  f7+&
      (1.d0-t) *   u      *   w     *  f8


 end subroutine interpoll3d_tri


 subroutine write_mod(model_to_write,iteration,name_data)

   implicit none
   real(kind=CUSTOM_REAL) model_to_write(nx,ny,nz)
   character(len=5) name_data
   character(len=100) name_file
   integer iteration

   write(name_file,'(a4,a5,a1,i5.5)') 'mod_',name_data(1:5),'.',iteration
   open(27,file=trim(name_file),access='direct',recl=CUSTOM_REAL*nx*ny*nz)
   write(27,rec=1) model_to_write
   close(27)

 end subroutine write_mod

 subroutine read_mod(model_to_write,iteration,name_data)

   implicit none
   real(kind=CUSTOM_REAL) model_to_write(nx,ny,nz)
   character(len=5) name_data
   character(len=100) name_file
   integer iteration

   write(name_file,'(a4,a5,a1,i5.5)') 'mod_',name_data(1:5),'.',iteration
   open(27,file=trim(name_file),access='direct',recl=CUSTOM_REAL*nx*ny*nz)
   read(27,rec=1) model_to_write
   close(27)

 end subroutine read_mod


 subroutine write_kl(model_to_write,i_source,iteration,name_data,nsrc)

    implicit none
    real(kind=CUSTOM_REAL) model_to_write(nx,ny,nz)
    character(len=5) name_data
    character(len=100) name_file
    integer i_source,iteration,nsrc

    !write(name_file,'(a4,a5,a1,i5.5,a1,i5.5)') 'grad_',name_data(1:5),'_',iteration,'_',i_source
    write(name_file,'(a4,a5,a1,i5.5)') 'grad',name_data(1:5),'.',iteration
    open(27,file=trim(name_file),access='direct',recl=CUSTOM_REAL*nx*ny*nz)
    write(27,rec=1) model_to_write(:,:,:)
    close(27)

 end subroutine write_kl


 subroutine read_kl(model_to_write,i_source,iteration,name_data)

   implicit none
   real(kind=CUSTOM_REAL) model_to_write(nx,ny,nz)
   character(len=5) name_data
   character(len=100) name_file
   integer i_source,iteration

   !write(name_file,'(a4,a5,a1,i5.5,a1,i5.5)') 'grad_',name_data(1:5),'_',iteration,'_',i_source
   write(name_file,'(a4,a5,a1,i5.5)') 'grad',name_data(1:5),'.',iteration
   open(27,file=trim(name_file),access='direct',recl=CUSTOM_REAL*nx*ny*nz)
   read(27,rec=1) model_to_write
   close(27)

 end subroutine read_kl







 subroutine  project_tomo_grid_4(dat_to_project1,dat_to_project2,dat_to_project3,dat_to_project4,ires)

   integer ispec,ires,ix,iy,iz
   integer i,j,k,igll,jgll,kgll
   integer imin,imax,jmin,jmax,kmin,kmax
   integer index_cellule,ierr,np
   double precision  ANGULAR_WIDTH_ETA_RAD,ANGULAR_WIDTH_XI_RAD,deg2rad
   double precision, parameter :: R_EARTH=6371000.d0
   double precision profondeur_courante
   double precision rx,ry,rz
   double precision x,y,z,x_centred,y_centred,z_centred
   double precision x0,x1,y0,y1,p0,p1,p,volume_cellule
   double precision volume1,volume2,volume3,volume4
   real(kind=CUSTOM_REAL)   dat_to_project1(NGLLX,NGLLY,NGLLZ,NSPEC_AB)
   real(kind=CUSTOM_REAL)   dat_to_project2(NGLLX,NGLLY,NGLLZ,NSPEC_AB)
   real(kind=CUSTOM_REAL)   dat_to_project3(NGLLX,NGLLY,NGLLZ,NSPEC_AB)
   real(kind=CUSTOM_REAL)   dat_to_project4(NGLLX,NGLLY,NGLLZ,NSPEC_AB)
   integer :: iii,jjj,kkk


   ! il faut definir dat (ie grad ou modele)
   !write(debug_file,'(a3,i4.4)') 'dbp',myrank
   !open(666,file=trim(debug_file))
   deg2rad = 3.141592653589793d0/180.d0

   ! initialization
   !! DEBUG
   !dat_to_project(:,:,:,:) = 10._CUSTOM_REAL
   !!
   data_tomo1(:,:,:)=0._CUSTOM_REAL
   data_tomo2(:,:,:)=0._CUSTOM_REAL
   data_tomo3(:,:,:)=0._CUSTOM_REAL
   data_tomo4(:,:,:)=0._CUSTOM_REAL

   vol_data(:,:,:)=0._CUSTOM_REAL
   indice(:,:,:)=0
   indice_glob(:,:,:)=0
   indice_int_glob(:,:,:,:)=0
   valeur_int_glob(:,:,:,:)=0._CUSTOM_REAL

   valeur_integration1(:,:,:,:) = 0._CUSTOM_REAL
   valeur_integration2(:,:,:,:) = 0._CUSTOM_REAL
   valeur_integration3(:,:,:,:) = 0._CUSTOM_REAL
   valeur_integration4(:,:,:,:) = 0._CUSTOM_REAL

   volume_integration(:,:,:,:) = 1._CUSTOM_REAL
   indice_integration(:,:,:,:) = 0
   indice_grille_i(:,:,:,:)=0
   indice_grille_j(:,:,:,:)=0
   indice_grille_k(:,:,:,:)=0

   xstore_tomo_grid(:,:,:)=0.d0
   ystore_tomo_grid(:,:,:)=0.d0
   zstore_tomo_grid(:,:,:)=0.d0
   profondeur_tomo_grid(:,:,:)=0.d0

!======== on se met dans la sphere cubique et on definit le pas des grilles

   ANGULAR_WIDTH_XI_RAD =  2.d0* dasin(dabs(xmin*1.00001)/R_EARTH)
   ANGULAR_WIDTH_ETA_RAD = 2.d0* dasin(dabs(ymin*1.00001)/R_EARTH)

   hx = ANGULAR_WIDTH_XI_RAD / (nx)
   hy = ANGULAR_WIDTH_ETA_RAD / (ny)
   hz = (r_earth-r_bot) / (nz)

   if (myrank==0) then
      !print *, 'Total number of points: ', np
      !print *, 'domain :'
      !print *, 'Xmin, Xmax :', xmin,xmax
      !print *, 'Ymin, Ymax :', ymin,ymax
      !print *, 'Zmin, Zmax :', zmin,zmax
      !print *, ' '
      !print *, 'angular width xi : ',ANGULAR_WIDTH_XI_RAD/deg2rad
      !print *, 'angular width eta : ',ANGULAR_WIDTH_ETA_RAD/deg2rad
      !print *, 'hx , hy, hz : ' , hx/deg2rad, hy/deg2rad, hz/1000.

     ! informations relatives a la grille
      open(10,file='grille.par')
      write(10,*)  nx,ny,nz
      write(10,*)  ANGULAR_WIDTH_XI_RAD, ANGULAR_WIDTH_ETA_RAD
      write(10,*)  r_earth,r_bot,hz
      write(10,*)
      close(10)
   endif

   if (myrank==0)   call cpu_time(start)
  !================================== 1/ interpolation grille tomographique  =========================================
    !write(666,*) 'nb el spectraux :' ,NSPEC_AB
   do ispec = 1, NSPEC_AB   !===== boucle element ispec

      xmin_local =  HUGEVAL
      xmax_local = -HUGEVAL
      ymin_local =  HUGEVAL
      ymax_local = -HUGEVAL
      zmin_local =  HUGEVAL
      zmax_local = -HUGEVAL

      !======================= Sommets de l'element ============
!1
      xstore_local(1)=xstore(ibool(1,1,1,ispec))
      ystore_local(1)=ystore(ibool(1,1,1,ispec))
      zstore_local(1)=zstore(ibool(1,1,1,ispec))
!2
      xstore_local(2)=xstore(ibool(NGLLX,1,1,ispec))
      ystore_local(2)=ystore(ibool(NGLLX,1,1,ispec))
      zstore_local(2)=zstore(ibool(NGLLX,1,1,ispec))
!3
      xstore_local(3)=xstore(ibool(NGLLX,NGLLY,1,ispec))
      ystore_local(3)=ystore(ibool(NGLLX,NGLLY,1,ispec))
      zstore_local(3)=zstore(ibool(NGLLX,NGLLY,1,ispec))
!4
      xstore_local(4)=xstore(ibool(1,NGLLY,1,ispec))
      ystore_local(4)=ystore(ibool(1,NGLLY,1,ispec))
      zstore_local(4)=zstore(ibool(1,NGLLY,1,ispec))
!5
      xstore_local(5)=xstore(ibool(1,1,NGLLZ,ispec))
      ystore_local(5)=ystore(ibool(1,1,NGLLZ,ispec))
      zstore_local(5)=zstore(ibool(1,1,NGLLZ,ispec))
!6
      xstore_local(6)=xstore(ibool(NGLLX,1,NGLLZ,ispec))
      ystore_local(6)=ystore(ibool(NGLLX,1,NGLLZ,ispec))
      zstore_local(6)=zstore(ibool(NGLLX,1,NGLLZ,ispec))
!7
      xstore_local(7)=xstore(ibool(NGLLX,NGLLY,NGLLZ,ispec))
      ystore_local(7)=ystore(ibool(NGLLX,NGLLY,NGLLZ,ispec))
      zstore_local(7)=zstore(ibool(NGLLX,NGLLY,NGLLZ,ispec))
!8
      xstore_local(8)=xstore(ibool(1,NGLLY,NGLLZ,ispec))
      ystore_local(8)=ystore(ibool(1,NGLLY,NGLLZ,ispec))
      zstore_local(8)=zstore(ibool(1,NGLLY,NGLLZ,ispec))
      !write(*,*) zstore_local
      !=============== on cherche le pave circonscrit a l'element en se basant sur les point GLL
      do iz=1,NGLLZ
         do iy=1,NGLLY
            do ix=1,NGLLX
               xmin_local=min(xmin_local,xstore(ibool(ix,iy,iz,ispec)))
               xmax_local=max(xmin_local,xstore(ibool(ix,iy,iz,ispec)))
               ymin_local=min(ymin_local,ystore(ibool(ix,iy,iz,ispec)))
               ymax_local=max(ymin_local,ystore(ibool(ix,iy,iz,ispec)))
               zmin_local=min(zmin_local,zstore(ibool(ix,iy,iz,ispec)))
               zmax_local=max(zmax_local,zstore(ibool(ix,iy,iz,ispec)))
            enddo
         enddo
      enddo

       ! =========== calcul du rayon min et max correspondant au pave
      call rayon_min_max(R_EARTH-zmax,xmin_local,xmax_local,ymin_local,ymax_local,zmin_local,zmax_local,pmin,pmax)

      ! indice sur grille fine d'integration (sphere cubique) !! ajout de +/-2 pour deborder sur l'element adjacent
      kmin = 1+ floor((pmin-r_bot)/(r_earth-r_bot)*(nnz-1) )  !- 2
      kmax = 1+ floor((pmax-r_bot)/(r_earth-r_bot)*(nnz-1) )  !+ 2

      call ksi_eta_min_max(R_EARTH-zmax,xmin_local,xmax_local,ymin_local,ymax_local,zmin_local,zmax_local,&
           ksimin,ksimax,etamin,etamax)

      ! indice sur grille fine d'integration (sphere cubique)
      imin = floor(1. + (NNX - 1) * (ksimin + 0.5d0 *  ANGULAR_WIDTH_XI_RAD)  / ANGULAR_WIDTH_XI_RAD) !- 2
      jmin = floor(1. + (NNY - 1) * (etamin + 0.5d0 * ANGULAR_WIDTH_ETA_RAD) / ANGULAR_WIDTH_ETA_RAD) !- 2

      ! indice sur grille fine d'integration (sphere cubique)
      imax = floor(1. + (NNX - 1) * (ksimax + 0.5d0 *  ANGULAR_WIDTH_XI_RAD)  / ANGULAR_WIDTH_XI_RAD) !+ 2
      jmax = floor(1. + (NNY - 1) * (etamax + 0.5d0 * ANGULAR_WIDTH_ETA_RAD) / ANGULAR_WIDTH_ETA_RAD) !+ 2
      !write(666,*) ksimax ,  ANGULAR_WIDTH_XI_RAD , ((ksimax + 0.5d0 * ANGULAR_WIDTH_XI_RAD) / ANGULAR_WIDTH_XI_RAD)

      imin=max(imin,1)
      imax=min(imax,nx)
      jmin=max(jmin,1)
      jmax=min(jmax,ny)
      kmin=max(kmin,1)
      kmax=min(kmax,nz)

      !================calcul de la valeur du noyau au centre de chaque cellule de la grille tomo

      do k=kmin,kmax

         z_centred = r_bot + (k-0.5d0)*hz  ! valeur au centre de la cellule (sphere cubique)
         profondeur_courante = R_EARTH  - z_centred
         !write(*,*) R_EARTH, z_centred, profondeur_courante
         do j=jmin,jmax

            ratio_eta = (dble(j)-0.5d0) / dble(NY)
            y_centred = 2.d0*ratio_eta-1
            y_centred = tan((ANGULAR_WIDTH_ETA_RAD/2.d0) * y_centred)


            do i=imin,imax

               index_cellule = i + (j-1)*(nnx-1) + (k-1)*(nnx-1)*(nny-1)

               ratio_xi = (dble(i)-0.5d0) / dble(NX)
               x_centred = 2.d0*ratio_xi-1
               x_centred = tan((ANGULAR_WIDTH_XI_RAD/2.d0) * x_centred)


               ! on repasse en cartesien pour travailler dans le repere SEM
               px = x_centred*z_centred
               py = y_centred*z_centred
               pz =  -(R_EARTH - z_centred/dsqrt(1.d0 + y_centred**2 + x_centred**2) - zmax)

               ! on cherche le xix, eta et gamma du point central qui est px,py,pz
               call Find_xix_eta_gamma(xstore_local,ystore_local,zstore_local,xi,eta,gamma,px,py,pz)

               ! on regarde si on est bien dans l'element ispec pour faire la correspondance SEM <=> tomo
               if (dabs(xi)<=1.05d0.and.dabs(eta)<=1.05d0.and.dabs(gamma)<=1.05d0) then


                  indice(i,j,k)=1 ! on a bien visite la cellule de la grille
                  xstore_tomo_grid(i,j,k)=px
                  ystore_tomo_grid(i,j,k)=py
                  zstore_tomo_grid(i,j,k)=pz
                  profondeur_tomo_grid(i,j,k)=profondeur_courante

               else

               endif

               ! == calcul des 8 sommets de la cellule d'integration =====================
               ! 1
               rx=1.d0;ry=1.d0;rz=1.d0

               pz = r_bot + (dble(k)-rz)*hz

               ratio_eta = (dble(j)-ry) / dble(NY)
               y = 2.d0*ratio_eta-1;y = tan((ANGULAR_WIDTH_ETA_RAD/2.d0) * y)

               ratio_xi = (dble(i)-rx) / dble(NX)
               x = 2.d0*ratio_xi-1;x = tan((ANGULAR_WIDTH_XI_RAD/2.d0) * x)

               px = x*pz; py = y*pz
               pz =  -(R_EARTH - pz/dsqrt(1.d0 + y*y + x*x) - zmax)



               xstore_local_int(1) = px
               ystore_local_int(1) = py
               zstore_local_int(1) = pz

            ! 2
               rx=0.d0;ry=1.d0;rz=1.d0

               pz = r_bot + (dble(k)-rz)*hz

               ratio_eta = (dble(j)-ry) / dble(NY)
               y = 2.d0*ratio_eta-1;y = tan((ANGULAR_WIDTH_ETA_RAD/2.d0) * y)

               ratio_xi = (dble(i)-rx) / dble(NX)
               x = 2.d0*ratio_xi-1;x = tan((ANGULAR_WIDTH_XI_RAD/2.d0) * x)

               px = x*pz; py = y*pz
               pz =  -(R_EARTH - pz/dsqrt(1.d0 + y*y + x*x) - zmax)



               xstore_local_int(2) = px
               ystore_local_int(2) = py
               zstore_local_int(2) = pz

               ! 3
               rx=0.d0;ry=0.d0;rz=1.d0

               pz = r_bot + (dble(k)-rz)*hz

               ratio_eta = (dble(j)-ry) / dble(NY)
               y = 2.d0*ratio_eta-1;y = tan((ANGULAR_WIDTH_ETA_RAD/2.d0) * y)

               ratio_xi = (dble(i)-rx) / dble(NX)
               x = 2.d0*ratio_xi-1;x = tan((ANGULAR_WIDTH_XI_RAD/2.d0) * x)

               px = x*pz; py = y*pz
               pz =  -(R_EARTH - pz/dsqrt(1.d0 + y*y + x*x) - zmax)



               xstore_local_int(3) = px
               ystore_local_int(3) = py
               zstore_local_int(3) = pz


               ! 4
               rx=1.d0;ry=0.d0;rz=1.d0

               pz = r_bot + (dble(k)-rz)*hz

               ratio_eta = (dble(j)-ry) / dble(NY)
               y = 2.d0*ratio_eta-1;y = tan((ANGULAR_WIDTH_ETA_RAD/2.d0) * y)

               ratio_xi = (dble(i)-rx) / dble(NX)
               x = 2.d0*ratio_xi-1;x = tan((ANGULAR_WIDTH_XI_RAD/2.d0) * x)

               px = x*pz; py = y*pz
               pz =  -(R_EARTH - pz/dsqrt(1.d0 + y*y + x*x) - zmax)



               xstore_local_int(4) = px
               ystore_local_int(4) = py
               zstore_local_int(4) = pz

               ! 5
               rx=1.d0;ry=1.d0;rz=0.d0

               pz = r_bot + (dble(k)-rz)*hz

               ratio_eta = (dble(j)-ry) / dble(NY)
               y = 2.d0*ratio_eta-1;y = tan((ANGULAR_WIDTH_ETA_RAD/2.d0) * y)

               ratio_xi = (dble(i)-rx) / dble(NX)
               x = 2.d0*ratio_xi-1;x = tan((ANGULAR_WIDTH_XI_RAD/2.d0) * x)

               px = x*pz; py = y*pz
               pz =  -(R_EARTH - pz/dsqrt(1.d0 + y*y + x*x) - zmax)



               xstore_local_int(5) = px
               ystore_local_int(5) = py
               zstore_local_int(5) = pz



               ! 6
               rx=0.d0;ry=1.d0;rz=0.d0

               pz = r_bot + (dble(k)-rz)*hz

               ratio_eta = (dble(j)-ry) / dble(NY)
               y = 2.d0*ratio_eta-1;y = tan((ANGULAR_WIDTH_ETA_RAD/2.d0) * y)

               ratio_xi = (dble(i)-rx) / dble(NX)
               x = 2.d0*ratio_xi-1;x = tan((ANGULAR_WIDTH_XI_RAD/2.d0) * x)

               px = x*pz; py = y*pz
               pz =  -(R_EARTH - pz/dsqrt(1.d0 + y*y + x*x) - zmax)



               xstore_local_int(6) = px
               ystore_local_int(6) = py
               zstore_local_int(6) = pz

               ! 7
               rx=0.d0;ry=0.d0;rz=0.d0

               pz = r_bot + (dble(k)-rz)*hz

               ratio_eta = (dble(j)-ry) / dble(NY)
               y = 2.d0*ratio_eta-1;y = tan((ANGULAR_WIDTH_ETA_RAD/2.d0) * y)

               ratio_xi = (dble(i)-rx) / dble(NX)
               x = 2.d0*ratio_xi-1;x = tan((ANGULAR_WIDTH_XI_RAD/2.d0) * x)

               px = x*pz; py = y*pz
               pz =  -(R_EARTH - pz/dsqrt(1.d0 + y*y + x*x) - zmax)



               xstore_local_int(7) = px
               ystore_local_int(7) = py
               zstore_local_int(7) = pz


               ! 8
               rx=1.d0;ry=0.d0;rz=0.d0

               pz = r_bot + (dble(k)-rz)*hz

               ratio_eta = (dble(j)-ry) / dble(NY)
               y = 2.d0*ratio_eta-1;y = tan((ANGULAR_WIDTH_ETA_RAD/2.d0) * y)

               ratio_xi = (dble(i)-rx) / dble(NX)
               x = 2.d0*ratio_xi-1;x = tan((ANGULAR_WIDTH_XI_RAD/2.d0) * x)

               px = x*pz; py = y*pz
               pz =  -(R_EARTH - pz/dsqrt(1.d0 + y*y + x*x) - zmax)



               xstore_local_int(8) = px
               ystore_local_int(8) = py
               zstore_local_int(8) = pz

               ! ========== calcul des points GLL
               call calcule_point_gll(xgll,ygll,zgll,xstore_local_int,ystore_local_int,zstore_local_int,shape3D)

               ! ========== interpolation de la fonciton dat aux points GLL
               do kgll = 1,NGLLZ
                  do jgll = 1,NGLLY
                     do igll = 1,NGLLX



                        px = xgll(igll,jgll,kgll)
                        py = ygll(igll,jgll,kgll)
                        pz = zgll(igll,jgll,kgll)

                        call Find_xix_eta_gamma(xstore_local,ystore_local,zstore_local,xi,eta,gamma,px,py,pz)

                        if (dabs(xi)<=1.05d0.and.dabs(eta)<=1.05d0.and.dabs(gamma)<=1.05d0) then

                           call interpolation_valeur(valeur,px,py,pz,dat_to_project1,xstore_local,ystore_local,zstore_local,&
                                ispec,NSPEC_AB,xigll,yigll,zigll)
                           valeur_integration1(igll,jgll,kgll,index_cellule)=valeur


                           call interpolation_valeur(valeur,px,py,pz,dat_to_project2,xstore_local,ystore_local,zstore_local,&
                                ispec,NSPEC_AB,xigll,yigll,zigll)
                           valeur_integration2(igll,jgll,kgll,index_cellule)=valeur

                           call interpolation_valeur(valeur,px,py,pz,dat_to_project3,xstore_local,ystore_local,zstore_local,&
                                ispec,NSPEC_AB,xigll,yigll,zigll)
                           valeur_integration3(igll,jgll,kgll,index_cellule)=valeur


                           call interpolation_valeur(valeur,px,py,pz,dat_to_project4,xstore_local,ystore_local,zstore_local,&
                                ispec,NSPEC_AB,xigll,yigll,zigll)
                           valeur_integration4(igll,jgll,kgll,index_cellule)=valeur

                           indice_integration(igll,jgll,kgll,index_cellule)=1

                        else


                        endif

                     enddo
                  enddo
               enddo


               ! ========================= on cherche les points GLL qui sont dans la cellule de la grille tomo
               call calcule_point_gll(xgll,ygll,zgll,xstore_local,ystore_local,zstore_local,shape3D)
               do kgll = 1,NGLLZ
                  do jgll = 1,NGLLY
                     do igll = 1,NGLLX



                        px = xgll(igll,jgll,kgll)
                        py = ygll(igll,jgll,kgll)
                        pz = zgll(igll,jgll,kgll)
                        call Find_xix_eta_gamma(xstore_local_int,ystore_local_int,zstore_local_int,xi,eta,gamma,px,py,pz)

                        if (dabs(xi)<=1.05d0.and.dabs(eta)<=1.05d0.and.dabs(gamma)<=1.05d0) then

                           indice_grille_i(igll,jgll,kgll,ispec) = i
                           indice_grille_i1(igll,jgll,kgll,ispec)= i

                           indice_grille_j(igll,jgll,kgll,ispec) = j
                           indice_grille_j1(igll,jgll,kgll,ispec)= j

                           indice_grille_k(igll,jgll,kgll,ispec) = k
                           indice_grille_k1(igll,jgll,kgll,ispec)= k

                           if (xi > 0. ) then
                              if (i<nx) indice_grille_i1(igll,jgll,kgll,ispec) = i+1
                           else
                              if (i>1)  indice_grille_i1(igll,jgll,kgll,ispec) = i-1
                           endif

                           if (eta > 0. ) then
                              if (j<ny) indice_grille_j1(igll,jgll,kgll,ispec) = j+1
                           else
                              if (j>1)  indice_grille_j1(igll,jgll,kgll,ispec) = j-1
                           endif

                           if (gamma > 0. ) then
                              if (k<nz) indice_grille_k1(igll,jgll,kgll,ispec) = k+1
                           else
                              if (k>1)  indice_grille_k1(igll,jgll,kgll,ispec) = k-1
                           endif

                        endif

                     enddo
                  enddo
               enddo

               ! ===========================================================


            enddo
         enddo
      enddo

   enddo !==== boucle ispec


   if (myrank==0) then
      call cpu_time(finish)
      write(*,*) 'time // proj 0: ', finish-start
      call cpu_time(start)
   endif
   !! === verif des correspondances (gll,ispec) -> tomo
   do ispec=1,nspec_ab
      do kgll = 1,NGLLZ
         do jgll = 1,NGLLY
            do igll = 1,NGLLX

               !! patch ad hoc pour effet de bords
               if (indice_grille_i(igll,jgll,kgll,ispec)==0 .and. igll==1) then
                  indice_grille_i(igll,jgll,kgll,ispec)=indice_grille_i(igll+1,jgll,kgll,ispec)
                  indice_grille_j(igll,jgll,kgll,ispec)=indice_grille_j(igll+1,jgll,kgll,ispec)
                  indice_grille_k(igll,jgll,kgll,ispec)=indice_grille_k(igll+1,jgll,kgll,ispec)
               endif

               if (indice_grille_i(igll,jgll,kgll,ispec)==0 .and. igll==5) then
                  indice_grille_i(igll,jgll,kgll,ispec)=indice_grille_i(igll-1,jgll,kgll,ispec)
                  indice_grille_j(igll,jgll,kgll,ispec)=indice_grille_j(igll-1,jgll,kgll,ispec)
                  indice_grille_k(igll,jgll,kgll,ispec)=indice_grille_k(igll-1,jgll,kgll,ispec)
               endif


               if (indice_grille_i(igll,jgll,kgll,ispec)==0 .and. jgll==1) then
                  indice_grille_i(igll,jgll,kgll,ispec)=indice_grille_i(igll,jgll+1,kgll,ispec)
                  indice_grille_j(igll,jgll,kgll,ispec)=indice_grille_j(igll,jgll+1,kgll,ispec)
                  indice_grille_k(igll,jgll,kgll,ispec)=indice_grille_k(igll,jgll+1,kgll,ispec)
               endif


               if (indice_grille_i(igll,jgll,kgll,ispec)==0 .and. jgll==5) then
                  indice_grille_i(igll,jgll,kgll,ispec)=indice_grille_i(igll,jgll-1,kgll,ispec)
                  indice_grille_j(igll,jgll,kgll,ispec)=indice_grille_j(igll,jgll-1,kgll,ispec)
                  indice_grille_k(igll,jgll,kgll,ispec)=indice_grille_k(igll,jgll-1,kgll,ispec)
               endif

               if (indice_grille_i(igll,jgll,kgll,ispec)==0 .and. kgll==1) then
                  indice_grille_i(igll,jgll,kgll,ispec)=indice_grille_i(igll,jgll,kgll+1,ispec)
                  indice_grille_j(igll,jgll,kgll,ispec)=indice_grille_j(igll,jgll,kgll+1,ispec)
                  indice_grille_k(igll,jgll,kgll,ispec)=indice_grille_k(igll,jgll,kgll+1,ispec)
               endif


               if (indice_grille_i(igll,jgll,kgll,ispec)==0 .and. kgll==5) then
                  indice_grille_i(igll,jgll,kgll,ispec)=indice_grille_i(igll,jgll,kgll-1,ispec)
                  indice_grille_j(igll,jgll,kgll,ispec)=indice_grille_j(igll,jgll,kgll-1,ispec)
                  indice_grille_k(igll,jgll,kgll,ispec)=indice_grille_k(igll,jgll,kgll-1,ispec)
               endif


               if (indice_grille_i(igll,jgll,kgll,ispec)==0 .and. jgll==1 .and. kgll==1) then
                  indice_grille_i(igll,jgll,kgll,ispec)=indice_grille_i(igll,jgll+1,kgll+1,ispec)
                  indice_grille_j(igll,jgll,kgll,ispec)=indice_grille_j(igll,jgll+1,kgll+1,ispec)
                  indice_grille_k(igll,jgll,kgll,ispec)=indice_grille_k(igll,jgll+1,kgll+1,ispec)
               endif


               if (indice_grille_i(igll,jgll,kgll,ispec)==0 .and. jgll==5 .and. kgll==5) then
                  indice_grille_i(igll,jgll,kgll,ispec)=indice_grille_i(igll,jgll-1,kgll-1,ispec)
                  indice_grille_j(igll,jgll,kgll,ispec)=indice_grille_j(igll,jgll-1,kgll-1,ispec)
                  indice_grille_k(igll,jgll,kgll,ispec)=indice_grille_k(igll,jgll-1,kgll-1,ispec)
               endif

               if (indice_grille_i(igll,jgll,kgll,ispec)==0 .and. igll==1 .and. jgll==1) then
                  indice_grille_i(igll,jgll,kgll,ispec)=indice_grille_i(igll+1,jgll+1,kgll,ispec)
                  indice_grille_j(igll,jgll,kgll,ispec)=indice_grille_j(igll+1,jgll+1,kgll,ispec)
                  indice_grille_k(igll,jgll,kgll,ispec)=indice_grille_k(igll+1,jgll+1,kgll,ispec)
               endif


               if (indice_grille_i(igll,jgll,kgll,ispec)==0 .and. igll==5 .and. jgll==5) then
                  indice_grille_i(igll,jgll,kgll,ispec)=indice_grille_i(igll-1,jgll-1,kgll,ispec)
                  indice_grille_j(igll,jgll,kgll,ispec)=indice_grille_j(igll-1,jgll-1,kgll,ispec)
                  indice_grille_k(igll,jgll,kgll,ispec)=indice_grille_k(igll-1,jgll-1,kgll,ispec)
               endif



            enddo
         enddo
      enddo
   enddo
   if (myrank==0) then
      call cpu_time(finish)
      write(*,*) 'time // proj : ', finish-start
   endif

   valeur_int_glob(:,:,:,:)=0._CUSTOM_REAL
   call mpi_reduce(valeur_integration1,valeur_int_glob,NGLLX*NGLLY*NGLLZ*(nnx-1)*(nny-1)*(nnz-1),CUSTOM_MPI_TYPE, MPI_SUM,0,MPI_COMM_WORLD,ier)
   valeur_integration1(:,:,:,:) = valeur_int_glob(:,:,:,:)

   valeur_int_glob(:,:,:,:)=0._CUSTOM_REAL
   call mpi_reduce(valeur_integration2,valeur_int_glob,NGLLX*NGLLY*NGLLZ*(nnx-1)*(nny-1)*(nnz-1),CUSTOM_MPI_TYPE, MPI_SUM,0,MPI_COMM_WORLD,ier)
   valeur_integration2(:,:,:,:) = valeur_int_glob(:,:,:,:)

   valeur_int_glob(:,:,:,:)=0._CUSTOM_REAL
   call mpi_reduce(valeur_integration3,valeur_int_glob,NGLLX*NGLLY*NGLLZ*(nnx-1)*(nny-1)*(nnz-1),CUSTOM_MPI_TYPE, MPI_SUM,0,MPI_COMM_WORLD,ier)
   valeur_integration3(:,:,:,:) = valeur_int_glob(:,:,:,:)

   valeur_int_glob(:,:,:,:)=0._CUSTOM_REAL
   call mpi_reduce(valeur_integration4,valeur_int_glob,NGLLX*NGLLY*NGLLZ*(nnx-1)*(nny-1)*(nnz-1),CUSTOM_MPI_TYPE, MPI_SUM,0,MPI_COMM_WORLD,ier)
   valeur_integration4(:,:,:,:) = valeur_int_glob(:,:,:,:)

   call mpi_reduce(indice_integration,indice_int_glob,NGLLX*NGLLY*NGLLZ*(nnx-1)*(nny-1)*(nnz-1),MPI_INTEGER, MPI_SUM,0,MPI_COMM_WORLD,ier)
   call mpi_reduce(indice,indice_glob,(nnx-1)*(nny-1)*(nnz-1),MPI_INTEGER ,MPI_SUM,0, MPI_COMM_WORLD, ierr)

    call  mpi_reduce(xstore_tomo_grid,work_array,(nnx-1)*(nny-1)*(nnz-1),MPI_DOUBLE_PRECISION ,MPI_SUM,0, MPI_COMM_WORLD, ierr)
    xstore_tomo_grid(:,:,:)=work_array(:,:,:)

    call  mpi_reduce(ystore_tomo_grid,work_array,(nnx-1)*(nny-1)*(nnz-1),MPI_DOUBLE_PRECISION ,MPI_SUM,0, MPI_COMM_WORLD, ierr)
    ystore_tomo_grid(:,:,:)=work_array(:,:,:)

    call  mpi_reduce(zstore_tomo_grid,work_array,(nnx-1)*(nny-1)*(nnz-1),MPI_DOUBLE_PRECISION ,MPI_SUM,0, MPI_COMM_WORLD, ierr)
    zstore_tomo_grid(:,:,:)=work_array(:,:,:)

    call  mpi_reduce( profondeur_tomo_grid,work_array,(nnx-1)*(nny-1)*(nnz-1),MPI_DOUBLE_PRECISION ,MPI_SUM,0, MPI_COMM_WORLD, ierr)
    profondeur_tomo_grid(:,:,:)=work_array(:,:,:)

    indice(:,:,:) = indice_glob(:,:,:)
    indice_integration(:,:,:,:) = indice_int_glob(:,:,:,:)


!================== 2/ projection sur grille tomo =====================================

    if (myrank==0) then  !======================== myrank=0
       call  cpu_time(start)
       index_cellule=0
       do k=1,nnz-1
          do j=1,nny-1
             do i=1,nnx-1
                index_cellule=index_cellule+1

                if (indice(i,j,k)/=0)  then
                   xstore_tomo_grid(i,j,k) =  xstore_tomo_grid(i,j,k) / real(indice(i,j,k),8)
                   ystore_tomo_grid(i,j,k) =  ystore_tomo_grid(i,j,k) / real(indice(i,j,k),8)
                   zstore_tomo_grid(i,j,k) =  zstore_tomo_grid(i,j,k) / real(indice(i,j,k),8)
                   profondeur_tomo_grid(i,j,k) = profondeur_tomo_grid(i,j,k) / real(indice(i,j,k),8)
                else
                   write(*,*) 'point oublie :'  ,i,j,k
                endif


                !! on reagarde si les points on ete visites et combien de fois?
                do kgll=1,NGLLZ
                   do jgll=1,NGLLY
                      do igll=1,NGLLX

                         if (indice_integration(igll,jgll,kgll,index_cellule) == 0) then  !! le point n'a pas ete visite

                            !! chercher le point le plus proche qui est a ete visite
                            !  premier patch : devrait marcher pour mon cas cartesien
                            if (kgll==1 .and. indice_integration(igll,jgll,kgll+1,index_cellule) > 0) then

                               valeur_integration1(igll,jgll,kgll,index_cellule) = valeur_integration1(igll,jgll,kgll+1,index_cellule) /&
                                    real(indice_integration(igll,jgll,kgll+1,index_cellule),CUSTOM_REAL)
                               valeur_integration2(igll,jgll,kgll,index_cellule) = valeur_integration2(igll,jgll,kgll+1,index_cellule) /&
                                    real(indice_integration(igll,jgll,kgll+1,index_cellule),CUSTOM_REAL)
                               valeur_integration3(igll,jgll,kgll,index_cellule) = valeur_integration3(igll,jgll,kgll+1,index_cellule) /&
                                    real(indice_integration(igll,jgll,kgll+1,index_cellule),CUSTOM_REAL)
                               valeur_integration4(igll,jgll,kgll,index_cellule) = valeur_integration4(igll,jgll,kgll+1,index_cellule) /&
                                    real(indice_integration(igll,jgll,kgll+1,index_cellule),CUSTOM_REAL)

                            endif
                         else  ! le point a ete visite "indice_integration" fois
!!$                            ! on divise par indice pour avoir la valeur moyenne car on a fait un MPI_SUM

                            valeur_integration1(igll,jgll,kgll,index_cellule) =  valeur_integration1(igll,jgll,kgll,index_cellule)/&
                                 real(indice_integration(igll,jgll,kgll,index_cellule),CUSTOM_REAL)
                            valeur_integration2(igll,jgll,kgll,index_cellule) =  valeur_integration2(igll,jgll,kgll,index_cellule)/&
                                 real(indice_integration(igll,jgll,kgll,index_cellule),CUSTOM_REAL)
                            valeur_integration3(igll,jgll,kgll,index_cellule) =  valeur_integration3(igll,jgll,kgll,index_cellule)/&
                                 real(indice_integration(igll,jgll,kgll,index_cellule),CUSTOM_REAL)
                            valeur_integration4(igll,jgll,kgll,index_cellule) =  valeur_integration4(igll,jgll,kgll,index_cellule)/&
                                 real(indice_integration(igll,jgll,kgll,index_cellule),CUSTOM_REAL)

                         endif

                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo


       !===== boucle grille tomo : verifier qu'on fait exactement la meme boucle que precedemment
       np=0
       do k=1,nnz-1  ! grille integration (sphere cubique)


          p0= r_bot + (k-1)*hz
          p1= r_bot + k*hz

          do j=2,nny ! grille integration (sphere cubique)


             ratio_eta = (dble(j-2)) / dble(NY)
             y0 = 2.d0*ratio_eta-1
             y0 = tan((ANGULAR_WIDTH_ETA_RAD/2.d0) * y0)

             ratio_eta = (dble(j-1)) / dble(NY)
             y1 = 2.d0*ratio_eta-1
             y1 = tan((ANGULAR_WIDTH_ETA_RAD/2.d0) * y1)


             do i=2,nnx ! grille integration (sphere cubique)

                index_cellule=i-1+(j-2)*(nnx-1)+(k-1)*(nnx-1)*(nny-1)

                ratio_xi = (dble(i-2)) / dble(NX)
                x0 = 2.d0*ratio_xi-1
                x0 = tan((ANGULAR_WIDTH_XI_RAD/2.d0) * x0)

                ratio_xi = (dble(i-1)) / dble(NX)
                x1 = 2.d0*ratio_xi-1
                x1 = tan((ANGULAR_WIDTH_XI_RAD/2.d0) * x1)

                !===== calcul des 8 sommets de la cellule d'integration dans le repere cartesien SEM

                ! 1
                x=x0;y=y0;p=p0
                xstore_local(1) = x*p
                ystore_local(1) = y*p
                zstore_local(1) =  -(R_EARTH - p/dsqrt(1.d0 + y*y + x*x) - zmax)

                ! 2
                x=x1;y=y0;p=p0
                xstore_local(2) = x*p
                ystore_local(2) = y*p
                zstore_local(2) =  -(R_EARTH - p/dsqrt(1.d0 + y*y + x*x) - zmax)

                ! 3
                x=x1;y=y1;p=p0
                xstore_local(3) = x*p
                ystore_local(3) = y*p
                zstore_local(3) =  -(R_EARTH - p/dsqrt(1.d0 + y*y + x*x) - zmax)

                ! 4
                x=x0;y=y1;p=p0
                xstore_local(4) = x*p
                ystore_local(4) = y*p
                zstore_local(4) =  -(R_EARTH - p/dsqrt(1.d0 + y*y + x*x) - zmax)

                ! 5
                x=x0;y=y0;p=p1
                xstore_local(5) = x*p
                ystore_local(5) = y*p
                zstore_local(5) =  -(R_EARTH - p/dsqrt(1.d0 + y*y + x*x) - zmax)

                ! 6
                x=x1;y=y0;p=p1
                xstore_local(6) = x*p
                ystore_local(6) = y*p
                zstore_local(6) =  -(R_EARTH - p/dsqrt(1.d0 + y*y + x*x) - zmax)

                ! 7
                x=x1;y=y1;p=p1
                xstore_local(7) = x*p
                ystore_local(7) = y*p
                zstore_local(7) =  -(R_EARTH - p/dsqrt(1.d0 + y*y + x*x) - zmax)


                ! 8
                x=x0;y=y1;p=p1
                xstore_local(8) = x*p
                ystore_local(8) = y*p
                zstore_local(8) =  -(R_EARTH - p/dsqrt(1.d0 + y*y + x*x) - zmax)

                ! === calcul de l'integrale de volume
                call compute_volume_int(volume1,xstore_local,ystore_local,zstore_local,wxgll,wygll,wzgll,&
                     shape3D,dershape3D,valeur_integration1,index_cellule,(nnx-1)*(nny-1)*(nnz-1))

                call compute_volume_int(volume2,xstore_local,ystore_local,zstore_local,wxgll,wygll,wzgll,&
                     shape3D,dershape3D,valeur_integration2,index_cellule,(nnx-1)*(nny-1)*(nnz-1))

                call compute_volume_int(volume3,xstore_local,ystore_local,zstore_local,wxgll,wygll,wzgll,&
                     shape3D,dershape3D,valeur_integration3,index_cellule,(nnx-1)*(nny-1)*(nnz-1))

                call compute_volume_int(volume4,xstore_local,ystore_local,zstore_local,wxgll,wygll,wzgll,&
                     shape3D,dershape3D,valeur_integration4,index_cellule,(nnx-1)*(nny-1)*(nnz-1))

                if (ires==1)  then
                   ! === si on veut faire juste la projection du modele, on ne doit pas normaliser par
                   ! === le volume de la cellule, on prend plutot la valeur moyenne.
                   call compute_volume_int(volume_cellule,xstore_local,ystore_local,zstore_local,wxgll,wygll,wzgll,&
                        shape3D,dershape3D,volume_integration,index_cellule,(nnx-1)*(nny-1)*(nnz-1))
                   volume1=volume1/volume_cellule
                   volume2=volume2/volume_cellule
                   volume3=volume3/volume_cellule
                   volume4=volume4/volume_cellule
                endif

                if (indice(i-1,j-1,k)==0) then
                   ! on a oublie le point i,j,k
                   np=np+1
                   !write(30,*) i,j,k
                else

                   data_tomo1(i-1,j-1,k) = volume1
                   data_tomo2(i-1,j-1,k) = volume2
                   data_tomo3(i-1,j-1,k) = volume3
                   data_tomo4(i-1,j-1,k) = volume4

                endif

             enddo
          enddo

       enddo
       call cpu_time(finish)
       write(*,*) 'time seq. ',finish-start
       !! DEBUG

       !!
    endif
    call mpi_bcast(data_tomo1,nx*ny*nz,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ier)
    call mpi_bcast(data_tomo2,nx*ny*nz,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ier)
    call mpi_bcast(data_tomo3,nx*ny*nz,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ier)
    call mpi_bcast(data_tomo4,nx*ny*nz,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ier)
    call mpi_bcast(profondeur_tomo_grid,nx*ny*nz,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)


  end subroutine project_tomo_grid_4









 subroutine  project_tomo_grid_2(dat_to_project1,dat_to_project2,ires)

   integer ispec,ires,ix,iy,iz
   integer i,j,k,igll,jgll,kgll
   integer imin,imax,jmin,jmax,kmin,kmax
   integer index_cellule,ierr,np
   double precision  ANGULAR_WIDTH_ETA_RAD,ANGULAR_WIDTH_XI_RAD,deg2rad
   double precision, parameter :: R_EARTH=6371000.d0
   double precision profondeur_courante
   double precision rx,ry,rz
   double precision x,y,z,x_centred,y_centred,z_centred
   double precision x0,x1,y0,y1,p0,p1,p,volume_cellule
   double precision volume1,volume2,volume3,volume4
   real(kind=CUSTOM_REAL)   dat_to_project1(NGLLX,NGLLY,NGLLZ,NSPEC_AB)
   real(kind=CUSTOM_REAL)   dat_to_project2(NGLLX,NGLLY,NGLLZ,NSPEC_AB)
   !real(kind=CUSTOM_REAL)   dat_to_project3(NGLLX,NGLLY,NGLLZ,NSPEC_AB)
   !real(kind=CUSTOM_REAL)   dat_to_project4(NGLLX,NGLLY,NGLLZ,NSPEC_AB)
   integer :: iii,jjj,kkk


   ! il faut definir dat (ie grad ou modele)
   !write(debug_file,'(a3,i4.4)') 'dbp',myrank
   !open(666,file=trim(debug_file))
   deg2rad = 3.141592653589793d0/180.d0

   ! initialization
   !! DEBUG
   !dat_to_project(:,:,:,:) = 10._CUSTOM_REAL
   !!
   data_tomo1(:,:,:)=0._CUSTOM_REAL
   data_tomo2(:,:,:)=0._CUSTOM_REAL
   data_tomo3(:,:,:)=0._CUSTOM_REAL
   data_tomo4(:,:,:)=0._CUSTOM_REAL

   vol_data(:,:,:)=0._CUSTOM_REAL
   indice(:,:,:)=0
   indice_glob(:,:,:)=0
   indice_int_glob(:,:,:,:)=0
   valeur_int_glob(:,:,:,:)=0._CUSTOM_REAL

   valeur_integration1(:,:,:,:) = 0._CUSTOM_REAL
   valeur_integration2(:,:,:,:) = 0._CUSTOM_REAL
   valeur_integration3(:,:,:,:) = 0._CUSTOM_REAL
   valeur_integration4(:,:,:,:) = 0._CUSTOM_REAL

   volume_integration(:,:,:,:) = 1._CUSTOM_REAL
   indice_integration(:,:,:,:) = 0
   indice_grille_i(:,:,:,:)=0
   indice_grille_j(:,:,:,:)=0
   indice_grille_k(:,:,:,:)=0

   xstore_tomo_grid(:,:,:)=0.d0
   ystore_tomo_grid(:,:,:)=0.d0
   zstore_tomo_grid(:,:,:)=0.d0
   profondeur_tomo_grid(:,:,:)=0.d0

!======== on se met dans la sphere cubique et on definit le pas des grilles

   ANGULAR_WIDTH_XI_RAD =  2.d0* dasin(dabs(xmin*1.00001)/R_EARTH)
   ANGULAR_WIDTH_ETA_RAD = 2.d0* dasin(dabs(ymin*1.00001)/R_EARTH)

   hx = ANGULAR_WIDTH_XI_RAD / (nx)
   hy = ANGULAR_WIDTH_ETA_RAD / (ny)
   hz = (r_earth-r_bot) / (nz)

   if (myrank==0) then
      !print *, 'Total number of points: ', np
      !print *, 'domain :'
      !print *, 'Xmin, Xmax :', xmin,xmax
      !print *, 'Ymin, Ymax :', ymin,ymax
      !print *, 'Zmin, Zmax :', zmin,zmax
      !print *, ' '
      !print *, 'angular width xi : ',ANGULAR_WIDTH_XI_RAD/deg2rad
      !print *, 'angular width eta : ',ANGULAR_WIDTH_ETA_RAD/deg2rad
      !print *, 'hx , hy, hz : ' , hx/deg2rad, hy/deg2rad, hz/1000.

     ! informations relatives a la grille
      open(10,file='grille.par')
      write(10,*)  nx,ny,nz
      write(10,*)  ANGULAR_WIDTH_XI_RAD, ANGULAR_WIDTH_ETA_RAD
      write(10,*)  r_earth,r_bot,hz
      write(10,*)
      close(10)
   endif

   if (myrank==0)   call cpu_time(start)
  !================================== 1/ interpolation grille tomographique  =========================================
    !write(666,*) 'nb el spectraux :' ,NSPEC_AB
   do ispec = 1, NSPEC_AB   !===== boucle element ispec

      xmin_local =  HUGEVAL
      xmax_local = -HUGEVAL
      ymin_local =  HUGEVAL
      ymax_local = -HUGEVAL
      zmin_local =  HUGEVAL
      zmax_local = -HUGEVAL

      !======================= Sommets de l'element ============
!1
      xstore_local(1)=xstore(ibool(1,1,1,ispec))
      ystore_local(1)=ystore(ibool(1,1,1,ispec))
      zstore_local(1)=zstore(ibool(1,1,1,ispec))
!2
      xstore_local(2)=xstore(ibool(NGLLX,1,1,ispec))
      ystore_local(2)=ystore(ibool(NGLLX,1,1,ispec))
      zstore_local(2)=zstore(ibool(NGLLX,1,1,ispec))
!3
      xstore_local(3)=xstore(ibool(NGLLX,NGLLY,1,ispec))
      ystore_local(3)=ystore(ibool(NGLLX,NGLLY,1,ispec))
      zstore_local(3)=zstore(ibool(NGLLX,NGLLY,1,ispec))
!4
      xstore_local(4)=xstore(ibool(1,NGLLY,1,ispec))
      ystore_local(4)=ystore(ibool(1,NGLLY,1,ispec))
      zstore_local(4)=zstore(ibool(1,NGLLY,1,ispec))
!5
      xstore_local(5)=xstore(ibool(1,1,NGLLZ,ispec))
      ystore_local(5)=ystore(ibool(1,1,NGLLZ,ispec))
      zstore_local(5)=zstore(ibool(1,1,NGLLZ,ispec))
!6
      xstore_local(6)=xstore(ibool(NGLLX,1,NGLLZ,ispec))
      ystore_local(6)=ystore(ibool(NGLLX,1,NGLLZ,ispec))
      zstore_local(6)=zstore(ibool(NGLLX,1,NGLLZ,ispec))
!7
      xstore_local(7)=xstore(ibool(NGLLX,NGLLY,NGLLZ,ispec))
      ystore_local(7)=ystore(ibool(NGLLX,NGLLY,NGLLZ,ispec))
      zstore_local(7)=zstore(ibool(NGLLX,NGLLY,NGLLZ,ispec))
!8
      xstore_local(8)=xstore(ibool(1,NGLLY,NGLLZ,ispec))
      ystore_local(8)=ystore(ibool(1,NGLLY,NGLLZ,ispec))
      zstore_local(8)=zstore(ibool(1,NGLLY,NGLLZ,ispec))
      !write(*,*) zstore_local
      !=============== on cherche le pave circonscrit a l'element en se basant sur les point GLL
      do iz=1,NGLLZ
         do iy=1,NGLLY
            do ix=1,NGLLX
               xmin_local=min(xmin_local,xstore(ibool(ix,iy,iz,ispec)))
               xmax_local=max(xmin_local,xstore(ibool(ix,iy,iz,ispec)))
               ymin_local=min(ymin_local,ystore(ibool(ix,iy,iz,ispec)))
               ymax_local=max(ymin_local,ystore(ibool(ix,iy,iz,ispec)))
               zmin_local=min(zmin_local,zstore(ibool(ix,iy,iz,ispec)))
               zmax_local=max(zmax_local,zstore(ibool(ix,iy,iz,ispec)))
            enddo
         enddo
      enddo

       ! =========== calcul du rayon min et max correspondant au pave
      call rayon_min_max(R_EARTH-zmax,xmin_local,xmax_local,ymin_local,ymax_local,zmin_local,zmax_local,pmin,pmax)

      ! indice sur grille fine d'integration (sphere cubique) !! ajout de +/-2 pour deborder sur l'element adjacent
      kmin = 1+ floor((pmin-r_bot)/(r_earth-r_bot)*(nnz-1) )  !- 2
      kmax = 1+ floor((pmax-r_bot)/(r_earth-r_bot)*(nnz-1) )  !+ 2

      call ksi_eta_min_max(R_EARTH-zmax,xmin_local,xmax_local,ymin_local,ymax_local,zmin_local,zmax_local,&
           ksimin,ksimax,etamin,etamax)

      ! indice sur grille fine d'integration (sphere cubique)
      imin = floor(1. + (NNX - 1) * (ksimin + 0.5d0 *  ANGULAR_WIDTH_XI_RAD)  / ANGULAR_WIDTH_XI_RAD) !- 2
      jmin = floor(1. + (NNY - 1) * (etamin + 0.5d0 * ANGULAR_WIDTH_ETA_RAD) / ANGULAR_WIDTH_ETA_RAD) !- 2

      ! indice sur grille fine d'integration (sphere cubique)
      imax = floor(1. + (NNX - 1) * (ksimax + 0.5d0 *  ANGULAR_WIDTH_XI_RAD)  / ANGULAR_WIDTH_XI_RAD) !+ 2
      jmax = floor(1. + (NNY - 1) * (etamax + 0.5d0 * ANGULAR_WIDTH_ETA_RAD) / ANGULAR_WIDTH_ETA_RAD) !+ 2
      !write(666,*) ksimax ,  ANGULAR_WIDTH_XI_RAD , ((ksimax + 0.5d0 * ANGULAR_WIDTH_XI_RAD) / ANGULAR_WIDTH_XI_RAD)

      imin=max(imin,1)
      imax=min(imax,nx)
      jmin=max(jmin,1)
      jmax=min(jmax,ny)
      kmin=max(kmin,1)
      kmax=min(kmax,nz)

      !================calcul de la valeur du noyau au centre de chaque cellule de la grille tomo

      do k=kmin,kmax

         z_centred = r_bot + (k-0.5d0)*hz  ! valeur au centre de la cellule (sphere cubique)
         profondeur_courante = R_EARTH  - z_centred
         !write(*,*) R_EARTH, z_centred, profondeur_courante
         do j=jmin,jmax

            ratio_eta = (dble(j)-0.5d0) / dble(NY)
            y_centred = 2.d0*ratio_eta-1
            y_centred = tan((ANGULAR_WIDTH_ETA_RAD/2.d0) * y_centred)


            do i=imin,imax

               index_cellule = i + (j-1)*(nnx-1) + (k-1)*(nnx-1)*(nny-1)

               ratio_xi = (dble(i)-0.5d0) / dble(NX)
               x_centred = 2.d0*ratio_xi-1
               x_centred = tan((ANGULAR_WIDTH_XI_RAD/2.d0) * x_centred)


               ! on repasse en cartesien pour travailler dans le repere SEM
               px = x_centred*z_centred
               py = y_centred*z_centred
               pz =  -(R_EARTH - z_centred/dsqrt(1.d0 + y_centred**2 + x_centred**2) - zmax)

               ! on cherche le xix, eta et gamma du point central qui est px,py,pz
               call Find_xix_eta_gamma(xstore_local,ystore_local,zstore_local,xi,eta,gamma,px,py,pz)

               ! on regarde si on est bien dans l'element ispec pour faire la correspondance SEM <=> tomo
               if (dabs(xi)<=1.05d0.and.dabs(eta)<=1.05d0.and.dabs(gamma)<=1.05d0) then


                  indice(i,j,k)=1 ! on a bien visite la cellule de la grille
                  xstore_tomo_grid(i,j,k)=px
                  ystore_tomo_grid(i,j,k)=py
                  zstore_tomo_grid(i,j,k)=pz
                  profondeur_tomo_grid(i,j,k)=profondeur_courante

               else

               endif

               ! == calcul des 8 sommets de la cellule d'integration =====================
               ! 1
               rx=1.d0;ry=1.d0;rz=1.d0

               pz = r_bot + (dble(k)-rz)*hz

               ratio_eta = (dble(j)-ry) / dble(NY)
               y = 2.d0*ratio_eta-1;y = tan((ANGULAR_WIDTH_ETA_RAD/2.d0) * y)

               ratio_xi = (dble(i)-rx) / dble(NX)
               x = 2.d0*ratio_xi-1;x = tan((ANGULAR_WIDTH_XI_RAD/2.d0) * x)

               px = x*pz; py = y*pz
               pz =  -(R_EARTH - pz/dsqrt(1.d0 + y*y + x*x) - zmax)



               xstore_local_int(1) = px
               ystore_local_int(1) = py
               zstore_local_int(1) = pz

               ! 2
               rx=0.d0;ry=1.d0;rz=1.d0

               pz = r_bot + (dble(k)-rz)*hz

               ratio_eta = (dble(j)-ry) / dble(NY)
               y = 2.d0*ratio_eta-1;y = tan((ANGULAR_WIDTH_ETA_RAD/2.d0) * y)

               ratio_xi = (dble(i)-rx) / dble(NX)
               x = 2.d0*ratio_xi-1;x = tan((ANGULAR_WIDTH_XI_RAD/2.d0) * x)

               px = x*pz; py = y*pz
               pz =  -(R_EARTH - pz/dsqrt(1.d0 + y*y + x*x) - zmax)



               xstore_local_int(2) = px
               ystore_local_int(2) = py
               zstore_local_int(2) = pz

               ! 3
               rx=0.d0;ry=0.d0;rz=1.d0

               pz = r_bot + (dble(k)-rz)*hz

               ratio_eta = (dble(j)-ry) / dble(NY)
               y = 2.d0*ratio_eta-1;y = tan((ANGULAR_WIDTH_ETA_RAD/2.d0) * y)

               ratio_xi = (dble(i)-rx) / dble(NX)
               x = 2.d0*ratio_xi-1;x = tan((ANGULAR_WIDTH_XI_RAD/2.d0) * x)

               px = x*pz; py = y*pz
               pz =  -(R_EARTH - pz/dsqrt(1.d0 + y*y + x*x) - zmax)



               xstore_local_int(3) = px
               ystore_local_int(3) = py
               zstore_local_int(3) = pz


               ! 4
               rx=1.d0;ry=0.d0;rz=1.d0

               pz = r_bot + (dble(k)-rz)*hz

               ratio_eta = (dble(j)-ry) / dble(NY)
               y = 2.d0*ratio_eta-1;y = tan((ANGULAR_WIDTH_ETA_RAD/2.d0) * y)

               ratio_xi = (dble(i)-rx) / dble(NX)
               x = 2.d0*ratio_xi-1;x = tan((ANGULAR_WIDTH_XI_RAD/2.d0) * x)

               px = x*pz; py = y*pz
               pz =  -(R_EARTH - pz/dsqrt(1.d0 + y*y + x*x) - zmax)



               xstore_local_int(4) = px
               ystore_local_int(4) = py
               zstore_local_int(4) = pz

               ! 5
               rx=1.d0;ry=1.d0;rz=0.d0

               pz = r_bot + (dble(k)-rz)*hz

               ratio_eta = (dble(j)-ry) / dble(NY)
               y = 2.d0*ratio_eta-1;y = tan((ANGULAR_WIDTH_ETA_RAD/2.d0) * y)

               ratio_xi = (dble(i)-rx) / dble(NX)
               x = 2.d0*ratio_xi-1;x = tan((ANGULAR_WIDTH_XI_RAD/2.d0) * x)

               px = x*pz; py = y*pz
               pz =  -(R_EARTH - pz/dsqrt(1.d0 + y*y + x*x) - zmax)



               xstore_local_int(5) = px
               ystore_local_int(5) = py
               zstore_local_int(5) = pz



               ! 6
               rx=0.d0;ry=1.d0;rz=0.d0

               pz = r_bot + (dble(k)-rz)*hz

               ratio_eta = (dble(j)-ry) / dble(NY)
               y = 2.d0*ratio_eta-1;y = tan((ANGULAR_WIDTH_ETA_RAD/2.d0) * y)

               ratio_xi = (dble(i)-rx) / dble(NX)
               x = 2.d0*ratio_xi-1;x = tan((ANGULAR_WIDTH_XI_RAD/2.d0) * x)

               px = x*pz; py = y*pz
               pz =  -(R_EARTH - pz/dsqrt(1.d0 + y*y + x*x) - zmax)



               xstore_local_int(6) = px
               ystore_local_int(6) = py
               zstore_local_int(6) = pz

               ! 7
               rx=0.d0;ry=0.d0;rz=0.d0

               pz = r_bot + (dble(k)-rz)*hz

               ratio_eta = (dble(j)-ry) / dble(NY)
               y = 2.d0*ratio_eta-1;y = tan((ANGULAR_WIDTH_ETA_RAD/2.d0) * y)

               ratio_xi = (dble(i)-rx) / dble(NX)
               x = 2.d0*ratio_xi-1;x = tan((ANGULAR_WIDTH_XI_RAD/2.d0) * x)

               px = x*pz; py = y*pz
               pz =  -(R_EARTH - pz/dsqrt(1.d0 + y*y + x*x) - zmax)



               xstore_local_int(7) = px
               ystore_local_int(7) = py
               zstore_local_int(7) = pz


               ! 8
               rx=1.d0;ry=0.d0;rz=0.d0

               pz = r_bot + (dble(k)-rz)*hz

               ratio_eta = (dble(j)-ry) / dble(NY)
               y = 2.d0*ratio_eta-1;y = tan((ANGULAR_WIDTH_ETA_RAD/2.d0) * y)

               ratio_xi = (dble(i)-rx) / dble(NX)
               x = 2.d0*ratio_xi-1;x = tan((ANGULAR_WIDTH_XI_RAD/2.d0) * x)

               px = x*pz; py = y*pz
               pz =  -(R_EARTH - pz/dsqrt(1.d0 + y*y + x*x) - zmax)



               xstore_local_int(8) = px
               ystore_local_int(8) = py
               zstore_local_int(8) = pz

               ! ========== calcul des points GLL
               call calcule_point_gll(xgll,ygll,zgll,xstore_local_int,ystore_local_int,zstore_local_int,shape3D)

               ! ========== interpolation de la fonciton dat aux points GLL
               do kgll = 1,NGLLZ
                  do jgll = 1,NGLLY
                     do igll = 1,NGLLX



                        px = xgll(igll,jgll,kgll)
                        py = ygll(igll,jgll,kgll)
                        pz = zgll(igll,jgll,kgll)

                        call Find_xix_eta_gamma(xstore_local,ystore_local,zstore_local,xi,eta,gamma,px,py,pz)

                        if (dabs(xi)<=1.05d0.and.dabs(eta)<=1.05d0.and.dabs(gamma)<=1.05d0) then

                           call interpolation_valeur(valeur,px,py,pz,dat_to_project1,xstore_local,ystore_local,zstore_local,&
                                ispec,NSPEC_AB,xigll,yigll,zigll)
                           valeur_integration1(igll,jgll,kgll,index_cellule)=valeur


                           call interpolation_valeur(valeur,px,py,pz,dat_to_project2,xstore_local,ystore_local,zstore_local,&
                                ispec,NSPEC_AB,xigll,yigll,zigll)
                           valeur_integration2(igll,jgll,kgll,index_cellule)=valeur

                           !call interpolation_valeur(valeur,px,py,pz,dat_to_project3,xstore_local,ystore_local,zstore_local,&
                           !     ispec,NSPEC_AB,xigll,yigll,zigll)
                           !valeur_integration3(igll,jgll,kgll,index_cellule)=valeur


                           !call interpolation_valeur(valeur,px,py,pz,dat_to_project4,xstore_local,ystore_local,zstore_local,&
                           !     ispec,NSPEC_AB,xigll,yigll,zigll)
                           !valeur_integration4(igll,jgll,kgll,index_cellule)=valeur

                           indice_integration(igll,jgll,kgll,index_cellule)=1

                        else


                        endif

                     enddo
                  enddo
               enddo


               ! ========================= on cherche les points GLL qui sont dans la cellule de la grille tomo
               call calcule_point_gll(xgll,ygll,zgll,xstore_local,ystore_local,zstore_local,shape3D)
               do kgll = 1,NGLLZ
                  do jgll = 1,NGLLY
                     do igll = 1,NGLLX



                        px = xgll(igll,jgll,kgll)
                        py = ygll(igll,jgll,kgll)
                        pz = zgll(igll,jgll,kgll)
                        call Find_xix_eta_gamma(xstore_local_int,ystore_local_int,zstore_local_int,xi,eta,gamma,px,py,pz)

                        if (dabs(xi)<=1.05d0.and.dabs(eta)<=1.05d0.and.dabs(gamma)<=1.05d0) then

                           indice_grille_i(igll,jgll,kgll,ispec) = i
                           indice_grille_i1(igll,jgll,kgll,ispec)= i

                           indice_grille_j(igll,jgll,kgll,ispec) = j
                           indice_grille_j1(igll,jgll,kgll,ispec)= j

                           indice_grille_k(igll,jgll,kgll,ispec) = k
                           indice_grille_k1(igll,jgll,kgll,ispec)= k

                           if (xi > 0. ) then
                              if (i<nx) indice_grille_i1(igll,jgll,kgll,ispec) = i+1
                           else
                              if (i>1)  indice_grille_i1(igll,jgll,kgll,ispec) = i-1
                           endif

                           if (eta > 0. ) then
                              if (j<ny) indice_grille_j1(igll,jgll,kgll,ispec) = j+1
                           else
                              if (j>1)  indice_grille_j1(igll,jgll,kgll,ispec) = j-1
                           endif

                           if (gamma > 0. ) then
                              if (k<nz) indice_grille_k1(igll,jgll,kgll,ispec) = k+1
                           else
                              if (k>1)  indice_grille_k1(igll,jgll,kgll,ispec) = k-1
                           endif

                        endif

                     enddo
                  enddo
               enddo

               ! ===========================================================


            enddo
         enddo
      enddo

   enddo !==== boucle ispec


   if (myrank==0) then
      call cpu_time(finish)
      write(*,*) 'time // proj 0: ', finish-start
      call cpu_time(start)
   endif
   !! === verif des correspondances (gll,ispec) -> tomo
   do ispec=1,nspec_ab
      do kgll = 1,NGLLZ
         do jgll = 1,NGLLY
            do igll = 1,NGLLX

               !! patch ad hoc pour effet de bords
               if (indice_grille_i(igll,jgll,kgll,ispec)==0 .and. igll==1) then
                  indice_grille_i(igll,jgll,kgll,ispec)=indice_grille_i(igll+1,jgll,kgll,ispec)
                  indice_grille_j(igll,jgll,kgll,ispec)=indice_grille_j(igll+1,jgll,kgll,ispec)
                  indice_grille_k(igll,jgll,kgll,ispec)=indice_grille_k(igll+1,jgll,kgll,ispec)
               endif

               if (indice_grille_i(igll,jgll,kgll,ispec)==0 .and. igll==5) then
                  indice_grille_i(igll,jgll,kgll,ispec)=indice_grille_i(igll-1,jgll,kgll,ispec)
                  indice_grille_j(igll,jgll,kgll,ispec)=indice_grille_j(igll-1,jgll,kgll,ispec)
                  indice_grille_k(igll,jgll,kgll,ispec)=indice_grille_k(igll-1,jgll,kgll,ispec)
               endif


               if (indice_grille_i(igll,jgll,kgll,ispec)==0 .and. jgll==1) then
                  indice_grille_i(igll,jgll,kgll,ispec)=indice_grille_i(igll,jgll+1,kgll,ispec)
                  indice_grille_j(igll,jgll,kgll,ispec)=indice_grille_j(igll,jgll+1,kgll,ispec)
                  indice_grille_k(igll,jgll,kgll,ispec)=indice_grille_k(igll,jgll+1,kgll,ispec)
               endif


               if (indice_grille_i(igll,jgll,kgll,ispec)==0 .and. jgll==5) then
                  indice_grille_i(igll,jgll,kgll,ispec)=indice_grille_i(igll,jgll-1,kgll,ispec)
                  indice_grille_j(igll,jgll,kgll,ispec)=indice_grille_j(igll,jgll-1,kgll,ispec)
                  indice_grille_k(igll,jgll,kgll,ispec)=indice_grille_k(igll,jgll-1,kgll,ispec)
               endif

               if (indice_grille_i(igll,jgll,kgll,ispec)==0 .and. kgll==1) then
                  indice_grille_i(igll,jgll,kgll,ispec)=indice_grille_i(igll,jgll,kgll+1,ispec)
                  indice_grille_j(igll,jgll,kgll,ispec)=indice_grille_j(igll,jgll,kgll+1,ispec)
                  indice_grille_k(igll,jgll,kgll,ispec)=indice_grille_k(igll,jgll,kgll+1,ispec)
               endif


               if (indice_grille_i(igll,jgll,kgll,ispec)==0 .and. kgll==5) then
                  indice_grille_i(igll,jgll,kgll,ispec)=indice_grille_i(igll,jgll,kgll-1,ispec)
                  indice_grille_j(igll,jgll,kgll,ispec)=indice_grille_j(igll,jgll,kgll-1,ispec)
                  indice_grille_k(igll,jgll,kgll,ispec)=indice_grille_k(igll,jgll,kgll-1,ispec)
               endif


               if (indice_grille_i(igll,jgll,kgll,ispec)==0 .and. jgll==1 .and. kgll==1) then
                  indice_grille_i(igll,jgll,kgll,ispec)=indice_grille_i(igll,jgll+1,kgll+1,ispec)
                  indice_grille_j(igll,jgll,kgll,ispec)=indice_grille_j(igll,jgll+1,kgll+1,ispec)
                  indice_grille_k(igll,jgll,kgll,ispec)=indice_grille_k(igll,jgll+1,kgll+1,ispec)
               endif


               if (indice_grille_i(igll,jgll,kgll,ispec)==0 .and. jgll==5 .and. kgll==5) then
                  indice_grille_i(igll,jgll,kgll,ispec)=indice_grille_i(igll,jgll-1,kgll-1,ispec)
                  indice_grille_j(igll,jgll,kgll,ispec)=indice_grille_j(igll,jgll-1,kgll-1,ispec)
                  indice_grille_k(igll,jgll,kgll,ispec)=indice_grille_k(igll,jgll-1,kgll-1,ispec)
               endif

               if (indice_grille_i(igll,jgll,kgll,ispec)==0 .and. igll==1 .and. jgll==1) then
                  indice_grille_i(igll,jgll,kgll,ispec)=indice_grille_i(igll+1,jgll+1,kgll,ispec)
                  indice_grille_j(igll,jgll,kgll,ispec)=indice_grille_j(igll+1,jgll+1,kgll,ispec)
                  indice_grille_k(igll,jgll,kgll,ispec)=indice_grille_k(igll+1,jgll+1,kgll,ispec)
               endif


               if (indice_grille_i(igll,jgll,kgll,ispec)==0 .and. igll==5 .and. jgll==5) then
                  indice_grille_i(igll,jgll,kgll,ispec)=indice_grille_i(igll-1,jgll-1,kgll,ispec)
                  indice_grille_j(igll,jgll,kgll,ispec)=indice_grille_j(igll-1,jgll-1,kgll,ispec)
                  indice_grille_k(igll,jgll,kgll,ispec)=indice_grille_k(igll-1,jgll-1,kgll,ispec)
               endif



            enddo
         enddo
      enddo
   enddo
   if (myrank==0) then
      call cpu_time(finish)
      write(*,*) 'time // proj : ', finish-start
   endif

   valeur_int_glob(:,:,:,:)=0._CUSTOM_REAL
   call mpi_reduce(valeur_integration1,valeur_int_glob,NGLLX*NGLLY*NGLLZ*(nnx-1)*(nny-1)*(nnz-1),CUSTOM_MPI_TYPE, MPI_SUM,0,MPI_COMM_WORLD,ier)
   valeur_integration1(:,:,:,:) = valeur_int_glob(:,:,:,:)

   valeur_int_glob(:,:,:,:)=0._CUSTOM_REAL
   call mpi_reduce(valeur_integration2,valeur_int_glob,NGLLX*NGLLY*NGLLZ*(nnx-1)*(nny-1)*(nnz-1),CUSTOM_MPI_TYPE, MPI_SUM,0,MPI_COMM_WORLD,ier)
   valeur_integration2(:,:,:,:) = valeur_int_glob(:,:,:,:)

   !valeur_int_glob(:,:,:,:)=0._CUSTOM_REAL
   !call mpi_reduce(valeur_integration3,valeur_int_glob,NGLLX*NGLLY*NGLLZ*(nnx-1)*(nny-1)*(nnz-1),CUSTOM_MPI_TYPE, MPI_SUM,0,MPI_COMM_WORLD,ier)
   !valeur_integration3(:,:,:,:) = valeur_int_glob(:,:,:,:)

   !valeur_int_glob(:,:,:,:)=0._CUSTOM_REAL
   !call mpi_reduce(valeur_integration4,valeur_int_glob,NGLLX*NGLLY*NGLLZ*(nnx-1)*(nny-1)*(nnz-1),CUSTOM_MPI_TYPE, MPI_SUM,0,MPI_COMM_WORLD,ier)
   !valeur_integration4(:,:,:,:) = valeur_int_glob(:,:,:,:)

   call mpi_reduce(indice_integration,indice_int_glob,NGLLX*NGLLY*NGLLZ*(nnx-1)*(nny-1)*(nnz-1),MPI_INTEGER, MPI_SUM,0,MPI_COMM_WORLD,ier)
   call mpi_reduce(indice,indice_glob,(nnx-1)*(nny-1)*(nnz-1),MPI_INTEGER ,MPI_SUM,0, MPI_COMM_WORLD, ierr)

    call  mpi_reduce(xstore_tomo_grid,work_array,(nnx-1)*(nny-1)*(nnz-1),MPI_DOUBLE_PRECISION ,MPI_SUM,0, MPI_COMM_WORLD, ierr)
    xstore_tomo_grid(:,:,:)=work_array(:,:,:)

    call  mpi_reduce(ystore_tomo_grid,work_array,(nnx-1)*(nny-1)*(nnz-1),MPI_DOUBLE_PRECISION ,MPI_SUM,0, MPI_COMM_WORLD, ierr)
    ystore_tomo_grid(:,:,:)=work_array(:,:,:)

    call  mpi_reduce(zstore_tomo_grid,work_array,(nnx-1)*(nny-1)*(nnz-1),MPI_DOUBLE_PRECISION ,MPI_SUM,0, MPI_COMM_WORLD, ierr)
    zstore_tomo_grid(:,:,:)=work_array(:,:,:)

    call  mpi_reduce( profondeur_tomo_grid,work_array,(nnx-1)*(nny-1)*(nnz-1),MPI_DOUBLE_PRECISION ,MPI_SUM,0, MPI_COMM_WORLD, ierr)
    profondeur_tomo_grid(:,:,:)=work_array(:,:,:)

    indice(:,:,:) = indice_glob(:,:,:)
    indice_integration(:,:,:,:) = indice_int_glob(:,:,:,:)


!================== 2/ projection sur grille tomo =====================================

    if (myrank==0) then  !======================== myrank=0
       call  cpu_time(start)
       index_cellule=0
       do k=1,nnz-1
          do j=1,nny-1
             do i=1,nnx-1
                index_cellule=index_cellule+1

                if (indice(i,j,k)/=0)  then
                   xstore_tomo_grid(i,j,k) =  xstore_tomo_grid(i,j,k) / real(indice(i,j,k),8)
                   ystore_tomo_grid(i,j,k) =  ystore_tomo_grid(i,j,k) / real(indice(i,j,k),8)
                   zstore_tomo_grid(i,j,k) =  zstore_tomo_grid(i,j,k) / real(indice(i,j,k),8)
                   profondeur_tomo_grid(i,j,k) = profondeur_tomo_grid(i,j,k) / real(indice(i,j,k),8)
                else
                   write(*,*) 'point oublie :'  ,i,j,k
                endif


                !! on reagarde si les points on ete visites et combien de fois?
                do kgll=1,NGLLZ
                   do jgll=1,NGLLY
                      do igll=1,NGLLX

                         if (indice_integration(igll,jgll,kgll,index_cellule) == 0) then  !! le point n'a pas ete visite

                            !! chercher le point le plus proche qui est a ete visite
                            !  premier patch : devrait marcher pour mon cas cartesien
                            if (kgll==1 .and. indice_integration(igll,jgll,kgll+1,index_cellule) > 0) then

                               valeur_integration1(igll,jgll,kgll,index_cellule) = valeur_integration1(igll,jgll,kgll+1,index_cellule) /&
                                    real(indice_integration(igll,jgll,kgll+1,index_cellule),CUSTOM_REAL)
                               valeur_integration2(igll,jgll,kgll,index_cellule) = valeur_integration2(igll,jgll,kgll+1,index_cellule) /&
                                    real(indice_integration(igll,jgll,kgll+1,index_cellule),CUSTOM_REAL)
                             !  valeur_integration3(igll,jgll,kgll,index_cellule) = valeur_integration3(igll,jgll,kgll+1,index_cellule) /&
                             !       real(indice_integration(igll,jgll,kgll+1,index_cellule),CUSTOM_REAL)
                             !  valeur_integration4(igll,jgll,kgll,index_cellule) = valeur_integration4(igll,jgll,kgll+1,index_cellule) /&
                             !       real(indice_integration(igll,jgll,kgll+1,index_cellule),CUSTOM_REAL)

                            endif
                         else  ! le point a ete visite "indice_integration" fois
!!$                            ! on divise par indice pour avoir la valeur moyenne car on a fait un MPI_SUM

                            valeur_integration1(igll,jgll,kgll,index_cellule) =  valeur_integration1(igll,jgll,kgll,index_cellule)/&
                                 real(indice_integration(igll,jgll,kgll,index_cellule),CUSTOM_REAL)
                            valeur_integration2(igll,jgll,kgll,index_cellule) =  valeur_integration2(igll,jgll,kgll,index_cellule)/&
                                 real(indice_integration(igll,jgll,kgll,index_cellule),CUSTOM_REAL)
                            !valeur_integration3(igll,jgll,kgll,index_cellule) =  valeur_integration3(igll,jgll,kgll,index_cellule)/&
                            !     real(indice_integration(igll,jgll,kgll,index_cellule),CUSTOM_REAL)
                            !valeur_integration4(igll,jgll,kgll,index_cellule) =  valeur_integration4(igll,jgll,kgll,index_cellule)/&
                            !     real(indice_integration(igll,jgll,kgll,index_cellule),CUSTOM_REAL)

                         endif

                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo


       !===== boucle grille tomo : verifier qu'on fait exactement la meme boucle que precedemment
       np=0
       do k=1,nnz-1  ! grille integration (sphere cubique)


          p0= r_bot + (k-1)*hz
          p1= r_bot + k*hz

          do j=2,nny ! grille integration (sphere cubique)


             ratio_eta = (dble(j-2)) / dble(NY)
             y0 = 2.d0*ratio_eta-1
             y0 = tan((ANGULAR_WIDTH_ETA_RAD/2.d0) * y0)

             ratio_eta = (dble(j-1)) / dble(NY)
             y1 = 2.d0*ratio_eta-1
             y1 = tan((ANGULAR_WIDTH_ETA_RAD/2.d0) * y1)


             do i=2,nnx ! grille integration (sphere cubique)

                index_cellule=i-1+(j-2)*(nnx-1)+(k-1)*(nnx-1)*(nny-1)

                ratio_xi = (dble(i-2)) / dble(NX)
                x0 = 2.d0*ratio_xi-1
                x0 = tan((ANGULAR_WIDTH_XI_RAD/2.d0) * x0)

                ratio_xi = (dble(i-1)) / dble(NX)
                x1 = 2.d0*ratio_xi-1
                x1 = tan((ANGULAR_WIDTH_XI_RAD/2.d0) * x1)

                !===== calcul des 8 sommets de la cellule d'integration dans le repere cartesien SEM

                ! 1
                x=x0;y=y0;p=p0
                xstore_local(1) = x*p
                ystore_local(1) = y*p
                zstore_local(1) =  -(R_EARTH - p/dsqrt(1.d0 + y*y + x*x) - zmax)

                ! 2
                x=x1;y=y0;p=p0
                xstore_local(2) = x*p
                ystore_local(2) = y*p
                zstore_local(2) =  -(R_EARTH - p/dsqrt(1.d0 + y*y + x*x) - zmax)

                ! 3
                x=x1;y=y1;p=p0
                xstore_local(3) = x*p
                ystore_local(3) = y*p
                zstore_local(3) =  -(R_EARTH - p/dsqrt(1.d0 + y*y + x*x) - zmax)

                ! 4
                x=x0;y=y1;p=p0
                xstore_local(4) = x*p
                ystore_local(4) = y*p
                zstore_local(4) =  -(R_EARTH - p/dsqrt(1.d0 + y*y + x*x) - zmax)

                ! 5
                x=x0;y=y0;p=p1
                xstore_local(5) = x*p
                ystore_local(5) = y*p
                zstore_local(5) =  -(R_EARTH - p/dsqrt(1.d0 + y*y + x*x) - zmax)

                ! 6
                x=x1;y=y0;p=p1
                xstore_local(6) = x*p
                ystore_local(6) = y*p
                zstore_local(6) =  -(R_EARTH - p/dsqrt(1.d0 + y*y + x*x) - zmax)

                ! 7
                x=x1;y=y1;p=p1
                xstore_local(7) = x*p
                ystore_local(7) = y*p
                zstore_local(7) =  -(R_EARTH - p/dsqrt(1.d0 + y*y + x*x) - zmax)


                ! 8
                x=x0;y=y1;p=p1
                xstore_local(8) = x*p
                ystore_local(8) = y*p
                zstore_local(8) =  -(R_EARTH - p/dsqrt(1.d0 + y*y + x*x) - zmax)

                ! === calcul de l'integrale de volume
                call compute_volume_int(volume1,xstore_local,ystore_local,zstore_local,wxgll,wygll,wzgll,&
                     shape3D,dershape3D,valeur_integration1,index_cellule,(nnx-1)*(nny-1)*(nnz-1))

                call compute_volume_int(volume2,xstore_local,ystore_local,zstore_local,wxgll,wygll,wzgll,&
                     shape3D,dershape3D,valeur_integration2,index_cellule,(nnx-1)*(nny-1)*(nnz-1))

                !call compute_volume_int(volume3,xstore_local,ystore_local,zstore_local,wxgll,wygll,wzgll,&
                !     shape3D,dershape3D,valeur_integration3,index_cellule,(nnx-1)*(nny-1)*(nnz-1))

                !call compute_volume_int(volume4,xstore_local,ystore_local,zstore_local,wxgll,wygll,wzgll,&
                !     shape3D,dershape3D,valeur_integration4,index_cellule,(nnx-1)*(nny-1)*(nnz-1))

                if (ires==1)  then
                   ! === si on veut faire juste la projection du modele, on ne doit pas normaliser par
                   ! === le volume de la cellule, on prend plutot la valeur moyenne.
                   call compute_volume_int(volume_cellule,xstore_local,ystore_local,zstore_local,wxgll,wygll,wzgll,&
                        shape3D,dershape3D,volume_integration,index_cellule,(nnx-1)*(nny-1)*(nnz-1))
                   volume1=volume1/volume_cellule
                   volume2=volume2/volume_cellule
                   !volume3=volume3/volume_cellule
                   !volume4=volume4/volume_cellule
                endif

                if (indice(i-1,j-1,k)==0) then
                   ! on a oublie le point i,j,k
                   np=np+1
                   !write(30,*) i,j,k
                else

                   data_tomo1(i-1,j-1,k) = volume1
                   data_tomo2(i-1,j-1,k) = volume2
                   !data_tomo3(i-1,j-1,k) = volume3
                   !data_tomo4(i-1,j-1,k) = volume4

                endif

             enddo
          enddo

       enddo
       call cpu_time(finish)
       write(*,*) 'time seq. ',finish-start

    endif
    call mpi_bcast(data_tomo1,nx*ny*nz,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ier)
    call mpi_bcast(data_tomo2,nx*ny*nz,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ier)
    !call mpi_bcast(data_tomo3,nx*ny*nz,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ier)
    !call mpi_bcast(data_tomo4,nx*ny*nz,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ier)
    call mpi_bcast(profondeur_tomo_grid,nx*ny*nz,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)


  end subroutine project_tomo_grid_2






  subroutine create_box_grid_projection()
    implicit none

    integer ispec,ix,iy,iz,igrid_local
    integer imin,imax,jmin,jmax,kmin,kmax
    integer i,j,k,igll,jgll,kgll
    integer index_grid,ierr
    double precision x_center,y_center,z_center
    integer ii,jj,kk


    indice_grille_i(:,:,:,:)=0
    indice_grille_j(:,:,:,:)=0
    indice_grille_k(:,:,:,:)=0
    indice_grille_i1(:,:,:,:)=0
    indice_grille_j1(:,:,:,:)=0
    indice_grille_k1(:,:,:,:)=0
    indice(:,:,:)=0
    indice_glob(:,:,:)=0

    xstore_tomo_grid(:,:,:)=0.d0
    ystore_tomo_grid(:,:,:)=0.d0
    zstore_tomo_grid(:,:,:)=0.d0
    profondeur_tomo_grid(:,:,:)=0.d0
    volume_integration(:,:,:,:)=1._CUSTOM_REAL

    hx = (xmax - xmin) / (nx)
    hy = (ymax - ymin) / (ny)
    hz = (zmax - zmin) / (nz)

    if (myrank==0) then
       write(*,*)
       write(*,*) xmin,xmax
       write(*,*) ymin,ymax
       write(*,*) zmin,zmax
       write(*,*) nx,ny,nz
       write(*,*) hx,hy,hz
       write(*,*)
    endif

    !! define interp tri grid (tomo -> SEM)

    do k=1,nz+1
       do j=1,ny+1
          do i=1,nx+1

             xstore_tomo_grid_interp_tri(i,j,k) = xmin + (i-1)*hx
             ystore_tomo_grid_interp_tri(i,j,k) = ymin + (j-1)*hy
             zstore_tomo_grid_interp_tri(i,j,k) = zmin + (k-1)*hz

             enddo
       enddo
    enddo

    do ispec = 1, NSPEC_AB   !===== boucle element ispec

       ni(ispec)=0

       xmin_local =  HUGEVAL
       xmax_local = -HUGEVAL
       ymin_local =  HUGEVAL
       ymax_local = -HUGEVAL
       zmin_local =  HUGEVAL
       zmax_local = -HUGEVAL

!======================= Sommets de l'element ============
!1
       xstore_local(1)=xstore(ibool(1,1,1,ispec))
       ystore_local(1)=ystore(ibool(1,1,1,ispec))
       zstore_local(1)=zstore(ibool(1,1,1,ispec))
!2
       xstore_local(2)=xstore(ibool(NGLLX,1,1,ispec))
       ystore_local(2)=ystore(ibool(NGLLX,1,1,ispec))
       zstore_local(2)=zstore(ibool(NGLLX,1,1,ispec))
!3
       xstore_local(3)=xstore(ibool(NGLLX,NGLLY,1,ispec))
       ystore_local(3)=ystore(ibool(NGLLX,NGLLY,1,ispec))
       zstore_local(3)=zstore(ibool(NGLLX,NGLLY,1,ispec))
!4
       xstore_local(4)=xstore(ibool(1,NGLLY,1,ispec))
       ystore_local(4)=ystore(ibool(1,NGLLY,1,ispec))
       zstore_local(4)=zstore(ibool(1,NGLLY,1,ispec))
!5
       xstore_local(5)=xstore(ibool(1,1,NGLLZ,ispec))
       ystore_local(5)=ystore(ibool(1,1,NGLLZ,ispec))
       zstore_local(5)=zstore(ibool(1,1,NGLLZ,ispec))
!6
       xstore_local(6)=xstore(ibool(NGLLX,1,NGLLZ,ispec))
       ystore_local(6)=ystore(ibool(NGLLX,1,NGLLZ,ispec))
       zstore_local(6)=zstore(ibool(NGLLX,1,NGLLZ,ispec))
!7
       xstore_local(7)=xstore(ibool(NGLLX,NGLLY,NGLLZ,ispec))
       ystore_local(7)=ystore(ibool(NGLLX,NGLLY,NGLLZ,ispec))
       zstore_local(7)=zstore(ibool(NGLLX,NGLLY,NGLLZ,ispec))
!8
       xstore_local(8)=xstore(ibool(1,NGLLY,NGLLZ,ispec))
       ystore_local(8)=ystore(ibool(1,NGLLY,NGLLZ,ispec))
       zstore_local(8)=zstore(ibool(1,NGLLY,NGLLZ,ispec))

!=============== on cherche le pave circonscrit a l'element en se basant sur les point GLL

       do iz=1,NGLLZ
          do iy=1,NGLLY
             do ix=1,NGLLX
                xmin_local=min(xmin_local,xstore(ibool(ix,iy,iz,ispec)))
                xmax_local=max(xmin_local,xstore(ibool(ix,iy,iz,ispec)))
                ymin_local=min(ymin_local,ystore(ibool(ix,iy,iz,ispec)))
                ymax_local=max(ymin_local,ystore(ibool(ix,iy,iz,ispec)))
                zmin_local=min(zmin_local,zstore(ibool(ix,iy,iz,ispec)))
                zmax_local=max(zmax_local,zstore(ibool(ix,iy,iz,ispec)))
             enddo
          enddo
       enddo


! =========== calcul des indices min et max de la grille tomo correspondant au pave

       kmin = MAX(1,floor((zmin_local  - zmin) / hz) +1)
       kmax=  MIN(NZ,floor((zmax_local  - zmin) / hz) +1)

       jmin=  MAX(1,floor((ymin_local  - ymin) / hy) +1)
       jmax=  MIN(NY,floor((ymax_local  - ymin) / hy) +1)

       imin=  MAX(1,floor((xmin_local - xmin) / hx) +1)
       imax=  MIN(NX,floor((xmax_local - xmin) / hx) +1)
!!$       if (myrank==0) then
!!$       write(*,*) ymin_local,  ymin + (jmin+0.5)*hy, ymin + (jmax+0.5)*hy,ymax_local
!!$    endif
!ymin_local,ymax_local,jmin,jmax,hy,ny!ispec,imin,imax,jmin,jmax,kmin,kmax
       do k = kmin, kmax
          z_center = zmin + (k-0.5)*hz
          do j = jmin, jmax
             y_center = ymin + (j-0.5)*hy
             !if (myrank==0) write(*,*) j,y_center,jmin,jmax
             do i = imin, imax
                x_center = xmin + (i-0.5)*hx

                index_grid = i + (j-1)*nx + (k-1)*nx*ny ! cell ijk

                ! 8 corners
                xstore_local_int(1)= xmin + (i-1)*hx
                ystore_local_int(1)= ymin + (j-1)*hy
                zstore_local_int(1)= zmin + (k-1)*hz

                xstore_local_int(2)= xmin + (i  )*hx
                ystore_local_int(2)= ymin + (j-1)*hy
                zstore_local_int(2)= zmin + (k-1)*hz

                xstore_local_int(3)= xmin + (i  )*hx
                ystore_local_int(3)= ymin + (j  )*hy
                zstore_local_int(3)= zmin + (k-1)*hz

                xstore_local_int(4)= xmin + (i-1)*hx
                ystore_local_int(4)= ymin + (j  )*hy
                zstore_local_int(4)= zmin + (k-1)*hz

                xstore_local_int(5)= xmin + (i-1)*hx
                ystore_local_int(5)= ymin + (j-1)*hy
                zstore_local_int(5)= zmin + (k  )*hz

                xstore_local_int(6)= xmin + (i  )*hx
                ystore_local_int(6)= ymin + (j-1)*hy
                zstore_local_int(6)= zmin + (k  )*hz

                xstore_local_int(7)= xmin + (i  )*hx
                ystore_local_int(7)= ymin + (j  )*hy
                zstore_local_int(7)= zmin + (k  )*hz

                xstore_local_int(8)= xmin + (i-1)*hx
                ystore_local_int(8)= ymin + (j  )*hy
                zstore_local_int(8)= zmin + (k  )*hz



                ! ========== compute GLL points
                ! compute xix, eta et gamma for the central point central x_center,y_center,z_center
                call Find_xix_eta_gamma(xstore_local,ystore_local,zstore_local,xi,eta,gamma,x_center,y_center,z_center)

                ! are we inside ispec?
                if (dabs(xi)<=1.05d0.and.dabs(eta)<=1.05d0.and.dabs(gamma)<=1.05d0) then
                   indice(i,j,k)=1 ! on a bien visite la cellule de la grille
                   xstore_tomo_grid(i,j,k) = x_center
                   ystore_tomo_grid(i,j,k) = y_center
                   zstore_tomo_grid(i,j,k) = z_center
                   profondeur_tomo_grid(i,j,k) = abs(zmax - z_center)
                   indice_spec(index_grid)=ispec !! correspondance cellule_tomo <-> elemements
                   ! la cellule index_cellule est cencee appartenir a ispec et tous ces pts gll aussi.
                   ni(ispec)=ni(ispec)+1
                   indice_grid(1,ni(ispec),ispec)=i
                   indice_grid(2,ni(ispec),ispec)=j
                   indice_grid(3,ni(ispec),ispec)=k

                   ! on calcule les points gll de la cellule :
!                   call calcule_point_gll(xgll,ygll,zgll,xstore_local_int,ystore_local_int,zstore_local_int,shape3D)
!                   do kgll = 1,NGLLZ
!                      do jgll = 1,NGLLY
!                         do igll = 1,NGLLX
                            !indice_integration(igll,jgll,kgll,index_grid)=1
                            !x_grid_gll(igll,jgll,kgll,index_grid) = xgll(igll,jgll,kgll)
                            !y_grid_gll(igll,jgll,kgll,index_grid) = ygll(igll,jgll,kgll)
                            !z_grid_gll(igll,jgll,kgll,index_grid) = zgll(igll,jgll,kgll)
!!$                            if (igll==3 .and. jgll==3 .and. kgll==3 .and. i==75 .and. j==78 ) then
!!$                               write(*,*)  myrank,xgll(igll,jgll,kgll) , ygll(igll,jgll,kgll), zgll(igll,jgll,kgll)
!!$                            endif
!                         enddo
!                      enddo
!                   enddo



!                else
                   !write(*,*)  xi,eta,gamma,x_center,y_center, z_center
!                endif

                   ! ========================= on cherche les points GLL qui sont dans la cellule de la grille tomo
                   call calcule_point_gll(xgll,ygll,zgll,xstore_local,ystore_local,zstore_local,shape3D)
                   do kgll = 1,NGLLZ
                      do jgll = 1,NGLLY
                         do igll = 1,NGLLX

                            ! on peut eliminer calcule_point_gll
                            ! et faire ca :
                            !px=xstore(ibool(igll,jgll,kgll,ispec))
                            !py=ystore(ibool(igll,jgll,kgll,ispec))
                            !pz=zstore(ibool(igll,jgll,kgll,ispec))

                            px = xgll(igll,jgll,kgll)
                            py = ygll(igll,jgll,kgll)
                            pz = zgll(igll,jgll,kgll)

                            indice_grille_i(igll,jgll,kgll,ispec) = floor((px - xmin) / hx) +1
                            indice_grille_j(igll,jgll,kgll,ispec) = floor((py - ymin) / hy) +1
                            indice_grille_k(igll,jgll,kgll,ispec) = floor((pz - zmin) / hz) +1

                            if (indice_grille_i(igll,jgll,kgll,ispec) < 1) indice_grille_i(igll,jgll,kgll,ispec)=1
                            if (indice_grille_j(igll,jgll,kgll,ispec) < 1) indice_grille_j(igll,jgll,kgll,ispec)=1
                            if (indice_grille_k(igll,jgll,kgll,ispec) < 1) indice_grille_k(igll,jgll,kgll,ispec)=1

                            indice_grille_i1(igll,jgll,kgll,ispec) = indice_grille_i(igll,jgll,kgll,ispec) + 1
                            indice_grille_j1(igll,jgll,kgll,ispec) = indice_grille_j(igll,jgll,kgll,ispec) + 1
                            indice_grille_k1(igll,jgll,kgll,ispec) = indice_grille_k(igll,jgll,kgll,ispec) + 1

                            if (indice_grille_i1(igll,jgll,kgll,ispec) > nx+1) indice_grille_i1(igll,jgll,kgll,ispec)=nx+1
                            if (indice_grille_j1(igll,jgll,kgll,ispec) > ny+1) indice_grille_j1(igll,jgll,kgll,ispec)=ny+1
                            if (indice_grille_k1(igll,jgll,kgll,ispec) > nz+1) indice_grille_k1(igll,jgll,kgll,ispec)=nz+1


                            if (indice_grille_i(igll,jgll,kgll,ispec) > nx) indice_grille_i(igll,jgll,kgll,ispec)=nx
                            if (indice_grille_j(igll,jgll,kgll,ispec) > ny) indice_grille_j(igll,jgll,kgll,ispec)=ny
                            if (indice_grille_k(igll,jgll,kgll,ispec) > nz) indice_grille_k(igll,jgll,kgll,ispec)=nz

!!$                            if (px >= xmax) then
!!$                               ii=indice_grille_i(igll,jgll,kgll,ispec)
!!$                               jj=indice_grille_j(igll,jgll,kgll,ispec)
!!$                               kk=indice_grille_k(igll,jgll,kgll,ispec)
!!$                               write(*,*) myrank,ii,jj,kk
!!$                               write(*,*) myrank, px, xstore_tomo_grid(ii,jj,kk)
!!$                               stop
!!$                            endif
!!$                            if (myrank==0) then
!!$                               write(*,*) px,py,pz
!!$                               ii=indice_grille_i(igll,jgll,kgll,ispec)
!!$                               jj=indice_grille_j(igll,jgll,kgll,ispec)
!!$                               kk=indice_grille_k(igll,jgll,kgll,ispec)
!!$                               WRITE(*,*) II,JJ,KK
!!$                               write(*,*) xstore_tomo_grid(ii,jj,kk),ystore_tomo_grid(ii,jj,kk),zstore_tomo_grid(ii,jj,kk)
!!$                               ii=indice_grille_i1(igll,jgll,kgll,ispec)
!!$                               jj=indice_grille_j1(igll,jgll,kgll,ispec)
!!$                               kk=indice_grille_k1(igll,jgll,kgll,ispec)
!!$                               WRITE(*,*) II,JJ,KK
!!$                               write(*,*) xstore_tomo_grid(ii,jj,kk),ystore_tomo_grid(ii,jj,kk),zstore_tomo_grid(ii,jj,kk)
!!$                               write(*,*) '-------------------------------------'
!!$                            endif
!!$                            call Find_xix_eta_gamma(xstore_local_int,ystore_local_int,zstore_local_int,xi,eta,gamma,px,py,pz)
!!$                            if (dabs(xi)<=1.05d0.and.dabs(eta)<=1.05d0.and.dabs(gamma)<=1.05d0) then
!!$
!!$                               indice_grille_i(igll,jgll,kgll,ispec) = i
!!$                               indice_grille_i1(igll,jgll,kgll,ispec)= i
!!$
!!$                               indice_grille_j(igll,jgll,kgll,ispec) = j
!!$                               indice_grille_j1(igll,jgll,kgll,ispec)= j
!!$
!!$                               indice_grille_k(igll,jgll,kgll,ispec) = k
!!$                               indice_grille_k1(igll,jgll,kgll,ispec)= k

                            ! je laisse ca mais je ne suis pas encore sur de son utilite
!!$                            if (xi > 0. ) then
!!$                               if (i<nx) indice_grille_i1(igll,jgll,kgll,ispec) = i+1
!!$                            else
!!$                               if (i>1)  indice_grille_i1(igll,jgll,kgll,ispec) = i-1
!!$                            endif
!!$
!!$                            if (eta > 0. ) then
!!$                               if (j<ny) indice_grille_j1(igll,jgll,kgll,ispec) = j+1
!!$                            else
!!$                               if (j>1)  indice_grille_j1(igll,jgll,kgll,ispec) = j-1
!!$                            endif
!!$
!!$                            if (gamma > 0. ) then
!!$                               if (k<nz) indice_grille_k1(igll,jgll,kgll,ispec) = k+1
!!$                            else
!!$                               if (k>1)  indice_grille_k1(igll,jgll,kgll,ispec) = k-1
!!$                            endif


!!$                            endif
                         enddo
                      enddo
                   enddo
                endif

             enddo
          enddo
       enddo



    enddo

    ! allocation des tableaux locaux
    ngrid_local=sum(ni(:))
    allocate(indice_integration(NGLLX,NGLLY,NGLLZ,ngrid_local))
    allocate(valeur_integration(NGLLX,NGLLY,NGLLZ,ngrid_local))
    allocate(volume_integration(NGLLX,NGLLY,NGLLZ,ngrid_local))
    allocate(x_grid_gll(NGLLX,NGLLY,NGLLZ,ngrid_local))
    allocate(y_grid_gll(NGLLX,NGLLY,NGLLZ,ngrid_local))
    allocate(z_grid_gll(NGLLX,NGLLY,NGLLZ,ngrid_local))

    index_grid=0
    do ispec=1,nspec_ab           ! loop on spectral elements
       do igrid_local=1,ni(ispec) ! loop on grid cells instide ispec
          index_grid=index_grid+1
          ! cell index
          i= indice_grid(1,igrid_local,ispec)
          j= indice_grid(2,igrid_local,ispec)
          k= indice_grid(3,igrid_local,ispec)


          !8 corners
          xstore_local_int(1)= xmin + (i-1)*hx
          ystore_local_int(1)= ymin + (j-1)*hy
          zstore_local_int(1)= zmin + (k-1)*hz

          xstore_local_int(2)= xmin + (i  )*hx
          ystore_local_int(2)= ymin + (j-1)*hy
          zstore_local_int(2)= zmin + (k-1)*hz

          xstore_local_int(3)= xmin + (i  )*hx
          ystore_local_int(3)= ymin + (j  )*hy
          zstore_local_int(3)= zmin + (k-1)*hz

          xstore_local_int(4)= xmin + (i-1)*hx
          ystore_local_int(4)= ymin + (j  )*hy
          zstore_local_int(4)= zmin + (k-1)*hz

          xstore_local_int(5)= xmin + (i-1)*hx
          ystore_local_int(5)= ymin + (j-1)*hy
          zstore_local_int(5)= zmin + (k  )*hz

          xstore_local_int(6)= xmin + (i  )*hx
          ystore_local_int(6)= ymin + (j-1)*hy
          zstore_local_int(6)= zmin + (k  )*hz

          xstore_local_int(7)= xmin + (i  )*hx
          ystore_local_int(7)= ymin + (j  )*hy
          zstore_local_int(7)= zmin + (k  )*hz

          xstore_local_int(8)= xmin + (i-1)*hx
          ystore_local_int(8)= ymin + (j  )*hy
          zstore_local_int(8)= zmin + (k  )*hz

          ! compute gll points of cell :
          call calcule_point_gll(xgll,ygll,zgll,xstore_local_int,ystore_local_int,zstore_local_int,shape3D)
          do kgll = 1,NGLLZ
             do jgll = 1,NGLLY
                do igll = 1,NGLLX
                   !index_grid = index_grid + 1
                   indice_integration(igll,jgll,kgll,index_grid) = 1
                   x_grid_gll(igll,jgll,kgll,index_grid) = xgll(igll,jgll,kgll)
                   y_grid_gll(igll,jgll,kgll,index_grid) = ygll(igll,jgll,kgll)
                   z_grid_gll(igll,jgll,kgll,index_grid) = zgll(igll,jgll,kgll)
                enddo
             enddo
          enddo

       enddo
    enddo

    ! mpi comm
    call mpi_reduce(indice,indice_glob,(nnx-1)*(nny-1)*(nnz-1),MPI_INTEGER ,MPI_SUM,0, MPI_COMM_WORLD, ierr)
    indice(:,:,:) = indice_glob(:,:,:)
    call mpi_bcast(indice,(nnx-1)*(nny-1)*(nnz-1),MPI_INTEGER,0, MPI_COMM_WORLD, ierr)

    call  mpi_reduce(xstore_tomo_grid,work_array,(nnx-1)*(nny-1)*(nnz-1),MPI_DOUBLE_PRECISION ,MPI_SUM,0, MPI_COMM_WORLD, ierr)
    xstore_tomo_grid(:,:,:)=work_array(:,:,:)


    call  mpi_reduce(ystore_tomo_grid,work_array,(nnx-1)*(nny-1)*(nnz-1),MPI_DOUBLE_PRECISION ,MPI_SUM,0, MPI_COMM_WORLD, ierr)
    ystore_tomo_grid(:,:,:)=work_array(:,:,:)


    call  mpi_reduce(zstore_tomo_grid,work_array,(nnx-1)*(nny-1)*(nnz-1),MPI_DOUBLE_PRECISION ,MPI_SUM,0, MPI_COMM_WORLD, ierr)
    zstore_tomo_grid(:,:,:)=work_array(:,:,:)


    call  mpi_reduce(profondeur_tomo_grid,work_array,(nnx-1)*(nny-1)*(nnz-1),MPI_DOUBLE_PRECISION ,MPI_SUM,0, MPI_COMM_WORLD, ierr)
    profondeur_tomo_grid(:,:,:)=work_array(:,:,:)


    ! ces tableaux deviennent locaux :
    !call mpi_reduce(indice_integration,indice_int_glob,NGLLX*NGLLY*NGLLZ*(nnx-1)*(nny-1)*(nnz-1),MPI_INTEGER, MPI_SUM,0,MPI_COMM_WORLD,ier)
    !indice_integration(:,:,:,:) = indice_int_glob(:,:,:,:)

    !call  mpi_reduce(x_grid_gll,wk_reduce,NGLLX*NGLLY*NGLLZ*(nnx-1)*(nny-1)*(nnz-1),MPI_DOUBLE_PRECISION ,MPI_SUM,0, MPI_COMM_WORLD, ierr)
    !x_grid_gll(:,:,:,:)=wk_reduce(:,:,:,:)
    !call mpi_bcast(x_grid_gll,NGLLX*NGLLY*NGLLZ*(nnx-1)*(nny-1)*(nnz-1),MPI_DOUBLE_PRECISION,0, MPI_COMM_WORLD, ierr)

    !call  mpi_reduce(y_grid_gll,wk_reduce,NGLLX*NGLLY*NGLLZ*(nnx-1)*(nny-1)*(nnz-1),MPI_DOUBLE_PRECISION ,MPI_SUM,0, MPI_COMM_WORLD, ierr)
    !y_grid_gll(:,:,:,:)=wk_reduce(:,:,:,:)
    !call mpi_bcast(y_grid_gll,NGLLX*NGLLY*NGLLZ*(nnx-1)*(nny-1)*(nnz-1),MPI_DOUBLE_PRECISION,0, MPI_COMM_WORLD, ierr)

    !call  mpi_reduce(z_grid_gll,wk_reduce,NGLLX*NGLLY*NGLLZ*(nnx-1)*(nny-1)*(nnz-1),MPI_DOUBLE_PRECISION ,MPI_SUM,0, MPI_COMM_WORLD, ierr)
    !z_grid_gll(:,:,:,:)=wk_reduce(:,:,:,:)
    !call mpi_bcast(z_grid_gll,NGLLX*NGLLY*NGLLZ*(nnx-1)*(nny-1)*(nnz-1),MPI_DOUBLE_PRECISION,0, MPI_COMM_WORLD, ierr)



    !write(*,*) 'MAX VALUE NI : ', myrank, maxval(ni)
  end subroutine create_box_grid_projection


  subroutine check_projection_grid()
    implicit none
    integer i,j,k
    integer igll,jgll,kgll
    integer ierr
    !integer index_cellule
    if (myrank==0) then
       !index_cellule=0
       do k=1,nnz-1
          do j=1,nny-1
             do i=1,nnx-1
                !index_cellule=index_cellule+1
                if (indice(i,j,k)==0) then
                   write(*,*) 'Warnning forget cell :',i,j,k,&
                        xstore_tomo_grid(i,j,k),ystore_tomo_grid(i,j,k),zstore_tomo_grid(i,j,k)
                else
                   xstore_tomo_grid(i,j,k) =  xstore_tomo_grid(i,j,k) / real(indice(i,j,k),8)
                   ystore_tomo_grid(i,j,k) =  ystore_tomo_grid(i,j,k) / real(indice(i,j,k),8)
                   zstore_tomo_grid(i,j,k) =  zstore_tomo_grid(i,j,k) / real(indice(i,j,k),8)
                   profondeur_tomo_grid(i,j,k) =  profondeur_tomo_grid(i,j,k) / real(indice(i,j,k),8)

                endif

!!$                do kgll=1,NGLLZ
!!$                   do jgll=1,NGLLY
!!$                      do igll=1,NGLLX
!!$
!!$                         if (indice_integration(igll,jgll,kgll,index_cellule)==0) then
!!$                            write(*,*) 'Warnning forget gll point ', igll,jgll,kgll,i,j,k
!!$                         else
!!$
!!$                            x_grid_gll(igll,jgll,kgll,index_cellule) = x_grid_gll(igll,jgll,kgll,index_cellule) &
!!$                                 / real(indice_integration(igll,jgll,kgll,index_cellule),8)
!!$
!!$                            y_grid_gll(igll,jgll,kgll,index_cellule) = y_grid_gll(igll,jgll,kgll,index_cellule) &
!!$                                 / real(indice_integration(igll,jgll,kgll,index_cellule),8)
!!$
!!$                            z_grid_gll(igll,jgll,kgll,index_cellule) = z_grid_gll(igll,jgll,kgll,index_cellule) &
!!$                                 / real(indice_integration(igll,jgll,kgll,index_cellule),8)
!!$
!!$
!!$
!!$
!!$                        endif
!!$
!!$                      enddo
!!$                   enddo
!!$                enddo

             enddo
          enddo
       enddo
    endif
    call mpi_bcast(xstore_tomo_grid,(nnx-1)*(nny-1)*(nnz-1),MPI_DOUBLE_PRECISION ,0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(ystore_tomo_grid,(nnx-1)*(nny-1)*(nnz-1),MPI_DOUBLE_PRECISION ,0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(zstore_tomo_grid,(nnx-1)*(nny-1)*(nnz-1),MPI_DOUBLE_PRECISION ,0, MPI_COMM_WORLD, ierr)
    call mpi_bcast(profondeur_tomo_grid,nx*ny*nz,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
    if (myrank == 0 ) write(*,*) 'tomo grid checked '
  end subroutine check_projection_grid



  subroutine interpole_one_field(dat_to_project)
    implicit none
    integer ispec,ires,i,j,k,icell,index_grid
    integer igll,jgll,kgll
    real(kind=CUSTOM_REAL) dat_to_project(NGLLX,NGLLY,NGLLZ,NSPEC_AB)

    valeur_integration(:,:,:,:)=0._CUSTOM_REAL
    indice_integration(:,:,:,:)=0

    !write(*,*) 'in interpole '
    index_grid=0
    do ispec=1,NSPEC_AB ! loop on spectral elements

       !======================= Sommets de l'element ============
!1
       xstore_local(1)=xstore(ibool(1,1,1,ispec))
       ystore_local(1)=ystore(ibool(1,1,1,ispec))
       zstore_local(1)=zstore(ibool(1,1,1,ispec))
!2
       xstore_local(2)=xstore(ibool(NGLLX,1,1,ispec))
       ystore_local(2)=ystore(ibool(NGLLX,1,1,ispec))
       zstore_local(2)=zstore(ibool(NGLLX,1,1,ispec))
!3
       xstore_local(3)=xstore(ibool(NGLLX,NGLLY,1,ispec))
       ystore_local(3)=ystore(ibool(NGLLX,NGLLY,1,ispec))
       zstore_local(3)=zstore(ibool(NGLLX,NGLLY,1,ispec))
!4
       xstore_local(4)=xstore(ibool(1,NGLLY,1,ispec))
       ystore_local(4)=ystore(ibool(1,NGLLY,1,ispec))
       zstore_local(4)=zstore(ibool(1,NGLLY,1,ispec))
!5
       xstore_local(5)=xstore(ibool(1,1,NGLLZ,ispec))
       ystore_local(5)=ystore(ibool(1,1,NGLLZ,ispec))
       zstore_local(5)=zstore(ibool(1,1,NGLLZ,ispec))
!6
       xstore_local(6)=xstore(ibool(NGLLX,1,NGLLZ,ispec))
       ystore_local(6)=ystore(ibool(NGLLX,1,NGLLZ,ispec))
       zstore_local(6)=zstore(ibool(NGLLX,1,NGLLZ,ispec))
!7
       xstore_local(7)=xstore(ibool(NGLLX,NGLLY,NGLLZ,ispec))
       ystore_local(7)=ystore(ibool(NGLLX,NGLLY,NGLLZ,ispec))
       zstore_local(7)=zstore(ibool(NGLLX,NGLLY,NGLLZ,ispec))
!8
       xstore_local(8)=xstore(ibool(1,NGLLY,NGLLZ,ispec))
       ystore_local(8)=ystore(ibool(1,NGLLY,NGLLZ,ispec))
       zstore_local(8)=zstore(ibool(1,NGLLY,NGLLZ,ispec))



       do icell=1,ni(ispec) ! loop on tomo cells
          index_grid=index_grid+1
          i=indice_grid(1,icell,ispec)
          j=indice_grid(2,icell,ispec)
          k=indice_grid(3,icell,ispec)

          !index_grid=i+nx*(j-1)+nx*ny*(k-1) ! il faut verifier que index_grid a vu toute la grille

          do kgll=1,NGLLZ
             do jgll=1,NGLLY
                do igll=1,NGLLX
                   !index_grid=index_grid+1
                   px=x_grid_gll(igll,jgll,kgll,index_grid)
                   py=y_grid_gll(igll,jgll,kgll,index_grid)
                   pz=z_grid_gll(igll,jgll,kgll,index_grid)

!!$                   if (igll==3 .and. jgll==3 .and. kgll==3 .and. i==75 .and. j==78 ) then
!!$
!!$                      write(*,*)'----------------------------------------------',myrank
!!$                      write(*,*) xstore_local(1),xstore_local(7)
!!$                      write(*,*) ystore_local(1),ystore_local(7)
!!$                      write(*,*) zstore_local(1),zstore_local(7)
!!$                      write(*,*)
!!$
!!$                      write(*,'(3f15.5)') px,py,pz
!!$                      write(*,*) valeur
!!$                   endif

!                   if (myrank==0) then
!                      write(*,*) '-> ',myrank,px,py,pz
!                      write(*,*) '-> ', igll,jgll,kgll,index_grid,i,j,k
!                      write(*,*)
!                   endif
                   call interpolation_valeur(valeur,px,py,pz,dat_to_project,xstore_local,ystore_local,zstore_local,&
                                                 ispec,NSPEC_AB,xigll,yigll,zigll)

                   valeur_integration(igll,jgll,kgll,index_grid)=valeur
                   indice_integration(igll,jgll,kgll,index_grid)=1
                   !if (i==1 .And. j==1) then
                   !    write(*,*) myrank, k,index_grid,valeur
                   !    write(*,*) dat_to_project(1,1,1,ispec)
                   !endif
                enddo
             enddo
          enddo


       enddo
    enddo

    ! mpi comm
    !valeur_int_glob(:,:,:,:)=0._CUSTOM_REAL
    !call mpi_reduce(valeur_integration,valeur_int_glob,NGLLX*NGLLY*NGLLZ*(nnx-1)*(nny-1)*(nnz-1),CUSTOM_MPI_TYPE, MPI_SUM,0,MPI_COMM_WORLD,ier)
    !valeur_integration(:,:,:,:)=valeur_int_glob(:,:,:,:)
    !call mpi_reduce(indice_integration,indice_int_glob,NGLLX*NGLLY*NGLLZ*(nnx-1)*(nny-1)*(nnz-1),MPI_INTEGER, MPI_SUM,0,MPI_COMM_WORLD,ier)
    !indice_integration(:,:,:,:)=indice_int_glob(:,:,:,:)
!!$    if (myrank==0) then
!!$       open(155,file='debug_vel.bin',access='direct',recl=4*nx*ny*nz)
!!$       write(155,rec=1)  valeur_integration(3,3,3,:)
!!$       close(155)
!!$       stop
!!$    endif
  end subroutine interpole_one_field


  subroutine project_one_field(ires)

    implicit none
    integer ires,iii
    integer ispec,icell,i,j,k,index_cellule
    integer igll,jgll,kgll
    volume_integration(:,:,:,:)=1._CUSTOM_REAL
    !valeur_integration(:,:,:,:)=0._CUSTOM_REAL
    data_tomo(:,:,:)=0.
    !if (myrank==0) then  !======================== myrank=0
    index_cellule=0
    do ispec=1,nspec_ab
       !write(*,*) myrank, ni(ispec)
       do icell=1,ni(ispec)
          index_cellule=index_cellule+1
          i=indice_grid(1,icell,ispec)
          j=indice_grid(2,icell,ispec)
          k=indice_grid(3,icell,ispec)
          !do k=1,nz
          !do j=1,ny
             !do i=1,nx
                !index_cellule=index_cellule+1


                !do kgll=1,NGLLZ
                   !do jgll=1,NGLLY
                      !do igll=1,NGLLX

                       !  if(indice_integration(igll,jgll,kgll,index_cellule)==0) then
                        !    write(*,*) 'Warning , not found gll : ', igll,jgll,kgll,i,j,k
                        ! else
                         !   valeur_integration(igll,jgll,kgll,index_cellule) = &
                         !        valeur_integration(igll,jgll,kgll,index_cellule) /&
                         !        real(indice_integration(igll,jgll,kgll,index_cellule),CUSTOM_REAL)
                         !endif

                      !enddo
                   !enddo
                !enddo

          xstore_local(1)=x_grid_gll(1,1,1,index_cellule)
          ystore_local(1)=y_grid_gll(1,1,1,index_cellule)
          zstore_local(1)=z_grid_gll(1,1,1,index_cellule)

          xstore_local(2)=x_grid_gll(NGLLX,1,1,index_cellule)
          ystore_local(2)=y_grid_gll(NGLLX,1,1,index_cellule)
          zstore_local(2)=z_grid_gll(NGLLX,1,1,index_cellule)

          xstore_local(3)=x_grid_gll(NGLLX,NGLLY,1,index_cellule)
          ystore_local(3)=y_grid_gll(NGLLX,NGLLY,1,index_cellule)
          zstore_local(3)=z_grid_gll(NGLLX,NGLLY,1,index_cellule)

          xstore_local(4)=x_grid_gll(1,NGLLY,1,index_cellule)
          ystore_local(4)=y_grid_gll(1,NGLLY,1,index_cellule)
          zstore_local(4)=z_grid_gll(1,NGLLY,1,index_cellule)

          xstore_local(5)=x_grid_gll(1,1,NGLLZ,index_cellule)
          ystore_local(5)=y_grid_gll(1,1,NGLLZ,index_cellule)
          zstore_local(5)=z_grid_gll(1,1,NGLLZ,index_cellule)

          xstore_local(6)=x_grid_gll(NGLLX,1,NGLLZ,index_cellule)
          ystore_local(6)=y_grid_gll(NGLLX,1,NGLLZ,index_cellule)
          zstore_local(6)=z_grid_gll(NGLLX,1,NGLLZ,index_cellule)

          xstore_local(7)=x_grid_gll(NGLLX,NGLLY,NGLLZ,index_cellule)
          ystore_local(7)=y_grid_gll(NGLLX,NGLLY,NGLLZ,index_cellule)
          zstore_local(7)=z_grid_gll(NGLLX,NGLLY,NGLLZ,index_cellule)

          xstore_local(8)=x_grid_gll(1,NGLLY,NGLLZ,index_cellule)
          ystore_local(8)=y_grid_gll(1,NGLLY,NGLLZ,index_cellule)
          zstore_local(8)=z_grid_gll(1,NGLLY,NGLLZ,index_cellule)

          ! volumic integral
          call compute_volume_int(volume,xstore_local,ystore_local,zstore_local,wxgll,wygll,wzgll,&
               shape3D,dershape3D,valeur_integration,index_cellule,ngrid_local)
!!$          if (myrank==1) then
!!$             write(*,*)
!!$             write(*,*) valeur_integration(:,:,:,index_cellule)
!!$             do iii=1,8
!!$                write(*,'(3f20.5)') xstore_local(iii),ystore_local(iii),zstore_local(iii)
!!$             enddo
!!$             !write(*,*) ystore_local(1),ystore_local(7)
!!$             !write(*,*) zstore_local(1),zstore_local(7)
!!$             write(*,*) index_cellule
!!$             write(*,*) volume
!!$          endif
          if (ires==1) then ! le volume de la cellule ne change pas, on peut le stocker
             call compute_volume_int(volume_cellule,xstore_local,ystore_local,zstore_local,wxgll,wygll,wzgll,&
                  shape3D,dershape3D,volume_integration,index_cellule,ngrid_local)
             volume=volume/volume_cellule ! divide by volume to retrieve the mean value
          endif

          data_tomo(i,j,k) = volume
          !if (i==1 .and. j==1) then
          !   write(*,*) myrank,index_cellule,volume
          !   write(*,*) xstore_local
          !   write(*,*) ystore_local
          !   write(*,*) zstore_local
          !   write(*,*) valeur_integration(:,:,:,index_cellule)
          !   write(*,*)
          !endif
       enddo
    enddo

    !endif
    vol_data(:,:,:)=0._CUSTOM_REAL
    call mpi_reduce(data_tomo,vol_data,(nnx-1)*(nny-1)*(nnz-1),CUSTOM_MPI_TYPE, MPI_SUM,0,MPI_COMM_WORLD,ier)
    if (myrank==0) then
       do k=1,nz
          do j=1,ny
             do i=1,nx
                data_tomo(i,j,k)=vol_data(i,j,k)/real(indice(i,j,k),CUSTOM_REAL)

             enddo
          enddo
         !write(*,*) k,data_tomo(1,1,k),vol_data(1,1,k),indice(1,1,k)
       enddo
    endif
    call mpi_bcast(data_tomo,nx*ny*nz,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ier)

  end subroutine project_one_field





 subroutine create_chunk_grid_projection()
   implicit none

   integer ispec,ix,iy,iz,igrid_local
   integer imin,imax,jmin,jmax,kmin,kmax
   integer i,j,k,igll,jgll,kgll,ijkgll
   integer ii,jj,kk
   integer index_grid,ierr
   double precision x_center,y_center,z_center
   double precision ANGULAR_WIDTH_ETA_RAD,ANGULAR_WIDTH_XI_RAD,deg2rad
   double precision rx,ry,rz
   double precision profondeur_courante
   double precision, parameter :: R_EARTH=6371000.d0
   double precision xc,yc,zc,pcx,pcy,pcz,pxgll,pygll,pzgll
   double precision ksi00,eta00,pmr00

   deg2rad = 3.141592653589793d0/180.d0

   indice_grille_i(:,:,:,:)=0
   indice_grille_j(:,:,:,:)=0
   indice_grille_k(:,:,:,:)=0
   indice_grille_i1(:,:,:,:)=0
   indice_grille_j1(:,:,:,:)=0
   indice_grille_k1(:,:,:,:)=0
   indice(:,:,:)=0
   indice_glob(:,:,:)=0

   xstore_tomo_grid(:,:,:)=0.d0
   ystore_tomo_grid(:,:,:)=0.d0
   zstore_tomo_grid(:,:,:)=0.d0
   profondeur_tomo_grid(:,:,:)=0.d0
   volume_integration(:,:,:,:)=1._CUSTOM_REAL


!======== on se met dans la sphere cubique et on definit le pas des grilles

   ANGULAR_WIDTH_XI_RAD =  2.d0* dasin(dabs(xmin*1.0000d0)/R_EARTH)
   ANGULAR_WIDTH_ETA_RAD = 2.d0* dasin(dabs(ymin*1.0000d0)/R_EARTH)

   hx = ANGULAR_WIDTH_XI_RAD / (nx)
   hy = ANGULAR_WIDTH_ETA_RAD / (ny)
   hz = (r_earth-r_bot) / (nz)


   if (myrank==0) then
      write(*,*)
      write(*,*) xmin,xmax
      write(*,*) ymin,ymax
      write(*,*) zmin,zmax
      write(*,*) nx,ny,nz
      write(*,*) hx,hy,hz
      write(*,*)
   endif
   !! define interp tri grid (tomo -> SEM)
   do k=1,nz+1
      do j=1,ny+1
         do i=1,nx+1
            rx=1.d0;ry=1.d0;rz=1.d0
            pcz = r_bot + (dble(k)-rz)*hz
            ratio_eta = (dble(j)-ry) / dble(NY)
            yc = 2.d0*ratio_eta-1;yc = tan((ANGULAR_WIDTH_ETA_RAD/2.d0) * yc)
            ratio_xi = (dble(i)-rx) / dble(NX)
            xc = 2.d0*ratio_xi-1;xc = tan((ANGULAR_WIDTH_XI_RAD/2.d0) * xc)
            pcx = xc*pcz; pcy = yc*pcz
            pcz =  -(R_EARTH - pcz/dsqrt(1.d0 + yc*yc + xc*xc) - zmax)
            xstore_tomo_grid_interp_tri(i,j,k)=pcx
            ystore_tomo_grid_interp_tri(i,j,k)=pcy
            zstore_tomo_grid_interp_tri(i,j,k)=pcz
         enddo
      enddo
   enddo

   do ispec = 1, NSPEC_AB   !===== boucle element ispec

      ni(ispec)=0

      xmin_local =  HUGEVAL
      xmax_local = -HUGEVAL
      ymin_local =  HUGEVAL
      ymax_local = -HUGEVAL
      zmin_local =  HUGEVAL
      zmax_local = -HUGEVAL

      !======================= Sommets de l'element ============
      !1
      xstore_local(1)=xstore(ibool(1,1,1,ispec))
      ystore_local(1)=ystore(ibool(1,1,1,ispec))
      zstore_local(1)=zstore(ibool(1,1,1,ispec))
      !2
      xstore_local(2)=xstore(ibool(NGLLX,1,1,ispec))
      ystore_local(2)=ystore(ibool(NGLLX,1,1,ispec))
      zstore_local(2)=zstore(ibool(NGLLX,1,1,ispec))
      !3
      xstore_local(3)=xstore(ibool(NGLLX,NGLLY,1,ispec))
      ystore_local(3)=ystore(ibool(NGLLX,NGLLY,1,ispec))
      zstore_local(3)=zstore(ibool(NGLLX,NGLLY,1,ispec))
      !4
      xstore_local(4)=xstore(ibool(1,NGLLY,1,ispec))
      ystore_local(4)=ystore(ibool(1,NGLLY,1,ispec))
      zstore_local(4)=zstore(ibool(1,NGLLY,1,ispec))
      !5
      xstore_local(5)=xstore(ibool(1,1,NGLLZ,ispec))
      ystore_local(5)=ystore(ibool(1,1,NGLLZ,ispec))
      zstore_local(5)=zstore(ibool(1,1,NGLLZ,ispec))
      !6
      xstore_local(6)=xstore(ibool(NGLLX,1,NGLLZ,ispec))
      ystore_local(6)=ystore(ibool(NGLLX,1,NGLLZ,ispec))
      zstore_local(6)=zstore(ibool(NGLLX,1,NGLLZ,ispec))
      !7
      xstore_local(7)=xstore(ibool(NGLLX,NGLLY,NGLLZ,ispec))
      ystore_local(7)=ystore(ibool(NGLLX,NGLLY,NGLLZ,ispec))
      zstore_local(7)=zstore(ibool(NGLLX,NGLLY,NGLLZ,ispec))
      !8
      xstore_local(8)=xstore(ibool(1,NGLLY,NGLLZ,ispec))
      ystore_local(8)=ystore(ibool(1,NGLLY,NGLLZ,ispec))
      zstore_local(8)=zstore(ibool(1,NGLLY,NGLLZ,ispec))

      !=============== on cherche le pave circonscrit a l'element en se basant sur les point GLL

      do iz=1,NGLLZ
         do iy=1,NGLLY
            do ix=1,NGLLX
               xmin_local=min(xmin_local,xstore(ibool(ix,iy,iz,ispec)))
               xmax_local=max(xmin_local,xstore(ibool(ix,iy,iz,ispec)))
               ymin_local=min(ymin_local,ystore(ibool(ix,iy,iz,ispec)))
               ymax_local=max(ymin_local,ystore(ibool(ix,iy,iz,ispec)))
               zmin_local=min(zmin_local,zstore(ibool(ix,iy,iz,ispec)))
               zmax_local=max(zmax_local,zstore(ibool(ix,iy,iz,ispec)))
            enddo
         enddo
      enddo

      ! =========== calcul du rayon min et max correspondant au pave
      call rayon_min_max(R_EARTH-zmax,xmin_local,xmax_local,ymin_local,ymax_local,zmin_local,zmax_local,pmin,pmax)

      ! =========== calcul des indices min et max de la grille tomo correspondant au pave
      ! indice sur grille fine d'integration (sphere cubique) !! ajout de +/-2 pour deborder sur l'element adjacent
      kmin = 1+ floor((pmin-r_bot)/(r_earth-r_bot)*(nnz-1) )  !- 2
      kmax = 1+ floor((pmax-r_bot)/(r_earth-r_bot)*(nnz-1) )  !+ 2

      call ksi_eta_min_max(R_EARTH-zmax,xmin_local,xmax_local,ymin_local,ymax_local,zmin_local,zmax_local,&
           ksimin,ksimax,etamin,etamax)

      ! indice sur grille fine d'integration (sphere cubique)
      !imin = floor(1. + (NNX - 1) * (ksimin + 0.5d0 *  ANGULAR_WIDTH_XI_RAD)  / ANGULAR_WIDTH_XI_RAD) !- 2
      !jmin = floor(1. + (NNY - 1) * (etamin + 0.5d0 * ANGULAR_WIDTH_ETA_RAD) / ANGULAR_WIDTH_ETA_RAD) !- 2

      imin = 1+ floor( 0.5*(NNX - 1) * (1. + 2.*ksimin / ANGULAR_WIDTH_XI_RAD ))
      jmin = 1+ floor( 0.5*(NNY - 1) * (1. + 2.*etamin / ANGULAR_WIDTH_ETA_RAD))

      ! indice sur grille fine d'integration (sphere cubique)
      !imax = floor(1. + (NNX - 1) * (ksimax + 0.5d0 *  ANGULAR_WIDTH_XI_RAD)  / ANGULAR_WIDTH_XI_RAD) !+ 2
      !jmax = floor(1. + (NNY - 1) * (etamax + 0.5d0 * ANGULAR_WIDTH_ETA_RAD) / ANGULAR_WIDTH_ETA_RAD) !+ 2

      imax = 1+ floor( 0.5*(NNX - 1) * (1. + 2.*ksimax / ANGULAR_WIDTH_XI_RAD ))
      jmax = 1+ floor( 0.5*(NNY - 1) * (1. + 2.*etamax / ANGULAR_WIDTH_ETA_RAD))

      imin=max(imin,1)
      imax=min(imax,nx)
      jmin=max(jmin,1)
      jmax=min(jmax,ny)
      kmin=max(kmin,1)
      kmax=min(kmax,nz)

      do k = kmin, kmax
         pz =  r_bot + (k-0.5d0)*hz  ! valeur au centre de la cellule (sphere cubique)
         profondeur_courante = R_EARTH  - pz

         do j = jmin, jmax
            ratio_eta = (dble(j)-0.5d0) / dble(NY)
            py = 2.d0*ratio_eta-1
            py = tan((ANGULAR_WIDTH_ETA_RAD/2.d0) * py)
            !if (myrank==0) write(*,*) ratio_eta,py
            do i = imin, imax
               ratio_xi = (dble(i)-0.5d0) / dble(NX)
               px = 2.d0*ratio_xi-1
               px = tan((ANGULAR_WIDTH_XI_RAD/2.d0) * px)
               !if (myrank==0) write(*,'(3f30.6)') px,py,pz
               ! on repasse en cartesien pour travailler dans le repere SEM
               x_center = px*pz
               y_center = py*pz
               z_center =  -(R_EARTH - pz/dsqrt(1.d0 + py**2 + px**2) - zmax)
               !if (myrank==0) write(*,*) x_center,y_center,z_center
               !if (myrank==0) &
               !write(*,'(a10,i4,3f30.6,3i4)')'------->:',ispec,x_center,y_center,z_center,&
               !i,j,k

               index_grid = i + (j-1)*nx + (k-1)*nx*ny ! cell ijk

               ! 8 corners
               rx=1.d0;ry=1.d0;rz=1.d0
               pcz = r_bot + (dble(k)-rz)*hz
               ratio_eta = (dble(j)-ry) / dble(NY)
               yc = 2.d0*ratio_eta-1;yc = tan((ANGULAR_WIDTH_ETA_RAD/2.d0) * yc)
               ratio_xi = (dble(i)-rx) / dble(NX)
               xc = 2.d0*ratio_xi-1;xc = tan((ANGULAR_WIDTH_XI_RAD/2.d0) * xc)
               pcx = xc*pcz; pcy = yc*pcz
               pcz =  -(R_EARTH - pcz/dsqrt(1.d0 + yc*yc + xc*xc) - zmax)
               xstore_local_int(1) = pcx
               ystore_local_int(1) = pcy
               zstore_local_int(1) = pcz
               !if (myrank==0) then
                  !write(*,*)
               !   write(*,*)  pcx,pcy,pcz
               !   write(*,*)
               !endif

               rx=0.d0;ry=1.d0;rz=1.d0
               pcz = r_bot + (dble(k)-rz)*hz
               ratio_eta = (dble(j)-ry) / dble(NY)
               yc = 2.d0*ratio_eta-1;yc = tan((ANGULAR_WIDTH_ETA_RAD/2.d0) * yc)
               ratio_xi = (dble(i)-rx) / dble(NX)
               xc = 2.d0*ratio_xi-1;xc = tan((ANGULAR_WIDTH_XI_RAD/2.d0) * xc)
               pcx = xc*pcz; pcy = yc*pcz
               pcz =  -(R_EARTH - pcz/dsqrt(1.d0 + yc*yc + xc*xc) - zmax)
               xstore_local_int(2) = pcx
               ystore_local_int(2) = pcy
               zstore_local_int(2) = pcz

               rx=0.d0;ry=0.d0;rz=1.d0
               pcz = r_bot + (dble(k)-rz)*hz
               ratio_eta = (dble(j)-ry) / dble(NY)
               yc = 2.d0*ratio_eta-1;yc = tan((ANGULAR_WIDTH_ETA_RAD/2.d0) * yc)
               ratio_xi = (dble(i)-rx) / dble(NX)
               xc = 2.d0*ratio_xi-1;xc = tan((ANGULAR_WIDTH_XI_RAD/2.d0) * xc)
               pcx = xc*pcz; pcy = yc*pcz
               pcz =  -(R_EARTH - pcz/dsqrt(1.d0 + yc*yc + xc*xc) - zmax)
               xstore_local_int(3) = pcx
               ystore_local_int(3) = pcy
               zstore_local_int(3) = pcz

               rx=1.d0;ry=0.d0;rz=1.d0
               pcz = r_bot + (dble(k)-rz)*hz
               ratio_eta = (dble(j)-ry) / dble(NY)
               yc = 2.d0*ratio_eta-1;yc = tan((ANGULAR_WIDTH_ETA_RAD/2.d0) * yc)
               ratio_xi = (dble(i)-rx) / dble(NX)
               xc = 2.d0*ratio_xi-1;xc = tan((ANGULAR_WIDTH_XI_RAD/2.d0) * xc)
               pcx = xc*pcz; pcy = yc*pcz
               pcz =  -(R_EARTH - pcz/dsqrt(1.d0 + yc*yc + xc*xc) - zmax)
               xstore_local_int(4) = pcx
               ystore_local_int(4) = pcy
               zstore_local_int(4) = pcz

               rx=1.d0;ry=1.d0;rz=0.d0
               pcz = r_bot + (dble(k)-rz)*hz
               ratio_eta = (dble(j)-ry) / dble(NY)
               yc = 2.d0*ratio_eta-1;yc = tan((ANGULAR_WIDTH_ETA_RAD/2.d0) * yc)
               ratio_xi = (dble(i)-rx) / dble(NX)
               xc = 2.d0*ratio_xi-1;xc = tan((ANGULAR_WIDTH_XI_RAD/2.d0) * xc)
               pcx = xc*pcz; pcy = yc*pcz
               pcz =  -(R_EARTH - pcz/dsqrt(1.d0 + yc*yc + xc*xc) - zmax)
               xstore_local_int(5) = pcx
               ystore_local_int(5) = pcy
               zstore_local_int(5) = pcz

               rx=0.d0;ry=1.d0;rz=0.d0
               pcz = r_bot + (dble(k)-rz)*hz
               ratio_eta = (dble(j)-ry) / dble(NY)
               yc = 2.d0*ratio_eta-1;yc = tan((ANGULAR_WIDTH_ETA_RAD/2.d0) * yc)
               ratio_xi = (dble(i)-rx) / dble(NX)
               xc = 2.d0*ratio_xi-1;xc = tan((ANGULAR_WIDTH_XI_RAD/2.d0) * xc)
               pcx = xc*pcz; pcy = yc*pcz
               pcz =  -(R_EARTH - pcz/dsqrt(1.d0 + yc*yc + xc*xc) - zmax)
               xstore_local_int(6) = pcx
               ystore_local_int(6) = pcy
               zstore_local_int(6) = pcz


               rx=0.d0;ry=0.d0;rz=0.d0
               pcz = r_bot + (dble(k)-rz)*hz
               ratio_eta = (dble(j)-ry) / dble(NY)
               yc = 2.d0*ratio_eta-1;yc = tan((ANGULAR_WIDTH_ETA_RAD/2.d0) * yc)
               ratio_xi = (dble(i)-rx) / dble(NX)
               xc = 2.d0*ratio_xi-1;xc = tan((ANGULAR_WIDTH_XI_RAD/2.d0) * xc)
               pcx = xc*pcz; pcy = yc*pcz
               pcz =  -(R_EARTH - pcz/dsqrt(1.d0 + yc*yc + xc*xc) - zmax)
               xstore_local_int(7) = pcx
               ystore_local_int(7) = pcy
               zstore_local_int(7) = pcz

               rx=1.d0;ry=0.d0;rz=0.d0
               pcz = r_bot + (dble(k)-rz)*hz
               ratio_eta = (dble(j)-ry) / dble(NY)
               yc = 2.d0*ratio_eta-1;yc = tan((ANGULAR_WIDTH_ETA_RAD/2.d0) * yc)
               ratio_xi = (dble(i)-rx) / dble(NX)
               xc = 2.d0*ratio_xi-1;xc = tan((ANGULAR_WIDTH_XI_RAD/2.d0) * xc)
               pcx = xc*pcz; pcy = yc*pcz
               pcz =  -(R_EARTH - pcz/dsqrt(1.d0 + yc*yc + xc*xc) - zmax)
               xstore_local_int(8) = pcx
               ystore_local_int(8) = pcy
               zstore_local_int(8) = pcz

               !if (myrank==0) then
               !    write(*,*) '-----------------'
               !    write(*,*) xstore_local_int
               !    write(*,*)
               !    write(*,*) ystore_local_int
               !    write(*,*)
               !   write(*,*) zstore_local_int
               !    write(*,*)
               !    write(*,*) x_center,y_center,z_center
               !endif
               ! ========== compute GLL points
               ! compute xix, eta et gamma for the central point central x_center,y_center,z_center
               !if (myrank==0) write(*,'(a10,i4,3f30.6)')'------->:',ispec,x_center,y_center,z_center

               call Find_xix_eta_gamma(xstore_local,ystore_local,zstore_local,xi,eta,gamma,x_center,y_center,z_center)
!!$               !! for debbug
!!$               xstore_tomo_grid(i,j,k) = x_center
!!$               ystore_tomo_grid(i,j,k) = y_center
!!$               zstore_tomo_grid(i,j,k) = z_center
               ! are we inside ispec?
               if (dabs(xi)<=1.05d0.and.dabs(eta)<=1.05d0.and.dabs(gamma)<=1.05d0) then
                  indice(i,j,k)=1 ! on a bien visite la cellule de la grille
                  xstore_tomo_grid(i,j,k) = x_center
                  ystore_tomo_grid(i,j,k) = y_center
                  zstore_tomo_grid(i,j,k) = z_center
                  profondeur_tomo_grid(i,j,k) = profondeur_courante !abs(zmax - z_center)
                  indice_spec(index_grid)=ispec !! correspondance cellule_tomo <-> elemements
                  ! la cellule index_cellule est cencee appartenir a ispec et tous ces pts gll aussi.
                  ni(ispec)=ni(ispec)+1
                  indice_grid(1,ni(ispec),ispec)=i
                  indice_grid(2,ni(ispec),ispec)=j
                  indice_grid(3,ni(ispec),ispec)=k

                  ! on calcule les points gll de la cellule :
!!$                  call calcule_point_gll(xgll,ygll,zgll,xstore_local_int,ystore_local_int,zstore_local_int,shape3D)
!!$                  do kgll = 1,NGLLZ
!!$                     do jgll = 1,NGLLY
!!$                        do igll = 1,NGLLX
!!$                           indice_integration(igll,jgll,kgll,index_grid)=1
!!$                           x_grid_gll(igll,jgll,kgll,index_grid) = xgll(igll,jgll,kgll)
!!$                           y_grid_gll(igll,jgll,kgll,index_grid) = ygll(igll,jgll,kgll)
!!$                           z_grid_gll(igll,jgll,kgll,index_grid) = zgll(igll,jgll,kgll)
!!$                            if (igll==3 .and. jgll==3 .and. kgll==3 .and. i==75 .and. j==78 ) then
!!$                               write(*,*)  myrank,xgll(igll,jgll,kgll) , ygll(igll,jgll,kgll), zgll(igll,jgll,kgll)
!!$                            endif
!!$                        enddo
!!$                     enddo
!!$                  enddo


                  !                else
                  !write(*,*)  xi,eta,gamma,x_center,y_center, z_center
                  !                endif

                  ! ========================= on cherche les points GLL de ispec qui sont dans la cellule de la grille tomo
                  call calcule_point_gll(xgll,ygll,zgll,xstore_local,ystore_local,zstore_local,shape3D)
                  !if (myrank==0) write(*,'(a10,i4,3f30.6)')'------->:',ispec,xstore_local(7),ystore_local(7),zstore_local(7)
                 ijkgll=0
                  do kgll = 1,NGLLZ
                     do jgll = 1,NGLLY
                        do igll = 1,NGLLX

                           ijkgll=ijkgll+1

                           pxgll = xgll(igll,jgll,kgll)
                           pygll = ygll(igll,jgll,kgll)
                           pzgll = zgll(igll,jgll,kgll)

                           !if (myrank==0) then
                              !write(*,'(a10,i4,3f30.6)') '------->:',myrank,px,py,pz
                           !   write(*,'(i4,3f30.6)') ijkgll,xgll(igll,jgll,kgll),ygll(igll,jgll,kgll),zgll(igll,jgll,kgll)

                           !endif

                           !call Find_xix_eta_gamma(xstore_local_int,ystore_local_int,zstore_local_int,xi,eta,gamma,pxgll,pygll,pzgll)
                           !if (dabs(xi)<=1.05d0.and.dabs(eta)<=1.05d0.and.dabs(gamma)<=1.05d0) then

                           !! inverser la formule du mapping de la sphere cubique
                           call inv_mapping_sph_cub(pxgll,pygll,pzgll,ksi00,eta00,R_EARTH-zmax)
                           call radius_comp(R_EARTH-zmax, pxgll,pygll,pzgll,pmr00)
                           ii = 1+ floor( 0.5*(NNX - 1) * (1. + 2.*ksi00 / ANGULAR_WIDTH_XI_RAD ))
                           jj = 1+ floor( 0.5*(NNY - 1) * (1. + 2.*eta00 / ANGULAR_WIDTH_ETA_RAD))
                           kk = 1+ floor((pmr00-r_bot)/(r_earth-r_bot)*(nnz-1) )

                           !write(*,*) ii,jj,kk

                           indice_grille_i(igll,jgll,kgll,ispec) = ii
                           indice_grille_i1(igll,jgll,kgll,ispec)= ii+1

                           indice_grille_j(igll,jgll,kgll,ispec) = jj
                           indice_grille_j1(igll,jgll,kgll,ispec)= jj+1

                           indice_grille_k(igll,jgll,kgll,ispec) = kk
                           indice_grille_k1(igll,jgll,kgll,ispec)= kk+1

                            if (indice_grille_i(igll,jgll,kgll,ispec) < 1) indice_grille_i(igll,jgll,kgll,ispec)=1
                            if (indice_grille_j(igll,jgll,kgll,ispec) < 1) indice_grille_j(igll,jgll,kgll,ispec)=1
                            if (indice_grille_k(igll,jgll,kgll,ispec) < 1) indice_grille_k(igll,jgll,kgll,ispec)=1

                            indice_grille_i1(igll,jgll,kgll,ispec) = indice_grille_i(igll,jgll,kgll,ispec) + 1
                            indice_grille_j1(igll,jgll,kgll,ispec) = indice_grille_j(igll,jgll,kgll,ispec) + 1
                            indice_grille_k1(igll,jgll,kgll,ispec) = indice_grille_k(igll,jgll,kgll,ispec) + 1

                            if (indice_grille_i1(igll,jgll,kgll,ispec) > nx+1) indice_grille_i1(igll,jgll,kgll,ispec)=nx+1
                            if (indice_grille_j1(igll,jgll,kgll,ispec) > ny+1) indice_grille_j1(igll,jgll,kgll,ispec)=ny+1
                            if (indice_grille_k1(igll,jgll,kgll,ispec) > nz+1) indice_grille_k1(igll,jgll,kgll,ispec)=nz+1


                            if (indice_grille_i(igll,jgll,kgll,ispec) > nx) indice_grille_i(igll,jgll,kgll,ispec)=nx
                            if (indice_grille_j(igll,jgll,kgll,ispec) > ny) indice_grille_j(igll,jgll,kgll,ispec)=ny
                            if (indice_grille_k(igll,jgll,kgll,ispec) > nz) indice_grille_k(igll,jgll,kgll,ispec)=nz

                            ! check debbug
                            !if (myrank==10) then
                            !   ii=indice_grille_i(igll,jgll,kgll,ispec)
                            !   jj=indice_grille_j(igll,jgll,kgll,ispec)
                            !   kk=indice_grille_k(igll,jgll,kgll,ispec)
                            !   write(*,*) ksi00*180/3.14,ANGULAR_WIDTH_XI_RAD*180/3.14
                            !   write(*,*) eta00*180/3.14,ANGULAR_WIDTH_ETA_RAD*180/3.14
                            !   write(*,*) pmr00,r_bot,pmr00-r_bot
                            !   write(*,*) ii,jj,kk
                            !   write(*,*) xstore_tomo_grid_interp_tri(ii,jj,kk),xstore_tomo_grid_interp_tri(ii+1,jj,kk)
                            !   write(*,*) xstore_tomo_grid_interp_tri(ii,jj+1,kk),xstore_tomo_grid_interp_tri(ii+1,jj+1,kk)
                            !   write(*,*) xstore_tomo_grid_interp_tri(ii,jj,kk+1),xstore_tomo_grid_interp_tri(ii+1,jj,kk+1)
                            !   write(*,*) xstore_tomo_grid_interp_tri(ii,jj+1,kk+1),xstore_tomo_grid_interp_tri(ii+1,jj+1,kk+1)
                            !   write(*,*) pxgll
                            !   write(*,*)
                            !   write(*,*) ystore_tomo_grid_interp_tri(ii,jj,kk),ystore_tomo_grid_interp_tri(ii,jj+1,kk)
                            !   write(*,*) ystore_tomo_grid_interp_tri(ii+1,jj,kk),ystore_tomo_grid_interp_tri(ii+1,jj+1,kk)
                            !   write(*,*) ystore_tomo_grid_interp_tri(ii,jj,kk+1),ystore_tomo_grid_interp_tri(ii,jj+1,kk+1)
                            !   write(*,*) ystore_tomo_grid_interp_tri(ii+1,jj,kk+1),ystore_tomo_grid_interp_tri(ii+1,jj+1,kk+1)
                            !   write(*,*) pygll
                            !   write(*,*)
                            !   write(*,*) zstore_tomo_grid_interp_tri(ii,jj,kk),zstore_tomo_grid_interp_tri(ii,jj+1,kk)
                            !   write(*,*) zstore_tomo_grid_interp_tri(ii+1,jj,kk),zstore_tomo_grid_interp_tri(ii+1,jj+1,kk)
                            !   write(*,*) zstore_tomo_grid_interp_tri(ii,jj,kk+1),zstore_tomo_grid_interp_tri(ii,jj+1,kk+1)
                            !   write(*,*) zstore_tomo_grid_interp_tri(ii+1,jj,kk+1),zstore_tomo_grid_interp_tri(ii+1,jj+1,kk+1)
                            !   write(*,*) pzgll

                            !   write(*,*) '-----'
                            !endif
                               ! je laisse ca mais je ne suis pas encore sur de son utilite
!!$                            if (xi > 0. ) then
!!$                               if (i<nx) indice_grille_i1(igll,jgll,kgll,ispec) = i+1
!!$                            else
!!$                               if (i>1)  indice_grille_i1(igll,jgll,kgll,ispec) = i-1
!!$                            endif
!!$
!!$                            if (eta > 0. ) then
!!$                               if (j<ny) indice_grille_j1(igll,jgll,kgll,ispec) = j+1
!!$                            else
!!$                               if (j>1)  indice_grille_j1(igll,jgll,kgll,ispec) = j-1
!!$                            endif
!!$
!!$                            if (gamma > 0. ) then
!!$                               if (k<nz) indice_grille_k1(igll,jgll,kgll,ispec) = k+1
!!$                            else
!!$                               if (k>1)  indice_grille_k1(igll,jgll,kgll,ispec) = k-1
!!$                            endif


                           !endif
                        enddo
                     enddo
                  enddo

               endif
            enddo
         enddo
      enddo

   enddo

   ! allocation des tableaux locaux
   ngrid_local=sum(ni(:))
   allocate(indice_integration(NGLLX,NGLLY,NGLLZ,ngrid_local))
   allocate(valeur_integration(NGLLX,NGLLY,NGLLZ,ngrid_local))
   allocate(volume_integration(NGLLX,NGLLY,NGLLZ,ngrid_local))
   allocate(x_grid_gll(NGLLX,NGLLY,NGLLZ,ngrid_local))
   allocate(y_grid_gll(NGLLX,NGLLY,NGLLZ,ngrid_local))
   allocate(z_grid_gll(NGLLX,NGLLY,NGLLZ,ngrid_local))


    index_grid=0
    do ispec=1,nspec_ab           ! loop on spectral elements
       do igrid_local=1,ni(ispec) ! loop on grid cells instide ispec
          index_grid=index_grid+1
          ! cell index
          i= indice_grid(1,igrid_local,ispec)
          j= indice_grid(2,igrid_local,ispec)
          k= indice_grid(3,igrid_local,ispec)


          ! 8 corners
          rx=1.d0;ry=1.d0;rz=1.d0
          pcz = r_bot + (dble(k)-rz)*hz
          ratio_eta = (dble(j)-ry) / dble(NY)
          yc = 2.d0*ratio_eta-1;yc = tan((ANGULAR_WIDTH_ETA_RAD/2.d0) * yc)
          ratio_xi = (dble(i)-rx) / dble(NX)
          xc = 2.d0*ratio_xi-1;xc = tan((ANGULAR_WIDTH_XI_RAD/2.d0) * xc)
          pcx = xc*pcz; pcy = yc*pcz
          pcz =  -(R_EARTH - pcz/dsqrt(1.d0 + yc*yc + xc*xc) - zmax)
          xstore_local_int(1) = pcx
          ystore_local_int(1) = pcy
          zstore_local_int(1) = pcz


          rx=0.d0;ry=1.d0;rz=1.d0
          pcz = r_bot + (dble(k)-rz)*hz
          ratio_eta = (dble(j)-ry) / dble(NY)
          yc = 2.d0*ratio_eta-1;yc = tan((ANGULAR_WIDTH_ETA_RAD/2.d0) * yc)
          ratio_xi = (dble(i)-rx) / dble(NX)
          xc = 2.d0*ratio_xi-1;xc = tan((ANGULAR_WIDTH_XI_RAD/2.d0) * xc)
          pcx = xc*pcz; pcy = yc*pcz
          pcz =  -(R_EARTH - pcz/dsqrt(1.d0 + yc*yc + xc*xc) - zmax)
          xstore_local_int(2) = pcx
          ystore_local_int(2) = pcy
          zstore_local_int(2) = pcz

          rx=0.d0;ry=0.d0;rz=1.d0
          pcz = r_bot + (dble(k)-rz)*hz
          ratio_eta = (dble(j)-ry) / dble(NY)
          yc = 2.d0*ratio_eta-1;yc = tan((ANGULAR_WIDTH_ETA_RAD/2.d0) * yc)
          ratio_xi = (dble(i)-rx) / dble(NX)
          xc = 2.d0*ratio_xi-1;xc = tan((ANGULAR_WIDTH_XI_RAD/2.d0) * xc)
          pcx = xc*pcz; pcy = yc*pcz
          pcz =  -(R_EARTH - pcz/dsqrt(1.d0 + yc*yc + xc*xc) - zmax)
          xstore_local_int(3) = pcx
          ystore_local_int(3) = pcy
          zstore_local_int(3) = pcz

          rx=1.d0;ry=0.d0;rz=1.d0
          pcz = r_bot + (dble(k)-rz)*hz
          ratio_eta = (dble(j)-ry) / dble(NY)
          yc = 2.d0*ratio_eta-1;yc = tan((ANGULAR_WIDTH_ETA_RAD/2.d0) * yc)
          ratio_xi = (dble(i)-rx) / dble(NX)
          xc = 2.d0*ratio_xi-1;xc = tan((ANGULAR_WIDTH_XI_RAD/2.d0) * xc)
          pcx = xc*pcz; pcy = yc*pcz
          pcz =  -(R_EARTH - pcz/dsqrt(1.d0 + yc*yc + xc*xc) - zmax)
          xstore_local_int(4) = pcx
          ystore_local_int(4) = pcy
          zstore_local_int(4) = pcz

          rx=1.d0;ry=1.d0;rz=0.d0
          pcz = r_bot + (dble(k)-rz)*hz
          ratio_eta = (dble(j)-ry) / dble(NY)
          yc = 2.d0*ratio_eta-1;yc = tan((ANGULAR_WIDTH_ETA_RAD/2.d0) * yc)
          ratio_xi = (dble(i)-rx) / dble(NX)
          xc = 2.d0*ratio_xi-1;xc = tan((ANGULAR_WIDTH_XI_RAD/2.d0) * xc)
          pcx = xc*pcz; pcy = yc*pcz
          pcz =  -(R_EARTH - pcz/dsqrt(1.d0 + yc*yc + xc*xc) - zmax)
          xstore_local_int(5) = pcx
          ystore_local_int(5) = pcy
          zstore_local_int(5) = pcz

          rx=0.d0;ry=1.d0;rz=0.d0
          pcz = r_bot + (dble(k)-rz)*hz
          ratio_eta = (dble(j)-ry) / dble(NY)
          yc = 2.d0*ratio_eta-1;yc = tan((ANGULAR_WIDTH_ETA_RAD/2.d0) * yc)
          ratio_xi = (dble(i)-rx) / dble(NX)
          xc = 2.d0*ratio_xi-1;xc = tan((ANGULAR_WIDTH_XI_RAD/2.d0) * xc)
          pcx = xc*pcz; pcy = yc*pcz
          pcz =  -(R_EARTH - pcz/dsqrt(1.d0 + yc*yc + xc*xc) - zmax)
          xstore_local_int(6) = pcx
          ystore_local_int(6) = pcy
          zstore_local_int(6) = pcz


          rx=0.d0;ry=0.d0;rz=0.d0
          pcz = r_bot + (dble(k)-rz)*hz
          ratio_eta = (dble(j)-ry) / dble(NY)
          yc = 2.d0*ratio_eta-1;yc = tan((ANGULAR_WIDTH_ETA_RAD/2.d0) * yc)
          ratio_xi = (dble(i)-rx) / dble(NX)
          xc = 2.d0*ratio_xi-1;xc = tan((ANGULAR_WIDTH_XI_RAD/2.d0) * xc)
          pcx = xc*pcz; pcy = yc*pcz
          pcz =  -(R_EARTH - pcz/dsqrt(1.d0 + yc*yc + xc*xc) - zmax)
          xstore_local_int(7) = pcx
          ystore_local_int(7) = pcy
          zstore_local_int(7) = pcz


          rx=1.d0;ry=0.d0;rz=0.d0
          pcz = r_bot + (dble(k)-rz)*hz
          ratio_eta = (dble(j)-ry) / dble(NY)
          yc = 2.d0*ratio_eta-1;yc = tan((ANGULAR_WIDTH_ETA_RAD/2.d0) * yc)
          ratio_xi = (dble(i)-rx) / dble(NX)
          xc = 2.d0*ratio_xi-1;xc = tan((ANGULAR_WIDTH_XI_RAD/2.d0) * xc)
          pcx = xc*pcz; pcy = yc*pcz
          pcz =  -(R_EARTH - pcz/dsqrt(1.d0 + yc*yc + xc*xc) - zmax)
          xstore_local_int(8) = pcx
          ystore_local_int(8) = pcy
          zstore_local_int(8) = pcz

          ! compute gll points of cell :
          call calcule_point_gll(xgll,ygll,zgll,xstore_local_int,ystore_local_int,zstore_local_int,shape3D)
          do kgll = 1,NGLLZ
             do jgll = 1,NGLLY
                do igll = 1,NGLLX
                   !index_grid = index_grid + 1
                   indice_integration(igll,jgll,kgll,index_grid) = 1
                   x_grid_gll(igll,jgll,kgll,index_grid) = xgll(igll,jgll,kgll)
                   y_grid_gll(igll,jgll,kgll,index_grid) = ygll(igll,jgll,kgll)
                   z_grid_gll(igll,jgll,kgll,index_grid) = zgll(igll,jgll,kgll)
                enddo
             enddo
          enddo

       enddo
    enddo

   ! mpi comm
   call mpi_reduce(indice,indice_glob,(nnx-1)*(nny-1)*(nnz-1),MPI_INTEGER ,MPI_SUM,0, MPI_COMM_WORLD, ierr)
   indice(:,:,:) = indice_glob(:,:,:)
   call mpi_bcast(indice,(nnx-1)*(nny-1)*(nnz-1),MPI_INTEGER,0, MPI_COMM_WORLD, ierr)

   call  mpi_reduce(xstore_tomo_grid,work_array,(nnx-1)*(nny-1)*(nnz-1),MPI_DOUBLE_PRECISION ,MPI_SUM,0, MPI_COMM_WORLD, ierr)
   xstore_tomo_grid(:,:,:)=work_array(:,:,:)

   call  mpi_reduce(ystore_tomo_grid,work_array,(nnx-1)*(nny-1)*(nnz-1),MPI_DOUBLE_PRECISION ,MPI_SUM,0, MPI_COMM_WORLD, ierr)
   ystore_tomo_grid(:,:,:)=work_array(:,:,:)

   call  mpi_reduce(zstore_tomo_grid,work_array,(nnx-1)*(nny-1)*(nnz-1),MPI_DOUBLE_PRECISION ,MPI_SUM,0, MPI_COMM_WORLD, ierr)
   zstore_tomo_grid(:,:,:)=work_array(:,:,:)

   call  mpi_reduce(profondeur_tomo_grid,work_array,(nnx-1)*(nny-1)*(nnz-1),MPI_DOUBLE_PRECISION ,MPI_SUM,0, MPI_COMM_WORLD, ierr)
   profondeur_tomo_grid(:,:,:)=work_array(:,:,:)
   call mpi_bcast(profondeur_tomo_grid,nx*ny*nz,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)

!!$   call mpi_reduce(indice_integration,indice_int_glob,NGLLX*NGLLY*NGLLZ*(nnx-1)*(nny-1)*(nnz-1),MPI_INTEGER, MPI_SUM,0,MPI_COMM_WORLD,ier)
!!$   indice_integration(:,:,:,:) = indice_int_glob(:,:,:,:)
!!$
!!$   call  mpi_reduce(x_grid_gll,wk_reduce,NGLLX*NGLLY*NGLLZ*(nnx-1)*(nny-1)*(nnz-1),MPI_DOUBLE_PRECISION ,MPI_SUM,0, MPI_COMM_WORLD, ierr)
!!$   x_grid_gll(:,:,:,:)=wk_reduce(:,:,:,:)
!!$   call mpi_bcast(x_grid_gll,NGLLX*NGLLY*NGLLZ*(nnx-1)*(nny-1)*(nnz-1),MPI_DOUBLE_PRECISION,0, MPI_COMM_WORLD, ierr)
!!$
!!$   call  mpi_reduce(y_grid_gll,wk_reduce,NGLLX*NGLLY*NGLLZ*(nnx-1)*(nny-1)*(nnz-1),MPI_DOUBLE_PRECISION ,MPI_SUM,0, MPI_COMM_WORLD, ierr)
!!$   y_grid_gll(:,:,:,:)=wk_reduce(:,:,:,:)
!!$   call mpi_bcast(y_grid_gll,NGLLX*NGLLY*NGLLZ*(nnx-1)*(nny-1)*(nnz-1),MPI_DOUBLE_PRECISION,0, MPI_COMM_WORLD, ierr)
!!$
!!$   call  mpi_reduce(z_grid_gll,wk_reduce,NGLLX*NGLLY*NGLLZ*(nnx-1)*(nny-1)*(nnz-1),MPI_DOUBLE_PRECISION ,MPI_SUM,0, MPI_COMM_WORLD, ierr)
!!$   z_grid_gll(:,:,:,:)=wk_reduce(:,:,:,:)
!!$   call mpi_bcast(z_grid_gll,NGLLX*NGLLY*NGLLZ*(nnx-1)*(nny-1)*(nnz-1),MPI_DOUBLE_PRECISION,0, MPI_COMM_WORLD, ierr)

   !write(*,*) 'MAX VALUE NI : ', myrank, maxval(ni)
 end subroutine create_chunk_grid_projection



end module project_tomo_grid_mod

