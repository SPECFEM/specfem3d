program re_format_outputs_files

  !! use the same numerber of MPI process than specfem !!!
  implicit none

  include 'mpif.h'
  INTEGER myrank,nbproc,ierr
  integer, dimension(MPI_STATUS_SIZE) :: statut
  integer, parameter :: etq=100
  character(len=500) output_veloc_name(3),output_stress_name(6),output_name

  character(len=500) output_displ_name(3),output_deriv_name(9)

  character(len=500) input_field_name,output_field_name,input_point_file
  character(len=500), allocatable :: working_axisem_dir(:)
  character(len=500) fichier,input_point_file_cart,meshdirectory
  character(len=500) incident_field,incident_field_tmp,tmp_file,prname,LOCAL_PATH,TRACT_PATH

  character(len=256)  :: simtwordtmp, simtvaluetmp, simutypevalue, line
  integer ioerr

  integer iti,a,nb_dump_samples,nbrec,irec,i,ntime,itime,n0,n00
  real, allocatable :: data_time(:,:),data_rec(:,:,:),data_rec0(:,:,:),data_tmp(:,:),data_tmp_to_send(:,:)
  real, allocatable :: ur(:,:),vr(:,:),tr(:,:), dur1(:,:),dur2(:,:),dur3(:,:)
  real nvx,nvy,nvz
  double precision, allocatable :: field_interp(:),field_to_write(:,:,:)
  double precision, allocatable :: field_to_write_all(:,:),buffer_to_write(:),buffer_to_send(:),buffer_to_recv(:)
  double precision, allocatable :: buffer_to_store_dd(:,:), buffer_to_store_v(:,:)
  double complex, allocatable :: zpad(:),padi(:)
  double precision, allocatable :: sinc_tab_displ_and_deriv(:),sinc_tab_veloc(:),sinc_tab_stress(:)
  real, allocatable :: xp(:),yp(:),zp(:)
  integer, allocatable :: igll_glob(:),jgll_glob(:),kgll_glob(:), inum_glob(:),iboun_gll(:)
  integer, allocatable :: nb_received(:),shift(:),nb_received_sv(:),shift_sv(:)
  integer, allocatable :: i_inf(:),i_sup(:),etiquette(:),indx_rec(:,:,:,:,:)
  double precision value_field
  real tt,dtt
  real frq_min,dt
  real tmin,tmax
  integer isim,nsim
  integer iunit, next_iunit
  integer itmin, itmax
  integer irank, nrec_to_store,irec0, nb_point
  double precision lat_src,lon_src,lat_mesh,lon_mesh,azi_rot
  integer, allocatable :: ivx(:),ivy(:),ivz(:), isxx(:),isyy(:),iszz(:),isxy(:),isxz(:),isyz(:), iux(:),iuy(:),iuz(:)
  integer, allocatable :: idu1d1(:),idu1d2(:),idu1d3(:),idu2d1(:),idu2d2(:),idu2d3(:),idu3d1(:),idu3d2(:),idu3d3(:)

  integer*8 iirec
  integer*8 icomp, nbproc0,ntime_interp
  integer ntime_to_store,it,iti0
  integer nb_rec_by_proc,nb_remain_proc,irecmin,irecmax
  integer itnew,ntnew, nSpecfem_proc
  integer ilayer, updown,ispec,code_face
  integer itime_decimate, ntime_decimate,decimate_time
  real current_time_step,current_time_step_half
  integer, allocatable :: flag_boundary_xmin(:),flag_boundary_xmax(:),flag_boundary_ymin(:)&
       ,flag_boundary_ymax(:),flag_boundary_zmin(:),flag_boundary_zmax(:)
  integer n2d_xmin,n2d_xmax,n2d_ymin,n2d_ymax,n2d_zmin,n2d_zmax,i_code_face,ispec2D
  integer max_spec,ispec_glob,max_spec_local,max_spec_glob,iproc,ispec_global
  integer, allocatable :: IndLoc2Glob(:,:),IndSpec_face(:,:)
!
  integer nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax,NSPEC2D_BOTTOM,nspec2D_top,NTIMESTEP
  integer MAX_MUN_ABS_BOUNDARY_FACES,IER,IFACE,IGLL,NGLLSQUARE,J,K,NUM_ABS_BOUNDARY_FACES
  integer KK,is
  integer, allocatable :: IND_REC2FACE(:,:,:),ABS_BOUNDARY_ISPEC(:), &
       ABS_BOUNDARY_IJK(:,:,:),ABS_BOUNDARY_JACOBIAN2DW(:,:),NREC_BY_PROC(:), &
       IND_FACE2REC(:,:,:),NUM_BOUNDARY_BY_PROC(:), &
       irec_glob(:,:)
  real, allocatable ::  abs_boundary_normal(:,:,:)

  logical :: recip_KH_integral
  integer Xk_force


  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nbproc,ierr)

  decimate_time=4
  NGLLSQUARE=25

  allocate(etiquette(nbproc))
  do irank=1,nbproc
     etiquette(irank)=10+irank
  enddo

  write(tmp_file,'(a18,i4.4)') 'trace_orig_new.txt',myrank
  open(399,file=trim(tmp_file))

  write(tmp_file,'(a18,i4.4)') 'trace_inte_new.txt',myrank
  open(499,file=trim(tmp_file))

  write(tmp_file,'(a18,i4.4)') 'trace_refo_new.txt',myrank
  open(599,file=trim(tmp_file))


!###################################################### SERIAL PORCESSING #############################################

  if (myrank == 0) then

     !-----------------------------------  GENERAL INPUT ----------------------------------

!!     open(unit=500, file='inparam_basic', status='old', action='read', iostat=ioerr)
!!
!!     do
!!       read(500, fmt='(a256)', iostat=ioerr) line
!!       if (ioerr < 0) exit
!!       if (len(trim(line)) < 1 .or. line(1:1) == '#') cycle
!!
!!       read(line,*) simtwordtmp, simtvaluetmp
!!       if (simtwordtmp == 'SIMULATION_TYPE') simutypevalue = simtvaluetmp
!!     enddo

!!     close(500)
!!     write(*,*) '--- Lecture simutypevalue OK ---', simutypevalue

     open(10,file='../expand_2D_3D.par')
     read(10,'(a)') input_point_file      !! meshfem3D bounday points (geographic)
     read(10,'(a)') input_point_file_cart !! meshfem3D boundary points (cartessian)
     read(10,*) nbproc0                   !! axisem MPI processes
     read(10,*) lat_src,lon_src           !! axisem source position
     read(10,*) lat_mesh,lon_mesh,azi_rot  !! mesh center position
     !! VM VM add azimuth rotation
     read(10,*) recip_KH_integral
     read(10,*) nsim                      !! AxiSEM simus
     read(10,*) Xk_force                  !! The Xk force source (X, Y, or Z) we simulate for KH and recip
     read(10,*) nSpecfem_proc             !! number of specfem procs
     read(10,'(a)') meshdirectory         !! mesfem3D results
     read(10,'(a)') LOCAL_PATH            !! specfem DATABASE PATH
     read(10,'(a)') TRACT_PATH
     close(10)

     !! count number of global elements
     open(90,file=trim(meshdirectory)//'/Numglob2loc_elmn.txt')
     max_spec_glob=0
     do
        read(90,*,end=89) ispec_glob,ispec,iproc
        max_spec_glob=max(max_spec_glob,ispec_glob)
     enddo
89  close(90)

     open(10,file='../'//trim(input_point_file_cart))
     read(10,*) nb_point

     allocate(xp(nb_point),yp(nb_point),zp(nb_point))
     allocate(inum_glob(nb_point))
     allocate(igll_glob(nb_point),jgll_glob(nb_point),kgll_glob(nb_point))
     allocate(indx_rec(5,5,5,max_spec_glob,6))  !!
     allocate(iboun_gll(nb_point))
     indx_rec=0

     do irec=1,nb_point
        read(10,*) xp(irec),yp(irec),zp(irec), &
             inum_glob(irec),igll_glob(irec),jgll_glob(irec),kgll_glob(irec), &
             iboun_gll(irec),ilayer,updown
        indx_rec(igll_glob(irec),jgll_glob(irec),kgll_glob(irec),inum_glob(irec),iboun_gll(irec))=irec
     enddo
     close(10)

     if (recip_KH_integral) nsim = 1 !! We simulate Xk forces one by one

     open(10,file='../reformat.par')
     read(10,*) frq_min
     read(10,*) tmin,tmax
     close(10)
     dtt=1./frq_min

     !---------------------- lecture des tables de correspondances Specfem glob and loc ---------

     open(90,file=trim(meshdirectory)//'/flags_boundary.txt')  ! table de correspondance ispec2D < - > ispec

     open(91,file=trim(meshdirectory)//'/Nb_ielm_faces.txt')   ! taille tableaux pour allocation memoire
     read(91,*) n2d_xmin
     read(91,*) n2d_xmax
     read(91,*) n2d_ymin
     read(91,*) n2d_ymax
     read(91,*) n2d_zmin
     read(91,*) n2d_zmax
     close(91)
     nspec2d_xmin=n2d_xmin
     nspec2d_xmax=n2d_xmax
     nspec2d_ymin=n2d_ymin
     nspec2d_ymax=n2d_ymax
     nspec2d_bottom=n2d_zmin
     nspec2D_top=n2d_zmax

     allocate(flag_boundary_xmin(n2d_xmin))
     allocate(flag_boundary_xmax(n2d_xmax))
     allocate(flag_boundary_ymin(n2d_ymin))
     allocate(flag_boundary_ymax(n2d_ymax))
     allocate(flag_boundary_zmin(n2d_zmin))
     allocate(flag_boundary_zmax(n2d_zmax))

     max_spec = 0

     do
        read(90,*,end=100) ispec,ispec2D,code_face
        if (code_face == 1) flag_boundary_xmin(ispec2D)=ispec
        if (code_face == 2) flag_boundary_xmax(ispec2D)=ispec
        if (code_face == 3) flag_boundary_ymin(ispec2D)=ispec
        if (code_face == 4) flag_boundary_ymax(ispec2D)=ispec
        if (code_face == 5) flag_boundary_zmin(ispec2D)=ispec
        if (code_face == 6) flag_boundary_zmax(ispec2D)=ispec
        max_spec = max(max_spec,ispec)
     enddo
100  close(90)

     allocate(IndSpec_face(max_spec,6))
     IndSpec_face(:,:) = 0
     open(90,file=trim(meshdirectory)//'/flags_boundary.txt')  ! table de correspondance ispec2D < - > ispec

     do
        read(90,*,end=101) ispec,ispec2D,code_face
        IndSpec_face(ispec,code_face)=ispec2D
     enddo
101  close(90)

     open(90,file=trim(meshdirectory)//'/Numglob2loc_elmn.txt')
     max_spec_local=0
     do
        read(90,*,end=102) ispec_glob,ispec,iproc
        max_spec_local=max(max_spec_local,ispec)
     enddo
102  close(90)

     allocate(IndLoc2Glob(max_spec_local,nbproc))
     IndLoc2Glob(:,:) = 0
     open(90,file=trim(meshdirectory)//'/Numglob2loc_elmn.txt')
     max_spec_local=0

     do
        read(90,*,end=103) ispec_glob,ispec,iproc
        IndLoc2Glob(ispec,iproc+1) = ispec_glob
     enddo
103  close(90)

     allocate(num_boundary_by_proc(nSpecfem_proc))
     max_mun_abs_boundary_faces=0
     do iproc = 1, nSpecfem_proc
        call create_name_database(prname,iproc-1,LOCAL_PATH)
        write(*,*) prname(1:len_trim(prname))//'absorb_dsm'
        open(27,file=prname(1:len_trim(prname))//'absorb_dsm',status='old', &
             action='read',form='unformatted',iostat=ier)
        if (ier /= 0)  write(*,*) 'error opening',prname(1:len_trim(prname))//'absorb_dsm'

        read(27) num_abs_boundary_faces
        num_boundary_by_proc(iproc)=num_abs_boundary_faces
        max_mun_abs_boundary_faces=max(num_abs_boundary_faces,max_mun_abs_boundary_faces)
        close(27)
     enddo

     allocate(ind_rec2face(2,nb_point, nSpecfem_proc))
     allocate(ind_face2rec(max_mun_abs_boundary_faces,25,nSpecfem_proc))
     allocate( nrec_by_proc( nSpecfem_proc))
     ind_face2rec=0
     ind_rec2face=0
     do iproc = 1, nSpecfem_proc
        call create_name_database(prname,iproc-1,LOCAL_PATH)
        write(*,*) prname(1:len_trim(prname))//'absorb_dsm'
        open(27,file=prname(1:len_trim(prname))//'absorb_dsm',status='old', &
             action='read',form='unformatted',iostat=ier)
        read(27) num_abs_boundary_faces
        write(*,*) num_abs_boundary_faces
        allocate(abs_boundary_ispec(num_abs_boundary_faces))
        allocate(abs_boundary_ijk(3,NGLLSQUARE,num_abs_boundary_faces))
        allocate(abs_boundary_jacobian2Dw(NGLLSQUARE,num_abs_boundary_faces))
        allocate(abs_boundary_normal(3,NGLLSQUARE,num_abs_boundary_faces))

        read(27) abs_boundary_ispec
        read(27) abs_boundary_ijk
        read(27) abs_boundary_jacobian2Dw
        read(27) abs_boundary_normal
        close(27)

        nrec_by_proc(iproc)=0

        do iface = 1, num_abs_boundary_faces


           ispec = abs_boundary_ispec(iface)
           ispec_glob=IndLoc2Glob(ispec,iproc)

           igll=13
           code_face=0
           if (abs_boundary_ijk(1,igll,iface) == 1) code_face=1
           if (abs_boundary_ijk(1,igll,iface) == 5) code_face=2
           if (abs_boundary_ijk(2,igll,iface) == 1) code_face=3
           if (abs_boundary_ijk(2,igll,iface) == 5) code_face=4
           if (abs_boundary_ijk(3,igll,iface) == 1) code_face=5
           if (abs_boundary_ijk(3,igll,iface) == 5) code_face=6
           if (code_face == 0) then
              write(*,*) 'wrong face '
              stop
           endif
           do igll = 1, NGLLSQUARE
              i = abs_boundary_ijk(1,igll,iface)
              j = abs_boundary_ijk(2,igll,iface)
              k = abs_boundary_ijk(3,igll,iface)
              !write(*,*) i,j,k,ispec_glob
              irec=indx_rec(i,j,k,ispec_glob,code_face)  !! pb: i,j,k,ispec_glob a plusieurs irec possible (ie arrete ou coin)
              if (irec == 0) then                !! il faut aussi indiquer la face qu'on regarde
                 write(*,*) 'point not found ',i,j,k,ispec_glob
                 write(*,*) 'ispec ', ispec
                 !irec=1
                 stop
              endif

              !           bijection
              !!  (irec,iproc) -> (iface,igll)
              !!  (iface,igll) -> (irec,iproc)
              !

              ind_rec2face(1,irec,iproc)=iface  !! un meme irec peut appartenir a 1 , 2 ou 3 faces
              ind_rec2face(2,irec,iproc)=igll   !! la il peut y avoir un pb ??? on peut se tromper de iface
              ind_face2rec(iface,igll,iproc)=irec
              nrec_by_proc(iproc) =  nrec_by_proc(iproc) + 1

           enddo


        enddo
        deallocate(abs_boundary_ispec)
        deallocate(abs_boundary_ijk)
        deallocate(abs_boundary_jacobian2Dw)
        deallocate(abs_boundary_normal)

     enddo


     ! -----  AxiSEM stuff ---------------------------------------------------------------------

     if (.not. recip_KH_integral) then

       allocate(working_axisem_dir(nsim))

       if (nsim == 1) then
         working_axisem_dir(1) = "./"
!!         if (simutypevalue == 'single') working_axisem_dir(1) = "./"
!!         if (simutypevalue == 'force')  working_axisem_dir(1) = "PX/"
       else if (nsim == 2) then
         working_axisem_dir(1) = "PZ/"
         working_axisem_dir(2) = "PX/"
       else
         working_axisem_dir(1) = "MZZ/"
         working_axisem_dir(2) = "MXX_P_MYY/"
         working_axisem_dir(3) = "MXZ_MYZ/"
         working_axisem_dir(4) = "MXY_MXX_M_MYY/"
       endif

       !open(10,file=trim( working_axisem_dir(1))//'Data/strain_info.dat0000')
       !read(10,*) nb_dump_samples
       !read(10,*) dt,i
       !close(10)
       open(10,file='info_for_specefm.txt')
       read(10,*) dt
       close(10)

       output_veloc_name(1)='velocityoutp_u1'
       output_veloc_name(2)='velocityoutp_u2'
       output_veloc_name(3)='velocityoutp_u3'

       output_stress_name(1)='stress_Sg11_out'
       output_stress_name(2)='stress_Sg22_out'
       output_stress_name(3)='stress_Sg33_out'
       output_stress_name(4)='stress_Sg12_out'
       output_stress_name(5)='stress_Sg13_out'
       output_stress_name(6)='stress_Sg23_out'

       iunit=6666
       allocate(ivx(nsim),ivy(nsim),ivz(nsim))
       allocate(isxx(nsim),isyy(nsim),iszz(nsim))
       allocate(isxy(nsim),isxz(nsim),isyz(nsim))

!!       write(*,*) working_axisem_dir(1), nsim, simutypevalue

       do isim=1,nsim

         ivx(isim)=next_iunit(iunit)
         ivy(isim)=next_iunit(iunit)
         ivz(isim)=next_iunit(iunit)
         isxx(isim)=next_iunit(iunit)
         isyy(isim)=next_iunit(iunit)
         iszz(isim)=next_iunit(iunit)
         isxy(isim)=next_iunit(iunit)
         isxz(isim)=next_iunit(iunit)
         isyz(isim)=next_iunit(iunit)

         write(fichier,'(a6,a15)') '/Data/',output_veloc_name(1)
         open(ivx(isim),file= trim(working_axisem_dir(isim))//trim(fichier), FORM="UNFORMATTED")
         write(fichier,'(a6,a15)') '/Data/',output_veloc_name(2)
         open(ivy(isim),file= trim(working_axisem_dir(isim))//trim(fichier), FORM="UNFORMATTED")
         write(fichier,'(a6,a15)') '/Data/',output_veloc_name(3)

         open(ivz(isim),file= trim(working_axisem_dir(isim))//trim(fichier), FORM="UNFORMATTED")
         write(fichier,'(a6,a15)') '/Data/',output_stress_name(1)
         open(isxx(isim),file= trim(working_axisem_dir(isim))//trim(fichier), FORM="UNFORMATTED")
         write(fichier,'(a6,a15)') '/Data/',output_stress_name(2)
         open(isyy(isim),file= trim(working_axisem_dir(isim))//trim(fichier), FORM="UNFORMATTED")
         write(fichier,'(a6,a15)') '/Data/',output_stress_name(3)
         open(iszz(isim),file= trim(working_axisem_dir(isim))//trim(fichier), FORM="UNFORMATTED")
         write(fichier,'(a6,a15)') '/Data/',output_stress_name(4)
         open(isxy(isim),file= trim(working_axisem_dir(isim))//trim(fichier), FORM="UNFORMATTED")
         write(fichier,'(a6,a15)') '/Data/',output_stress_name(5)
         open(isxz(isim),file= trim(working_axisem_dir(isim))//trim(fichier), FORM="UNFORMATTED")
         write(fichier,'(a6,a15)') '/Data/',output_stress_name(6)
         open(isyz(isim),file= trim(working_axisem_dir(isim))//trim(fichier), FORM="UNFORMATTED")

         write(*,*) 'openning ', trim(working_axisem_dir(isim))//trim(fichier)

       enddo

       do isim=1,nsim
         read(ivx(isim))  nbrec,ntime
         read(ivy(isim))  nbrec,ntime
         read(ivz(isim))  nbrec,ntime
         read(isxx(isim)) nbrec,ntime
         read(isyy(isim)) nbrec,ntime
         read(iszz(isim)) nbrec,ntime
         read(isxy(isim)) nbrec,ntime
         read(isxz(isim)) nbrec,ntime
         read(isyz(isim)) nbrec,ntime
       enddo

       write(*,*) ' time step ', dtt
       write(*,*) 'READING OK',ntime,nbrec

     else !! CD CD for KH

       nsim = 1 !! We simulate Xk forces one by one
       write(*,*) 'With reciprocity and KH integral, nsim have always to be 1'
       allocate(working_axisem_dir(nsim))

       working_axisem_dir(1) = "./"
!!       if (simutypevalue == 'single') working_axisem_dir(1) = "./"
!!       if (simutypevalue == 'force')  working_axisem_dir(1) = "PX/"

       open(10,file='info_for_specefm.txt')
       read(10,*) dt
       close(10)

       output_displ_name(1)='displ_out_u1'
       output_displ_name(2)='displ_out_u2'
       output_displ_name(3)='displ_out_u3'

       output_deriv_name(1)='deriv_out_du1d1' !! CD CD add this
       output_deriv_name(2)='deriv_out_du1d2'
       output_deriv_name(3)='deriv_out_du1d3'
       output_deriv_name(4)='deriv_out_du2d1'
       output_deriv_name(5)='deriv_out_du2d2'
       output_deriv_name(6)='deriv_out_du2d3'
       output_deriv_name(7)='deriv_out_du3d1'
       output_deriv_name(8)='deriv_out_du3d2'
       output_deriv_name(9)='deriv_out_du3d3'

       iunit=6666

       allocate(iux(nsim),iuy(nsim),iuz(nsim))

       allocate(idu1d1(nsim),idu1d2(nsim),idu1d3(nsim))
       allocate(idu2d1(nsim),idu2d2(nsim),idu2d3(nsim))
       allocate(idu3d1(nsim),idu3d2(nsim),idu3d3(nsim))

!!       write(*,*) working_axisem_dir(1), nsim, simutypevalue

       do isim=1,nsim

         iux(isim)=next_iunit(iunit)
         iuy(isim)=next_iunit(iunit)
         iuz(isim)=next_iunit(iunit)

         idu1d1(isim)=next_iunit(iunit)
         idu1d2(isim)=next_iunit(iunit)
         idu1d3(isim)=next_iunit(iunit)
         idu2d1(isim)=next_iunit(iunit)
         idu2d2(isim)=next_iunit(iunit)
         idu2d3(isim)=next_iunit(iunit)
         idu3d1(isim)=next_iunit(iunit)
         idu3d2(isim)=next_iunit(iunit)
         idu3d3(isim)=next_iunit(iunit)


         write(fichier,'(a15)') output_displ_name(1)
         open(iux(isim),file= trim(working_axisem_dir(isim))//trim(fichier), FORM="UNFORMATTED")
         write(fichier,'(a15)') output_displ_name(2)
         open(iuy(isim),file= trim(working_axisem_dir(isim))//trim(fichier), FORM="UNFORMATTED")
         write(fichier,'(a15)') output_displ_name(3)
         open(iuz(isim),file= trim(working_axisem_dir(isim))//trim(fichier), FORM="UNFORMATTED")

         write(fichier,'(a15)') output_deriv_name(1)
         open(idu1d1(isim),file= trim(working_axisem_dir(isim))//trim(fichier), FORM="UNFORMATTED")
         write(fichier,'(a15)') output_deriv_name(2)
         open(idu1d2(isim),file= trim(working_axisem_dir(isim))//trim(fichier), FORM="UNFORMATTED")
         write(fichier,'(a15)') output_deriv_name(3)
         open(idu1d3(isim),file= trim(working_axisem_dir(isim))//trim(fichier), FORM="UNFORMATTED")
         write(fichier,'(a15)') output_deriv_name(4)
         open(idu2d1(isim),file= trim(working_axisem_dir(isim))//trim(fichier), FORM="UNFORMATTED")
         write(fichier,'(a15)') output_deriv_name(5)
         open(idu2d2(isim),file= trim(working_axisem_dir(isim))//trim(fichier), FORM="UNFORMATTED")
         write(fichier,'(a15)') output_deriv_name(6)
         open(idu2d3(isim),file= trim(working_axisem_dir(isim))//trim(fichier), FORM="UNFORMATTED")
         write(fichier,'(a15)') output_deriv_name(7)
         open(idu3d1(isim),file= trim(working_axisem_dir(isim))//trim(fichier), FORM="UNFORMATTED")
         write(fichier,'(a15)') output_deriv_name(8)
         open(idu3d2(isim),file= trim(working_axisem_dir(isim))//trim(fichier), FORM="UNFORMATTED")
         write(fichier,'(a15)') output_deriv_name(9)
         open(idu3d3(isim),file= trim(working_axisem_dir(isim))//trim(fichier), FORM="UNFORMATTED")

         write(*,*) 'openning ', trim(working_axisem_dir(isim))//trim(fichier)

       enddo

       do isim=1,nsim

         read(iux(isim))  nbrec,ntime
         read(iuy(isim))  nbrec,ntime
         read(iuz(isim))  nbrec,ntime

         read(idu1d1(isim)) nbrec,ntime
         read(idu1d2(isim)) nbrec,ntime
         read(idu1d3(isim)) nbrec,ntime
         read(idu2d1(isim)) nbrec,ntime
         read(idu2d2(isim)) nbrec,ntime
         read(idu2d3(isim)) nbrec,ntime
         read(idu3d1(isim)) nbrec,ntime
         read(idu3d2(isim)) nbrec,ntime
         read(idu3d3(isim)) nbrec,ntime

       enddo

       write(*,*) ' time step ', dtt
       write(*,*) 'READING OK',ntime,nbrec

     endif

  endif  !! if (myrank == 0)

!###########################################  END SERIAL PROCESS ##################


  call mpi_bcast(ntime,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(nbrec,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(nSpecfem_proc,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(MAX_MUN_ABS_BOUNDARY_FACES,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  if (myrank > 0) allocate(nrec_by_proc(nSpecfem_proc))
  call mpi_bcast(nrec_by_proc,nSpecfem_proc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

  call mpi_bcast(nb_point,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  if (myrank > 0) allocate(ind_rec2face(2,nb_point, nSpecfem_proc))
  call mpi_bcast(ind_rec2face,2*nb_point*nSpecfem_proc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

  call mpi_bcast(nsim,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(n0,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(n00,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(ntime_interp,2,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(dtt,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(dt,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(frq_min,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(tmin,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(tmax,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)

  call mpi_bcast(LOCAL_PATH,500,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(TRACT_PATH,500,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)

  call mpi_bcast(recip_KH_integral,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)

  ! not useful ??
  if (myrank > 0) then
    allocate(working_axisem_dir(nsim))

    if (nsim == 1) then
      working_axisem_dir(1) = "./"
!!      if (simutypevalue == 'single') working_axisem_dir(1) = "./"
!!      if (simutypevalue == 'force')  working_axisem_dir(1) = "PX/"
    else if (nsim == 2) then
      working_axisem_dir(1) = "PZ/"
      working_axisem_dir(2) = "PX/"
    else
      working_axisem_dir(1) = "MZZ/"
      working_axisem_dir(2) = "MXX_P_MYY/"
      working_axisem_dir(3) = "MXZ_MYZ/"
      working_axisem_dir(4) = "MXY_MXX_M_MYY/"
    endif
  endif

  itmin=tmin / dtt
  itmax=tmax / dtt

  if (MOD(itmax - itmin + 1,2) > 0) ITMIN=ITMIN+1
  ntime_interp = itmax - itmin + 1

  write(*,*) ' dt ', dt
  write(*,*) ' itmin ',itmin
  write(*,*) ' itmax ', itmax
  write(*,*) ' nb ',itmax - itmin + 1
  write(*,*) ' ntime_interp ',ntime_interp

  if (.not. recip_KH_integral) then

    allocate(data_rec(9,ntime,nrec_by_proc(myrank+1)))
    allocate(data_rec0(1,nbrec,9),data_tmp(nbrec,9))
    allocate(data_tmp_to_send(9,nbrec))

  else

    allocate(data_rec(12,ntime,nrec_by_proc(myrank+1)))
    allocate(data_rec0(1,nbrec,12),data_tmp(nbrec,12))
    allocate(data_tmp_to_send(12,nbrec))

  endif

  call create_name_database(prname,myrank,LOCAL_PATH)
  write(*,*) prname(1:len_trim(prname))//'absorb_dsm'
  open(27,file=prname(1:len_trim(prname))//'absorb_dsm',status='old', &
       action='read',form='unformatted',iostat=ier)
  read(27) num_abs_boundary_faces
  write(*,*) num_abs_boundary_faces
  allocate(abs_boundary_ispec(num_abs_boundary_faces))
  allocate(abs_boundary_ijk(3,NGLLSQUARE,num_abs_boundary_faces))
  allocate(abs_boundary_jacobian2Dw(NGLLSQUARE,num_abs_boundary_faces))
  allocate(abs_boundary_normal(3,NGLLSQUARE,num_abs_boundary_faces))

  read(27) abs_boundary_ispec
  read(27) abs_boundary_ijk
  read(27) abs_boundary_jacobian2Dw
  read(27) abs_boundary_normal
  close(27)


! ################################ reading and scatter the data  ################################

  if (.not. recip_KH_integral) then

!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Case of classic coupling !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!

     allocate(irec_glob(NGLLSQUARE*MAX_MUN_ABS_BOUNDARY_FACES,nSpecfem_proc))
     allocate(vr(3,nrec_by_proc(myrank+1)),tr(3, nrec_by_proc(myrank+1)))

     do itime=1,ntime

        if (myrank == 0) then
           if (mod(itime,100) == 0) write(*,*) 'reading ', itime, ' / ',ntime
           data_rec0=0.

           do isim=1,nsim

              read(ivx(isim))  data_tmp(:,1)
              read(ivy(isim))  data_tmp(:,2)
              read(ivz(isim))  data_tmp(:,3)
              read(isxx(isim)) data_tmp(:,4)
              read(isyy(isim)) data_tmp(:,5)
              read(iszz(isim)) data_tmp(:,6)
              read(isxy(isim)) data_tmp(:,7)
              read(isxz(isim)) data_tmp(:,8)
              read(isyz(isim)) data_tmp(:,9)

              data_rec0(1,:,:)= data_rec0(1,:,:)+data_tmp(:,:)

           enddo
           !write(*,*) itime,data_rec0(1,100,3)
        endif


        !!scatter the data into the rigth MPI partition
        do icomp =1, 9

           if (myrank == 0) then

              irank=0
              kk=0
              do iface=1,num_boundary_by_proc(irank+1)
                 do igll=1,NGLLSQUARE
                    kk=kk+1
                    irec=ind_face2rec(iface,igll,irank+1)
                    irec_glob(kk,irank+1)=irec
                    if ( irec == 0) then
                       write(*,*) ' irec ', irec,iface,igll,irank+1
                       stop
                    endif
                    data_rec(icomp,itime,kk)=data_rec0(1,irec,icomp)
                 enddo
              enddo

              !if (icomp==3) write(*,*) itime,data_rec(icomp,itime,100),data_rec0(1,100,icomp)

              do irank=1,nbproc-1
                 kk=0
                 do iface=1,num_boundary_by_proc(irank+1)
                    do igll=1,NGLLSQUARE
                       kk=kk+1
                       irec=ind_face2rec(iface,igll,irank+1)
                       irec_glob(kk,irank+1)=irec
                       if ( irec == 0) then
                          write(*,*) ' irec ', irec,iface,igll,irank+1
                          stop
                       endif
                       data_tmp_to_send(icomp,kk) = data_rec0(1,irec,icomp)
                    enddo
                 enddo
                 call mpi_send(data_tmp_to_send(icomp,1:kk),nrec_by_proc(irank+1),MPI_REAL,irank,etq,MPI_COMM_WORLD,ierr)
              enddo

           else
              call mpi_recv(data_rec(icomp,itime,:),nrec_by_proc(myrank+1),MPI_REAL,0,etq,MPI_COMM_WORLD,statut,ierr)
!!$           if (myrank==1) write(*,*) itime,data_rec(3,itime,100)
           endif

        enddo  !! icomp
     enddo     !! itime
     call mpi_bcast(irec_glob,NGLLSQUARE*MAX_MUN_ABS_BOUNDARY_FACES*nSpecfem_proc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!!!


     if (myrank == 0) then
        do isim=1,nsim
           close(ivx(isim))
           close(ivy(isim))
           close(ivz(isim))
           close(isxx(isim))
           close(isyy(isim))
           close(iszz(isim))
           close(isxy(isim))
           close(isxz(isim))
           close(isyz(isim))
        enddo
     endif

     !############################ MPI INTERPOLATION ########################


     nrec_to_store = nrec_by_proc(myrank+1)
     allocate(buffer_to_store_v(9,nrec_to_store))
     allocate(sinc_tab_veloc(ntime),sinc_tab_stress(ntime))
     irecmin=1
     irecmax=nrec_to_store
     ntime_to_store = (itmax - itmin + 1) / 2


     call create_name_database(prname,myrank,TRACT_PATH)
     write(*,*) TRACT_PATH,prname,myrank

     open(28,file=prname(1:len_trim(prname))//'sol_axisem',status='unknown', &
          action='write',form='unformatted',iostat=ier)
     if (ier /= 0) write(*,*) 'error opening', prname(1:len_trim(prname))//'sol_axisem'

     if (nrec_to_store >= 100) then
        do i=1,ntime
           write(399,*) (i-1)*dt,data_rec(3,i,100)
        enddo
     endif

     close(399)

     ntnew=ntime_interp

     if (myrank == 0) write(*,*) 'time step ', dtt,ntnew

     tt=-dtt
     do itnew=1,ntnew
        !write(*,*) itnew,ntnew
        tt=tt+dtt
        current_time_step = tmin + (itnew-1)*dtt !+ 0.5*dtt
        !current_time_step_half = current_time_step - 0.5*dtt

        sinc_tab_veloc(:)      = 0.d0
        buffer_to_store_v(:,:) = 0.d0
        !sinc_tab_stress(:)=0.d0

        call compute_sinc(sinc_tab_veloc,current_time_step,ntime,dt)
        !call compute_sinc(sinc_tab_stress,current_time_step,ntime,dt)

        if (myrank == 0 .and. mod(itnew,100) == 0) write(*,*) myrank, tt,itnew,ntnew

        do irec0=irecmin, irecmax

           do is = 1, ntime

              buffer_to_store_v(1,irec0)=buffer_to_store_v(1,irec0) + sinc_tab_veloc(is) * dble(data_rec(1,is,irec0))
              buffer_to_store_v(2,irec0)=buffer_to_store_v(2,irec0) + sinc_tab_veloc(is) * dble(data_rec(2,is,irec0))
              buffer_to_store_v(3,irec0)=buffer_to_store_v(3,irec0) + sinc_tab_veloc(is) * dble(data_rec(3,is,irec0))
              buffer_to_store_v(4,irec0)=buffer_to_store_v(4,irec0) + sinc_tab_veloc(is) * dble(data_rec(4,is,irec0))
              buffer_to_store_v(5,irec0)=buffer_to_store_v(5,irec0) + sinc_tab_veloc(is) * dble(data_rec(5,is,irec0))
              buffer_to_store_v(6,irec0)=buffer_to_store_v(6,irec0) + sinc_tab_veloc(is) * dble(data_rec(6,is,irec0))
              buffer_to_store_v(7,irec0)=buffer_to_store_v(7,irec0) + sinc_tab_veloc(is) * dble(data_rec(7,is,irec0))
              buffer_to_store_v(8,irec0)=buffer_to_store_v(8,irec0) + sinc_tab_veloc(is) * dble(data_rec(8,is,irec0))
              buffer_to_store_v(9,irec0)=buffer_to_store_v(9,irec0) + sinc_tab_veloc(is) * dble(data_rec(9,is,irec0))

           enddo

           !call interpol_sinc(data_rec(:,irec0,1),buffer_to_store(irec0,1),current_time_step,ntime,dt,sinc_tab_veloc)
           !call interpol_sinc(data_rec(:,irec0,2),buffer_to_store(irec0,2),current_time_step,ntime,dt,sinc_tab_veloc)
           !call interpol_sinc(data_rec(:,irec0,3),buffer_to_store(irec0,3),current_time_step,ntime,dt,sinc_tab_veloc)

           if (irec0 == 100) write(499,*)  current_time_step, buffer_to_store_v(3,irec0)

           !call interpol_sinc(data_rec(:,irec0,4),buffer_to_store(irec0,4),current_time_step,ntime,dt,sinc_tab_stress!)
           !call interpol_sinc(data_rec(:,irec0,5),buffer_to_store(irec0,5),current_time_step,ntime,dt,sinc_tab_stress)
           !call interpol_sinc(data_rec(:,irec0,6),buffer_to_store(irec0,6),current_time_step,ntime,dt,sinc_tab_stress)
           !call interpol_sinc(data_rec(:,irec0,7),buffer_to_store(irec0,7),current_time_step,ntime,dt,sinc_tab_stress)
           !call interpol_sinc(data_rec(:,irec0,8),buffer_to_store(irec0,8),current_time_step,ntime,dt,sinc_tab_stress)
           !call interpol_sinc(data_rec(:,irec0,9),buffer_to_store(irec0,9),current_time_step,ntime,dt,sinc_tab_stress)

           ! compute traction

           irec =irec_glob(irec0,myrank+1)

           if (irec == 0 ) then
              write(*,*) myrank , irec0, irec
              stop
           endif

           iface=ind_rec2face(1,irec,myrank+1)
           igll =ind_rec2face(2,irec,myrank+1)

           nvx=abs_boundary_normal(1,igll,iface)
           nvy=abs_boundary_normal(2,igll,iface)
           nvz=abs_boundary_normal(3,igll,iface)

           !!
           if (myrank == 0) then

           endif
           !!

           vr(1,irec0) = buffer_to_store_v(1,irec0)
           vr(2,irec0) = buffer_to_store_v(2,irec0)
           vr(3,irec0) = buffer_to_store_v(3,irec0)

           tr(1,irec0) = buffer_to_store_v(4,irec0)*nvx + buffer_to_store_v(7,irec0) * nvy + buffer_to_store_v(8,irec0) * nvz
           tr(2,irec0) = buffer_to_store_v(7,irec0)*nvx + buffer_to_store_v(5,irec0) * nvy + buffer_to_store_v(9,irec0) * nvz
           tr(3,irec0) = buffer_to_store_v(8,irec0)*nvx + buffer_to_store_v(9,irec0) * nvy + buffer_to_store_v(6,irec0) * nvz

        enddo

        write(28) vr,tr

     enddo

     close(28)
     close(499)

     write(*,*) 'nbrec ', myrank, irecmax-irecmin+1, ntime_interp
     call MPI_Barrier(MPI_COMM_WORLD,ierr)
     write(*,*) myrank,ierr
     call MPI_FINALIZE(ierr)
     write(*,*) myrank,ierr
     stop

!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Case of recip & KH !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!

  else

    allocate(irec_glob(NGLLSQUARE*MAX_MUN_ABS_BOUNDARY_FACES,nSpecfem_proc))
    allocate(ur(3,nrec_by_proc(myrank+1)))
    allocate(dur1(3,nrec_by_proc(myrank+1)), dur2(3, nrec_by_proc(myrank+1)), dur3(3, nrec_by_proc(myrank+1)))

    do itime=1,ntime

      if (myrank == 0) then

        if (mod(itime,100) == 0) write(*,*) 'reading ', itime, ' / ', ntime

        data_rec0 = 0.

        do isim=1,nsim

          read(iux(isim)) data_tmp(:,1)
          read(iuy(isim)) data_tmp(:,2)
          read(iuz(isim)) data_tmp(:,3)

          read(idu1d1(isim)) data_tmp(:,4)
          read(idu1d2(isim)) data_tmp(:,5)
          read(idu1d3(isim)) data_tmp(:,6)
          read(idu2d1(isim)) data_tmp(:,7)
          read(idu2d2(isim)) data_tmp(:,8)
          read(idu2d3(isim)) data_tmp(:,9)
          read(idu3d1(isim)) data_tmp(:,10)
          read(idu3d2(isim)) data_tmp(:,11)
          read(idu3d3(isim)) data_tmp(:,12)

          data_rec0(1,:,:)= data_rec0(1,:,:)+data_tmp(:,:)

        enddo
      endif

      !! scatter the data into the rigth MPI partition
      do icomp =1, 12

        if (myrank == 0) then

          irank = 0
          kk    = 0

          do iface=1,num_boundary_by_proc(irank+1)
            do igll=1,NGLLSQUARE

              kk                    = kk+1
              irec                  = ind_face2rec(iface,igll,irank+1)
              irec_glob(kk,irank+1) = irec

              if ( irec == 0) then
                write(*,*) ' irec ', irec,iface,igll,irank+1
                stop
              endif

              data_rec(icomp,itime,kk)=data_rec0(1,irec,icomp)

            enddo
          enddo

          do irank=1,nbproc-1

            kk = 0

            do iface=1,num_boundary_by_proc(irank+1)
              do igll=1,NGLLSQUARE

                kk                    = kk+1
                irec                  = ind_face2rec(iface,igll,irank+1)
                irec_glob(kk,irank+1) = irec

                if ( irec == 0) then
                  write(*,*) ' irec ', irec,iface,igll,irank+1
                  stop
                endif

                data_tmp_to_send(icomp,kk) = data_rec0(1,irec,icomp)
              enddo
            enddo

            call mpi_send(data_tmp_to_send(icomp,1:kk),nrec_by_proc(irank+1),MPI_REAL,irank,etq,MPI_COMM_WORLD,ierr)

          enddo

        else

          call mpi_recv(data_rec(icomp,itime,:),nrec_by_proc(myrank+1),MPI_REAL,0,etq,MPI_COMM_WORLD,statut,ierr)

        endif

      enddo  !! icomp
    enddo    !! itime

    call mpi_bcast(irec_glob,NGLLSQUARE*MAX_MUN_ABS_BOUNDARY_FACES*nSpecfem_proc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

    if (myrank == 0) then
      do isim=1,nsim

        close(iux(isim))
        close(iuy(isim))
        close(iuz(isim))

        close(idu1d1(isim))
        close(idu1d2(isim))
        close(idu1d3(isim))
        close(idu2d1(isim))
        close(idu2d2(isim))
        close(idu2d3(isim))
        close(idu3d1(isim))
        close(idu3d2(isim))
        close(idu3d3(isim))

      enddo
    endif

    !######################### MPI INTERPOLATION #########################

    nrec_to_store = nrec_by_proc(myrank+1)

    allocate(buffer_to_store_dd(12,nrec_to_store))
    allocate(sinc_tab_displ_and_deriv(ntime))

    irecmin        = 1
    irecmax        = nrec_to_store
    ntime_to_store = (itmax - itmin + 1) / 2

    call create_name_database(prname,myrank,TRACT_PATH)
    write(*,*) TRACT_PATH,prname,myrank

    open(29,file=prname(1:len_trim(prname))//'axisem_displ_for_int_KH',status='unknown', &
            action='write',form='unformatted',iostat=ier)
    if (ier /= 0) write(*,*) 'error opening', prname(1:len_trim(prname))//'axisem_displ_for_int_KH'

    open(30,file=prname(1:len_trim(prname))//'axisem_Du_for_int_KH',status='unknown', &
            action='write',form='unformatted',iostat=ier)
    if (ier /= 0) write(*,*) 'error opening', prname(1:len_trim(prname))//'axisem_Du_for_int_KH'

    if (nrec_to_store >= 100) then
      do i=1,ntime
        write(399,*) (i-1)*dt,data_rec(3,i,100)
      enddo
    endif

    close(399)

    ntnew = ntime_interp

    if (myrank == 0) write(*,*) 'time step ', dtt,ntnew

    tt = -dtt
    do itnew=1,ntnew

      tt = tt + dtt
      current_time_step = tmin + (itnew-1)*dtt !+ 0.5*dtt

      sinc_tab_displ_and_deriv(:) = 0.d0
      buffer_to_store_dd(:,:) = 0.d0

      call compute_sinc(sinc_tab_displ_and_deriv,current_time_step,ntime,dt)

      if (myrank == 0 .and. mod(itnew,100) == 0) write(*,*) myrank, tt,itnew,ntnew

      do irec0=irecmin, irecmax

        do is = 1, ntime

          buffer_to_store_dd(1,irec0)  = buffer_to_store_dd(1,irec0)  + sinc_tab_displ_and_deriv(is)*dble(data_rec(1,is,irec0))
          buffer_to_store_dd(2,irec0)  = buffer_to_store_dd(2,irec0)  + sinc_tab_displ_and_deriv(is)*dble(data_rec(2,is,irec0))
          buffer_to_store_dd(3,irec0)  = buffer_to_store_dd(3,irec0)  + sinc_tab_displ_and_deriv(is)*dble(data_rec(3,is,irec0))

          buffer_to_store_dd(4,irec0)  = buffer_to_store_dd(4,irec0)  + sinc_tab_displ_and_deriv(is)*dble(data_rec(4,is,irec0))
          buffer_to_store_dd(5,irec0)  = buffer_to_store_dd(5,irec0)  + sinc_tab_displ_and_deriv(is)*dble(data_rec(5,is,irec0))
          buffer_to_store_dd(6,irec0)  = buffer_to_store_dd(6,irec0)  + sinc_tab_displ_and_deriv(is)*dble(data_rec(6,is,irec0))
          buffer_to_store_dd(7,irec0)  = buffer_to_store_dd(7,irec0)  + sinc_tab_displ_and_deriv(is)*dble(data_rec(7,is,irec0))
          buffer_to_store_dd(8,irec0)  = buffer_to_store_dd(8,irec0)  + sinc_tab_displ_and_deriv(is)*dble(data_rec(8,is,irec0))
          buffer_to_store_dd(9,irec0)  = buffer_to_store_dd(9,irec0)  + sinc_tab_displ_and_deriv(is)*dble(data_rec(9,is,irec0))
          buffer_to_store_dd(10,irec0) = buffer_to_store_dd(10,irec0) + sinc_tab_displ_and_deriv(is)*dble(data_rec(10,is,irec0))
          buffer_to_store_dd(11,irec0) = buffer_to_store_dd(11,irec0) + sinc_tab_displ_and_deriv(is)*dble(data_rec(11,is,irec0))
          buffer_to_store_dd(12,irec0) = buffer_to_store_dd(12,irec0) + sinc_tab_displ_and_deriv(is)*dble(data_rec(12,is,irec0))

        enddo

        if (irec0 == 100) write(499,*)  current_time_step, buffer_to_store_dd(3,irec0)

        irec =irec_glob(irec0,myrank+1)

        if (irec == 0 ) then
          write(*,*) myrank , irec0, irec
          stop
        endif

        ur(1,irec0) = buffer_to_store_dd(1,irec0)
        ur(2,irec0) = buffer_to_store_dd(2,irec0)
        ur(3,irec0) = buffer_to_store_dd(3,irec0)

        dur1(1,irec0) = buffer_to_store_dd(4,irec0)
        dur1(2,irec0) = buffer_to_store_dd(5,irec0)
        dur1(3,irec0) = buffer_to_store_dd(6,irec0)

        dur2(1,irec0) = buffer_to_store_dd(7,irec0)
        dur2(2,irec0) = buffer_to_store_dd(8,irec0)
        dur2(3,irec0) = buffer_to_store_dd(9,irec0)

        dur3(1,irec0) = buffer_to_store_dd(10,irec0)
        dur3(2,irec0) = buffer_to_store_dd(11,irec0)
        dur3(3,irec0) = buffer_to_store_dd(12,irec0)


      enddo

      write(29) ur

      write(30) dur1, dur2, dur3

    enddo

    close(29)
    close(30)

    close(499)

    write(*,*) 'nbrec ', myrank, irecmax-irecmin+1, ntime_interp
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    write(*,*) myrank,ierr
    call MPI_FINALIZE(ierr)
    write(*,*) myrank,ierr
    stop

  endif !! if (.not. recip_KH_integral)

end program re_format_outputs_files



!!$subroutine interpolate(padi,zero_pad,dat_int,dat,n00,n0,nt,n,dt)
!!$  implicit none
!!$  integer m,nt,n,n0,n00,i
!!$  real dat(n)
!!$  double complex  zero_pad(n0),padi(n00)
!!$  double precision dat_int(nt),value_debb
!!$  real dt
!!$  real, parameter :: pi=3.141592653589793
!!$
!!$  dat_int=0.
!!$  zero_pad=0.
!!$  padi=0.
!!$  value_debb=dat(1)
!!$  padi(1:n)=dat(:)
!!$  !write(399,*) '------------'
!!$  !write(399,*) padi(:)
!!$  call fft(padi,n00,1)
!!$  m=n00 / 2
!!$  zero_pad(1:m)=padi(1:m)
!!$  !write(399,*) zero_pad(1),zero_pad(m/2),zero_pad(m), zero_pad(m+1),zero_pad(n0)
!!$  call fft(zero_pad,n0,-1)
!!$  !write(399,*) zero_pad(1),zero_pad(m/2),zero_pad(m), zero_pad(m+1),zero_pad(n0)
!!$  dat_int(:)=real(zero_pad(1:nt)) /real(m)
!!$
!!$   !do i=1,nt
!!$   !write(499,*) dat_int(i)
!!$   !enddo
!!$   !do i=1,n
!!$   ! write(399,*) dat(i)
!!$   !enddo
!!$  !! debug
!!$  !dat_int(:)=value_debb
!!$  ! divide by m works but I don't know why.
!!$  ! * 2.*pi * (dt)
!!$  !write(*,*) pi, dt , 2.*pi / (dt)
!!$
!!$end subroutine interpolate


!-----------------------------------------

  subroutine fft(dat, nn, isign)
    !use deconvolution, only: CUSTOM_REAL
    implicit none
    integer, parameter :: CUSTOM_REAL=8
    integer :: nn, isign
    real(kind=CUSTOM_REAL) :: dat(2*nn)
    real(kind=CUSTOM_REAL) :: wr, wi, wpr, wpi, wtemp, theta
    real(kind=CUSTOM_REAL) :: tempr, tempi
    integer :: n,i,j,m,mmax,istep

    n = 2*nn
    j = 1
    do i = 1,n,2
       if (j > i) then
          tempr = dat(j)
          tempi = dat(j+1)
          dat(j) = dat(i)
          dat(j+1) = dat(i+1)
          dat(i) = tempr
          dat(i+1) = tempi
       endif
       m = n/2

       do while ((m >= 2) .and. (j > m))
          j = j-m
          m = m/2
       enddo

       j = j+m

    enddo

    mmax = 2

    do while (n > mmax)
       istep = 2*mmax
       theta = 6.28318530717959d0/(isign*mmax)
       wpr = -2.d0*dsin(0.5d0*theta)**2
       wpi = dsin(theta)
       wr = 1.d0
       wi = 0.d0
       do m = 1, mmax, 2
          do i = m, n, istep
             j = i+mmax
             tempr = wr*dat(j) -wi*dat(j+1)
             tempi = wr*dat(j+1)+wi*dat(j)
             dat(j) = dat(i) - tempr
             dat(j+1) = dat(i+1) -tempi
             dat(i) = dat(i) + tempr
             dat(i+1) = dat(i+1) + tempi
          enddo
          wtemp=wr
          wr=wr*wpr-wi*wpi+wr
          wi=wi*wpr+wtemp*wpi+wi
       enddo
       mmax = istep
    enddo
    return
  end subroutine fft



   function next_iunit(i)
     implicit none
     integer i,next_iunit
     i=i+1
     next_iunit=i
   end function next_iunit
!========================================================================================
  subroutine compute_sinc(sinc_tab,tt,n,dt)
     implicit none
     integer i,n
     real(kind=8) sinc_tab(n)
     real tt,dt
     double precision fi,mysinc,t_sinc

     do i=1,n
        t_sinc = (dble(tt) - dble(i-1)*dble(dt))/dble(dt)
        sinc_tab(i) =  mysinc(t_sinc)
     enddo

  end subroutine compute_sinc
!========================================================================================

   subroutine interpol_sinc(f,fi,tt,n,dt,sinc_tab)
     implicit none
     integer i,n
     real f(n),tt,dt
     double precision fi,mysinc,t_sinc,sinc_tab(n)

     fi=0.d0
     !write(399,*) ' ------ '

     do i=1,n
        !t_sinc = (dble(tt) - dble(i-1)*dble(dt))/dble(dt)

        fi = fi + dble(f(i)) * sinc_tab(i) !* mysinc(t_sinc)

     enddo
     !write(399,*) tt,fi,f(1)

   end subroutine interpol_sinc

!!$   !================================================================================
!!$   ! Tukey tapering windows
!!$   !--------------------------------------------------
!!$   ! N     : number of samples
!!$   ! alpha : percentage of signal to taper
!!$   ! tuk   : tapered window
!!$   function tuckeywin(N,alpha) result(tuk)
!!$
!!$    integer(kind=8), intent(in)   :: N
!!$    real(kind=4), intent(in)      :: alpha
!!$
!!$    integer(kind=8) :: i
!!$    real(kind=4), parameter :: pipi=3.141592653589793
!!$    real(kind=4), dimension(N) :: tuk
!!$
!!$    !*** Central part
!!$    tuk(:) = 1.
!!$
!!$    !*** Left part
!!$    do i=0,int(0.5*alpha*(N-1))
!!$       tuk(i+1) = 0.5*(1+cos(pipi*(2.*i/(alpha*(N-1.))-1.)))
!!$    enddo
!!$
!!$    !*** Right part
!!$    do i=int((N-1)*(1-alpha/2.)),N-1
!!$       tuk(i+1) = 0.5*(1+cos(pipi*(2.*i/(alpha*(N-1.))-(2./alpha)+1.)))
!!$    enddo
!!$
!!$  end function tuckeywin
  !--------------------------------------------------------------------------------
  !================================================================================
  ! Sinc functions
  real(kind=8) function mysinc(x)

    real(kind=8) :: x
    real(kind=8), parameter :: pipi=3.141592653589793

    if (abs(x) >= 1d-13) then
       mysinc = sin(pipi*x)/(pipi*x)
    else
       mysinc = 1.d0
    endif

  end function mysinc
  !-----------------------------
!=====================================================================

subroutine create_name_database(prname,iproc,LOCAL_PATH)

! create the name of the database for the mesher and the solver

  implicit none

  integer iproc

! name of the database file
  character(len=500) prname,procname,LOCAL_PATH,clean_LOCAL_PATH

! create the name for the database of the current slide and region
  write(procname,"('/proc',i6.6,'_')") iproc

! suppress white spaces if any
  clean_LOCAL_PATH = adjustl(LOCAL_PATH)

! create full name with path
  prname = clean_LOCAL_PATH(1:len_trim(clean_LOCAL_PATH)) // procname

end subroutine create_name_database

