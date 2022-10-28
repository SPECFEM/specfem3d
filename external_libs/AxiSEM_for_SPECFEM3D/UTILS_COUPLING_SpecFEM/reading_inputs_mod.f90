  module reading_inputs

    contains

!!------------------------------------------------------------------------------------------------------------!!

      subroutine read_inputs(isim)

      use global_parameters
      use rotation_matrix

      real(kind=CUSTOM_REAL) srclon,srclat,srccolat
      real(kind=CUSTOM_REAL) meshlon,meshlat,meshcolat
      real(kind=CUSTOM_REAL) x,y,z
      integer i,isim
      integer ispec,igll,jgll,kgll,iface,ilayer,updown

      open(10,file='../expand_2D_3D.par')
      read(10,'(a)') input_point_file
      read(10,'(a)') input_point_file_cart
      read(10,*) nbproc
      read(10,*) lat_src,lon_src
      read(10,*) lat_mesh,lon_mesh,azi_rot !! VM VM add azimuth rotation
      read(10,*) recip_KH_integral
      close(10)

      !! hardcoded the input file name
      input_displ_name(1)='displacement_us' !! CD CD add this
      input_displ_name(2)='displacement_up'
      input_displ_name(3)='displacement_uz'

      input_deriv_name(1)='pgrad_cyl_Dus_1' !! CD CD add this
      input_deriv_name(2)='pgrad_cyl_Dus_2'
      input_deriv_name(3)='pgrad_cyl_Dus_3'
      input_deriv_name(4)='pgrad_cyl_Dup_1'
      input_deriv_name(5)='pgrad_cyl_Dup_2'
      input_deriv_name(6)='pgrad_cyl_Dup_3'
      input_deriv_name(7)='pgrad_cyl_Duz_1'
      input_deriv_name(8)='pgrad_cyl_Duz_2'
      input_deriv_name(9)='pgrad_cyl_Duz_3'

      input_veloc_name(1)='velocityfiel_us'
      input_veloc_name(2)='velocityfiel_up'
      input_veloc_name(3)='velocityfiel_uz'

      input_stress_name(1)='stress_Sg11_sol'
      input_stress_name(2)='stress_Sg22_sol'
      input_stress_name(3)='stress_Sg33_sol'
      input_stress_name(4)='stress_Sg12_sol'
      input_stress_name(5)='stress_Sg13_sol'
      input_stress_name(6)='stress_Sg23_sol'

      !! hardcoded the output file name
      output_displ_name(1)='displ_out_u1' !! CD CD add this
      output_displ_name(2)='displ_out_u2'
      output_displ_name(3)='displ_out_u3'

      output_veloc_name(1)='velocityoutp_u1'
      output_veloc_name(2)='velocityoutp_u2'
      output_veloc_name(3)='velocityoutp_u3'

      output_stress_name(1)='stress_Sg11_out'
      output_stress_name(2)='stress_Sg22_out'
      output_stress_name(3)='stress_Sg33_out'
      output_stress_name(4)='stress_Sg12_out'
      output_stress_name(5)='stress_Sg13_out'
      output_stress_name(6)='stress_Sg23_out'

      output_deriv_name(1)='deriv_out_du1d1' !! CD CD add this
      output_deriv_name(2)='deriv_out_du1d2'
      output_deriv_name(3)='deriv_out_du1d3'
      output_deriv_name(4)='deriv_out_du2d1'
      output_deriv_name(5)='deriv_out_du2d2'
      output_deriv_name(6)='deriv_out_du2d3'
      output_deriv_name(7)='deriv_out_du3d1'
      output_deriv_name(8)='deriv_out_du3d2'
      output_deriv_name(9)='deriv_out_du3d3'


      open(10,file=trim(working_axisem_dir)//trim(simdir(isim))//'/nb_rec_to_read.par')
      read(10,*) ntime
      close(10)

      srclon = lon_src * pi / 180.
      srccolat =  (90. - lat_src) * pi / 180.
      call def_rot_matrix(srccolat,srclon,rot_mat,trans_rot_mat)

      meshlon = lon_mesh * pi / 180.
      meshcolat =  (90. - lat_mesh) * pi / 180.
      call def_rot_matrix_DG(meshcolat,meshlon,rot_mat_mesh,trans_rot_mat_mesh)

      !! VM VM add rot azi
      azi_rot = azi_rot * pi / 180.
      call def_rot_azi_chunk(rot_azi_chunk,trans_rot_azi_chunk,azi_rot)

      !! permutation matrix
      mat(1,1)=0.
      mat(1,2)=1.
      mat(1,3)=0.
      mat(2,1)=1.
      mat(2,2)=0.
      mat(2,3)=0.
      mat(3,1)=0.
      mat(3,2)=0.
      mat(3,3)=1.

      tmat(1,1)=mat(1,1)
      tmat(1,2)=mat(2,1)
      tmat(1,3)=mat(3,1)
      tmat(2,1)=mat(1,2)
      tmat(2,2)=mat(2,2)
      tmat(2,3)=mat(3,2)
      tmat(3,1)=mat(1,3)
      tmat(3,2)=mat(2,3)
      tmat(3,3)=mat(3,3)
      !! permutation matrix

      open(10,file=trim(working_axisem_dir)//trim(simdir(isim))//'../'//trim(input_point_file))
      read(10,*) nbrec
      open(11,file=trim(working_axisem_dir)//trim(simdir(isim))//'../'//trim(input_point_file_cart))
      read(11,*) nbrec

      allocate(up(nbrec))
      allocate(reciever_geogr(3,nbrec),reciever_sph(3,nbrec),reciever_cyl(3,nbrec),reciever_interp_value(nbrec))
      allocate(data_rec(nbrec,3),stress_rec(nbrec,6),stress_to_write(nbrec,6),strain_rec(nbrec,6),deriv_rec(nbrec,9))
      allocate(data_tmpKH_rec(nbrec,3,nsim), deriv_tmpKH_rec(nbrec,9,nsim))

      stress_rec=0.
      data_rec=0.
      allocate(f1(nbrec),f2(nbrec),phi(nbrec))

      do i=1,nbrec !! radius, latitude, longitude
        read(10,*) reciever_geogr(1,i),reciever_geogr(2,i),reciever_geogr(3,i)
        read(11,*) x,y,z,ispec,igll,jgll,kgll,iface,ilayer,updown
        if (kgll == 5) then
          up(i) = .false.
        else
          up(i) = .true.
        endif
        reciever_sph(1,i) = reciever_geogr(1,i)*1000.
        reciever_sph(2,i) = (90. - reciever_geogr(2,i)) * pi / 180.  ! Read latitudes instead of colatitudes
        reciever_sph(3,i) = reciever_geogr(3,i)  * pi / 180.

        call rotate_box(reciever_sph(1,i), reciever_sph(2,i), reciever_sph(3,i),trans_rot_mat)

        reciever_cyl(1,i) = reciever_sph(1,i)  * sin( reciever_sph(2,i))
        reciever_cyl(2,i) = reciever_sph(3,i)
        reciever_cyl(3,i) = reciever_sph(1,i)  * cos( reciever_sph(2,i))

        phi(i) = reciever_cyl(2,i)
      enddo

      close(10)
      close(11)

      end subroutine read_inputs

!
!--------------------!
!

      subroutine read_info_simu(nsim_to_send)

      use global_parameters, only: Mij,nsim,simdir,src_type,ntime,working_axisem_dir, &
                                   dtt,magnitude,recip_KH_integral,Xk_force

      character(len=256)  :: keyword, keyvalue, line
      character(len=100),allocatable, dimension(:) :: bkgrndmodel, rot_rec
      character(len=7),allocatable, dimension(:)   :: stf_type

      integer iostat,ioerr,isim,nsim_to_send
      integer, allocatable, dimension(:) :: ishift_deltat, ishift_seisdt, ishift_straindt
      integer, allocatable, dimension(:) :: ibeg, iend
      integer, allocatable, dimension(:) :: nt, nt_strain, nrec_tmp, nt_seis_tmp
      real, allocatable, dimension(:)    :: dt_tmp, period_tmp, dt_seis_tmp
      real, allocatable, dimension(:)    :: dt_strain, dt_snap, nt_snap, &
                                            srccolat_tmp, srclon_tmp, src_depth_tmp
      real, allocatable, dimension(:)    :: shift_fact_tmp
      !real, allocatable, dimension(:)     :: magnitude

      integer :: i_param_post = 500
      integer    i, nsimtmp
      logical    use_netcdf
      real       tshift
      character(len=12)  :: src_file_type
      character(len=100) :: junk

      !! CAUTION : this parameter is HARDCODED, if you change (add line for example) expand_2D_3D.par = =>
      !!           you have to CHANGE ALSO this parameter
      integer, parameter :: linestoskip_expandpar = 5

      working_axisem_dir='./'
      if (.not. recip_KH_integral) &
          open(unit=i_param_post, file='inparam_basic', status='old', action='read', iostat=ioerr)

      open(20,file='../expand_2D_3D.par')
      do i = 1, linestoskip_expandpar
        read(20,*)
      enddo
      read(20,*) recip_KH_integral
      read(20,*) nsimtmp
      read(20,*) Xk_force
      close(20)

      if (.not. recip_KH_integral) then

        do
          read(i_param_post, fmt='(a256)', iostat=ioerr) line
          if (ioerr < 0) exit
          if (len(trim(line)) < 1 .or. line(1:1) == '#') cycle

          read(line,*) keyword, keyvalue

          select case(trim(keyword))
          case('SIMULATION_TYPE')

            if (keyvalue == 'single') then
              nsim = 1
              if (nsimtmp /= nsim) stop 'NSIM INCORRECTLY CONFIGURED in expand_2D_3D.par'
              allocate(simdir(nsim))
              simdir(1) = "./"

            else if (keyvalue == 'moment') then
              nsim = 4
              if (nsimtmp /= nsim) stop 'NSIM INCORRECTLY CONFIGURED in expand_2D_3D.par'
              allocate(simdir(nsim))
              simdir(1) = "MZZ/"
              simdir(2) = "MXX_P_MYY/"
              simdir(3) = "MXZ_MYZ/"
              simdir(4) = "MXY_MXX_M_MYY/"

              !! read CMTSOLUTION
              write(*,*)'  reading CMTSOLUTION file....'
              open(unit=20000,file='CMTSOLUTION',POSITION='REWIND',status='old')
              read(20000,*) junk
              read(20000,*) junk
              read(20000,*) junk
              read(20000,*) junk
              read(20000,*) junk
              read(20000,*) junk
              read(20000,*) junk
              read(20000,*) junk, Mij(1) !Mrr
              read(20000,*) junk, Mij(2) !Mtt
              read(20000,*) junk, Mij(3) !Mpp
              read(20000,*) junk, Mij(4) !Mrt
              read(20000,*) junk, Mij(5) !Mrp
              read(20000,*) junk, Mij(6) !Mtp
              close(20000)

              Mij = Mij / 1.E7 ! CMTSOLUTION given in dyn-cm

            else if (keyvalue == 'force') then

              nsim = 2
              if (nsimtmp /= nsim) stop 'NSIM INCORRECTLY CONFIGURED in expand_2D_3D.par'
              allocate(simdir(nsim))
              simdir(1) = "PZ/" !! Vertical force
              simdir(2) = "PX/" !! Horizontal force


            endif
          end select
        enddo

      else

        nsim = 3
        if (nsimtmp /= nsim) stop 'NSIM INCORRECTLY CONFIGURED in expand_2D_3D.par'
        allocate(simdir(nsim))
        simdir(1) = "ThetaF/"
        simdir(2) = "PhiF/"
        simdir(3) = "VertF/"

      endif

      nsim_to_send=nsim
      !open(10,file=trim(working_axisem_dir)//trim(simdir(1))//'nb_rec_to_read.par')
      !read(10,*) ntime
      !close(10)

      close(i_param_post)

      tshift = 0.

      allocate(bkgrndmodel(nsim), stf_type(nsim))
      allocate(src_type(nsim,2))
      allocate(dt_tmp(nsim), period_tmp(nsim), magnitude(nsim), dt_seis_tmp(nsim))
      allocate(dt_strain(nsim), dt_snap(nsim))
      allocate(rot_rec(nsim))
      allocate(nt(nsim), nrec_tmp(nsim), nt_seis_tmp(nsim), nt_strain(nsim), nt_snap(nsim))
      allocate(ibeg(nsim), iend(nsim), srccolat_tmp(nsim), srclon_tmp(nsim), src_depth_tmp(nsim))
      allocate(shift_fact_tmp(nsim))
      allocate(ishift_deltat(nsim), ishift_seisdt(nsim), ishift_straindt(nsim))

      do isim = 1,nsim

        open(unit=99,file=trim(simdir(isim))//'/simulation.info')
        read(99,*) bkgrndmodel(isim)
        read(99,*) dt_tmp(isim)
        read(99,*) nt(isim)
        read(99,*) src_type(isim,1)
        read(99,*) src_type(isim,2)
        read(99,*) stf_type(isim)
        read(99,*) src_file_type
        read(99,*) period_tmp(isim)
        read(99,*) src_depth_tmp(isim)
        read(99,*) srccolat_tmp(isim)
        read(99,*) srclon_tmp(isim)
        read(99,*) magnitude(isim)
        read(99,*) nrec_tmp(isim)
        read(99,*) nt_seis_tmp(isim)
        read(99,*) dt_seis_tmp(isim)
        read(99,*) nt_strain(isim)
        read(99,*) dt_strain(isim)
        read(99,*) nt_snap(isim)
        read(99,*) dt_snap(isim)
        read(99,*) rot_rec(isim)
        read(99,*) ibeg(isim)
        read(99,*) iend(isim)
        read(99,*) shift_fact_tmp(isim)
        read(99,*) ishift_deltat(isim)
        read(99,*) ishift_seisdt(isim)
        read(99,*) ishift_straindt(isim)
        read(99,*)
        read(99,*)
        read(99,*) use_netcdf
        close(99)

      enddo

      dtt= dt_strain(1)
      open(666,file='info_for_specefm.txt')
      write(666,*) dt_strain(1)
      close(666)

      end subroutine read_info_simu

!!------------------------------------------------------------------------------------------------------------!!

  end module reading_inputs
