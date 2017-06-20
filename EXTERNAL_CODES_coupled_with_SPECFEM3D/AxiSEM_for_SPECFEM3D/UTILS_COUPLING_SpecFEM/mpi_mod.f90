module mpi_mod


  use global_parameters
  use rotation_matrix
  !use reading_field
  include 'mpif.h'
  INTEGER myrank,nbrank,ierr

contains

  subroutine init_mpi
    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nbrank,ierr)

  end subroutine init_mpi

  subroutine finalize_mpi
   call MPI_FINALIZE(ierr)
  end subroutine finalize_mpi

  subroutine barrier_mpi()
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  end subroutine barrier_mpi

  subroutine alloc_all_mpi()

   integer :: ier

   if (myrank == 0 ) then
      write(*,*) 'INISDE ALLOC MPI '
      write(*,*) nsim, ntime, nbrec, nbproc, ibeg, iend, nel, rot_mat, trans_rot_mat, rot_mat_mesh, trans_rot_mat_mesh
   endif

   call mpi_bcast(nsim,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
   call mpi_bcast(ntime,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
   call mpi_bcast(nbrec,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
   call mpi_bcast(nbproc,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
   call mpi_bcast(ibeg,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
   call mpi_bcast(iend,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
   call mpi_bcast(nel,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
   call mpi_bcast(rot_mat,9,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
   call mpi_bcast(trans_rot_mat,9,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
   call mpi_bcast(rot_mat_mesh,9,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
   call mpi_bcast(trans_rot_mat_mesh,9,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
   call mpi_bcast(rot_azi_chunk,9,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
   call mpi_bcast(trans_rot_azi_chunk,9,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)


   if (myrank > 0) then
      allocate(reciever_geogr(3,nbrec),reciever_sph(3,nbrec),reciever_cyl(3,nbrec),reciever_interp_value(nbrec))
      allocate(data_rec(nbrec,3),stress_rec(nbrec,6),stress_to_write(nbrec,6),strain_rec(nbrec,6),deriv_rec(nbrec,9))

      allocate(data_tmpKH_rec(nbrec,3,nsim), deriv_tmpKH_rec(nbrec,9,nsim))

      allocate(f1(nbrec),f2(nbrec),phi(nbrec))
      allocate(scoor(ibeg:iend,ibeg:iend,nel),zcoor(ibeg:iend,ibeg:iend,nel))
      allocate(data_read(ibeg:iend,ibeg:iend,nel))
      allocate(xi_rec(nbrec),eta_rec(nbrec))
      allocate(magnitude(nbrec))
      allocate(rec2elm(nbrec))
      allocate(src_type(nsim,2))
      allocate(depth_ele(nel))
   endif
   allocate(stress_reduce(nbrec,6),deriv_reduce(nbrec,9),data_reduce(nbrec,3))


  end subroutine alloc_all_mpi


  subroutine bcast_all_mpi


    ! single
    call mpi_bcast(data_rec,3*nbrec,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(stress_to_write,6*nbrec,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(stress_rec,6*nbrec,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(strain_rec,6*nbrec,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(f1,nbrec,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(f2,nbrec,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(phi,nbrec,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(Mij,6,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(magnitude,nsim,MPI_REAL,0,MPI_COMM_WORLD,ierr)

    call mpi_bcast(lat_mesh,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(lon_mesh,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)

    ! double
    call mpi_bcast(scoor,(iend-ibeg+1)*(iend-ibeg+1)*nel,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(zcoor,(iend-ibeg+1)*(iend-ibeg+1)*nel,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(depth_ele,nel,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(reciever_geogr,3*nbrec,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(reciever_sph,3*nbrec,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(reciever_cyl,3*nbrec,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(xi_rec,nbrec,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(eta_rec,nbrec,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

    ! integer
    call mpi_bcast(rec2elm,nbrec,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(Xk_force,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr) !! CD CD

    ! character
    call mpi_bcast(src_type,10*2*nsim,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)

    ! logical
    call mpi_bcast(recip_KH_integral,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)

  end subroutine bcast_all_mpi


  subroutine distrib_mpi
    integer nb_by_pc,left

    nb_by_pc=(nbrec/nbrank)
    left=nbrec-nb_by_pc*nbrank

    if (myrank < left) then
      irecmin = (myrank*(nb_by_pc+1)) + 1
      irecmax = irecmin + nb_by_pc
    else
      irecmin = myrank*(nb_by_pc) + left + 1
      irecmax = irecmin + nb_by_pc - 1
    endif

  end subroutine distrib_mpi

  subroutine reduce_mpi_veloc()
    !write(*,*) myrank, irecmin,irecmax
    !if (myrank > 12) data_rec(irecmin:irecmax,:)=1.
    !data_rec(irecmin:irecmax,:)=myrank
    call mpi_reduce(data_rec,data_reduce,nbrec*3,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    data_rec(:,:)=data_reduce(:,:)

  end subroutine reduce_mpi_veloc

  subroutine reduce_mpi_deriv()

    call mpi_reduce(deriv_rec,deriv_reduce,nbrec*9,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    deriv_rec(:,:)=deriv_reduce(:,:)

  end subroutine reduce_mpi_deriv

  subroutine reduce_mpi_stress()

    call mpi_reduce(stress_rec,stress_reduce,nbrec*6,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    stress_rec(:,:)=stress_reduce(:,:)

  end subroutine reduce_mpi_stress

end module mpi_mod
