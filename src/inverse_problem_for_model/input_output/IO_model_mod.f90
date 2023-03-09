!=====================================================================
!
!                          S p e c f e m 3 D
!                          -----------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                              CNRS, France
!                       and Princeton University, USA
!                 (there are currently many more authors!)
!                           (c) October 2017
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
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

!
!
! import or export external model from specfem mesh or FD grid
!
!

module IO_model

  !! IMPORT VARIABLES FROM SPECFEM -------------------------------------------------------------------------------------------------
  use shared_parameters, only: NUMBER_OF_SIMULTANEOUS_RUNS,  ANISOTROPY, ACOUSTIC_SIMULATION, ELASTIC_SIMULATION, LOCAL_PATH

  use specfem_par, only: CUSTOM_REAL,  NGLLX, NGLLY, NGLLZ, NSPEC_AB,  myrank, mygroup, &
                         rhostore, mustore, kappastore, FOUR_THIRDS, PI

  use specfem_par_elastic, only: rho_vp, rho_vs,  NSPEC_ANISO, &
                                  c11store,c12store,c13store,c14store,c15store,c16store, &
                                  c22store,c23store,c24store,c25store,c26store,c33store, &
                                  c34store,c35store,c36store,c44store,c45store,c46store, &
                                  c55store,c56store,c66store

  !---------------------------------------------------------------------------------------------------------------------------------
  use inverse_problem_par
  use mesh_tools

  implicit none

  !---- locals --
  integer, private :: ier

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!--------------------------------------------------------------------------------------------------------------------
!> read input fd model
!--------------------------------------------------------------------------------------------------------------------
  subroutine ReadInputFDmodel(inversion_param)
    type(inver),                                                intent(inout) :: inversion_param
    character(len=MAX_STRING_LEN)                                             :: fd_grid_model
    if (ANISOTROPY) then

       if (DEBUG_MODE) then
          write(IIDD, *)
          write(IIDD, *) 'READING INPUT FD MODEL FROM DISK FOR C_IJKL'
          write(IIDD, *)
          write(IIDD, *) 'Family parameter:',  trim(inversion_param%parameter_family_name)
          write(IIDD, *)
       endif

       !! HRDCODED NAME FOR NOW
       fd_grid_model='fd_grid.txt'
       call import_FD_model_ANISO(fd_grid_model, inversion_param)

    else

       if (DEBUG_MODE) then
          write(IIDD, *)
          write(IIDD, *) ' READING INPUT FD MODEL FROM DISK FOR RHO VP VS '
          write(IIDD, *)
          write(IIDD, *) ' Family Parameter : ',  trim(inversion_param%parameter_family_name)
          write(IIDD, *)
       endif

       if (ELASTIC_SIMULATION) then
          !! use a simple nearest neighbor interpolation
          !! HARDCODED NAME FOR NOW
          fd_grid_model='fd_grid.txt'
          call import_FD_model_Elastic_ISO(fd_grid_model)
       endif

       if (ACOUSTIC_SIMULATION) then
          !! use a simple nearest neighbor interpolation
          !! HARDCODED NAME FOR NOW
          fd_grid_model='fd_grid.txt'
          call import_FD_model_ACOUSTIC(fd_grid_model)
       endif

    endif

    ! create mass matrix (duplicate subroutine because of conflicts in
    ! module from generate_databases and specfem
    ! eg : rmass is defined in module create_regions_mesh_ext_par
    !      and is defined in module specfem_par_elastic
    ! to avoid that kind of conflicts we need to define
    ! each array once in a module and all components
    ! of package must use the same module
    call create_mass_matrices_Stacey_duplication_routine()

  end subroutine ReadInputFDmodel

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!--------------------------------------------------------------------------------------------------------------------
! read sem model already decomposed in NROC files
!--------------------------------------------------------------------------------------------------------------------
  subroutine ReadInputSEMmodel(inversion_param)

    type(inver),                                                intent(in)    :: inversion_param
    ! local
    character(len=256)                                                        :: path_file,name_file
    real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable                   :: wks_model_rho, wks_model_vp, wks_model_vs
    real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable                   :: wks_model_cij

    !! rwriting SEM model stored in mesh and in (rho,vp,vs) for isotropic case
    !! and in cijkl for anisotropic case.

    if (ANISOTROPY) then

       if (DEBUG_MODE) then
          write(IIDD, *)
          write(IIDD, *) ' READING INPUT MODEL FROM DISK IN : model_cij_input '
          write(IIDD, *)
          write(IIDD, *) ' Family Parameter : ',  trim(inversion_param%parameter_family_name)
          write(IIDD, *)
       endif

       allocate(wks_model_cij(NGLLX,NGLLY,NGLLZ,NSPEC_ANISO), stat=ier)
       if (ier /= 0) call exit_MPI_without_rank('error allocating array 298')
       if (ier /= 0) call exit_MPI(myrank,"error allocation wks_model_cij in ReadInputSEMmodel subroutine, IO_model_mod")

       path_file = trim(LOCAL_PATH) // '/proc'

       if (mygroup <= 0) then !! only the first group read model and need to bcast at all other
          write(name_file,'(i6.6,a20)') myrank, '_model_cij_input.bin'
          open(888,file=trim(path_file),form='unformatted')
          read(888)  c11store
          read(888)  c12store
          read(888)  c13store
          read(888)  c14store
          read(888)  c15store
          read(888)  c16store
          read(888)  c22store
          read(888)  c23store
          read(888)  c24store
          read(888)  c25store
          read(888)  c26store
          read(888)  c33store
          read(888)  c34store
          read(888)  c35store
          read(888)  c36store
          read(888)  c44store
          read(888)  c45store
          read(888)  c46store
          read(888)  c55store
          read(888)  c56store
          read(888)  c66store
          !! add density
          read(888)  rhostore
          !! for stacey (those arrays are always defined in specfem)
          read(888)  rho_vp
          read(888)  rho_vs
          close(888)
       endif

       !! bcast to others groups
       if (NUMBER_OF_SIMULTANEOUS_RUNS > 1) then


          call bcast_all_cr_for_database(c11store(1,1,1,1), NGLLX*NGLLY*NGLLZ*NSPEC_ANISO)
          call bcast_all_cr_for_database(c12store(1,1,1,1), NGLLX*NGLLY*NGLLZ*NSPEC_ANISO)
          call bcast_all_cr_for_database(c13store(1,1,1,1), NGLLX*NGLLY*NGLLZ*NSPEC_ANISO)
          call bcast_all_cr_for_database(c14store(1,1,1,1), NGLLX*NGLLY*NGLLZ*NSPEC_ANISO)
          call bcast_all_cr_for_database(c15store(1,1,1,1), NGLLX*NGLLY*NGLLZ*NSPEC_ANISO)
          call bcast_all_cr_for_database(c16store(1,1,1,1), NGLLX*NGLLY*NGLLZ*NSPEC_ANISO)
          call bcast_all_cr_for_database(c22store(1,1,1,1), NGLLX*NGLLY*NGLLZ*NSPEC_ANISO)
          call bcast_all_cr_for_database(c23store(1,1,1,1), NGLLX*NGLLY*NGLLZ*NSPEC_ANISO)
          call bcast_all_cr_for_database(c24store(1,1,1,1), NGLLX*NGLLY*NGLLZ*NSPEC_ANISO)
          call bcast_all_cr_for_database(c25store(1,1,1,1), NGLLX*NGLLY*NGLLZ*NSPEC_ANISO)
          call bcast_all_cr_for_database(c26store(1,1,1,1), NGLLX*NGLLY*NGLLZ*NSPEC_ANISO)
          call bcast_all_cr_for_database(c33store(1,1,1,1), NGLLX*NGLLY*NGLLZ*NSPEC_ANISO)
          call bcast_all_cr_for_database(c34store(1,1,1,1), NGLLX*NGLLY*NGLLZ*NSPEC_ANISO)
          call bcast_all_cr_for_database(c35store(1,1,1,1), NGLLX*NGLLY*NGLLZ*NSPEC_ANISO)
          call bcast_all_cr_for_database(c36store(1,1,1,1), NGLLX*NGLLY*NGLLZ*NSPEC_ANISO)
          call bcast_all_cr_for_database(c44store(1,1,1,1), NGLLX*NGLLY*NGLLZ*NSPEC_ANISO)
          call bcast_all_cr_for_database(c45store(1,1,1,1), NGLLX*NGLLY*NGLLZ*NSPEC_ANISO)
          call bcast_all_cr_for_database(c46store(1,1,1,1), NGLLX*NGLLY*NGLLZ*NSPEC_ANISO)
          call bcast_all_cr_for_database(c55store(1,1,1,1), NGLLX*NGLLY*NGLLZ*NSPEC_ANISO)
          call bcast_all_cr_for_database(c56store(1,1,1,1), NGLLX*NGLLY*NGLLZ*NSPEC_ANISO)
          call bcast_all_cr_for_database(c66store(1,1,1,1), NGLLX*NGLLY*NGLLZ*NSPEC_ANISO)
          call bcast_all_cr_for_database(rhostore(1,1,1,1), NGLLX*NGLLY*NGLLZ*NSPEC_AB)
          call bcast_all_cr_for_database(rho_vp(1,1,1,1), NGLLX*NGLLY*NGLLZ*NSPEC_AB)
          call bcast_all_cr_for_database(rho_vs(1,1,1,1), NGLLX*NGLLY*NGLLZ*NSPEC_AB)

       endif

       deallocate(wks_model_cij)

   else

      if (DEBUG_MODE) then
         write(IIDD, *)
         write(IIDD, *) ' READING INPUT MODEL FROM DISK IN : model_vp_input, model_vs_input, model_rh_input'
         write(IIDD, *)
         write(IIDD, *) ' Family Parameter : ',  trim(inversion_param%parameter_family_name)
         write(IIDD, *)
      endif

       allocate(wks_model_rho(NGLLX,NGLLY,NGLLZ,NSPEC_AB), stat=ier)
       if (ier /= 0) call exit_MPI_without_rank('error allocating array 299')
       if (ier /= 0) call exit_MPI(myrank,"error allocation wks_model_rho in ReadInputSEMmodel subroutine, IO_model_mod")

       allocate(wks_model_vp(NGLLX,NGLLY,NGLLZ,NSPEC_AB), stat=ier)
       if (ier /= 0) call exit_MPI_without_rank('error allocating array 300')
       if (ier /= 0) call exit_MPI(myrank,"error allocation wks_model_vp in ReadInputSEMmodel subroutine, IO_model_mod")

       allocate(wks_model_vs(NGLLX,NGLLY,NGLLZ,NSPEC_AB), stat=ier)
       if (ier /= 0) call exit_MPI_without_rank('error allocating array 301')
       if (ier /= 0) call exit_MPI(myrank,"error allocation wks_model_vs in ReadInputSEMmodel subroutine, IO_model_mod")

       if (mygroup <= 0) then !! only the first group read model and need to bcast at all other
          !! read input model
          path_file = trim(LOCAL_PATH) // '/proc'
          write(name_file,'(i6.6,a19)') myrank, '_model_vp_input.bin'
          path_file=(trim(path_file))//trim(name_file)
          open(888,file=trim(path_file),form='unformatted')
          read(888) wks_model_vp
          close(888)

          path_file = trim(LOCAL_PATH) // '/proc'
          write(name_file,'(i6.6,a19)') myrank, '_model_vs_input.bin'
          path_file=(trim(path_file))//trim(name_file)
          open(888,file=trim(path_file),form='unformatted')
          read(888) wks_model_vs
          close(888)

          path_file = trim(LOCAL_PATH) // '/proc'
          write(name_file,'(i6.6,a19)') myrank, '_model_rh_input.bin'
          path_file=(trim(path_file))//trim(name_file)
          open(888,file=trim(path_file),form='unformatted')
          read(888) wks_model_rho
          close(888)

          !! put input model in specfem database
          ! new model
          if (ELASTIC_SIMULATION) then
             rho_vs(:,:,:,:)  = wks_model_rho(:,:,:,:) *  wks_model_vs(:,:,:,:)
             rho_vp(:,:,:,:)  = wks_model_rho(:,:,:,:) *  wks_model_vp(:,:,:,:)
             mustore(:,:,:,:) = wks_model_rho(:,:,:,:) * (wks_model_vs(:,:,:,:)**2)
          endif
          kappastore(:,:,:,:) = wks_model_rho(:,:,:,:) * ( (wks_model_vp(:,:,:,:)**2) -&
               FOUR_THIRDS*(wks_model_vs(:,:,:,:)**2))
          rhostore(:,:,:,:) = wks_model_rho(:,:,:,:)

       endif

       !! bcast to others groups
       if (NUMBER_OF_SIMULTANEOUS_RUNS > 1) then

          call bcast_all_cr_for_database(kappastore(1,1,1,1), NGLLX*NGLLY*NGLLZ*NSPEC_AB)
          call bcast_all_cr_for_database(rhostore(1,1,1,1), NGLLX*NGLLY*NGLLZ*NSPEC_AB)

          if (ELASTIC_SIMULATION) then
             call bcast_all_cr_for_database(mustore(1,1,1,1), NGLLX*NGLLY*NGLLZ*NSPEC_AB)
             call bcast_all_cr_for_database(rho_vs(1,1,1,1), NGLLX*NGLLY*NGLLZ*NSPEC_AB)
             call bcast_all_cr_for_database(rho_vp(1,1,1,1), NGLLX*NGLLY*NGLLZ*NSPEC_AB)
          endif

       endif

       deallocate(wks_model_rho, wks_model_vp, wks_model_vs)

    endif

    !! define the mass matrix that corresponds on model just read
    call create_mass_matrices_Stacey_duplication_routine()

  end subroutine ReadInputSEMmodel

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!--------------------------------------------------------------------------------------------------------------------
! read sem model already decomposed in NROC files
!--------------------------------------------------------------------------------------------------------------------
  subroutine ReadInputSEMpriormodel(inversion_param)

    type(inver),                                                intent(inout) :: inversion_param
    ! local
    character(len=256)                                                        :: path_file,name_file
    real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable                   :: wks_model
    !real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable                   :: wks_model_cij

    !! rwriting SEM model stored in mesh and in (rho,vp,vs) for isotropic case
    !! and in cijkl for anisotropic case.

    if (ANISOTROPY) then

       write(*,*) " Aniso Not implemented yet for prior model "
       stop

    else

       allocate(wks_model(NGLLX,NGLLY,NGLLZ,NSPEC_AB), stat=ier)
       if (ier /= 0) call exit_MPI_without_rank('error allocating array 302')
       if (ier /= 0) call exit_MPI(myrank,"error allocation wks_model in ReadInputSEMpriormodel subroutine, IO_model_mod")
       wks_model(:,:,:,:) = 0._CUSTOM_REAL

       allocate(inversion_param%prior_model(NGLLX,NGLLY,NGLLZ,NSPEC_AB,inversion_param%NinvPar), stat=ier)
       if (ier /= 0) call exit_MPI_without_rank('error allocating array 303')
       if (ier /= 0) call exit_MPI(myrank,"error allocation wks_model_vp in ReadInputSEMmodel subroutine, IO_model_mod")
       inversion_param%prior_model(:,:,:,:,:) = 0._CUSTOM_REAL

       if (mygroup <= 0) then !! only the first group read model and need to bcast at all other
          !! read input model
          path_file = trim(LOCAL_PATH) // '/proc'
          write(name_file,'(i6.6,a19)') myrank, '_model_vp_prior.bin'
          path_file=(trim(path_file))//trim(name_file)
          open(888,file=trim(path_file),form='unformatted')
          read(888) wks_model
          close(888)
          inversion_param%prior_model(:,:,:,:,1) = wks_model(:,:,:,:)

          path_file = trim(LOCAL_PATH) // '/proc'
          write(name_file,'(i6.6,a19)') myrank, '_model_vs_prior.bin'
          path_file=(trim(path_file))//trim(name_file)
          open(888,file=trim(path_file),form='unformatted')
          read(888) wks_model
          close(888)
          inversion_param%prior_model(:,:,:,:,2) = wks_model(:,:,:,:)

          path_file = trim(LOCAL_PATH) // '/proc'
          write(name_file,'(i6.6,a19)') myrank, '_model_rh_prior.bin'
          path_file=(trim(path_file))//trim(name_file)
          open(888,file=trim(path_file),form='unformatted')
          read(888) wks_model
          close(888)
          inversion_param%prior_model(:,:,:,:,3) = wks_model(:,:,:,:)

       endif

       !! bcast to others groups
       if (NUMBER_OF_SIMULTANEOUS_RUNS > 1) then
          call bcast_all_cr_for_database(inversion_param%prior_model(1,1,1,1,1), &
                                         NGLLX*NGLLY*NGLLZ*NSPEC_AB*inversion_param%NinvPar)
       endif

       deallocate(wks_model)

    endif

  end subroutine ReadInputSEMpriormodel

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!--------------------------------------------------------------------------------------------------------------------
! write sem model decomposed in NROC files
!--------------------------------------------------------------------------------------------------------------------
  subroutine WriteOutputSEMmodel(inversion_param)

    type(inver),                                                intent(in)    :: inversion_param
    ! local
    character(len=256)                                                        :: path_file,name_file
    real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable                   :: wks_model_rho, wks_model_vp, wks_model_vs
    real(kind=CUSTOM_REAL)                                                    :: vp_min,vp_max,vs_min,vs_max,rho_min,rho_max
    real(kind=CUSTOM_REAL)                                                    :: vp_min_glob,vp_max_glob,vs_min_glob,vs_max_glob, &
                                                                                 rho_min_glob,rho_max_glob
    integer                                                                   :: ispec, ifrq

    ! current frequency iteration
    ifrq = iter_frq

    !! reading SEM model stored in mesh and in (rho,vp,vs) for isotropic case
    !! and in cijkl for anisotropic case.

    if (ANISOTROPY) then
       ! user output
       if (myrank == 0) then
          write(INVERSE_LOG_FILE,*)
          write(INVERSE_LOG_FILE,*) '******************************************************'
          write(INVERSE_LOG_FILE,*) '*** WRITING OUPUT MODEL ON DISK : model_cij_output ***'
          write(INVERSE_LOG_FILE,*) '******************************************************'
          write(INVERSE_LOG_FILE,*)
          write(INVERSE_LOG_FILE,*) ' Family Parameter : ',  trim(inversion_param%parameter_family_name)
          write(INVERSE_LOG_FILE,*)
          if (inversion_param%only_forward) then
            write(INVERSE_LOG_FILE,*) ' frequency band   : ',ifrq
          else
            write(INVERSE_LOG_FILE,*) ' frequency band   : ',ifrq,' out of ',inversion_param%Nifrq
          endif
          write(INVERSE_LOG_FILE,*)
          call flush_iunit(INVERSE_LOG_FILE)
       endif

       ! stores model
       if (trim(inversion_param%parameter_family_name) == "VTI") then

          call write_vti_sem_model(ifrq)

       else

          if (mygroup <= 0) then !! only the first group write the model and need to bcast at all other
             path_file = trim(LOCAL_PATH) // '/proc'
             write(name_file,'(i6.6,a5,i6.6,a21)') myrank,'_freq', ifrq , '_model_cij_output.bin'
             path_file=(trim(path_file))//trim(name_file)
             open(888,file=trim(path_file),form='unformatted')
             !! 21 coefs elastic tensor
             write(888)  c11store
             write(888)  c12store
             write(888)  c13store
             write(888)  c14store
             write(888)  c15store
             write(888)  c16store
             write(888)  c22store
             write(888)  c23store
             write(888)  c24store
             write(888)  c25store
             write(888)  c26store
             write(888)  c33store
             write(888)  c34store
             write(888)  c35store
             write(888)  c36store
             write(888)  c44store
             write(888)  c45store
             write(888)  c46store
             write(888)  c55store
             write(888)  c56store
             write(888)  c66store
             !! add density
             write(888)  rhostore
             !! for stacey (those arrays are always defined in specfem)
             write(888)  rho_vp
             write(888)  rho_vs
             close(888)

          endif
       endif

    else
       ! isotropic models
       allocate(wks_model_rho(NGLLX,NGLLY,NGLLZ,NSPEC_AB), stat=ier)
       if (ier /= 0) call exit_MPI_without_rank('error allocating array 304')
       if (ier /= 0) call exit_MPI(myrank,"error allocation wks_model_rho in ReadInputSEMmodel subroutine, IO_model_mod")

       allocate(wks_model_vp(NGLLX,NGLLY,NGLLZ,NSPEC_AB), stat=ier)
       if (ier /= 0) call exit_MPI_without_rank('error allocating array 305')
       if (ier /= 0) call exit_MPI(myrank,"error allocation wks_model_vp in ReadInputSEMmodel subroutine, IO_model_mod")

       allocate(wks_model_vs(NGLLX,NGLLY,NGLLZ,NSPEC_AB), stat=ier)
       if (ier /= 0) call exit_MPI_without_rank('error allocating array 306')
       if (ier /= 0) call exit_MPI(myrank,"error allocation wks_model_vs in ReadInputSEMmodel subroutine, IO_model_mod")

       wks_model_rho(:,:,:,:) = 0._CUSTOM_REAL
       wks_model_vp(:,:,:,:) = 0._CUSTOM_REAL
       wks_model_vs(:,:,:,:) = 0._CUSTOM_REAL

       ! get model from specfem database
       if (ELASTIC_SIMULATION) then
          do ispec=1, NSPEC_AB  !! update just the elastic elements
             if (ispec_is_elastic(ispec)) then
                wks_model_rho(:,:,:,ispec) =  rho_vs(:,:,:,ispec) * rho_vs(:,:,:,ispec) / mustore(:,:,:,ispec)
                wks_model_vp(:,:,:,ispec) = (kappastore(:,:,:,ispec) + (4./3.) * mustore(:,:,:,ispec) ) / rho_vp(:,:,:,ispec)
                wks_model_vs(:,:,:,ispec) = mustore(:,:,:,ispec) /  rho_vs(:,:,:,ispec)
             endif
          enddo
       endif

       if (ACOUSTIC_SIMULATION) then
          do ispec=1, NSPEC_AB  !! update just the acoustic elements
             if (ispec_is_acoustic(ispec)) then
                wks_model_rho(:,:,:,ispec) = rhostore(:,:,:,ispec)
                wks_model_vp(:,:,:,ispec) = sqrt(kappastore(:,:,:,ispec)/rhostore(:,:,:,ispec))
                !wks_model_vs(:,:,:,ispec) = 0._CUSTOM_REAL  !! put vs=0 only in acoustic elements
             endif
          enddo
       endif

       ! model stats
       rho_min = minval(wks_model_rho(:,:,:,:))
       rho_max = maxval(wks_model_rho(:,:,:,:))
       call min_all_cr(rho_min,rho_min_glob)
       call max_all_cr(rho_max,rho_max_glob)

       vp_min = minval(wks_model_vp(:,:,:,:))
       vp_max = maxval(wks_model_vp(:,:,:,:))
       call min_all_cr(vp_min,vp_min_glob)
       call max_all_cr(vp_max,vp_max_glob)

       if (ELASTIC_SIMULATION) then
         vs_min = minval(wks_model_vs(:,:,:,:))
         vs_max = maxval(wks_model_vs(:,:,:,:))
         call min_all_cr(vs_min,vs_min_glob)
         call max_all_cr(vs_max,vs_max_glob)
       endif

       ! user output
       if (myrank == 0) then
          write(INVERSE_LOG_FILE,*)
          write(INVERSE_LOG_FILE,*) &
               '*****************************************************************************************'
          write(INVERSE_LOG_FILE,*) &
               '*** WRITING OUTPUT MODEL ON DISK : model_vp_output, model_vs_output, model_rho_output ***'
          write(INVERSE_LOG_FILE,*) &
               '*****************************************************************************************'
          write(INVERSE_LOG_FILE,*)
          write(INVERSE_LOG_FILE,*) ' Family Parameter : ',  trim(inversion_param%parameter_family_name)
          write(INVERSE_LOG_FILE,*)
          if (inversion_param%only_forward) then
            write(INVERSE_LOG_FILE,*) ' frequency band   : ',ifrq
          else
            write(INVERSE_LOG_FILE,*) ' frequency band   : ',ifrq,' out of ',inversion_param%Nifrq
          endif
          write(INVERSE_LOG_FILE,*)
          if (ACOUSTIC_SIMULATION .and. ELASTIC_SIMULATION) then
            write(INVERSE_LOG_FILE,*) ' Model: acoustic/elastic'
          else if (ACOUSTIC_SIMULATION) then
            write(INVERSE_LOG_FILE,*) ' Model: acoustic'
          else if (ELASTIC_SIMULATION) then
            write(INVERSE_LOG_FILE,*) ' Model: elastic'
          endif
          ! stats
          write(INVERSE_LOG_FILE,*) '        P velocity min,max = ',vp_min_glob,vp_max_glob
          if (ELASTIC_SIMULATION) then
            write(INVERSE_LOG_FILE,*) '        S velocity min,max = ',vs_min_glob,vs_max_glob
          endif
          write(INVERSE_LOG_FILE,*) '        density    min,max = ',rho_min_glob,rho_max_glob
          write(INVERSE_LOG_FILE,*)
          call flush_iunit(INVERSE_LOG_FILE)
       endif

       if (mygroup <= 0) then !! only the first group write the model and need to bcast at all other
          !! read input model
          path_file = trim(LOCAL_PATH) // '/proc'
          write(name_file,'(i6.6,a5,i6.6,a20)') myrank,'_freq', ifrq , '_model_vp_output.bin'
          path_file=(trim(path_file))//trim(name_file)
          open(888,file=trim(path_file),form='unformatted')
          write(888) wks_model_vp
          close(888)

          path_file = trim(LOCAL_PATH) // '/proc'
          write(name_file,'(i6.6,a5,i6.6,a20)') myrank,'_freq', ifrq , '_model_vs_output.bin'
          path_file=(trim(path_file))//trim(name_file)
          open(888,file=trim(path_file),form='unformatted')
          write(888) wks_model_vs
          close(888)

          path_file = trim(LOCAL_PATH) // '/proc'
          write(name_file,'(i6.6,a5,i6.6,a21)') myrank,'_freq', ifrq , '_model_rho_output.bin'
          path_file=(trim(path_file))//trim(name_file)
          open(888,file=trim(path_file),form='unformatted')
          write(888) wks_model_rho
          close(888)

       endif

       deallocate(wks_model_rho, wks_model_vp, wks_model_vs)

    endif

  end subroutine WriteOutputSEMmodel
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!--------------------------------------------------------------------------------------------------------------------
! write fields on disk
!--------------------------------------------------------------------------------------------------------------------
  subroutine  DumpArraysMeshSpecfem(model, gradient, descent, inversion_param)

    type(inver),                                                intent(in)    :: inversion_param
    real(kind=CUSTOM_REAL), dimension(:,:,:,:,:),               intent(in)    :: model, gradient, descent
    ! local
    character(len=MAX_LEN_STRING)                                             :: prefix_name
    integer                                                                   :: global_iter

    ! combined index for current frequency stage and inversion iteration
    global_iter = iter_frq * 10000  + iter_inverse

    if (inversion_param%dump_model_at_each_iteration) then
       prefix_name='model_'
       if (mygroup <= 0) call DumpModelArray(model, inversion_param, global_iter, prefix_name)
    endif

    if (inversion_param%dump_gradient_at_each_iteration) then
       prefix_name='gradient'
       if (mygroup <= 0) call DumpArray(gradient, inversion_param, global_iter, prefix_name)
    endif

    if (inversion_param%dump_descent_direction_at_each_iteration) then
       prefix_name='descent_direction'
       if (mygroup <= 0) call DumpArray(descent, inversion_param, global_iter, prefix_name)
    endif

    if (inversion_param%dump_image_at_each_iteration) then
      call DumpImage(inversion_param, global_iter)
    endif

  end subroutine DumpArraysMeshSpecfem
  !------------------------------------------------------------------------------------------------------------
  subroutine DumpImage(inversion_param,current_iter)

    use specfem_par, only: it
    use image_PNM_par, only: PNM_IMAGE_TYPE,PNM_COLOR_PALETTE

    implicit none
    type(inver),                                                intent(in)    :: inversion_param
    integer,                                                    intent(in)    :: current_iter
    integer                                                                   :: i,it_bak

    ! backup specfem it
    it_bak = it

    ! to name output files
    it = current_iter

    ! sets colormap (1==red-blue,2==rainbow,3==tomo,4==terrain)
    PNM_COLOR_PALETTE = 4

    do i = 1,inversion_param%NinvPar
      select case(trim(inversion_param%param_inv_name(i)))
      case ('vp')
        PNM_IMAGE_TYPE = 5
        call write_PNM_create_image()
      case ('vs')
        PNM_IMAGE_TYPE = 6
        call write_PNM_create_image()
      case ('rho')
        PNM_IMAGE_TYPE = 7
        call write_PNM_create_image()
      end select
    enddo

    ! resets specfem it
    it = it_bak

  end subroutine DumpImage
  !------------------------------------------------------------------------------------------------------------
  subroutine DumpArray(field, inversion_param, current_iter, current_name)

    real(kind=CUSTOM_REAL), dimension(:,:,:,:,:),               intent(in)    :: field
    integer,                                                    intent(in)    :: current_iter
    type(inver),                                                intent(in)    :: inversion_param
    character(len=MAX_LEN_STRING),                              intent(in)    :: current_name
    character(len=MAX_LEN_STRING)                                             :: file_name, file_prefix, file_sufix
    integer                                                                   :: i

    do i = 1, inversion_param%NinvPar

       file_prefix = trim(LOCAL_PATH)
       write(file_sufix,'(a5,i5.5,a4)') 'iter_',current_iter,'.bin'
       file_sufix = trim(inversion_param%param_inv_name(i))//'_'//trim(file_sufix)
       write(file_name,'(a5,i6.6)') '/proc', myrank

       file_name = trim(file_prefix)//trim(file_name)//'_'//trim(current_name)//trim(file_sufix)

       if (DEBUG_MODE) then
          write(IIDD,*) 'Dump : ', trim(file_name)
       endif

       open(888,file=trim(file_name),form='unformatted')
       write(888) field(:,:,:,:,i)
       close(888)

    enddo

  end subroutine DumpArray
  !------------------------------------------------------------------------------------------------------------
  subroutine DumpModelArray(field, inversion_param, current_iter, current_name)

    real(kind=CUSTOM_REAL), dimension(:,:,:,:,:),               intent(in)    :: field
    integer,                                                    intent(in)    :: current_iter
    type(inver),                                                intent(in)    :: inversion_param
    character(len=MAX_LEN_STRING),                              intent(in)    :: current_name
    character(len=MAX_LEN_STRING)                                             :: file_name, file_prefix, file_sufix
    integer                                                                   :: i

    do i = 1, inversion_param%NinvPar

      file_prefix = trim(LOCAL_PATH)
      write(file_sufix,'(a5,i5.5,a4)') 'iter_',current_iter,'.bin'
      file_sufix = trim(inversion_param%param_inv_name(i))//'_'//trim(file_sufix)
      write(file_name,'(a5,i6.6)') '/proc', myrank

      file_name = trim(file_prefix)//trim(file_name)//'_'//trim(current_name)//trim(file_sufix)

      if (DEBUG_MODE) then
        write(IIDD,*) 'Dump : ', trim(file_name)
      endif

      open(888,file=trim(file_name),form='unformatted')

      select case(inversion_param%parameter_metric)
      case(0)  !! directly the parameter P
        write(888) field(:,:,:,:,i)
      case(1)  !! P / Pref
        write(888) field(:,:,:,:,i)
      case(2)  !! log(P)
        write(888) exp(field(:,:,:,:,i))
      case(3) !! log(P/Pref)
        write(888) exp(field(:,:,:,:,i))
      end select

      close(888)

    enddo

  end subroutine DumpModelArray


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!--------------------------------------------------------------------------------------------------------------------
!>  ! inport acoustic model from FD grid  store in current rho kappa and  mass matrix
!--------------------------------------------------------------------------------------------------------------------
  subroutine import_FD_model_ACOUSTIC(fd_grid)

    use specfem_par, only: myrank, xstore, ystore, zstore,  kappastore, rhostore, ibool, &
         NGLLX, NGLLY, NGLLZ, NSPEC_AB, MAX_STRING_LEN

    implicit none

    character(len=MAX_STRING_LEN),            intent(in)  :: fd_grid
    character(len=MAX_STRING_LEN)                         :: rho_file, vp_file
    real(kind=CUSTOM_REAL)                                :: ox_fd, oy_fd, oz_fd
    real(kind=CUSTOM_REAL)                                :: hx_fd, hy_fd, hz_fd
    integer                                               :: nx_fd, ny_fd, nz_fd
    real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: vp_fd, rho_fd
    integer                                               :: i, j, k, ispec, iglob
    integer                                               :: ii, jj, kk
    real(kind=CUSTOM_REAL)                                :: xmin, xmax, ymin, ymax, zmin, zmax
    real(kind=CUSTOM_REAL)                                :: xmin_glob, xmax_glob
    real(kind=CUSTOM_REAL)                                :: ymin_glob, ymax_glob
    real(kind=CUSTOM_REAL)                                :: zmin_glob, zmax_glob

    !! READ FD MODEL :: TODO only the group 0 must read this and then bcast to other
    if (myrank == 0) then
       open(4444, file=trim(fd_grid))
       read(4444,*) ox_fd, oy_fd, oz_fd
       read(4444,*) nx_fd, ny_fd, nz_fd
       read(4444,*) hx_fd, hy_fd, hz_fd
       read(4444,'(a)') rho_file
       read(4444,'(a)') vp_file
       close(4444)

       allocate(vp_fd(nx_fd,ny_fd,nz_fd), rho_fd(nx_fd,ny_fd,nz_fd),stat=ier)
       if (ier /= 0) call exit_MPI_without_rank('error allocating array 307')

       open(4444,file=trim(rho_file),access='direct',recl=CUSTOM_REAL*nx_fd*ny_fd*nz_fd)
       read(4444,rec=1) rho_fd
       close(4444)

       open(4444,file=trim(vp_file),access='direct',recl=CUSTOM_REAL*nx_fd*ny_fd*nz_fd)
       read(4444,rec=1) vp_fd
       close(4444)

    endif

    call bcast_all_singlei(nx_fd)
    call bcast_all_singlei(ny_fd)
    call bcast_all_singlei(nz_fd)

    call bcast_all_singlecr(ox_fd)
    call bcast_all_singlecr(oy_fd)
    call bcast_all_singlecr(oz_fd)

    call bcast_all_singlecr(hx_fd)
    call bcast_all_singlecr(hy_fd)
    call bcast_all_singlecr(hz_fd)

    if (myrank > 0) then
      allocate(vp_fd(nx_fd,ny_fd,nz_fd), rho_fd(nx_fd,ny_fd,nz_fd),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 308')
    endif

    call bcast_all_cr(rho_fd,nx_fd*ny_fd*nz_fd)
    call bcast_all_cr(vp_fd,nx_fd*ny_fd*nz_fd)

    !! check if fd model is bigger than SEM model ------------------------
    ! get mesh properties
    xmin=minval(xstore)
    xmax=maxval(xstore)
    ymin=minval(ystore)
    ymax=maxval(ystore)
    zmin=minval(zstore)
    zmax=maxval(zstore)

    call min_all_all_cr(xmin,xmin_glob)
    call max_all_all_cr(xmax,xmax_glob)

    call min_all_all_cr(ymin,ymin_glob)
    call max_all_cr(ymax,ymax_glob)

    call min_all_all_cr(zmin,zmin_glob)
    call max_all_all_cr(zmax,zmax_glob)


    if (ox_fd > xmin_glob ) stop " Error in fd input model : x origin  begin after sem model  "
    if (oy_fd > ymin_glob ) stop " Error in fd input model : y origin  begin after sem model  "
    if (oz_fd > zmin_glob ) stop " Error in fd input model : z origin  begin after sem model  "

    if (ox_fd + (nx_fd-1)*hx_fd < xmax_glob) stop " Error in fd input model : x end before sem model "
    if (oy_fd + (ny_fd-1)*hy_fd < ymax_glob) stop " Error in fd input model : y end before sem model "
    if (oz_fd + (nz_fd-1)*hz_fd < zmax_glob) stop " Error in fd input model : z end before sem model "


    !! PROJECT FD MODEL IN SEM GRID
    do ispec=1,NSPEC_AB
       do  k=1,NGLLZ
          do j=1,NGLLY
             do i=1,NGLLX
                iglob=ibool(i,j,k,ispec)

                !! nearest neighbor
                ii = 1+ floor( (xstore(iglob) - ox_fd)/hx_fd)
                jj = 1+ floor( (ystore(iglob) - oy_fd)/hy_fd)
                kk = 1+ floor( (zstore(iglob) - oz_fd)/hz_fd)

                rhostore(i,j,k,ispec) = rho_fd(ii,jj,kk)
                kappastore(i,j,k,ispec) = rho_fd(ii,jj,kk) * vp_fd(ii,jj,kk)**2

             enddo
          enddo
       enddo
    enddo

    deallocate(vp_fd, rho_fd)

  end subroutine import_FD_model_ACOUSTIC

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!--------------------------------------------------------------------------------------------------------------------
!  ! inport model from FD grid  store in current kappa, mu, mass matrix
!--------------------------------------------------------------------------------------------------------------------
  subroutine import_FD_model_Elastic_ISO(fd_grid)

    use specfem_par, only: myrank, xstore, ystore, zstore,  kappastore, mustore, rhostore, ibool, &
         NGLLX, NGLLY, NGLLZ, NSPEC_AB, FOUR_THIRDS, MAX_STRING_LEN
    use specfem_par_elastic, only: rho_vp, rho_vs

    implicit none

    character(len=MAX_STRING_LEN),            intent(in)  :: fd_grid
    character(len=MAX_STRING_LEN)                         :: rho_file, vp_file, vs_file
    real(kind=CUSTOM_REAL)                                :: ox_fd, oy_fd, oz_fd
    real(kind=CUSTOM_REAL)                                :: hx_fd, hy_fd, hz_fd
    integer                                               :: nx_fd, ny_fd, nz_fd
    real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: vp_fd, vs_fd, rho_fd
    integer                                               :: i, j, k, ispec, iglob
    integer                                               :: ii, jj, kk
    real(kind=CUSTOM_REAL)                                :: xmin, xmax, ymin, ymax, zmin, zmax
    real(kind=CUSTOM_REAL)                                :: xmin_glob, xmax_glob
    real(kind=CUSTOM_REAL)                                :: ymin_glob, ymax_glob
    real(kind=CUSTOM_REAL)                                :: zmin_glob, zmax_glob
    real(kind=CUSTOM_REAL)                                :: xp, yp, zp
    real(kind=CUSTOM_REAL)                                :: rh_interp, vp_interp, vs_interp

    !! READ FD MODEL :: TODO only the group 0 must read this and then bcast to other
    if (myrank == 0) then
       open(4444, file=trim(fd_grid))
       read(4444,*) ox_fd, oy_fd, oz_fd
       read(4444,*) nx_fd, ny_fd, nz_fd
       read(4444,*) hx_fd, hy_fd, hz_fd
       read(4444,'(a)') rho_file
       read(4444,'(a)') vp_file
       read(4444,'(a)') vs_file
       close(4444)

       allocate(vp_fd(nx_fd,ny_fd,nz_fd), vs_fd(nx_fd,ny_fd,nz_fd), rho_fd(nx_fd,ny_fd,nz_fd),stat=ier)
       if (ier /= 0) call exit_MPI_without_rank('error allocating array 309')

       open(4444,file=trim(rho_file),access='direct',recl=CUSTOM_REAL*nx_fd*ny_fd*nz_fd)
       read(4444,rec=1) rho_fd
       close(4444)

       open(4444,file=trim(vp_file),access='direct',recl=CUSTOM_REAL*nx_fd*ny_fd*nz_fd)
       read(4444,rec=1) vp_fd
       close(4444)

       open(4444,file=trim(vs_file),access='direct',recl=CUSTOM_REAL*nx_fd*ny_fd*nz_fd)
       read(4444,rec=1) vs_fd
       close(4444)

    endif

    call bcast_all_singlei(nx_fd)
    call bcast_all_singlei(ny_fd)
    call bcast_all_singlei(nz_fd)

    call bcast_all_singlecr(ox_fd)
    call bcast_all_singlecr(oy_fd)
    call bcast_all_singlecr(oz_fd)

    call bcast_all_singlecr(hx_fd)
    call bcast_all_singlecr(hy_fd)
    call bcast_all_singlecr(hz_fd)

    if (myrank > 0) then
      allocate(vp_fd(nx_fd,ny_fd,nz_fd), vs_fd(nx_fd,ny_fd,nz_fd), rho_fd(nx_fd,ny_fd,nz_fd),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 310')
    endif

    call bcast_all_cr(rho_fd,nx_fd*ny_fd*nz_fd)
    call bcast_all_cr(vp_fd,nx_fd*ny_fd*nz_fd)
    call bcast_all_cr(vs_fd,nx_fd*ny_fd*nz_fd)

    !! check if fd model is bigger than SEM model ------------------------
    ! get mesh properties
    xmin=minval(xstore)
    xmax=maxval(xstore)
    ymin=minval(ystore)
    ymax=maxval(ystore)
    zmin=minval(zstore)
    zmax=maxval(zstore)

    call min_all_all_cr(xmin,xmin_glob)
    call max_all_all_cr(xmax,xmax_glob)

    call min_all_all_cr(ymin,ymin_glob)
    call max_all_all_cr(ymax,ymax_glob)

    call min_all_all_cr(zmin,zmin_glob)
    call max_all_all_cr(zmax,zmax_glob)


    if (ox_fd > xmin_glob ) stop " Error in fd input model : x origin  begin after sem model  "
    if (oy_fd > ymin_glob ) stop " Error in fd input model : y origin  begin after sem model  "
    if (oz_fd > zmin_glob ) stop " Error in fd input model : z origin  begin after sem model  "

    if (ox_fd + (nx_fd-1)*hx_fd < xmax_glob) stop " Error in fd input model : x end before sem model "
    if (oy_fd + (ny_fd-1)*hy_fd < ymax_glob) stop " Error in fd input model : y end before sem model "
    if (oz_fd + (nz_fd-1)*hz_fd < zmax_glob) stop " Error in fd input model : z end before sem model "


    !! PROJECT FD MODEL IN SEM GRID
    do ispec=1,NSPEC_AB
       do  k=1,NGLLZ
          do j=1,NGLLY
             do i=1,NGLLX
                iglob=ibool(i,j,k,ispec)

                !! nearest neighbor
                ii = 1+ floor( (xstore(iglob) - ox_fd)/hx_fd)
                jj = 1+ floor( (ystore(iglob) - oy_fd)/hy_fd)
                kk = 1+ floor( (zstore(iglob) - oz_fd)/hz_fd)

                xp=xstore(iglob)
                yp=ystore(iglob)
                zp=zstore(iglob)

                !! trilinear interpolation
                call Get_value_by_trilinear_interp(rh_interp, xp, yp, zp, rho_fd, &
                     nx_fd, ny_fd, nz_fd, ox_fd, oy_fd, oz_fd, hx_fd, hy_fd, hz_fd)
                call Get_value_by_trilinear_interp(vp_interp, xp, yp, zp, vp_fd, &
                     nx_fd, ny_fd, nz_fd, ox_fd, oy_fd, oz_fd, hx_fd, hy_fd, hz_fd)
                call Get_value_by_trilinear_interp(vs_interp, xp, yp, zp, vs_fd, &
                     nx_fd, ny_fd, nz_fd, ox_fd, oy_fd, oz_fd, hx_fd, hy_fd, hz_fd)

                rhostore(i,j,k,ispec) = rh_interp
                kappastore(i,j,k,ispec) = rh_interp * ( (vp_interp**2) -&
                     FOUR_THIRDS*(vs_interp**2))
                mustore(i,j,k,ispec) = rh_interp * (vs_interp**2)
                rho_vs(i,j,k,ispec) = rh_interp * vs_interp
                rho_vp(i,j,k,ispec) = rh_interp * vp_interp

             enddo
          enddo
       enddo
    enddo

    deallocate(vp_fd, vs_fd, rho_fd)

  end subroutine import_FD_model_Elastic_ISO

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!--------------------------------------------------------------------------------------------------------------------
!  ! inport model from FD grid  store in current anisotropic specfem model
!--------------------------------------------------------------------------------------------------------------------
  subroutine import_FD_model_ANISO(fd_grid, inversion_param)

    use constants, only: MAX_STRING_LEN
    use specfem_par, only: myrank, xstore, ystore, zstore, rhostore, ibool, &
         NGLLX, NGLLY, NGLLZ, NSPEC_AB, MAX_STRING_LEN
    !use specfem_par_elastic, only: rho_vp, rho_vs  here i commented because they are already called in head of module
    use interpolation_mod, only: trilin_interp

    implicit none
    type(inver),                           intent(inout)    :: inversion_param
    character(len=MAX_STRING_LEN),            intent(in)    :: fd_grid
    character(len=MAX_STRING_LEN)                           :: model_file, type_model
    real(kind=CUSTOM_REAL)                                  :: ox_fd, oy_fd, oz_fd
    real(kind=CUSTOM_REAL)                                  :: hx_fd, hy_fd, hz_fd
    integer                                                 :: nx_fd, ny_fd, nz_fd
    integer                                                 :: nb_model_to_read
    real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: model_fd
    integer                                                 :: i, j, k, ispec, iglob, ipar
    integer                                                 :: ii, jj, kk
    real(kind=CUSTOM_REAL)                                  :: xmin, xmax, ymin, ymax, zmin, zmax
    real(kind=CUSTOM_REAL)                                  :: xmin_glob, xmax_glob
    real(kind=CUSTOM_REAL)                                  :: ymin_glob, ymax_glob
    real(kind=CUSTOM_REAL)                                  :: zmin_glob, zmax_glob
    real(kind=CUSTOM_REAL)                                  :: rho, vp, vs, gm, ep, de, lambda, mu, kappa, val
    real(kind=CUSTOM_REAL)                                  :: phi, theta !! fast axis direction TrISO
    real(kind=CUSTOM_REAL)                                  :: CIJ(6,6)
    real(kind=CUSTOM_REAL)                                  :: c11, c22, c33
    real(kind=CUSTOM_REAL)                                  :: c44, c55, c66
    real(kind=CUSTOM_REAL)                                  :: c12, c13, c23
    real(kind=CUSTOM_REAL), dimension(22)                   :: cval
    real(kind=CUSTOM_REAL), dimension(:),  allocatable      :: xcrd_fd, ycrd_fd, zcrd_fd

    !! READ FD MODEL :: TODO only the group 0 must read this and then bcast to other
    if (myrank == 0) then
       open(4444, file=trim(fd_grid))
       read(4444,*) ox_fd, oy_fd, oz_fd
       read(4444,*) nx_fd, ny_fd, nz_fd
       read(4444,*) hx_fd, hy_fd, hz_fd
       read(4444,*) nb_model_to_read
       read(4444, '(a)' ) type_model
       read(4444,'(a)') model_file
       close(4444)

       allocate(model_fd(nx_fd,ny_fd,nz_fd, nb_model_to_read),stat=ier)
       if (ier /= 0) call exit_MPI_without_rank('error allocating array 311')

       open(4444,file=trim(model_file),access='direct',recl=CUSTOM_REAL*nx_fd*ny_fd*nz_fd)
       do i=1,nb_model_to_read
          read(4444,rec=i) model_fd(:,:,:,i)
       enddo
       close(4444)
    endif

    call bcast_all_singlei(nx_fd)
    call bcast_all_singlei(ny_fd)
    call bcast_all_singlei(nz_fd)
    call bcast_all_singlei(nb_model_to_read)

    call bcast_all_singlecr(ox_fd)
    call bcast_all_singlecr(oy_fd)
    call bcast_all_singlecr(oz_fd)

    call bcast_all_singlecr(hx_fd)
    call bcast_all_singlecr(hy_fd)
    call bcast_all_singlecr(hz_fd)

    call bcast_all_ch_array(type_model,1,MAX_STRING_LEN)

    if (myrank > 0) then
      allocate(model_fd(nx_fd, ny_fd, nz_fd, nb_model_to_read),stat=ier)
      if (ier /= 0) call exit_MPI_without_rank('error allocating array 312')
    endif
    call bcast_all_cr(model_fd,nx_fd*ny_fd*nz_fd*nb_model_to_read)

    if (myrank == 0) write(INVERSE_LOG_FILE, *) ' MPI Bcast full elastic tensor, size : ',nx_fd*ny_fd*nz_fd*nb_model_to_read

    !! check if fd model is bigger than SEM model ------------------------
    ! get mesh properties
    xmin=minval(xstore)
    xmax=maxval(xstore)
    ymin=minval(ystore)
    ymax=maxval(ystore)
    zmin=minval(zstore)
    zmax=maxval(zstore)

    call min_all_all_cr(xmin,xmin_glob)
    call max_all_all_cr(xmax,xmax_glob)

    call min_all_all_cr(ymin,ymin_glob)
    call max_all_all_cr(ymax,ymax_glob)

    call min_all_all_cr(zmin,zmin_glob)
    call max_all_all_cr(zmax,zmax_glob)

    if (ox_fd > xmin_glob ) stop " Error in fd input model : x origin  begin after sem model  "
    if (oy_fd > ymin_glob ) stop " Error in fd input model : y origin  begin after sem model  "
    if (oz_fd > zmin_glob ) stop " Error in fd input model : z origin  begin after sem model  "

    if (ox_fd + (nx_fd-1)*hx_fd < xmax_glob) stop " Error in fd input model : x end before sem model "
    if (oy_fd + (ny_fd-1)*hy_fd < ymax_glob) stop " Error in fd input model : y end before sem model "
    if (oz_fd + (nz_fd-1)*hz_fd < zmax_glob) stop " Error in fd input model : z end before sem model "


    select case(trim(adjustl(type_model)))

    case('VTI') ! rho, vp, vs, epsillon, delta, gamma (oil-industry like)

       if (myrank == 0) then
          write(INVERSE_LOG_FILE,*)
          write(INVERSE_LOG_FILE,*) '*********************************************'
          write(INVERSE_LOG_FILE,*) '***         READING FD VTI MODEL          ***'
          write(INVERSE_LOG_FILE,*) '*********************************************'
          write(INVERSE_LOG_FILE,*)
          write(INVERSE_LOG_FILE,*) ' Param family :', trim(inversion_param%parameter_family_name)
          write(INVERSE_LOG_FILE,*)
       endif

       !! PROJECT FD MODEL IN SEM GRID
       do ispec=1,NSPEC_AB
          do  k=1,NGLLZ
             do j=1,NGLLY
                do i=1,NGLLX
                   iglob=ibool(i,j,k,ispec)

                   !! nearest neighbor
                   ii = 1+ floor( (xstore(iglob) - ox_fd)/hx_fd)
                   jj = 1+ floor( (ystore(iglob) - oy_fd)/hy_fd)
                   kk = 1+ floor( (zstore(iglob) - oz_fd)/hz_fd)

                   rho= model_fd(ii,jj,kk,1)
                   vp = model_fd(ii,jj,kk,2)
                   vs = model_fd(ii,jj,kk,3)
                   ep = model_fd(ii,jj,kk,4)
                   gm = model_fd(ii,jj,kk,5)
                   de = model_fd(ii,jj,kk,6)

                   !
                   rhostore(i,j,k,ispec) = rho
                   rho_vs(i,j,k,ispec)   = rho * vs
                   rho_vp(i,j,k,ispec)   = rho * vp

                   c11store(i,j,k,ispec) = (1._CUSTOM_REAL+2._CUSTOM_REAL*ep)* rho * vp * vp
                   c12store(i,j,k,ispec) = -2*(2*gm+1)*vs*vs*rho + (1._CUSTOM_REAL+2._CUSTOM_REAL*ep)* rho * vp * vp
                   c13store(i,j,k,ispec) =  rho * (sqrt( (vp**2-vs**2)*((1+2*de)*vp**2-vs**2)) - vs**2)
                   c14store(i,j,k,ispec) =  0._CUSTOM_REAL
                   c15store(i,j,k,ispec) =  0._CUSTOM_REAL
                   c16store(i,j,k,ispec) =  0._CUSTOM_REAL
                   c22store(i,j,k,ispec) = (1._CUSTOM_REAL+2._CUSTOM_REAL*ep)* rho * vp * vp
                   c23store(i,j,k,ispec) = rho * ( sqrt( (vp**2-vs**2) * (  (1+2*de)*vp**2-vs**2) ) - vs**2 )
                   c24store(i,j,k,ispec) =  0._CUSTOM_REAL
                   c25store(i,j,k,ispec) =  0._CUSTOM_REAL
                   c26store(i,j,k,ispec) =  0._CUSTOM_REAL
                   c33store(i,j,k,ispec) = rho *vp *vp
                   c34store(i,j,k,ispec) =  0._CUSTOM_REAL
                   c35store(i,j,k,ispec) =  0._CUSTOM_REAL
                   c36store(i,j,k,ispec) =  0._CUSTOM_REAL
                   c44store(i,j,k,ispec) = rho*vs*vs
                   c45store(i,j,k,ispec) =  0._CUSTOM_REAL
                   c46store(i,j,k,ispec) =  0._CUSTOM_REAL
                   c55store(i,j,k,ispec) = rho*vs*vs
                   c56store(i,j,k,ispec) =  0._CUSTOM_REAL
                   c66store(i,j,k,ispec) = (2*gm+1)*vs*vs*rho

                enddo
             enddo
          enddo
       enddo

    case ('TRISO') !! tansverse isotropic

       !! PROJECT FD MODEL IN SEM GRID
       do ispec=1,NSPEC_AB
          do  k=1,NGLLZ
             do j=1,NGLLY
                do i=1,NGLLX
                   iglob=ibool(i,j,k,ispec)

                   !! nearest neighbor
                   ii = 1+ floor( (xstore(iglob) - ox_fd)/hx_fd)
                   jj = 1+ floor( (ystore(iglob) - oy_fd)/hy_fd)
                   kk = 1+ floor( (zstore(iglob) - oz_fd)/hz_fd)


                   CIJ(:,:)=0.
                   rho   = model_fd(ii,jj,kk,1)
                   vp    = model_fd(ii,jj,kk,2)
                   vs    = model_fd(ii,jj,kk,3)
                   ep    = model_fd(ii,jj,kk,4)
                   gm    = model_fd(ii,jj,kk,5)
                   de    = model_fd(ii,jj,kk,6)
                   phi   = model_fd(ii,jj,kk,7)
                   theta = model_fd(ii,jj,kk,8)

                   ! Set rhostore and Stacey's arrays
                   rhostore(i,j,k,ispec) = rho
                   rho_vs(i,j,k,ispec)   = rho * vs
                   rho_vp(i,j,k,ispec)   = rho * vp

                   ! define tensor in aniso axis reference
                   CIJ(1,1)=rho*vp*vp  !! will be overwrite
                   CIJ(4,4)=rho*vs*vs
                   CIJ(1,3)=CIJ(1,1)-2.*CIJ(4,4)+(de-ep)*CIJ(1,1)+2.*gm*CIJ(4,4)
                   CIJ(2,3)=CIJ(1,3)
                   CIJ(3,3)=CIJ(1,1)*(1.-ep)
                   CIJ(1,1)=CIJ(1,1)*(1.+ep)
                   CIJ(2,2)=CIJ(1,1)
                   CIJ(6,6)=CIJ(4,4)*(1.+gm)
                   CIJ(4,4)=CIJ(4,4)*(1.-gm)
                   CIJ(5,5)=CIJ(4,4)
                   CIJ(1,2)=CIJ(1,1)-2.*CIJ(6,6)

                   !! rotate and store
                   call RotateTransIso(phi, theta, CIJ)

                   !! store CIJ
                   c11store(i,j,k,ispec) = CIJ(1,1)
                   c12store(i,j,k,ispec) = CIJ(1,2)
                   c13store(i,j,k,ispec) = CIJ(1,3)
                   c14store(i,j,k,ispec) = CIJ(1,4)
                   c15store(i,j,k,ispec) = CIJ(1,5)
                   c16store(i,j,k,ispec) = CIJ(1,6)
                   c22store(i,j,k,ispec) = CIJ(2,2)
                   c23store(i,j,k,ispec) = CIJ(2,3)
                   c24store(i,j,k,ispec) = CIJ(2,4)
                   c25store(i,j,k,ispec) = CIJ(2,5)
                   c26store(i,j,k,ispec) = CIJ(2,6)
                   c33store(i,j,k,ispec) = CIJ(3,3)
                   c34store(i,j,k,ispec) = CIJ(3,4)
                   c35store(i,j,k,ispec) = CIJ(3,5)
                   c36store(i,j,k,ispec) = CIJ(3,6)
                   c44store(i,j,k,ispec) = CIJ(4,4)
                   c45store(i,j,k,ispec) = CIJ(4,5)
                   c46store(i,j,k,ispec) = CIJ(4,6)
                   c55store(i,j,k,ispec) = CIJ(5,5)
                   c56store(i,j,k,ispec) = CIJ(5,6)
                   c66store(i,j,k,ispec) = CIJ(6,6)

                enddo
             enddo
          enddo

       enddo

    case ('FULL') !! fully anisotropic elastic tensor

       !! DEFINE TOMO GRID LOOKUP TABLES
       if (.not. allocated(xcrd_fd)) then
         allocate(xcrd_fd(nx_fd),stat=ier)
         if (ier /= 0) call exit_MPI_without_rank('error allocating array 313')
       endif
       if (.not. allocated(ycrd_fd)) then
         allocate(ycrd_fd(ny_fd),stat=ier)
         if (ier /= 0) call exit_MPI_without_rank('error allocating array 314')
       endif
       if (.not. allocated(zcrd_fd)) then
         allocate(zcrd_fd(nz_fd),stat=ier)
         if (ier /= 0) call exit_MPI_without_rank('error allocating array 315')
       endif
       do i=1,nx_fd
          xcrd_fd(i) = ox_fd + hx_fd * real(i-1)
       enddo
       do i=1,ny_fd
          ycrd_fd(i) = oy_fd + hy_fd * real(i-1)
       enddo
       do i=1,nz_fd
          zcrd_fd(i) = oz_fd + hz_fd * real(i-1)
       enddo

       !! PROJECT FD MODEL IN SEM GRID
       do ispec=1,NSPEC_AB
          do  k=1,NGLLZ
             do j=1,NGLLY
                do i=1,NGLLX
                   iglob=ibool(i,j,k,ispec)

                   !! nearest neighbor
                   !ii = 1+ floor( (xstore(iglob) - ox_fd)/hx_fd)
                   !jj = 1+ floor( (ystore(iglob) - oy_fd)/hy_fd)
                   !kk = 1+ floor( (zstore(iglob) - oz_fd)/hz_fd)

                   !! trilinear (not optimized here but can be done if needed,
                   !            (to do so, split trilinear routine such that coefficients
                   !             are computed only once per gll)
                   do ipar=1,22

                      call trilin_interp(xstore(iglob), ystore(iglob), zstore(iglob), &
                                               xcrd_fd,       ycrd_fd,       zcrd_fd, &
                                                 nx_fd,         ny_fd,         nz_fd, &
                                                 model_fd(:,:,:,ipar),           val)
                      ! Get tensors components
                      cval(ipar) = val

                   enddo

                   ! store c_ij
                   c11store(i,j,k,ispec) = cval(1)
                   c12store(i,j,k,ispec) = cval(2)
                   c13store(i,j,k,ispec) = cval(3)
                   c14store(i,j,k,ispec) = cval(4)
                   c15store(i,j,k,ispec) = cval(5)
                   c16store(i,j,k,ispec) = cval(6)
                   c22store(i,j,k,ispec) = cval(7)
                   c23store(i,j,k,ispec) = cval(8)
                   c24store(i,j,k,ispec) = cval(9)
                   c25store(i,j,k,ispec) = cval(10)
                   c26store(i,j,k,ispec) = cval(11)
                   c33store(i,j,k,ispec) = cval(12)
                   c34store(i,j,k,ispec) = cval(13)
                   c35store(i,j,k,ispec) = cval(14)
                   c36store(i,j,k,ispec) = cval(15)
                   c44store(i,j,k,ispec) = cval(16)
                   c45store(i,j,k,ispec) = cval(17)
                   c46store(i,j,k,ispec) = cval(18)
                   c55store(i,j,k,ispec) = cval(19)
                   c56store(i,j,k,ispec) = cval(20)
                   c66store(i,j,k,ispec) = cval(21)

                   !* Useless but easier to read and to use after
                   c11 =  cval(1);     c22 =  cval(7);    c33 = cval(12);
                   c12 =  cval(2);     c13 =  cval(3);    c23 = cval(8);
                   c44 = cval(16);     c55 = cval(19);    c66 = cval(21);
                   rho = cval(22);

                   !* Store rho
                   rhostore(i,j,k,ispec) = rho

                   !* Stacey (need to determine isotropic part of tensor to determine vp and vs)
                   !  I should look at Browaeys and Chevrot (2004) or Norris (2006)
                   !  For now, use Fedorov (1968) approximation
                   kappa  = ((c11 + c22 + c33) + 2._CUSTOM_REAL * (c12 + c13 + c23)) / 9._CUSTOM_REAL
                   mu     =  ( 2._CUSTOM_REAL * (c11 + c22 + c33 - c12 - c23 - c13) &
                             + 6._CUSTOM_REAL * (c44 + c55 + c66) ) / 30._CUSTOM_REAL

                   lambda = kappa - (2._CUSTOM_REAL * mu /3._CUSTOM_REAL)
                   vp     = sqrt((lambda + 2._CUSTOM_REAL * mu) / rho)
                   vs     = sqrt(mu / rho)

                   kappastore(i,j,k,ispec) = kappa
                   mustore(i,j,k,ispec)    = mu

                   !* Store stacey
                   rho_vp(i,j,k,ispec)   = rho * vp
                   rho_vs(i,j,k,ispec)   = rho * vs

                enddo
             enddo
          enddo

       enddo

    case default

       write(*,*) 'ERROR : anisotropic model ', trim(adjustl(type_model)), ' Not yet implemented '
       stop

    end select

    !! free memory
    if (allocated(model_fd)) deallocate(model_fd)

  end subroutine import_FD_model_ANISO

!#################################################################################################################################
!> Rotate transverse isotropic elastic tensor from fast axis direction reference to specfem reference
!################################################################################################################################
  subroutine RotateTransIso(phi, theta, CIJ)

    real(kind=CUSTOM_REAL),              intent(in)    :: phi, theta
    real(kind=CUSTOM_REAL),              intent(inout) :: CIJ(6,6)

    real(kind=CUSTOM_REAL)                             :: a(3,3), m(6,6), C1(6,6), C2(6,6)
    real(kind=CUSTOM_REAL)                             :: cos_phi, sin_phi
    real(kind=CUSTOM_REAL)                             :: cos_theta, sin_theta
    integer                                            :: i1,i2,j1,j2

    cos_phi = cos(phi * PI / 180.)
    sin_phi = sin(phi * PI / 180.)
    cos_theta = cos(theta * PI / 180.)
    sin_theta = sin(theta * PI / 180.)

    a(1,1) = cos_theta*cos_phi
    a(1,2) = -sin_phi
    a(1,3) = sin_theta*cos_phi
    a(2,1) = cos_theta*sin_phi
    a(2,2) = cos_phi
    a(2,3) = sin_theta*sin_phi
    a(3,1) = -sin_theta
    a(3,2) = 0.d0
    a(3,3) = cos_theta

    M(1,1) = a(1,1)**2
    M(1,2) = a(1,2)**2
    M(1,3) = a(1,3)**2
    M(1,4) = 2.d0*a(1,2)*a(1,3)
    M(1,5) = 2.d0*a(1,1)*a(1,3)
    M(1,6) = 2.d0*a(1,1)*a(1,2)
    M(2,1) = a(2,1)**2
    M(2,2) = a(2,2)**2
    M(2,3) = a(2,3)**2
    M(2,4) = 2.d0*a(2,2)*a(2,3)
    M(2,5) = 2.d0*a(2,1)*a(2,3)
    M(2,6) = 2.d0*a(2,1)*a(2,2)
    M(3,1) = a(3,1)**2
    M(3,2) = a(3,2)**2
    M(3,3) = a(3,3)**2
    M(3,4) = 2.d0*a(3,2)*a(3,3)
    M(3,5) = 2.d0*a(3,1)*a(3,3)
    M(3,6) = 2.d0*a(3,1)*a(3,2)
    M(4,1) = a(2,1)*a(3,1)
    M(4,2) = a(2,2)*a(3,2)
    M(4,3) = a(2,3)*a(3,3)
    M(4,4) = a(2,2)*a(3,3)+a(2,3)*a(3,2)
    M(4,5) = a(2,1)*a(3,3)+a(2,3)*a(3,1)
    M(4,6) = a(2,1)*a(3,2)+a(2,2)*a(3,1)
    M(5,1) = a(1,1)*a(3,1)
    M(5,2) = a(1,2)*a(3,2)
    M(5,3) = a(1,3)*a(3,3)
    M(5,4) = a(1,2)*a(3,3)+a(1,3)*a(3,2)
    M(5,5) = a(1,1)*a(3,3)+a(1,3)*a(3,1)
    M(5,6) = a(1,1)*a(3,2)+a(1,2)*a(3,1)
    M(6,1) = a(1,1)*a(2,1)
    M(6,2) = a(1,2)*a(2,2)
    M(6,3) = a(1,3)*a(2,3)
    M(6,4) = a(1,2)*a(2,3)+a(1,3)*a(2,2)
    M(6,5) = a(1,1)*a(2,3)+a(1,3)*a(2,1)
    M(6,6) = a(1,1)*a(2,2)+a(1,2)*a(2,1)

    C1(:,:)=CIJ(:,:)
    do i1 = 1,6
       do i2 = i1,6
          C2(i1,i2) = 0.d0
          do j1 = 1,6
             do j2 = j1,6
                if (j1 /= j2) then
                   C2(i1,i2) = C2(i1,i2)+M(i1,j1)*M(i2,j2)*C1(j1,j2) &
                        +M(i1,j2)*M(i2,j1)*C1(j1,j2)
                else
                   C2(i1,i2) = C2(i1,i2)+M(i1,j1)*M(i2,j2)*C1(j1,j2)
                endif
             enddo
          enddo
       enddo
    enddo
    CIJ(:,:)=C2(:,:)

  end subroutine RotateTransIso


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!--------------------------------------------------------------------------------------------------------------------
! write sem model decomposed in NROC files for VTI  (Thomsen's parameters)
!--------------------------------------------------------------------------------------------------------------------

  subroutine write_vti_sem_model(ifrq)
    integer,                                                       intent(in) :: ifrq
    character(len=256)                                                        :: path_file,name_file
    real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable                   :: wks_model_rho, wks_model_vp, wks_model_vs
    real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable                   :: wks_model_ep, wks_model_de, wks_model_ga

    allocate(wks_model_rho(NGLLX,NGLLY,NGLLZ,NSPEC_AB), stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 316')
    if (ier /= 0) call exit_MPI(myrank,"error allocation wks_model_rho in ReadInputSEMmodel subroutine, IO_model_mod")

    allocate(wks_model_vp(NGLLX,NGLLY,NGLLZ,NSPEC_AB), stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 317')
    if (ier /= 0) call exit_MPI(myrank,"error allocation wks_model_vp in ReadInputSEMmodel subroutine, IO_model_mod")

    allocate(wks_model_vs(NGLLX,NGLLY,NGLLZ,NSPEC_AB), stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 318')
    if (ier /= 0) call exit_MPI(myrank,"error allocation wks_model_vs in ReadInputSEMmodel subroutine, IO_model_mod")

    allocate(wks_model_ep(NGLLX,NGLLY,NGLLZ,NSPEC_AB), stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 319')
    if (ier /= 0) call exit_MPI(myrank,"error allocation wks_model_ep in ReadInputSEMmodel subroutine, IO_model_mod")

    allocate(wks_model_de(NGLLX,NGLLY,NGLLZ,NSPEC_AB), stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 320')
    if (ier /= 0) call exit_MPI(myrank,"error allocation wks_model_ge in ReadInputSEMmodel subroutine, IO_model_mod")

    allocate(wks_model_ga(NGLLX,NGLLY,NGLLZ,NSPEC_AB), stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 321')
    if (ier /= 0) call exit_MPI(myrank,"error allocation wks_model_da in ReadInputSEMmodel subroutine, IO_model_mod")

    wks_model_rho(:,:,:,:) = rhostore(:,:,:,:)
    wks_model_vp(:,:,:,:) = sqrt( c33store(:,:,:,:) / rhostore(:,:,:,:))
    wks_model_vs(:,:,:,:) = sqrt( c44store(:,:,:,:) / rhostore(:,:,:,:))
    wks_model_ep(:,:,:,:) = (c11store(:,:,:,:) - c33store(:,:,:,:)) / (2.*c33store(:,:,:,:))
    wks_model_ga(:,:,:,:) = (c66store(:,:,:,:) - c44store(:,:,:,:)) / (2.*c44store(:,:,:,:))
    wks_model_de(:,:,:,:) = 0.5 * ( (c13store(:,:,:,:) + c44store(:,:,:,:) )**2 - (c33store(:,:,:,:) - c44store(:,:,:,:))**2 )&
         / (c33store(:,:,:,:)*(c33store(:,:,:,:) - c44store(:,:,:,:)))

    !! write output  model
    path_file = trim(LOCAL_PATH) // '/proc'
    write(name_file,'(i6.6,a5,i6.6,a20)') myrank,'_freq', ifrq , '_model_vp_output.bin'
    path_file=(trim(path_file))//trim(name_file)
    open(888,file=trim(path_file),form='unformatted')
    write(888) wks_model_vp
    close(888)

    path_file = trim(LOCAL_PATH) // '/proc'
    write(name_file,'(i6.6,a5,i6.6,a20)') myrank,'_freq', ifrq , '_model_vs_output.bin'
    path_file=(trim(path_file))//trim(name_file)
    open(888,file=trim(path_file),form='unformatted')
    write(888) wks_model_vs
    close(888)

    path_file = trim(LOCAL_PATH) // '/proc'
    write(name_file,'(i6.6,a5,i6.6,a21)') myrank,'_freq', ifrq , '_model_rho_output.bin'
    path_file=(trim(path_file))//trim(name_file)
    open(888,file=trim(path_file),form='unformatted')
    write(888) wks_model_rho
    close(888)

    !! write output  model
    path_file = trim(LOCAL_PATH) // '/proc'
    write(name_file,'(i6.6,a5,i6.6,a20)') myrank,'_freq', ifrq , '_model_ep_output.bin'
    path_file=(trim(path_file))//trim(name_file)
    open(888,file=trim(path_file),form='unformatted')
    write(888) wks_model_ep
    close(888)

    path_file = trim(LOCAL_PATH) // '/proc'
    write(name_file,'(i6.6,a5,i6.6,a20)') myrank,'_freq', ifrq , '_model_ga_output.bin'
    path_file=(trim(path_file))//trim(name_file)
    open(888,file=trim(path_file),form='unformatted')
    write(888) wks_model_ga
    close(888)

    path_file = trim(LOCAL_PATH) // '/proc'
    write(name_file,'(i6.6,a5,i6.6,a20)') myrank,'_freq', ifrq , '_model_de_output.bin'
    path_file=(trim(path_file))//trim(name_file)
    open(888,file=trim(path_file),form='unformatted')
    write(888) wks_model_de
    close(888)

    deallocate(wks_model_rho, wks_model_vp, wks_model_vs,wks_model_ep, wks_model_de, wks_model_ga)

  end subroutine write_vti_sem_model

end module IO_model
