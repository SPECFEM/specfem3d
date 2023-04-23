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

! This module implements kinematic faults: prescribed spatio-temporal slip history
!
! Authors:
! Percy Galvez, Jean-Paul Ampuero, Javier Ruiz, Surendra Somala

module fault_solver_kinematic

  use fault_solver_common

  implicit none

  private

  !Number of time steps defined by the user : NTOUT
  integer, save :: NTOUT,NSNAP

  logical, save :: SIMULATION_TYPE_KIN = .false.

  public :: BC_KINFLT_init, BC_KINFLT_set_all, SIMULATION_TYPE_KIN

contains

!=====================================================================
! BC_KINFLT_init initializes kinematic faults
!
! prname        fault database is read from file prname_fault_db.bin

  subroutine BC_KINFLT_init(prname)

  use specfem_par, only: nt => NSTEP, DTglobal => DT

  use constants, only: myrank,IIN_PAR,IIN_BIN

  implicit none

  character(len=MAX_STRING_LEN), intent(in) :: prname ! 'proc***'

  ! local parameters
  type(bc_dynandkinflt_type),pointer :: bc
  real(kind=CUSTOM_REAL) :: dt_real
  integer :: iflt,ier,dummy_idfault
  integer :: nbfaults,nbfaults_bin
  integer :: size_Kelvin_Voigt
  integer :: rupture_type
  character(len=MAX_STRING_LEN) :: filename

  real(kind=CUSTOM_REAL) :: DUMMY
  NAMELIST / BEGIN_FAULT / dummy_idfault

  ! initializes
  dummy_idfault = 0

  ! note: all processes will open this file
  open(unit=IIN_PAR,file=IN_DATA_FILES(1:len_trim(IN_DATA_FILES))//'Par_file_faults',status='old',iostat=ier)

  ! checks if file exists
  if (ier /= 0) then
    if (myrank == 0) write(IMAIN,*) '  no kinematic faults'
    close(IIN_PAR)
    ! all done
    return
  endif

  ! number of faults
  read(IIN_PAR,*) nbfaults
  if (nbfaults == 0) then
    !if (myrank == 0) write(IMAIN,*) 'No faults found in file DATA/Par_file_faults'
    ! nothing to do
    return
  endif

  ! checks if binary fault file was created
  filename = prname(1:len_trim(prname))//'fault_db.bin'
  open(unit=IIN_BIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
  if (ier /= 0) then
    write(IMAIN,*) 'Fatal error: file ',trim(filename),' not found. Abort'
    call exit_MPI(myrank,'Error: fault file fault_db.bin not found')
  endif

  ! Reading etas of each fault
  ! Skip reading viscosity eta of each fault, will be done with binary file
  do iflt = 1,nbfaults
    read(IIN_PAR,*) ! eta
  enddo

  ! fault rupture type: 1 = dyn 2 = kin
  read(IIN_PAR,*) rupture_type

  ! fault simulation type == 2 for kinematic rupture simulation
  ! checks if anything to do
  if (rupture_type /= 2) then
    close(IIN_BIN)
    close(IIN_PAR)
    ! all done
    return
  endif

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  incorporating kinematic rupture simulation'
    write(IMAIN,*) '  found ', nbfaults, ' fault(s) in file DATA/Par_file_faults'
  endif

  ! sets kinematic rupture flag
  SIMULATION_TYPE_KIN = .true.

  ! reads parameters:
  !   NTOUT : Number of time steps
  !   NTSNAP: time interation of snapshots
  !   V_HEALING (-1 : Healing off)
  !   V_RUPT
  read(IIN_PAR,*) NTOUT
  read(IIN_PAR,*) NSNAP
  read(IIN_PAR,*) DUMMY
  read(IIN_PAR,*) DUMMY

  ! from binary fault file
  read(IIN_BIN) nbfaults_bin ! should be the same as in IIN_PAR

  ! checks
  if (nbfaults /= nbfaults_bin) then
    print *,'Error: number of faults ',nbfaults,' in Par_file_faults should match number stored in fault_db.bin: ',nbfaults_bin
    print *,'Please check setup and rerun simulation...'
    call exit_MPI(myrank,'Error invalid number of faults in Par_file_faults')
  endif

  ! saves number of faults
  Nfaults = nbfaults
  allocate( faults(nbfaults) ,stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1993')

  dt_real = real(DTglobal)

  ! reads in fault parameters
  do iflt = 1,nbfaults
    ! example lines:
    ! &BEGIN_FAULT /
    ! &STRESS_TENSOR Sigma=0e0,0e0,0e0,0e0,0e0,0e0/
    ! ..
    read(IIN_PAR,nml=BEGIN_FAULT,end=100)

    ! initializes fault
    bc => faults(iflt)
    call init_one_fault(bc,IIN_BIN,IIN_PAR,dt_real,nt,iflt)
  enddo

  ! close files
  close(IIN_BIN)
  close(IIN_PAR)

  ! reads Kelvin-Voigt parameters
  filename = prname(1:len_trim(prname))//'Kelvin_voigt_eta.bin'
  open(unit=IIN_BIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
  if (ier /= 0) then
    write(IMAIN,*) 'Fatal error: file ',trim(filename),' not found. Abort'
    call exit_MPI(myrank,'Error opening file Kelvin_voigt_eta.bin')
  endif

  ! reads in values
  read(IIN_BIN) size_Kelvin_Voigt
  if (size_Kelvin_Voigt > 0) then
    allocate(Kelvin_Voigt_eta(size_Kelvin_Voigt),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1362')
    read(IIN_BIN) Kelvin_Voigt_eta
  endif
  close(IIN_BIN)
  ! sets flag if this process has damping on fault elements
  if (allocated(Kelvin_Voigt_eta)) then
    USE_KELVIN_VOIGT_DAMPING = .true.
  else
    USE_KELVIN_VOIGT_DAMPING = .false.
  endif

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  number of time steps                 = ',NTOUT
    write(IMAIN,*) '  time interaction of snapshots NTSNAP = ',NSNAP
    call flush_IMAIN()
  endif

  return

100 if (myrank == 0)write(IMAIN,*) 'Fatal error: did not find BEGIN_FAULT input block in file DATA/Par_file_faults. Abort.'
    call exit_MPI(myrank,'Error: BEGIN_FAULT entry missing')

  end subroutine BC_KINFLT_init


!---------------------------------------------------------------------

  subroutine init_one_fault(bc,IIN_BIN,IIN_PAR,dt_real,NT,iflt)

  implicit none

  type(bc_dynandkinflt_type), intent(inout) :: bc
  integer, intent(in)                 :: IIN_BIN,IIN_PAR
  integer, intent(in)                 :: NT,iflt
  real(kind=CUSTOM_REAL), intent(in)  :: dt_real

  ! local parameters
  real(kind=CUSTOM_REAL) :: kindt
  integer :: ier,recordlength

  NAMELIST / KINPAR / kindt

  ! reads in fault_db binary file and initializes fault arrays
  call initialize_fault(bc,IIN_BIN)

  ! sets up initial fault state
  if (bc%nspec > 0) then
    allocate(bc%T(3,bc%nglob),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1994')
    allocate(bc%D(3,bc%nglob),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1995')
    allocate(bc%V(3,bc%nglob),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1996')

    bc%T(:,:) = 0.0_CUSTOM_REAL
    bc%D(:,:) = 0.0_CUSTOM_REAL
    bc%V(:,:) = 0.0_CUSTOM_REAL

    ! time interval between two loaded slip rates
    read(IIN_PAR,nml=KINPAR)

    bc%kin_dt = kindt
    bc%kin_it = 0

    ! Always have in memory the slip-rate model at two times, t1 and t2,
    ! spatially interpolated in the spectral element grid
    allocate(bc%v_kin_t1(2,bc%nglob),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1997')
    allocate(bc%v_kin_t2(2,bc%nglob),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1998')
    bc%v_kin_t1(:,:) = 0.0_CUSTOM_REAL
    bc%v_kin_t2(:,:) = 0.0_CUSTOM_REAL
  else
    ! dummy allocations (for subroutine arguments)
    allocate(bc%T(3,1), &
             bc%D(3,1), &
             bc%V(3,1))
  endif

  ! output dataT structure
  recordlength = 7  ! number of quantities to store for each record
  call init_dataT(bc%dataT,bc%coord,bc%nglob,NT,dt_real,recordlength,iflt)

  call init_dataXZ(bc%dataXZ,bc)

  end subroutine init_one_fault


!=====================================================================
! adds boundary term Bt to Force array for each fault.
!
  subroutine BC_KINFLT_set_all(F,Vel,Dis)

  implicit none
  real(kind=CUSTOM_REAL), dimension(:,:), intent(in) :: Vel,Dis
  real(kind=CUSTOM_REAL), dimension(:,:), intent(inout) :: F

  ! local parameters
  integer :: iflt

  ! checks if anything to do
  if (Nfaults == 0) return

  ! loops over faults
  do iflt = 1,Nfaults
    ! note: this routine should be called by all processes, regardless if they contain no fault elements,
    !       for managing MPI calls and file outputs
    call BC_KINFLT_set_single(faults(iflt),F,Vel,Dis,iflt)
  enddo

  end subroutine BC_KINFLT_set_all

!---------------------------------------------------------------------
!
!NOTE: On non-split nodes at fault edges, dD=dV=dA=0 but bc%T is corrupted.
!      That does not affect computations: the net contribution of B*T is =0.
!      However, the output T in these nodes should be ignored.
!      It is =0 if the user sets bc%V=0 there in the input slip rates.
!
  subroutine BC_KINFLT_set_single(bc,MxA,V,D,iflt)

  use specfem_par, only: it,NSTEP,myrank

  implicit none
  real(kind=CUSTOM_REAL), intent(inout) :: MxA(:,:)
  type(bc_dynandkinflt_type), intent(inout) :: bc
  real(kind=CUSTOM_REAL), intent(in) :: V(:,:),D(:,:)
  integer,intent(in) :: iflt
  integer :: it_kin,itime
  real(kind=CUSTOM_REAL), dimension(3,bc%nglob) :: T
  real(kind=CUSTOM_REAL), dimension(3,bc%nglob) :: dD,dV,dA,dV_free
  real(kind=CUSTOM_REAL) :: t1,t2
  real(kind=CUSTOM_REAL) :: half_dt,timeval

  if (bc%nspec > 0) then !Surendra : for parallel faults

    half_dt = 0.5e0_CUSTOM_REAL*bc%dt

    ! get predicted values
    dD = get_jump(bc,D) ! dD_predictor
    dV = get_jump(bc,V) ! dV_predictor
    dA = get_weighted_jump(bc,MxA) ! dA_free

    ! rotate to fault frame (tangent,normal)
    ! component 3 is normal to the fault
    dD = rotate(bc,dD,1)
    dV = rotate(bc,dV,1)
    dA = rotate(bc,dA,1)

    ! Time marching
    timeval = it*bc%dt
    ! Slip_rate step "it_kin"
    it_kin = bc%kin_it*nint(bc%kin_dt/bc%dt)
    ! (nint : Fortran round (nearest whole number) ,
    !  if nint(a)=0.5 then "a" get upper bound )

    ! Loading the next slip_rate one ahead it.
    ! This is done in case bc%kin_dt
    ! if (it_kin == it) it_kin=it_kin+1 !


    !NOTE : it and it_kin is being used due to integers are exact numbers.
    if (it > it_kin) then

      print *, 'it :', it
      print *, 'it_kin :', it_kin

      bc%kin_it = bc%kin_it +1
      bc%v_kin_t1 = bc%v_kin_t2
      print *, 'loading v_kin_t2'
      !Temporal : just for snapshots file names kin_dt=0.1 , dt=0.0001
      !snapshot(100=itime).. : itime=kin_it*(kin_dt/dt)
      itime = bc%kin_it*nint(bc%kin_dt/bc%dt)
      call load_vslip_snapshots(bc%dataXZ,itime,iflt,myrank)
!     loading slip rates
      bc%v_kin_t2(1,:)=bc%dataXZ%v1
      bc%v_kin_t2(2,:)=bc%dataXZ%v2

      !linear interpolation in time between t1 and t2
      !REMARK , bc%kin_dt is the delta "t" between two snapshots.

    endif

    t1 = (bc%kin_it-1) * bc%kin_dt
    t2 = bc%kin_it * bc%kin_dt

    ! Kinematic velocity_rate
    ! bc%V : Imposed a priori and read from slip rate snapshots (from time reversal)
    !        Linear interpolation between consecutive kinematic time steps.
    !        V will be given at each time step.
    bc%V(1,:) = ( (t2 - timeval)*bc%v_kin_t1(1,:) + (timeval - t1)*bc%v_kin_t2(1,:) )/ bc%kin_dt
    bc%V(2,:) = ( (t2 - timeval)*bc%v_kin_t1(2,:) + (timeval - t1)*bc%v_kin_t2(2,:) )/ bc%kin_dt

    !dV_free = dV_predictor + (dt/2)*dA_free
    dV_free(1,:) = dV(1,:) + half_dt*dA(1,:)
    dV_free(2,:) = dV(2,:) + half_dt*dA(2,:)
    dV_free(3,:) = dV(3,:) + half_dt*dA(3,:)

    ! T = Z*( dV_free - V) , V known apriori as input.
    ! CONVENTION : T(ibulk1) = T = -T(ibulk2)
    T(1,:) = bc%Z * ( dV_free(1,:) - bc%V(1,:) )
    T(2,:) = bc%Z * ( dV_free(2,:) - bc%V(2,:) )
    T(3,:) = bc%Z * ( dV_free(3,:) )

    ! Save tractions
    bc%T = T

    ! Update slip in fault frame
    bc%D = dD

    ! Rotate tractions back to (x,y,z) frame
    T = rotate(bc,T,-1)

    ! Add boundary term B*T to M*a
    call add_BT(bc,MxA,T)

    !-- intermediate storage of outputs --
    call store_dataT(bc%dataT,bc%D,bc%V,bc%T,it)

    !-- OUTPUTS --
    ! write dataT every NTOUT time steps or at the end of simulation
    if (mod(it,NTOUT) == 0 .or. it == NSTEP) call SCEC_write_dataT(bc%dataT)
    ! write dataXZ every NSNAP time steps
    if (mod(it,NSNAP) == 0) call write_dataXZ(bc%dataXZ,it,iflt)

  endif

  end subroutine BC_KINFLT_set_single

!===============================================================

  subroutine init_dataXZ(dataXZ,bc)

  implicit none
  type(dataXZ_type), intent(inout) :: dataXZ
  type(bc_dynandkinflt_type) :: bc

  integer :: ier

  if (bc%nglob > 0) then
    dataXZ%d1 => bc%d(1,:)
    dataXZ%d2 => bc%d(2,:)
    dataXZ%v1 => bc%v(1,:)
    dataXZ%v2 => bc%v(2,:)
    dataXZ%t1 => bc%t(1,:)
    dataXZ%t2 => bc%t(2,:)
    dataXZ%t3 => bc%t(3,:)

    allocate(dataXZ%xcoord(bc%nglob),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1999')
    allocate(dataXZ%ycoord(bc%nglob),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2000')
    allocate(dataXZ%zcoord(bc%nglob),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 2001')
  else
    ! dummy allocations (for subroutine arguments)
    allocate(dataXZ%xcoord(1), &
             dataXZ%ycoord(1), &
             dataXZ%zcoord(1))
  endif

  end subroutine init_dataXZ

!===============================================================

  subroutine write_dataXZ(dataXZ,itime,iflt)

  use specfem_par, only: myrank,OUTPUT_FILES

  implicit none
  type(dataXZ_type), intent(in) :: dataXZ
  integer, intent(in) :: itime,iflt

  ! local parameters
  integer :: ier
  character(len=MAX_STRING_LEN) :: filename
  integer, parameter :: IOUT_SN = 171

  write(filename,"(a,'/Proc',I0,'Snapshot',I0,'_F',I0,'.bin')") trim(OUTPUT_FILES),myrank,itime,iflt

  open(unit=IOUT_SN, file=trim(filename), status='unknown', form='unformatted',action='write',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening Snapshot file: ',trim(filename)
    stop 'Error opening Snapshot file'
  endif

  write(IOUT_SN) dataXZ%xcoord
  write(IOUT_SN) dataXZ%ycoord
  write(IOUT_SN) dataXZ%zcoord
  write(IOUT_SN) dataXZ%d1
  write(IOUT_SN) dataXZ%d2
  write(IOUT_SN) dataXZ%v1
  write(IOUT_SN) dataXZ%v2
  write(IOUT_SN) dataXZ%t1
  write(IOUT_SN) dataXZ%t2
  write(IOUT_SN) dataXZ%t3
  close(IOUT_SN)

  end subroutine write_dataXZ


!---------------------------------------------------------------
!LOAD_VSLIP_SNAPSHOTS(v,dataXZ,itime,coord,npoin,nglob,iflt)
!Loading slip velocity from snapshots.
!   INPUT  itime : iteration time
!          coord : Receivers coordinates
!          npoin : number of Receivers.
!          nglob : number of GLL points along the fault.
!          dataXZ : slip rate .
!          iflt : number of faults.

!   OUTPUT v : slip rate on receivers.
  subroutine load_vslip_snapshots(dataXZ,itime,iflt,myrank)

  use constants, only: IIN_FLT

  implicit none

  integer, intent(in) :: itime,iflt
  type(dataXZ_type), intent(inout) :: dataXZ
  character(len=70) :: filename
  integer :: ier,myrank

  write(filename,"('../INPUT_FILES/Proc',I0,'Snapshot',I0,'_F',I0,'.bin')") myrank,itime,iflt
  print *, trim(filename)

  !open(unit=IIN_FLT, file= trim(filename), status='old', form='formatted', action='read',iostat=ier)
  open(unit=IIN_FLT, file= trim(filename), status='old', form='unformatted', action='read',iostat=ier)

!  COMPILERS WRITE BINARY OUTPUTS IN DIFFERENT FORMATS!
  !open(unit=IIN_FLT, file= trim(filename), status='old', form='unformatted', action='read',iostat=ier)
  !if ( ier /= 0 ) stop 'Snapshots have been found'

  if (ier == 0) then
    !read(IIN_FLT,"(5F24.15)") dataXZ%xcoord,dataXZ%ycoord,dataXZ%zcoord,dataXZ%v1,dataXZ%v2
    write(IMAIN,*) 'Load vslip file for kinematic rupture simulation!'
    !write(IMAIN,*)   max(abs(dataXZ

    read(IIN_FLT)   dataXZ%xcoord
    read(IIN_FLT)   dataXZ%ycoord
    read(IIN_FLT)   dataXZ%zcoord
    read(IIN_FLT)   dataXZ%v1
    read(IIN_FLT)   dataXZ%v2
    close(IIN_FLT)
  else
    ! if file not found, set slip velocity to zero
    dataXZ%v1 = 0e0_CUSTOM_REAL
    dataXZ%v2 = 0e0_CUSTOM_REAL
  endif

  end subroutine load_vslip_snapshots


end module fault_solver_kinematic
