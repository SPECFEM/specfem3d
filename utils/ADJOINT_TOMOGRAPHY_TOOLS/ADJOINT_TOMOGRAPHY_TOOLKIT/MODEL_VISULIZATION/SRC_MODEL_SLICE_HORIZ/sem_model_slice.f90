program sem_model_slice
implicit none

include 'mpif.h'
include '../../SHARE_FILES/HEADER_FILES/constants.h'
include '../../SHARE_FILES/HEADER_FILES/values_from_mesher.h'
include '../../SHARE_FILES/HEADER_FILES/precision.h'

integer,parameter:: NUM_NODES=99  ! for recent mesher 100 processors
integer,parameter:: iregion=1    ! for region one
integer,parameter:: NMAXPTS=10000000
real(kind=CUSTOM_REAL):: R_CUT_RANGE=0.00785d0  ! dr=0.00785 ~ 50 km depth

integer::iproc,ipt,npts
character(len=150):: xyz_infile,topo_dir,model_dir,filename,gmt_outfile
character(len=256):: prname_topo, prname_file
character(len=256):: topo_file, data_file
real(kind=CUSTOM_REAL),dimension(NMAXPTS):: x,y,z
real(kind=CUSTOM_REAL),dimension(NMAXPTS):: xfound,yfound,zfound,vfound
real(kind=CUSTOM_REAL),dimension(NMAXPTS)::distmin,dist, distall, v, vall
integer,dimension(NMAXPTS):: ispec_found,ix_found,iy_found,iz_found
real(kind=CUSTOM_REAL),dimension(NGLOB_CRUST_MANTLE):: xstore, ystore, zstore
integer,dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE):: ibool
real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE):: vstore
integer:: ier,sizeprocs,myrank,ios
integer:: i,j,k,ispec,iglob
real(kind=CUSTOM_REAL):: r,theta,phi,lat,lon,dep, xmesh,ymesh,zmesh
real(kind=CUSTOM_REAL),dimension(2,NMAXPTS):: in,out
real(kind=CUSTOM_REAL):: R_CUT, r1

call MPI_INIT(ier)
call MPI_COMM_SIZE(MPI_COMM_WORLD,sizeprocs,ier)
call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ier)

! read input file
call getarg(1,xyz_infile)
call getarg(2,topo_dir)
call getarg(3,model_dir)
call getarg(4,filename)
call getarg(5,gmt_outfile)

! read interpolate points
if ( myrank == 0) then
        write(*,*) "INPUT FILE:", trim(xyz_infile)
        write(*,*) "TOPOLOGY FILE:", trim(topo_dir)
        write(*,*) "MODEL DIR:", trim(model_dir)
        write(*,*) "VALUE NAME:", trim(filename)
        write(*,*) "OUTPUT:", trim(gmt_outfile)
endif


open(1001,file=trim(xyz_infile),status='old',iostat=ios)
i=0
do while (1 == 1)
        i=i+1
        read(1001,*,iostat=ios) xmesh,ymesh,zmesh
        R_CUT=sqrt(xmesh**2+ymesh**2+zmesh**2)
        if (ios /= 0) exit
        x(i)=xmesh
        y(i)=ymesh
        z(i)=zmesh
enddo
close(1001)
npts=i-1

!endif
!call MPI_BCAST(npts,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
!call MPI_BCAST(x,NMAXPTS,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ier)
!call MPI_BCAST(y,NMAXPTS,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ier)
!call MPI_BCAST(z,NMAXPTS,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ier)


if ( myrank == 0 ) then
        write(*,*) 'Total number of points = ', npts
        if (npts > NMAXPTS .or. npts < 0) call exit_MPI(myrank,'Npts error...')
endif


write(prname_topo,'(a,i6.6,a,i1,a)') trim(topo_dir)//'/proc',myrank,'_reg',iregion,'_'
write(prname_file,'(a,i6.6,a,i1,a)') trim(model_dir)//'/proc',myrank,'_reg',iregion,'_'

! read value file
data_file = trim(prname_file) // trim(filename) // '.bin'
open(unit = 27,file = trim(data_file),status='old',action='read', iostat = ios,form ='unformatted')
if (ios /= 0) call exit_MPI(myrank,'Error reading value file')
read(27) vstore(:,:,:,1:NSPEC_CRUST_MANTLE)
close(27)

if (myrank == 0) write(*,*) 'DONE READING',trim(data_file)

! read topology
topo_file = trim(prname_topo) // 'solver_data_2' // '.bin'
open(unit = 28,file = trim(topo_file),status='old',action='read', iostat = ios, form='unformatted')
if (ios /= 0) call exit_MPI(myrank,'Error reading topology file')
xstore(:) = 0.0
ystore(:) = 0.0
zstore(:) = 0.0
ibool(:,:,:,:) = -1
read(28) xstore(1:NGLOB_CRUST_MANTLE)
read(28) ystore(1:NGLOB_CRUST_MANTLE)
read(28) zstore(1:NGLOB_CRUST_MANTLE)
read(28) ibool(:,:,:,1:NSPEC_CRUST_MANTLE)
close(28)

if (myrank == 0) write(*,*) 'DONE READING',trim(topo_file)

distmin(1:npts)=HUGEVAL

do ispec=1,NSPEC_CRUST_MANTLE

   if (myrank == 0) write(*,*) 'ispec=',ispec

   do k = 1,NGLLZ
      do j = 1,NGLLY
         do i = 1,NGLLX
            iglob=ibool(i,j,k,ispec)


            r1=sqrt((xstore(iglob))**2+(ystore(iglob))**2+(zstore(iglob))**2)

            if ( abs(r1-R_CUT) < R_CUT_RANGE ) then
                dist(1:npts) = dsqrt((x(1:npts)-dble(xstore(iglob)))**2 &
                                +(y(1:npts)-dble(ystore(iglob)))**2 &
                                +(z(1:npts)-dble(zstore(iglob)))**2)

                do ipt=1,npts
                        if (dist(ipt) < distmin(ipt)) then
                                distmin(ipt)=dist(ipt)
                                ispec_found(ipt)=ispec
                                ix_found(ipt)=i
                                iy_found(ipt)=j
                                iz_found(ipt)=k
                                vfound(ipt)=vstore(i,j,k,ispec)
                        endif
                enddo ! ipt loop
            endif ! radius within in a range

         enddo ! i loop
       enddo  ! j loop
    enddo ! k loop
enddo  ! ispec loop

call MPI_BARRIER(MPI_COMM_WORLD,ier)
if (myrank == 0) print *,'Done looping over global points'

do i = 1,npts
        in(1,i)=distmin(i)
        in(2,i)=myrank
enddo
call MPI_REDUCE(in,out,npts,CUSTOM_MPI_2REAL,MPI_MINLOC,0,MPI_COMM_WORLD,ier)

call MPI_BCAST(out,2*npts,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ier)

v(1:npts)=0
dist(1:npts)=0.

do i=1,npts
        if (myrank == nint(out(2,i))) then
                v(i)=vfound(i)
                dist(i)=distmin(i)
        endif
enddo

call MPI_REDUCE(v,vall,npts,CUSTOM_MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD,ier)
call MPI_REDUCE(dist,distall,npts,CUSTOM_MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD,ier)

if (myrank == 0) then
        open(1002, file=gmt_outfile,status='unknown')
        do i = 1,npts
               xmesh=x(i); ymesh=y(i); zmesh=z(i)
               call xyz_2_rthetaphi(xmesh,ymesh,zmesh,r,theta,phi)
               lat=90.0 - theta*180.0/PI
               lon=phi*180.0/PI
               dep=(1.0-r)*R_EARTH_KM
!               write(1002,*) lon,lat,dep,vall(i),distall(i)
               write(1002,*) lon,lat,r,vall(i),distall(i)
        enddo
        close(1002)
endif

202 FORMAT (5(F12.9,2X))
call MPI_FINALIZE(ier)
end program sem_model_slice
