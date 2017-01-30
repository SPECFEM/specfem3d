
  program check_the_CPML_flag

! Dimitri Komatitsch, CNRS Marseille, France, January 2017

  implicit none

! flags for the seven CPML regions
  integer, parameter :: CPML_X_ONLY = 1
  integer, parameter :: CPML_Y_ONLY = 2
  integer, parameter :: CPML_Z_ONLY = 3
  integer, parameter :: CPML_XY_ONLY = 4
  integer, parameter :: CPML_XZ_ONLY = 5
  integer, parameter :: CPML_YZ_ONLY = 6
  integer, parameter :: CPML_XYZ = 7

  integer :: NCPML,idummy,i

  integer, dimension(:), allocatable :: imaterialCPML

open(unit=23,file='absorbing_cpml_file',status='old',action='read')

read(23,*) NCPML

print *, 'reading ',NCPML,' CPML elements'

allocate(imaterialCPML(NCPML))

do i = 1,NCPML
  read(23,*) idummy,imaterialCPML(i)
enddo

close(23)

print *,'minval and maxval of CPML flags read (should be 1 and 7 in most cases) = ',minval(imaterialCPML),maxval(imaterialCPML)

print *,'number of purely X face PML elements = ',count(imaterialCPML == CPML_X_ONLY)
print *,'number of purely Y face PML elements = ',count(imaterialCPML == CPML_Y_ONLY)
print *,'number of purely Z face PML elements = ',count(imaterialCPML == CPML_Z_ONLY)

print *,'number of XY edge PML elements = ',count(imaterialCPML == CPML_XY_ONLY)
print *,'number of XZ edge PML elements = ',count(imaterialCPML == CPML_XZ_ONLY)
print *,'number of YZ edge PML elements = ',count(imaterialCPML == CPML_YZ_ONLY)

print *,'number of XYZ corner PML elements = ',count(imaterialCPML == CPML_XYZ)

end program check_the_CPML_flag

