
 program rescalenonlinear

! determine GMT colorbar mapping

 implicit none

 double precision, parameter :: POWERVAL = 0.25d0

 double precision a

 a = 0.1d0
 write(*,*) a,' gives ',a**POWERVAL

 a = 0.2d0
 write(*,*) a,' gives ',a**POWERVAL

 a = 0.3d0
 write(*,*) a,' gives ',a**POWERVAL

 a = 0.4d0
 write(*,*) a,' gives ',a**POWERVAL

 a = 0.5d0
 write(*,*) a,' gives ',a**POWERVAL

 a = 0.6d0
 write(*,*) a,' gives ',a**POWERVAL

 a = 0.7d0
 write(*,*) a,' gives ',a**POWERVAL

 a = 0.8d0
 write(*,*) a,' gives ',a**POWERVAL

 a = 0.9d0
 write(*,*) a,' gives ',a**POWERVAL

 end

