subroutine pinputTra(outputDir,psvmodel,modelname,stationsinf,tlen,imin,imax,r0min, &
  r0max,r0delta,r0lat,r0lon,itranslat,mt,dt,f0,f1,myrank)
  implicit none
  character(120) tmpfile
  character(120) :: dummy,outputDir,psvmodel,modelname,stationsinf
  real(kind(0d0)) :: tlen,r0min,r0max,r0delta,tb,te,f0,f1
  real(kind(0d0)) :: r0lat,r0lon,dt
  real(kind(0d0)) :: mt(6)
  integer :: i0,i1,imin,imax,itranslat,myrank


  write(tmpfile,'(a10,i5.5)') 'TmpWrkFile',myrank
  open(unit=1, file=tmpfile,status='unknown')
  open(10,file='inputIASP.infTra')
100 continue
  read(10,110) dummy
110 format(a120)
  if (dummy(1:1) == '#') goto 100
  if (dummy(1:3) == 'end') goto 120
  write(1,110) dummy
  goto 100
120 continue
  close(1)
  close(10)

  open(unit=1,file=tmpfile,status='unknown')
  read(1,110) outputDir
  read(1,110) psvmodel
  read(1,110) modelname
  read(1,110) stationsinf
  outputDir=trim(outputDir)
  psvmodel=trim(psvmodel)
  modelname=trim(modelname)
  stationsinf=trim(stationsinf)
  read(1,*) tlen
  read(1,*) r0min,r0lat,r0lon
  r0min = 6371.d0 -r0min ! because in this version we write the source DEPTH
  r0max=r0min
  r0delta=20.d0
  read(1,*) imin,imax
  read(1,*) itranslat
  ! MT reading for TraFFT
  read(1,*) mt(1),mt(2),mt(3),mt(4),mt(5),mt(6)
  read(1,110) dummy
  read(1,*) i0,i1
  read(1,*) dt,tb,te
  read(1,*) f0,f1
  close(1)
  dt = 1./dt

  close(1)
end subroutine pinputTra
