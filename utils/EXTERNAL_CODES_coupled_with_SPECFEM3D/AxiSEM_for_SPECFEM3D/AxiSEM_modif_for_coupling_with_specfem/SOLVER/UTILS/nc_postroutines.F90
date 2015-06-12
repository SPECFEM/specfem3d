!
!    Copyright 2013, Tarje Nissen-Meyer, Alexandre Fournier, Martin van Driel
!                    Simon Stähler, Kasra Hosseini, Stefanie Hempel
!
!    This file is part of AxiSEM.
!    It is distributed from the webpage <http://www.axisem.info>
!
!    AxiSEM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    AxiSEM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with AxiSEM.  If not, see <http://www.gnu.org/licenses/>.
!

module nc_postroutines

#ifdef unc
  use netcdf
#endif

  implicit none
  integer, parameter         :: sp = selected_real_kind(6, 37)
  type ncparamtype
      integer, dimension(4)  :: id
      integer, dimension(4)  :: seis_grpid
      integer, dimension(4)  :: disp_varid, recnam_varid
      integer, dimension(4)  :: phi_varid, theta_varid, thetar_varid
  end type

contains

!-----------------------------------------------------------------------------------------
subroutine nc_open(ncid, nsim, simdir)
  implicit none
  type(ncparamtype), intent(inout)            :: ncid
  integer, intent(in)                         :: nsim
  character(len=100), allocatable, intent(in) :: simdir(:)
  character(len=256)                          :: fnam
  integer                                     :: isim
    
#ifdef unc
  do isim = 1, nsim
      fnam = trim(simdir(isim))//'/Data/axisem_output.nc4'
      write(6,*) 'Opening forward file ', trim(fnam)
      call check( nf90_open ( path=trim(fnam), mode=NF90_NOWRITE, ncid=ncid%id(isim)))
      call check( nf90_inq_grp_ncid(ncid%id(isim), "Seismograms", ncid%seis_grpid(isim)) )
      call check( nf90_inq_varid(ncid%seis_grpid(isim), 'displacement', &
                                 ncid%disp_varid(isim)) )
      call check( nf90_inq_varid(ncid%seis_grpid(isim), 'receiver_name', &
                                 ncid%recnam_varid(isim)))
      call check( nf90_inq_varid(ncid%seis_grpid(isim), 'phi', ncid%phi_varid(isim)))
      call check( nf90_inq_varid(ncid%seis_grpid(isim), 'theta', &
                                 ncid%theta_varid(isim)))
      call check( nf90_inq_varid(ncid%seis_grpid(isim), 'theta_requested', &
                                 ncid%thetar_varid(isim)))
  end do
#endif


end subroutine nc_open
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine nc_close(ncid)
  implicit none
  type(ncparamtype), intent(inout)            :: ncid

#ifdef unc
  integer                                     :: state, isim
  
  do isim = 1,4
      state = nf90_close(ncid%id(isim))
  end do
#endif

end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine nc_read_recnames(ncid, nrec, recname, theta, thetar, phi)
  type(ncparamtype), intent(in)               :: ncid
  integer, intent(in)                         :: nrec
  character(len=40), intent(out)              :: recname(nrec)
  real(kind=sp), dimension(nrec), intent(out) :: theta, phi, thetar
  integer                                     :: irec

#ifdef unc
  do irec = 1, nrec
      call check( nf90_get_var(ncid%seis_grpid(1), ncid%recnam_varid(1), &
                               start = (/irec, 1/), &
                               count = (/1, 40/), & 
                               values = recname(irec)) )
  end do
  call check( nf90_get_var(ncid%seis_grpid(1), ncid%phi_varid(1), &
                           values = phi) )
  call check( nf90_get_var(ncid%seis_grpid(1), ncid%theta_varid(1), &
                           values  = theta) )
  call check( nf90_get_var(ncid%seis_grpid(1), ncid%thetar_varid(1), &
                           values  = thetar) )
  thetar = 90. - thetar
  theta = 90. - theta

#else
  ! To stop Ifort complaining that intent(out)-variables are not set
  phi     = 0.0
  theta   = 0.0
  thetar  = 0.0
  recname = ''

#endif
end subroutine nc_read_recnames
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine nc_read_seismograms(ncid, nt, nrec, nsim, seis_snglcomp)
  implicit none
  type(ncparamtype), intent(in)  :: ncid
  integer, intent(in)            :: nt, nrec, nsim
  real, intent(out)              :: seis_snglcomp(nt,3,nrec,nsim)
  integer                        :: isim

#ifdef unc
  do isim = 1, nsim
      call check( nf90_get_var( ncid%seis_grpid(isim), ncid%disp_varid(isim), &
                                start=(/1, 1, 1/), &
                                count = (/nt, 3, nrec/), &
                                values = seis_snglcomp(:,:,:,isim)) )
  end do
#else
  ! To stop Ifort complaining that intent(out)-variables are not set
  seis_snglcomp = 0.0

#endif
end subroutine nc_read_seismograms
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine check(status)
! Translates netcdf error codes into error messages
  implicit none
  integer, intent ( in) :: status

#ifdef unc
  if(status /= nf90_noerr) then 
    print *, trim(nf90_strerror(status))
    stop 1
  end if
#endif

end subroutine check  
!-----------------------------------------------------------------------------------------

end module nc_postroutines
