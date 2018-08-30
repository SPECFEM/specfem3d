!
!    Copyright 2013, Tarje Nissen-Meyer, Alexandre Fournier, Martin van Driel
!                    Simon Stahler, Kasra Hosseini, Stefanie Hempel
!
!    This file is part of AxiSEM.
!    It is distributed from the webpage < http://www.axisem.info>
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
!    along with AxiSEM.  If not, see < http://www.gnu.org/licenses/>.
!

!=========================================================================================
!> Contains some routines to ease interaction with netcdf files
module nc_helpers

#ifdef enable_netcdf
    use netcdf
#endif
    use data_io, only: verbose, ncid_out
    use data_proc, only: mynum
    use global_parameters

    implicit none
    private

    public :: check
    public :: putvar_real1d
    public :: putvar_real2d
    public :: putvar_real3d
    public :: getgrpid
    public :: getvarid
    public :: nc_write_att_int
    public :: nc_write_att_real
    public :: nc_write_att_dble
    public :: nc_write_att_char

contains

!-----------------------------------------------------------------------------------------
!> Translates NetCDF error code into readable message
subroutine check(status)
    integer, intent ( in) :: status ! < Error code
#ifdef enable_netcdf
    if (status /= nf90_noerr) then
        print *, trim(nf90_strerror(status))
        call abort()
    endif
#endif
end subroutine check
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine getvarid(ncid, name, varid)
    integer, intent(in)          :: ncid
    character(len=*), intent(in) :: name
    integer, intent(out)         :: varid
#ifdef enable_netcdf
    integer                      :: status

    status = nf90_inq_varid( ncid  = ncid, &
                             name  = name, &
                             varid = varid )
    if (status /= NF90_NOERR) then
        write(*,100) mynum, trim(name), ncid
        stop
    else if (verbose > 2) then
        write(*,101) trim(name), ncid, varid
        call flush(6)
    endif
100 format('ERROR: CPU ', I4, ' could not find variable: ''', A, ''' in NCID', I7)
101 format('    Variable ''', A, ''' found in NCID', I7, ', has ID:', I7)
#else
    varid = 0
#endif
end subroutine getvarid
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine getgrpid(ncid, name, grpid)
    integer, intent(in)          :: ncid
    character(len=*), intent(in) :: name
    integer, intent(out)         :: grpid
#ifdef enable_netcdf
    integer                      :: status

    status = nf90_inq_ncid( ncid     = ncid, &
                            name     = name, &
                            grp_ncid = grpid )
    if (status /= NF90_NOERR) then
        write(*,100) mynum, trim(name), ncid
        stop
    else if (verbose > 2) then
        write(*,101) trim(name), ncid, grpid
        call flush(6)
    endif
100 format('ERROR: CPU ', I4, ' could not find group: ''', A, ''' in NCID', I7)
101 format('    Group ''', A, ''' found in NCID', I7, ', has ID:', I7)
#else
    grpid = 0
#endif
end subroutine getgrpid
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine putvar_real1d(ncid, varid, values, start, count)
! < Help interpret the inane NetCDF error messages
   integer, intent(in)          :: ncid, varid, start, count
   real, intent(in)             :: values(:)

#ifdef enable_netcdf
   integer                      :: xtype, ndims, status, dimsize
   integer                      :: dimid(10)
   character(len=nf90_max_name) :: varname, dimname


   status = nf90_inquire_variable(ncid  = ncid, &
                                  varid = varid, &
                                  name  = varname )

   if (status /= NF90_NOERR) then
       write(*,99) mynum, varid, ncid
       print *, trim(nf90_strerror(status))
       stop
   endif

   if (size(values) /= count) then
       write(*,100) mynum, trim(varname), varid, ncid, size(values), count
       stop
   endif

   status = nf90_put_var(ncid   = ncid, &
                         varid  = varid, &
                         values = values, &
                         start  = [start], &
                         count  = [count] )


   if (status /= NF90_NOERR) then
       status = nf90_inquire_variable(ncid  =  ncid, &
                                      varid = varid, &
                                      name  = varname, &
                                      ndims = ndims)
       if (ndims /= 1) then
           write(*,101) mynum, trim(varname), varid, ncid, ndims
           print *, trim(nf90_strerror(status))
           stop
       endif
       status = nf90_inquire_variable(ncid   = ncid, &
                                      varid  = varid, &
                                      name   = varname, &
                                      xtype  = xtype, &
                                      ndims  = ndims, &
                                      dimids = dimid  )

       status = nf90_inquire_dimension(ncid  = ncid, &
                                       dimid = dimid(1), &
                                       name  = dimname, &
                                       len   = dimsize )
       if (start + count - 1 > dimsize) then
           write(*,102) mynum, trim(varname), varid, ncid, start, count, dimsize, trim(dimname)
           print *, trim(nf90_strerror(status))
           stop
       endif

       write(*,103) mynum, trim(varname), varid, ncid, start, count, dimsize, trim(dimname)
       print *, trim(nf90_strerror(status))
       stop

   else if (verbose > 2) then
       write(*,200) mynum, real(count) * 4. / 1048576., ncid, varid
       call flush(6)
   endif

99  format('ERROR: CPU ', I4, ' could not find 1D variable: ',I7,' in NCID', I7)
100 format('ERROR: CPU ', I4, ' could not write 1D variable: ''', A, '''(',I7,') in NCID', I7, / &
           '       was given ', I10, ' values, but ''count'' is ', I10)
101 format('ERROR: CPU ', I4, ' could not write 1D variable: ''', A, '''(',I7,') in NCID', I7, / &
           '       Variable has ', I2,' dimensions instead of one')
102 format('ERROR: CPU ', I4, ' could not write 1D variable: ''', A, '''(',I7,') in NCID', I7, / &
           '       start (', I10, ') + count(', I10, ') is larger than size (', I10,')',    / &
           '       of dimension ', A)
103 format('ERROR: CPU ', I4, ' could not write 1D variable: ''', A, '''(',I7,') in NCID', I7, / &
           '       start:   ', I10, / &
           '       count:   ', I10, / &
           '       dimsize: ', I10, / &
           '       dimname: ', A)
200 format('    Proc ', I4, ': Wrote', F10.3, ' MB into 1D variable in NCID', I7, ', with ID:', I7)
#endif
end subroutine putvar_real1d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine putvar_real2d(ncid, varid, values, start, count)
! < Help interpret the inane NetCDF error messages
   integer, intent(in)          :: ncid, varid
   integer, intent(in)          :: start(2), count(2)
   real, intent(in)             :: values(:,:)

#ifdef enable_netcdf
   integer                      :: xtype, ndims, status, dimsize, idim
   integer                      :: dimid(10)
   character(len=nf90_max_name) :: varname, dimname


   status = nf90_inquire_variable(ncid  = ncid, &
                                  varid = varid, &
                                  name  = varname )

   if (status /= NF90_NOERR) then
       write(*,99) mynum, varid, ncid
       print *, trim(nf90_strerror(status))
       stop
   endif
   ! Check if variable size is consistent with values of 'count'
   do idim = 1, 2
       if (size(values,idim) /= count(idim)) then
           write(*,100) mynum, trim(varname), varid, ncid, idim, size(values, idim), count(idim)
           stop
       endif
   enddo

   ! Write data to file
   status = nf90_put_var(ncid   = ncid, &
                         varid  = varid, &
                         values = values, &
                         start  = start, &
                         count  = count )


   ! If an error has occurred, try to find a reason
   if (status /= NF90_NOERR) then
       status = nf90_inquire_variable(ncid  =  ncid, &
                                      varid = varid, &
                                      name  = varname, &
                                      ndims = ndims)

       ! Check whether variable in NetCDF file has more or less than three dimensions
       if (ndims /= 2) then
           write(*,101) mynum, trim(varname), varid, ncid, ndims
           print *, trim(nf90_strerror(status))
           stop
       endif

       ! Check whether dimension sizes are compatible with amount of data written
       status = nf90_inquire_variable(ncid   = ncid, &
                                      varid  = varid, &
                                      name   = varname, &
                                      xtype  = xtype, &
                                      ndims  = ndims, &
                                      dimids = dimid  )

       do idim = 1, 2
           status = nf90_inquire_dimension(ncid  = ncid, &
                                           dimid = dimid(idim), &
                                           name  = dimname, &
                                           len   = dimsize )
           if (start(idim) + count(idim) - 1 > dimsize) then
               write(*,102) mynum, trim(varname), varid, ncid, start(idim), count(idim), &
                            dimsize, trim(dimname), idim
               print *, trim(nf90_strerror(status))
               stop
           endif

           ! Otherwise just dump as much information as possible and stop
           write(*,103) mynum, trim(varname), varid, ncid, start(idim), count(idim), &
                        dimsize, trim(dimname)
           print *, trim(nf90_strerror(status))

       enddo

       stop

   else if (verbose > 2) then
       ! Everything okay
       write(*,200) mynum, real(product(count)) * 4. / 1048576., ncid, varid
       call flush(6)
   endif

99  format('ERROR: CPU ', I4, ' could not find 2D variable: ',I7,' in NCID', I7)
100 format('ERROR: CPU ', I4, ' could not write 2D variable: ''', A, '''(',I7,') in NCID', I7, / &
           '       dimension ', I1,' was given ', I10, ' values, but ''count'' is ', I10)
101 format('ERROR: CPU ', I4, ' could not write 2D variable: ''', A, '''(',I7,') in NCID', I7, / &
           '       Variable has ', I2,' dimensions instead of two')
102 format('ERROR: CPU ', I4, ' could not write 2D variable: ''', A, '''(',I7,') in NCID', I7, / &
           '       start (', I10, ') + count(', I10, ') is larger than size (', I10,')',    / &
           '       of dimension ', A, ' (', I1, ')')
103 format('ERROR: CPU ', I4, ' could not write 2D variable: ''', A, '''(',I7,') in NCID', I7, / &
           '       start:   ', I10, / &
           '       count:   ', I10, / &
           '       dimsize: ', I10, / &
           '       dimname: ', A)
200 format('    Proc ', I4, ': Wrote', F10.3, ' MB into 2D variable in NCID', I7, ', with ID:', I7)
#endif
end subroutine putvar_real2d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine putvar_real3d(ncid, varid, values, start, count)
! < Help interpret the inane NetCDF error messages
   integer, intent(in)          :: ncid, varid
   integer, intent(in)          :: start(3), count(3)
   real, intent(in)             :: values(:,:,:)

#ifdef enable_netcdf
   integer                      :: xtype, ndims, status, dimsize, idim
   integer                      :: dimid(10)
   character(len=nf90_max_name) :: varname, dimname


   status = nf90_inquire_variable(ncid  = ncid, &
                                  varid = varid, &
                                  name  = varname )

   if (status /= NF90_NOERR) then
       write(*,99) mynum, varid, ncid
       print *, trim(nf90_strerror(status))
       stop
   endif
   ! Check if variable size is consistent with values of 'count'
   do idim = 1, 3
       if (size(values,idim) /= count(idim)) then
           write(*,100) mynum, trim(varname), varid, ncid, idim, size(values, idim), count(idim)
           stop
       endif
   enddo

   ! Write data to file
   status = nf90_put_var(ncid   = ncid, &
                         varid  = varid, &
                         values = values, &
                         start  = start, &
                         count  = count )


   ! If an error has occurred, try to find a reason
   if (status /= NF90_NOERR) then
       status = nf90_inquire_variable(ncid  =  ncid, &
                                      varid = varid, &
                                      name  = varname, &
                                      ndims = ndims)

       ! Check whether variable in NetCDF file has more or less than three dimensions
       if (ndims /= 3) then
           write(*,101) mynum, trim(varname), varid, ncid, ndims
           print *, trim(nf90_strerror(status))
           stop
       endif

       ! Check whether dimension sizes are compatible with amount of data written
       status = nf90_inquire_variable(ncid   = ncid, &
                                      varid  = varid, &
                                      name   = varname, &
                                      xtype  = xtype, &
                                      ndims  = ndims, &
                                      dimids = dimid  )

       do idim = 1, 3
           status = nf90_inquire_dimension(ncid  = ncid, &
                                           dimid = dimid(idim), &
                                           name  = dimname, &
                                           len   = dimsize )
           if (start(idim) + count(idim) - 1 > dimsize) then
               write(*,102) mynum, trim(varname), varid, ncid, start(idim), count(idim), &
                            dimsize, trim(dimname), idim
               print *, trim(nf90_strerror(status))
               stop
           endif

           ! Otherwise just dump as much information as possible and stop
           write(*,103) mynum, trim(varname), varid, ncid, start(idim), count(idim), &
                        dimsize, trim(dimname)
           print *, trim(nf90_strerror(status))

       enddo

       stop

   else if (verbose > 2) then
       ! Everything okay
       write(*,200) mynum, real(product(count)) * 4. / 1048576., ncid, varid
       call flush(6)
   endif

99  format('ERROR: CPU ', I4, ' could not find 3D variable: ',I7,' in NCID', I7)
100 format('ERROR: CPU ', I4, ' could not write 3D variable: ''', A, '''(',I7,') in NCID', I7, / &
           '       dimension ', I1,' was given ', I10, ' values, but ''count'' is ', I10)
101 format('ERROR: CPU ', I4, ' could not write 3D variable: ''', A, '''(',I7,') in NCID', I7, / &
           '       Variable has ', I2,' dimensions instead of three')
102 format('ERROR: CPU ', I4, ' could not write 3D variable: ''', A, '''(',I7,') in NCID', I7, / &
           '       start (', I10, ') + count(', I10, ') is larger than size (', I10,')',    / &
           '       of dimension ', A, ' (', I1, ')')
103 format('ERROR: CPU ', I4, ' could not write 3D variable: ''', A, '''(',I7,') in NCID', I7, / &
           '       start:   ', I10, / &
           '       count:   ', I10, / &
           '       dimsize: ', I10, / &
           '       dimname: ', A)
200 format('    Proc ', I4, ': Wrote', F10.3, ' MB into 3D variable in NCID', I7, ', with ID:', I7)
#endif
end subroutine putvar_real3d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Write NetCDF attribute of type Character
subroutine nc_write_att_char(attribute_value, attribute_name)
    character(len=*), intent(in)    :: attribute_name, attribute_value

#ifdef enable_netcdf
    call check( nf90_put_att(ncid_out, NF90_GLOBAL, attribute_name, attribute_value) )
#endif
end subroutine nc_write_att_char
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Write NetCDF attribute of type Real
subroutine nc_write_att_real(attribute_value, attribute_name)
  character(len=*),  intent(in)      :: attribute_name
  real(sp), intent(in)               :: attribute_value

#ifdef enable_netcdf
  call check( nf90_put_att(ncid_out, NF90_GLOBAL, attribute_name, attribute_value) )
#endif
end subroutine nc_write_att_real
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Write NetCDF attribute of type Double
subroutine nc_write_att_dble(attribute_value, attribute_name)
  character(len=*),  intent(in)      :: attribute_name
  real(dp), intent(in)               :: attribute_value

#ifdef enable_netcdf
  call check( nf90_put_att(ncid_out, NF90_GLOBAL, attribute_name, attribute_value) )
#endif
end subroutine nc_write_att_dble
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Write NetCDF attribute of type Integer
subroutine nc_write_att_int(attribute_value, attribute_name)
  character(len=*),  intent(in)     :: attribute_name
  integer, intent(in)               :: attribute_value

#ifdef enable_netcdf
  call check( nf90_put_att(ncid_out, NF90_GLOBAL, attribute_name, attribute_value) )
#endif
end subroutine nc_write_att_int
!-----------------------------------------------------------------------------------------

end module nc_helpers
!=========================================================================================
