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

! read values from parameter file, ignoring white lines and comments

  subroutine read_value_integer_mesh(iunit,ignore_junk,value_to_read, name, ier)

  use constants, only: MAX_STRING_LEN

  implicit none

  logical :: ignore_junk
  integer :: iunit
  integer :: value_to_read
  integer :: ier
  character(len=*) :: name
  character(len=MAX_STRING_LEN) :: string_read

  call unused_string(name)
  ier = 0

  call read_next_line(iunit,ignore_junk,string_read,ier)
  if (ier /= 0) return

  read(string_read,*,iostat=ier) value_to_read

  end subroutine read_value_integer_mesh

!--------------------

  subroutine read_value_dble_precision_mesh(iunit,ignore_junk,value_to_read, name, ier)

  use constants, only: MAX_STRING_LEN

  implicit none

  logical :: ignore_junk
  integer :: iunit
  double precision :: value_to_read
  integer :: ier
  character(len=*) :: name
  character(len=MAX_STRING_LEN) :: string_read

  call unused_string(name)
  ier = 0

  call read_next_line(iunit,ignore_junk,string_read,ier)
  if (ier /= 0) return

  read(string_read,*,iostat=ier) value_to_read

  end subroutine read_value_dble_precision_mesh

!--------------------

  subroutine read_value_logical_mesh(iunit,ignore_junk,value_to_read, name, ier)

  use constants, only: MAX_STRING_LEN

  implicit none

  logical :: ignore_junk
  logical :: value_to_read
  integer :: iunit
  integer :: ier
  character(len=*) :: name
  character(len=MAX_STRING_LEN) :: string_read

  call unused_string(name)
  ier = 0

  call read_next_line(iunit,ignore_junk,string_read,ier)
  if (ier /= 0) return

  read(string_read,*,iostat=ier) value_to_read

  end subroutine read_value_logical_mesh

!--------------------

  subroutine read_value_string_mesh(iunit,ignore_junk,value_to_read, name, ier)

  use constants, only: MAX_STRING_LEN

  implicit none

  logical :: ignore_junk
  integer :: iunit
  character(len=*) :: value_to_read
  integer :: ier
  character(len=*) :: name
  character(len=MAX_STRING_LEN) :: string_read

  call unused_string(name)
  ier = 0

  call read_next_line(iunit,ignore_junk,string_read,ier)
  if (ier /= 0) return

  value_to_read = string_read

  end subroutine read_value_string_mesh

!--------------------


  subroutine read_value_doubling_integer_mesh(iunit,ignore_junk,value_to_read, name, ier)

  use constants, only: MAX_STRING_LEN

  implicit none

  logical :: ignore_junk
  integer :: iunit
  integer :: value_to_read
  integer :: ier
  character(len=*) :: name
  ! local parameters
  character(len=MAX_STRING_LEN) :: string_read
  character(len=MAX_STRING_LEN) :: value_read

  integer :: index_equal_sign

  string_read = ''
  value_read = ''

  call read_next_line(iunit,ignore_junk,string_read,ier)
  if (ier /= 0) return

  ! checks if line contains name string
  if (index(string_read,trim(name)) > 0) then

    ! suppress leading junk (up to the first equal sign, included)
    index_equal_sign = index(string_read,'=')
    if (index_equal_sign <= 1 .or. index_equal_sign == len_trim(string_read)) then
      stop 'incorrect syntax detected in Mesh_Par_file'
    else
      value_read(1:(len_trim(string_read)-index_equal_sign)) = string_read((index_equal_sign + 1):len_trim(string_read))

      ! suppress leading and trailing white spaces again, if any, after having suppressed the leading junk
      string_read = adjustl(value_read)
      string_read = string_read(1:len_trim(string_read))

      read(string_read,*,iostat=ier) value_to_read
    endif
  else
    ! returns an error
    ier = 1
    return
  endif

  end subroutine read_value_doubling_integer_mesh

!--------------------

  subroutine read_value_doubling_skip_mesh(iunit,ignore_junk, name, ier)

  use constants, only: MAX_STRING_LEN

  implicit none

  logical :: ignore_junk
  integer :: iunit
  integer :: ier
  character(len=*) :: name
  ! local parameters
  character(len=MAX_STRING_LEN) :: string_read

  call read_next_line(iunit,ignore_junk,string_read,ier)
  if (ier /= 0) return

  ! debug
  !print *,'skip line: ',trim(string_read),' index: ',index(string_read,trim(name))

  ! skip more lines
  ! reads next line if current line contains name string
  do while (index(string_read,trim(name)) > 0)
    call read_next_line(iunit,ignore_junk,string_read,ier)
    if (ier /= 0) return
    ! debug
    !print *,'skip line: ',trim(string_read),' index: ',index(string_read,trim(name))
  enddo

  ! now, line already contains new parameter
  ! reverts back file pointer to previous line for the next read statements
  backspace(iunit)

  end subroutine read_value_doubling_skip_mesh

!
!------------------------------------------------------------------------------
!

! interface file

  subroutine read_interface_parameters(iunit,SUPPRESS_UTM_PROJECTION,interface_file, &
                                       npx_interface,npy_interface, &
                                       orig_x_interface,orig_y_interface, &
                                       spacing_x_interface,spacing_y_interface,ier)

  use constants, only: MAX_STRING_LEN,DONT_IGNORE_JUNK

  implicit none

  logical :: SUPPRESS_UTM_PROJECTION
  integer :: iunit
  integer :: ier
  integer :: npx_interface,npy_interface
  double precision :: orig_x_interface,orig_y_interface
  double precision :: spacing_x_interface,spacing_y_interface
  character(len=MAX_STRING_LEN) :: interface_file
  character(len=MAX_STRING_LEN) :: string_read

  ier = 0
  call read_next_line(iunit,DONT_IGNORE_JUNK,string_read,ier)
  if (ier /= 0) return

  read(string_read,*,iostat=ier) SUPPRESS_UTM_PROJECTION,npx_interface,npy_interface, &
             orig_x_interface,orig_y_interface,spacing_x_interface,spacing_y_interface

  call read_value_string_mesh(iunit,DONT_IGNORE_JUNK,interface_file,'INTERFACE_FILE',ier)

  end subroutine read_interface_parameters

!--------------------

! material parameter list

  subroutine read_material_parameters(iunit,material_properties,material_properties_undef,imat,NMATERIALS,ier)

  use constants, only: MAX_STRING_LEN,DONT_IGNORE_JUNK,IDOMAIN_ACOUSTIC,IDOMAIN_ELASTIC,IDOMAIN_POROELASTIC
  use constants_meshfem, only: NUMBER_OF_MATERIAL_PROPERTIES
  implicit none

  integer,intent(in) :: iunit,imat,NMATERIALS
  double precision, dimension(NMATERIALS,NUMBER_OF_MATERIAL_PROPERTIES),intent(inout):: material_properties
  character(len=MAX_STRING_LEN), dimension(NMATERIALS,3),intent(inout) :: material_properties_undef
  integer,intent(out) :: ier

  ! local parameters
  integer :: mat_id
  double precision :: rho,vp,vs,Q_Kappa,Q_mu,anisotropy_flag
  double precision :: kxx,kxy,kxz,kyy,kyz,kzz
  double precision :: phi,tort,rho_s,rho_f,kappa_s,kappa_f,mu_fr,kappa_fr,eta
  integer :: domain_id
  integer :: i,counter,idummy
  character(len=MAX_STRING_LEN) :: string_read
  character(len=128) :: undef_keyword,undef_domain,undef_tomofile
  logical :: new_number
  logical, external :: is_numeric,is_digit_alpha

  ! initializes
  material_properties(imat,:) = 0.d0
  material_properties_undef(imat,:) = ""
  ier = 0

  rho = 0.d0
  vp = 0.d0
  vs = 0.d0
  Q_mu = 0.d0
  Q_Kappa = 0.d0
  anisotropy_flag = 0

  rho_s = 0.d0
  rho_f = 0.d0
  kappa_s = 0.d0
  kappa_f = 0.d0
  kappa_fr = 0.d0
  mu_fr = 0.d0
  eta = 0.d0
  phi = 0.d0
  tort = 0.d0
  kxx = 0.d0
  kxy = 0.d0
  kxz = 0.d0
  kyy = 0.d0
  kyz = 0.d0
  kzz = 0.d0

  mat_id = 0

  call read_next_line(iunit,DONT_IGNORE_JUNK,string_read,ier)
  if (ier /= 0) return

  ! counts numbers on line
  ! example:  acoustic/elastic
  !   1 1020.0  1500.0 0.0  9999.0 9999.0  0  1
  !   -> 8 numbers
  ! example: poroelastic
  !   1           2500.d0 1020.d0 0.4 2.0  1d-11 0.0  0.0 1d-11 0.0 1d-11 16.0554d9 2.295d9  10.0d9    0.0  9.63342d9  3
  !   -> 17 numbers
  ! example: tomographic
  !   -1 tomography elastic tomography_model.xyz 0 2
  !   -> 3 numbers (note that tomography_model_1.xyz 0 2 would count 4)
  counter = 0
  new_number = .true.
  do i = 1,len_trim(string_read)
    ! checks if character is a number
    if (is_digit_alpha(string_read(i:i))) then
      ! new numbers must start with a numeric
      if (new_number .and. is_numeric(string_read(i:i))) then
        ! new number starts
        counter = counter + 1
        new_number = .false.
      endif
    else
      ! a space, tab, .., between numbers reset flag
      if (.not. new_number) new_number = .true.
    endif
  enddo
  ! debug
  !print *,'debug: material parameters ',counter,' items: string ***',trim(string_read),'***'

  domain_id = 0
  ier = 0
  select case (counter)
  case (8)
    ! default line
    ! format: #material_id  #rho  #vp  #vs  #Q_Kappa  #Q_mu  #anisotropy_flag  #domain_id
    read(string_read,*,iostat=ier) mat_id,rho,vp,vs,Q_Kappa,Q_mu,anisotropy_flag,domain_id
    if (domain_id /= IDOMAIN_ACOUSTIC .and. domain_id /= IDOMAIN_ELASTIC) &
      stop 'Error material parameters for acoustic/elastic domains must have domain_id == 1 or 2'

  case (17)
    ! poroelastic line
    ! format:
    ! #material_id #rho_s #rho_f #phi #tort #kxx #kxy #kxz #kyy #kyz #kzz #kappa_s  #kappa_f #kappa_fr #eta #mu_fr    #domain_id
    ! 1           2500.d0 1020.d0 0.4 2.0  1d-11 0.0  0.0 1d-11 0.0 1d-11 16.0554d9 2.295d9  10.0d9    0.0  9.63342d9  3
    read(string_read,*,iostat=ier) mat_id,rho_s,rho_f,phi,tort,kxx,kxy,kxz,kyy,kyz,kzz, &
                                   kappa_s,kappa_f,kappa_fr,eta,mu_fr,domain_id
    if (domain_id /= IDOMAIN_POROELASTIC) &
      stop 'Error material parameters for poroelastic domains must have domain_id == 3'

  case (3,4)
    ! tomography line
    ! format: #material_id #type-keyword #domain-name #tomo-filename #tomo_id #domain_id
    ! example: -1 tomography elastic tomography_model.xyz 0 2
    !       or -1 tomography elastic tomography_model_1.xyz 0 2
    read(string_read,*,iostat=ier) mat_id,undef_keyword,undef_domain,undef_tomofile,idummy,domain_id
    ! sets dummy values (will be overimposed in xgenerate_databases when tomography values are taken)
    if (domain_id == IDOMAIN_ACOUSTIC) then
      ! acoustic (zero vs)
      rho = 1.d0
      vp = 1.d0
      vs = 0.d0
    else if (domain_id == IDOMAIN_ELASTIC) then
      ! elastic
      rho = 1.d0
      vp = 1.d0
      vs = 1.d0
    else
      ! poroelastic tomography models not supported yet
      print *,'Error: tomography material only supports acoustic or elastic domains. Please check domain_id value.'
      print *,'   line: ',trim(string_read)
      stop 'Error in input tomography material'
    endif
    ! checks with consistency between domain-name and domain_id
    if ( (domain_id == IDOMAIN_ACOUSTIC .and. trim(undef_domain) /= 'acoustic') .or. &
         (domain_id == IDOMAIN_ELASTIC .and. trim(undef_domain) /= 'elastic') .or. &
         (domain_id == IDOMAIN_POROELASTIC .and. trim(undef_domain) /= 'poroelastic') ) then
      print *,'Error: tomography material has inconsistent domain_id and domain-name. Please check domain given:'
      print *,'  line: ',trim(string_read)
      print *,'  domain name: ',trim(undef_domain),' and domain id: ',domain_id
      stop 'Error in tomography material domain definition'
    endif

  case default
    print *,'Error: material parameter line format not recognized:'
    print *,'  line: ',trim(string_read)
    print *,'  counted ',counter,' separate numbers on this line.'
    print *
    ! sets a read error
    ier = 1
  end select

  ! check
  if (ier /= 0) then
    print *,'Error while reading your input Mesh_Par_file in routine read_material_parameters()'
    print *
    print *,'For acoustic/elastic materials, we recently changed the input format from:'
    print *,'  mat_id,rho,vp,vs,Q_mu,anisotropy_flag,domain_id'
    print *,'to:'
    print *,'  mat_id,rho,vp,vs,Q_Kappa,Q_mu,anisotropy_flag,domain_id'
    print *,'in order to add support for Q_Kappa.'
    print *
    print *,'For poroelastic materials, we use the input format:'
    print *,'  mat_id,rho_s,rho_f,phi,tort,kxx,kxy,kxz,kyy,kyz,kzz,kappa_s,kappa_f,kappa_fr,eta,mu_fr,domain_id'
    print *
    print *,'It is likely that your input file still uses the old convention and needs to be updated.'
    print *,'If you do not know what value to add for Q_Kappa, add 9999., i.e negligible Q_Kappa attenuation'
    print *,'and then your results will be unchanged compared to older versions of the code.'
    print *
    stop 'Error in input Mesh_Par_file in routine read_material_parameters()'
  endif

  ! checks domain
  if (domain_id /= IDOMAIN_ACOUSTIC .and. domain_id /= IDOMAIN_ELASTIC .and. domain_id /= IDOMAIN_POROELASTIC) then
    print *,'Error: material has invalid domain_id, must be either 1 == acoustic,2 == elastic or 3 == poroelastic'
    print *,'  line : ',trim(string_read)
    stop 'Invalid domain_id in material parameters'
  endif

  ! stores material
  select case(domain_id)
  case (IDOMAIN_ACOUSTIC,IDOMAIN_ELASTIC)
    ! acoustic/elastic material properties
    ! input format: #material_id  #rho  #vp  #vs  #Q_Kappa  #Q_mu  #anisotropy_flag  #domain_id
    material_properties(imat,1) = rho
    material_properties(imat,2) = vp
    material_properties(imat,3) = vs
    material_properties(imat,4) = Q_Kappa
    material_properties(imat,5) = Q_mu
    material_properties(imat,6) = anisotropy_flag
    material_properties(imat,7) = domain_id
    material_properties(imat,8) = mat_id
  case (IDOMAIN_POROELASTIC)
    ! poroelastic
    ! input format: #material_id #rho_s #rho_f #phi #tort #kxx #kxy #kxz #kyy #kyz #kzz ..
    !                                                        .. #kappa_s  #kappa_f #kappa_fr #eta #mu_fr #domain_id
    ! output format:
    !   rhos,rhof,phi,tort,eta,0,material_domain_id,kxx,kxy,kxz,kyy,kyz,kzz,kappas,kappaf,kappafr,mufr
    material_properties(imat,1) = rho_s      ! rho
    material_properties(imat,2) = rho_f      ! vp
    material_properties(imat,3) = phi        ! vs
    material_properties(imat,4) = tort       ! Q_Kappa - not implemented yet
    material_properties(imat,5) = eta        ! Q_mu - not implemented yet
    material_properties(imat,6) = 0          ! anisotropy_flag - not implemented yet
    material_properties(imat,7) = domain_id  ! must be position 7
    material_properties(imat,8) = mat_id     ! must be position 8
    material_properties(imat,9) = kxx
    material_properties(imat,10) = kxy
    material_properties(imat,11) = kxz
    material_properties(imat,12) = kyy
    material_properties(imat,13) = kyz
    material_properties(imat,14) = kzz
    material_properties(imat,15) = kappa_s
    material_properties(imat,16) = kappa_f
    material_properties(imat,17) = kappa_fr
    material_properties(imat,18) = mu_fr      ! anisotropy_flag - not implemented yet
  end select

  ! undefined/tomography materials
  if (mat_id < 0) then
    ! undefined material:
    ! input format: #material_id #type-keyword #domain-name #tomo-filename #tomo_id #domain_id
    !      example: -1 tomography elastic tomography_model.xyz 0 2
    material_properties_undef(imat,1) = undef_keyword
    material_properties_undef(imat,2) = undef_domain
    material_properties_undef(imat,3) = undef_tomofile
  endif

  end subroutine read_material_parameters

!--------------------

  logical function is_numeric(char)

  ! returns .true. if input character is a number

  implicit none
  character(len=1), intent(in) :: char

  is_numeric = .false.

  if (index('0123456789', char) /= 0) then
    is_numeric = .true.
  endif

  end function


  logical function is_digit_alpha(char)

  ! returns .true. if input character is a number or a '.' or a 'd' or a 'e' or a '-'
  ! examples: 5.0, 1.d0, 7.e-12

  implicit none
  character(len=1), intent(in) :: char

  is_digit_alpha = .false.

  if (index('0123456789.de-', char) /= 0) then
    is_digit_alpha = .true.
  endif

  end function


!--------------------

! region parameter list

  subroutine read_region_parameters(iunit,ix_beg_region,ix_end_region,iy_beg_region,iy_end_region, &
                                    iz_beg_region,iz_end_region,imaterial_number,ier)

  use constants, only: MAX_STRING_LEN,DONT_IGNORE_JUNK

  implicit none

  integer,intent(in) :: iunit
  integer,intent(out) :: ix_beg_region,ix_end_region,iy_beg_region,iy_end_region
  integer,intent(out) :: iz_beg_region,iz_end_region,imaterial_number
  integer,intent(out) :: ier

  ! local parameters
  character(len=MAX_STRING_LEN) :: string_read

  ier = 0
  call read_next_line(iunit,DONT_IGNORE_JUNK,string_read,ier)
  if (ier /= 0) return

  ! format: #NEX_XI_BEGIN  #NEX_XI_END  #NEX_ETA_BEGIN  #NEX_ETA_END  #NZ_BEGIN #NZ_END  #material_id
  read(string_read,*,iostat=ier) ix_beg_region,ix_end_region,iy_beg_region,iy_end_region, &
                                 iz_beg_region,iz_end_region,imaterial_number

  end subroutine read_region_parameters

!--------------------

  subroutine read_next_line(iunit,suppress_junk,string_read, ier)

  use constants, only: MAX_STRING_LEN

  implicit none

  integer :: iunit,ier
  logical :: suppress_junk
  character(len=MAX_STRING_LEN) :: string_read
  integer :: index_equal_sign

  ier = 0
  do
    read(unit=iunit,fmt="(a)",iostat=ier) string_read
    if (ier /= 0) then
      print *,'Error reading parameter file Mesh_Par_file: missing next line. Please check file...'
      stop 'Error while reading parameter file Mesh_Par_file'
    endif

! suppress leading white spaces, if any
    string_read = adjustl(string_read)

! suppress trailing carriage return (ASCII code 13) if any (e.g. if input text file coming from Windows/DOS)
    if (index(string_read,achar(13)) > 0) string_read = string_read(1:index(string_read,achar(13))-1)

! exit loop when we find the first line that is not a comment or a white line
    if (len_trim(string_read) == 0) cycle
    if (string_read(1:1) /= '#') exit
  enddo

! suppress trailing white spaces, if any
  string_read = string_read(1:len_trim(string_read))

! suppress trailing comments, if any
  if (index(string_read,'#') > 0) string_read = string_read(1:index(string_read,'#')-1)

  if (suppress_junk) then
! suppress leading junk (up to the first equal sign, included)
    index_equal_sign = index(string_read,'=')
    if (index_equal_sign <= 1 .or. index_equal_sign == len_trim(string_read)) then
      print *,'Error reading Mesh_Par_file line: ',trim(string_read)
      stop 'Error incorrect syntax detected in Mesh_Par_file'
    else
      string_read = string_read((index_equal_sign + 1):len_trim(string_read))
    endif
  endif

! suppress leading and trailing white spaces again, if any, after having suppressed the leading junk
  string_read = adjustl(string_read)
  string_read = string_read(1:len_trim(string_read))

  end subroutine read_next_line

!--------------------

  subroutine open_parameter_file_mesh(filename)

  use constants, only: IIN,MAX_STRING_LEN

  implicit none

  character(len=MAX_STRING_LEN),intent(in) :: filename

  ! local parameters
  integer :: ier

  open(unit=IIN,file=trim(filename),status='old',action='read',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening file: ',trim(filename)
    print *
    print *,'Please check your file path and run-directory.'
    stop 'Error opening Mesh_Par_file'
  endif

  end subroutine open_parameter_file_mesh

!--------------------

  subroutine close_parameter_file_mesh

  use constants, only: IIN

  close(IIN)

  end subroutine close_parameter_file_mesh

!--------------------

! dummy subroutine to avoid warnings about variable not used in other subroutines
  subroutine unused_string(s)

  character(len=*) s

  if (len(s) == 1) continue

  end subroutine unused_string

