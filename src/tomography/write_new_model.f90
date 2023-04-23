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


subroutine write_new_model_iso()

! file output for new model

  use tomography_model_iso

  implicit none
  real(kind=CUSTOM_REAL) :: min_vp,min_vs,max_vp,max_vs,min_rho,max_rho
  character(len=MAX_STRING_LEN) :: m_file, fname

  ! user output
  if (myrank == 0) print *,'writing out new model...'

  ! vp model
  call max_all_cr(maxval(model_vp_new),max_vp)
  call min_all_cr(minval(model_vp_new),min_vp)

  fname = 'vp_new'
  write(m_file,'(a,i6.6,a)') trim(OUTPUT_MODEL_DIR)//'proc',myrank,trim(REG)//trim(fname)//'.bin'
  if (myrank == 0) print *,'  ',trim(OUTPUT_MODEL_DIR)//'proc**'//trim(REG)//trim(fname)//'.bin'

  open(IOUT,file=trim(m_file),form='unformatted',action='write')
  write(IOUT) model_vp_new
  close(IOUT)

  ! vs model
  call max_all_cr(maxval(model_vs_new),max_vs)
  call min_all_cr(minval(model_vs_new),min_vs)

  fname = 'vs_new'
  write(m_file,'(a,i6.6,a)') trim(OUTPUT_MODEL_DIR)//'proc',myrank,trim(REG)//trim(fname)//'.bin'
  if (myrank == 0) print *,'  ',trim(OUTPUT_MODEL_DIR)//'proc**'//trim(REG)//trim(fname)//'.bin'

  open(IOUT,file=trim(m_file),form='unformatted',action='write')
  write(IOUT) model_vs_new
  close(IOUT)

  ! rho model
  call max_all_cr(maxval(model_rho_new),max_rho)
  call min_all_cr(minval(model_rho_new),min_rho)

  fname = 'rho_new'
  write(m_file,'(a,i6.6,a)') trim(OUTPUT_MODEL_DIR)//'proc',myrank,trim(REG)//trim(fname)//'.bin'
  if (myrank == 0) print *,'  ',trim(OUTPUT_MODEL_DIR)//'proc**'//trim(REG)//trim(fname)//'.bin'

  open(IOUT,file=trim(m_file),form='unformatted',action='write')
  write(IOUT) model_rho_new
  close(IOUT)

  if (myrank == 0) then
    print *
    print *,'new models:'
    print *,'  vp min/max : ',min_vp,max_vp
    print *,'  vs min/max : ',min_vs,max_vs
    print *,'  rho min/max: ',min_rho,max_rho
    print *
  endif
  call synchronize_all()

  if (PRINT_STATISTICS_FILES .and. myrank == 0) then
    open(IOUT,file=trim(OUTPUT_STATISTICS_DIR)//'statistics_vs_vp_rho_new_minmax',status='unknown')
    write(IOUT,*) '#min_vs #max_vs #min_vp #max_vp #min_rho #max_rho'
    write(IOUT,'(6e24.12)') min_vs,max_vs,min_vp,max_vp,min_rho,max_rho
    close(IOUT)
  endif

end subroutine write_new_model_iso


!
!-------------------------------------------------------------------------------------------------
!

subroutine write_new_model_tiso()

! file output for TI new model

  use tomography_model_tiso
  implicit none
  real(kind=CUSTOM_REAL) :: min_vpv,min_vph,min_vsv,min_vsh, &
    max_vpv,max_vph,max_vsv,max_vsh,min_eta,max_eta, min_rho,max_rho
  character(len=MAX_STRING_LEN) :: m_file, fname

  ! user output
  if (myrank == 0) print *,'writing out new model...'

  ! vpv model
  call max_all_cr(maxval(model_vpv_new),max_vpv)
  call min_all_cr(minval(model_vpv_new),min_vpv)

  fname = 'vpv_new'
  write(m_file,'(a,i6.6,a)') trim(OUTPUT_MODEL_DIR)//'proc',myrank,trim(REG)//trim(fname)//'.bin'
  if (myrank == 0) print *,'  ',trim(OUTPUT_MODEL_DIR)//'proc**'//trim(REG)//trim(fname)//'.bin'

  open(IOUT,file=trim(m_file),form='unformatted',action='write')
  write(IOUT) model_vpv_new
  close(IOUT)

  ! vph model
  call max_all_cr(maxval(model_vph_new),max_vph)
  call min_all_cr(minval(model_vph_new),min_vph)

  fname = 'vph_new'
  write(m_file,'(a,i6.6,a)') trim(OUTPUT_MODEL_DIR)//'proc',myrank,trim(REG)//trim(fname)//'.bin'
  if (myrank == 0) print *,'  ',trim(OUTPUT_MODEL_DIR)//'proc**'//trim(REG)//trim(fname)//'.bin'

  open(IOUT,file=trim(m_file),form='unformatted',action='write')
  write(IOUT) model_vph_new
  close(IOUT)

  ! vsv model
  call max_all_cr(maxval(model_vsv_new),max_vsv)
  call min_all_cr(minval(model_vsv_new),min_vsv)

  fname = 'vsv_new'
  write(m_file,'(a,i6.6,a)') trim(OUTPUT_MODEL_DIR)//'proc',myrank,trim(REG)//trim(fname)//'.bin'
  if (myrank == 0) print *,'  ',trim(OUTPUT_MODEL_DIR)//'proc**'//trim(REG)//trim(fname)//'.bin'

  open(IOUT,file=trim(m_file),form='unformatted',action='write')
  write(IOUT) model_vsv_new
  close(IOUT)

  ! vsh model
  call max_all_cr(maxval(model_vsh_new),max_vsh)
  call min_all_cr(minval(model_vsh_new),min_vsh)

  fname = 'vsh_new'
  write(m_file,'(a,i6.6,a)') trim(OUTPUT_MODEL_DIR)//'proc',myrank,trim(REG)//trim(fname)//'.bin'
  if (myrank == 0) print *,'  ',trim(OUTPUT_MODEL_DIR)//'proc**'//trim(REG)//trim(fname)//'.bin'

  open(IOUT,file=trim(m_file),form='unformatted',action='write')
  write(IOUT) model_vsh_new
  close(IOUT)

  ! eta model
  call max_all_cr(maxval(model_eta_new),max_eta)
  call min_all_cr(minval(model_eta_new),min_eta)

  fname = 'eta_new'
  write(m_file,'(a,i6.6,a)') trim(OUTPUT_MODEL_DIR)//'proc',myrank,trim(REG)//trim(fname)//'.bin'
  if (myrank == 0) print *,'  ',trim(OUTPUT_MODEL_DIR)//'proc**'//trim(REG)//trim(fname)//'.bin'

  open(IOUT,file=trim(m_file),form='unformatted',action='write')
  write(IOUT) model_eta_new
  close(IOUT)

  ! rho model
  call max_all_cr(maxval(model_rho_new),max_rho)
  call min_all_cr(minval(model_rho_new),min_rho)

  fname = 'rho_new'
  write(m_file,'(a,i6.6,a)') trim(OUTPUT_MODEL_DIR)//'proc',myrank,trim(REG)//trim(fname)//'.bin'
  if (myrank == 0) print *,'  ',trim(OUTPUT_MODEL_DIR)//'proc**'//trim(REG)//trim(fname)//'.bin'

  open(IOUT,file=trim(m_file),form='unformatted',action='write')
  write(IOUT) model_rho_new
  close(IOUT)

  ! user output
  if (myrank == 0) then
    print *
    print *,'new models:'
    print *,'  vpv min/max: ',min_vpv,max_vpv
    print *,'  vph min/max: ',min_vph,max_vph
    print *,'  vsv min/max: ',min_vsv,max_vsv
    print *,'  vsh min/max: ',min_vsh,max_vsh
    print *,'  eta min/max: ',min_eta,max_eta
    print *,'  rho min/max: ',min_rho,max_rho
    print *
  endif
  call synchronize_all()

  if (PRINT_STATISTICS_FILES .and. myrank == 0) then
    open(IOUT,file=trim(OUTPUT_STATISTICS_DIR)//'statistics_vs_vp_rho_new_minmax',status='unknown')
    write(IOUT,*) '#min_vsv #max_vsv #min_vsh #max_vsh #min_vpv #max_vpv #min_vph #max_vph ' &
               // '#min_eta #max_eta #min_rho #max_rho'
    write(IOUT,'(12e24.12)') min_vsv,max_vsv,min_vsh,max_vsh,min_vpv,max_vpv,min_vph,max_vph, &
                             min_eta,max_eta,min_rho,max_rho
    close(IOUT)
  endif

end subroutine write_new_model_tiso


