!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  1 . 4
!               ---------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!         (c) California Institute of Technology September 2006
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
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
! United States and French Government Sponsorship Acknowledged.

  subroutine save_adjoint_kernels()

  use specfem_par
  use specfem_par_elastic
  use specfem_par_acoustic
  
  implicit none
  integer:: ispec,i,j,k,iglob
  
  ! finalizes calculation of rhop, beta, alpha kernels
  do ispec = 1, NSPEC_AB
  
    ! elastic simulations
    if( ispec_is_elastic(ispec) ) then
  
      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX
            iglob = ibool(i,j,k,ispec)
            
            ! isotropic adjoint kernels (see e.g. Tromp et al. 2005)
            
            ! density kernel
            ! multiplies with rho
            rho_kl(i,j,k,ispec) = - rho_vs(i,j,k,ispec)**2 / mustore(i,j,k,ispec) * rho_kl(i,j,k,ispec) 
            
            ! shear modulus kernel
            mu_kl(i,j,k,ispec) = - mustore(i,j,k,ispec) * mu_kl(i,j,k,ispec)
            
            ! bulk modulus kernel
            kappa_kl(i,j,k,ispec) = - kappastore(i,j,k,ispec) * kappa_kl(i,j,k,ispec)
            
            ! density prime kernel
            rhop_kl(i,j,k,ispec) = rho_kl(i,j,k,ispec) + kappa_kl(i,j,k,ispec) + mu_kl(i,j,k,ispec)
            
            ! vs kernel
            beta_kl(i,j,k,ispec) = 2._CUSTOM_REAL * (mu_kl(i,j,k,ispec) &
                  - 4._CUSTOM_REAL * mustore(i,j,k,ispec) &
                    / (3._CUSTOM_REAL * kappastore(i,j,k,ispec)) * kappa_kl(i,j,k,ispec))
                  
            ! vp kernel
            alpha_kl(i,j,k,ispec) = 2._CUSTOM_REAL * (1._CUSTOM_REAL &
                  + 4._CUSTOM_REAL * mustore(i,j,k,ispec) &
                    / (3._CUSTOM_REAL * kappastore(i,j,k,ispec))) * kappa_kl(i,j,k,ispec)
          enddo
        enddo
      enddo

    endif ! elastic

    ! acoustic simulations
    if( ispec_is_acoustic(ispec) ) then
  
      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX
            ! rho prime kernel
            rhop_ac_kl(i,j,k,ispec) = rho_ac_kl(i,j,k,ispec) + kappa_ac_kl(i,j,k,ispec)
            
            ! vp kernel
            alpha_ac_kl(i,j,k,ispec) = TWO *  kappa_ac_kl(i,j,k,ispec)
          enddo
        enddo
      enddo

    endif ! acoustic

    
  enddo

  ! save kernels to binary files  
  if( ELASTIC_SIMULATION ) then
    open(unit=27,file=prname(1:len_trim(prname))//'rho_kernel.bin',status='unknown',form='unformatted')
    write(27) rho_kl
    close(27)
    open(unit=27,file=prname(1:len_trim(prname))//'mu_kernel.bin',status='unknown',form='unformatted')
    write(27) mu_kl
    close(27)
    open(unit=27,file=prname(1:len_trim(prname))//'kappa_kernel.bin',status='unknown',form='unformatted')
    write(27) kappa_kl
    close(27)
    open(unit=27,file=prname(1:len_trim(prname))//'rhop_kernel.bin',status='unknown',form='unformatted')
    write(27) rhop_kl
    close(27)
    open(unit=27,file=prname(1:len_trim(prname))//'beta_kernel.bin',status='unknown',form='unformatted')
    write(27) beta_kl
    close(27)
    open(unit=27,file=prname(1:len_trim(prname))//'alpha_kernel.bin',status='unknown',form='unformatted')
    write(27) alpha_kl
    close(27)

    if (SAVE_MOHO_MESH) then
      open(unit=27,file=prname(1:len_trim(prname))//'moho_kernel.bin',status='unknown',form='unformatted')
      write(27) moho_kl
      close(27)
    endif

  endif


  ! save kernels to binary files  
  if( ACOUSTIC_SIMULATION ) then
    open(unit=27,file=prname(1:len_trim(prname))//'rho_acoustic_kernel.bin',status='unknown',form='unformatted')
    write(27) rho_ac_kl
    close(27)
    open(unit=27,file=prname(1:len_trim(prname))//'kappa_acoustic_kernel.bin',status='unknown',form='unformatted')
    write(27) kappa_ac_kl
    close(27)
    open(unit=27,file=prname(1:len_trim(prname))//'rho_prime_acoustic_kernel.bin',status='unknown',form='unformatted')
    write(27) rhop_ac_kl
    close(27)
    open(unit=27,file=prname(1:len_trim(prname))//'alpha_acoustic_kernel.bin',status='unknown',form='unformatted')
    write(27) alpha_ac_kl
    close(27)

  endif

  end subroutine save_adjoint_kernels
