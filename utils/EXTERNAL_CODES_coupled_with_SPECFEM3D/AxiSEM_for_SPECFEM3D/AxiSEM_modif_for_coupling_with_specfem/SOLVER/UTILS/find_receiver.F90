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

!===========================
program find_receiver
!===========================

use nc_routines
implicit none
real :: Mij(6),period,srcdepth,srclat,srclon,time_shift,reclon,reclat,deg2rad,rad2deg,rjunk
real :: period_sim,shift_fact,magnitude_sim,src_depth_sim,dt,mij_phi(3,4),ph_src,th_src,rot_mat(3,3)
real :: theta_over_dtheta,dtheta_rec,dtheta_off,argu,lam_deg,sumw
real, allocatable :: th_rec(:),ph_rec(:),th_rec_rot(:),ph_rec_rot(:),t(:),epidist(:),deltatheta(:)
real, allocatable :: seis(:,:),seis_snglcomp(:,:,:),seistmp(:,:),w(:)
character(len=30) :: junk,eventname,src_header(12)
character(len=3) :: components
character(len=200) :: dirname,interp_method
character(len=30),  allocatable :: recname(:)
integer :: nrec,irec,nt,i,dirind,nr_per_deg,ind1,ind2,ind_rec,ind_rec2,it,ijunk
integer :: nrec_sim,ishift_deltat,ishift_seisdt,num_interp,count_neg,p
integer, dimension(4) :: ncid_in, nc_disp_varid
character(len=4) :: ind_recchar,appidur
character(len=34) :: stf_type,colat,model
character(len=1),dimension(3) :: reccomp
logical, allocatable :: mask_min(:)
logical :: ljunk,interpolate_seis,even_spaced_theta,theta_discrete_smaller
logical :: usenetcdf
real, parameter :: pi = 3.14159265

! >>>>> Static parameters <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
interpolate_seis=.false.
num_interp = 2
! choose from spline (linear),sinc,idw (inverse distance weighting)
interp_method='sinc'
p=3

if (trim(interp_method)=='spline') num_interp=2
if (.not. interpolate_seis) interp_method='closest'

! >>>>> Receiver information <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
open(unit=99,file='receiver_info.dat')
  read(99,*)nrec,components,colat,model
  call def_rec_comp(components,reccomp)
  allocate(recname(nrec),ph_rec(nrec),th_rec(nrec))
  deg2rad=pi/180. ;  rad2deg=180./pi
  do irec=1,nrec
     read(99,*)recname(irec),reclat,reclon
     th_rec(irec)=reclat*deg2rad; ph_rec(irec)=reclon*deg2rad
  enddo
close(99)
if (colat/='colat') th_rec=pi/2.-th_rec

! >>>>> Source information <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
open(unit=20000,file='CMTSOLUTION',POSITION='REWIND',status='old')
  read(20000,*)stf_type,src_header(1:12)
  read(20000,*)junk,junk,eventname
  read(20000,*)junk,junk,time_shift
  read(20000,*)junk,junk,period
  read(20000,*)junk,srclat
  read(20000,*)junk,srclon
  read(20000,*)junk,srcdepth
  read(20000,*)junk,Mij(1) !Mrr
  read(20000,*)junk,Mij(2) !Mtt
  read(20000,*)junk,Mij(3) !Mpp
  read(20000,*)junk,Mij(4) !Mrt
  read(20000,*)junk,Mij(5) !Mrp
  read(20000,*)junk,Mij(6) !Mtp
close(20000)
Mij=Mij/1.E7  ! CMTSOLUTION given in dyn-cm
Mij=Mij/1.E20 ! simulated sources are for 1E20.
write(6,*)'maximum Mij:',maxval(Mij)

! >>>>> Find closest depth and model <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
call depths(srcdepth,model,dirname)

! >>>>> load simulation info <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
open(unit=99,file=trim(dirname)//'MZZ/simulation.info')
  read(99,*)junk
  read(99,*)rjunk
  read(99,*)ijunk
  read(99,*)junk
  read(99,*)junk
  read(99,*)junk
  read(99,*)junk
  read(99,*)period_sim
  read(99,*)src_depth_sim
  read(99,*)magnitude_sim
  read(99,*)nrec_sim
  read(99,*)nt
  read(99,*)dt
  read(99,*)ljunk
  read(99,*)ijunk
  read(99,*)rjunk
  read(99,*)ijunk
  read(99,*)rjunk
  read(99,*)junk
  read(99,*)ijunk
  read(99,*)ijunk
  read(99,*)shift_fact
  read(99,*)ishift_deltat
  read(99,*)ishift_seisdt
  read(99,*)ijunk
  read(99,*)junk
  read(99,*)dtheta_rec
  read(99,*)usenetcdf
close(99)
allocate(t(nt),seis(nt,3),seis_snglcomp(nt,3,4))
do it=1,nt
   t(it)=dt*real(it)
enddo
shift_fact = shift_fact+1.5*period
write(6,*)'shift factor:',shift_fact
call define_io_appendix(appidur,int(period))

! >>>>> Rotate & extract epicentral distance <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
th_src=(90.-srclat)*deg2rad; ph_src=srclon*deg2rad
allocate(th_rec_rot(nrec),ph_rec_rot(nrec))
call rotate_src_rec(th_src,ph_src,nrec,th_rec,ph_rec,th_rec_rot,ph_rec_rot,rot_mat)

even_spaced_theta=.false.
if (dtheta_rec>0.) even_spaced_theta=.true.
if (.not. even_spaced_theta) interpolate_seis=.false.

! plot original source and receiver locations in google earth kml file with link to seismogram images
  call save_google_earth_kml(th_src,ph_src,src_depth_sim,Mij,period,th_rec,ph_rec,interp_method, &
                             appidur,th_rec_rot,ph_rec_rot,reccomp,nrec,recname)

! >>>>> prepare receiver index + neighbors <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
if (.not. even_spaced_theta) then ! brute force search based on monotonous guesses. 
   allocate(epidist(nrec_sim),deltatheta(nrec_sim))
   open(unit=98,file=trim(dirname)//'/MZZ/Data/receiver_pts.dat')
   do irec=1,nrec_sim
      read(98,*)epidist(irec),rjunk,ijunk
      if (irec>1) deltatheta(irec-1)=epidist(irec)-epidist(irec-1)
   enddo
   deltatheta(nrec_sim)=deltatheta(nrec_sim-1)
   close(98)
   
   if (nrec_sim>360) then  ! require 2 points per interval
      nr_per_deg = int(real(nrec_sim)/180.)
   elseif (nrec_sim>180) then
      nr_per_deg = int(real(nrec_sim)/90.)
   elseif (nrec_sim>90) then
      nr_per_deg = int(real(nrec_sim)/45.)
   elseif (nrec_sim>45) then 
      nr_per_deg = int(real(nrec_sim)/22.5)
   elseif (nrec_sim>23) then 
      nr_per_deg = int(real(nrec_sim)/11.5)
   endif
   write(6,*)'Min/max epidist [deg]:',minval(epidist),maxval(epidist),nr_per_deg
   write(6,*)'Min/max delta theta [deg]:',minval(deltatheta),maxval(deltatheta)
endif
if (interpolate_seis) allocate(mask_min(nrec_sim))

if (period==0.) then 
   lam_deg = period_sim*5.8/111.195
else
   lam_deg = period*5.8/111.195
endif
write(6,*)'period,P-wavelength: [deg]',lam_deg

if (interpolate_seis) allocate(seistmp(nt,3))
allocate(w(num_interp))

!**************** Open Netcdf output files ******************
if (usenetcdf) call nc_open(dirname,ncid_in,nc_disp_varid)
!************************************************************


!--------------
do irec=1,nrec
!--------------
   theta_discrete_smaller = .false.

! >>>>> find receiver index + neighbors <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   if (.not. even_spaced_theta) then ! brute force minimization - slower
      ind1 = floor(th_rec_rot(irec)*rad2deg)*nr_per_deg 
      ind2 = max(min(nrec_sim,ceiling(th_rec_rot(irec)*rad2deg)*nr_per_deg),ind1+1)
      ind_rec = ind1+minloc(abs(th_rec_rot(irec)*rad2deg - epidist(ind1:ind2)),1)
!      write(6,*)'INDICES:',ind1,ind2,ind_rec
      dtheta_off = th_rec_rot(irec)*rad2deg-epidist(ind_rec)
      write(6,*)'Receiver #, desired/offered th [deg]:',irec,th_rec_rot(irec)*rad2deg,epidist(ind_rec)

   else  ! even-spaced theta... find via distance difference - faster (?)
      ind_rec = floor(th_rec_rot(irec)*rad2deg/dtheta_rec)+1
      if (abs(th_rec_rot(irec)*rad2deg-real(ind_rec-1)*dtheta_rec)>dtheta_rec/2.) ind_rec = ind_rec+1 
      dtheta_off = th_rec_rot(irec)*rad2deg-(ind_rec-1)*dtheta_rec
      write(6,*)'desired/offered th [deg], lam-frac diff :',&
            th_rec_rot(irec)*rad2deg,ind_rec,(ind_rec-1)*dtheta_rec,abs(dtheta_off)/lam_deg
   endif
   if (dtheta_off>0.) theta_discrete_smaller = .true.
   theta_over_dtheta = th_rec_rot(irec)*rad2deg/abs(dtheta_rec)
!   write(6,*)'theta/dtheta:',theta_over_dtheta

! >>>>> calculate moment tensor and azimuth prefactors/radiation patterns <<<<<<<
   call compute_radiation_prefactor(Mij,ph_rec_rot(irec),mij_phi)

! >>>>> load seismograms <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   if (.not. interpolate_seis) then ! choosing closest location
      if (usenetcdf) then
	      call nc_read_seis(ncid_in, nc_disp_varid, ind_rec, nt, seis_snglcomp)
	  else
          call read_seis(dirname,ind_rec,nt,seis_snglcomp)
      end if
      call sum_individual_wavefields(seis,seis_snglcomp,nt,mij_phi)

   else ! interpolation using num_interp locations
      count_neg = max(1,ind_rec-num_interp/2)
!      if (theta_discrete_smaller) count_neg = count_neg -1
      write(6,*)irec,'num_interp,count neg,ind_rec:',num_interp,count_neg,ind_rec
      seis=0.; seistmp=0.
      if (trim(interp_method)=='idw') then 
         do i=1,num_interp
            w(i) = 1./(th_rec_rot(irec) - real(count_neg+i-1)*dtheta_rec)**p
         enddo
         sumw=sum(w)
      endif
      do i=1,num_interp
         write(6,*)'loading seismogram for interpolation at theta=',real(count_neg+i-1)*dtheta_rec
         if (usenetcdf) then
           call nc_read_seis(ncid_in, nc_disp_varid, count_neg+i, nt, seis_snglcomp)
         else
           call read_seis(dirname,count_neg+i,nt,seis_snglcomp)
         end if  
         call sum_individual_wavefields(seistmp,seis_snglcomp,nt,mij_phi)

         if (trim(interp_method)=='sinc') then 
            argu = pi*(theta_over_dtheta-real(count_neg+i))
            seis = seis + sin(argu)/argu*seistmp

         elseif (interp_method=='idw') then  ! inverse distance weighting
            seis = seis + w(i)*seistmp

         elseif (interp_method=='spline') then  ! linear splines
            seistmp = -(seistmp-seis)/dtheta_rec ! slope
            seis = seistmp* th_rec_rot(irec)*rad2deg + seis - seistmp*real(count_neg+i)*dtheta_rec
         endif

      enddo
      if (interp_method=='idw') seis=seis/sumw
   endif

! >>>>> sum, filter, rotate <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   if (period>0.) then 
      call convolve_with_stf(period,dt,nt,'gauss_0',seis,seis_snglcomp(:,:,1))
      seis = seis_snglcomp(:,:,1)
   endif
   call rotate_receiver_comp(components,th_src,ph_src,&
        th_rec_rot(irec),ph_rec_rot(irec),th_rec_rot(irec),ph_rec_rot(irec),nt,rot_mat,seis)

! >>>>> save to file <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   open(unit=20,file=trim(recname(irec))//'_'//trim(interp_method)//'_filtperiod'//appidur//'_'//trim(reccomp(1))//'.ascii')
   open(unit=21,file=trim(recname(irec))//'_'//trim(interp_method)//'_filtperiod'//appidur//'_'//trim(reccomp(2))//'.ascii')
   open(unit=22,file=trim(recname(irec))//'_'//trim(interp_method)//'_filtperiod'//appidur//'_'//trim(reccomp(3))//'.ascii')
   do it=1,nt
      write(20,*)t(it)-shift_fact,seis(it,1)
      write(21,*)t(it)-shift_fact,seis(it,2)
      write(22,*)t(it)-shift_fact,seis(it,3)
   enddo
   close(20); close(21); close(22)
!--------------
enddo ! nrec
!--------------
if (usenetcdf) call nc_close(ncid_in)

write(6,*)' .... DONE.'

!===========================
end program find_receiver
!===========================

!################################################################################################

!--------------------------------------------------------------------
subroutine def_rec_comp(rec_comp_sys,reccomp)

implicit none
character(len=3), intent(in)  :: rec_comp_sys
character(len=1), intent(out) :: reccomp(3)

  if (rec_comp_sys=='sph') then
     reccomp(1)='th'; reccomp(2)='ph'; reccomp(3)='r'
  elseif (rec_comp_sys=='enz') then
     reccomp(1)='N'; reccomp(2)='E'; reccomp(3)='Z'
  elseif (rec_comp_sys=='cyl') then
     reccomp(1)='s'; reccomp(2)='ph'; reccomp(3)='z'
  elseif (rec_comp_sys=='xyz') then
     reccomp(1)='x'; reccomp(2)='y'; reccomp(3)='z'
  elseif (rec_comp_sys=='src') then
     reccomp(1)='R'; reccomp(2)='T'; reccomp(3)='Z'
  endif
  
end subroutine def_rec_comp
!--------------------------------------------------------------------

!--------------------------------------------------------------------
subroutine read_seis(dirname,ind_rec,nt,seis_snglcomp)

implicit none
integer, intent(in)            :: nt,ind_rec
character(len=200), intent(in) :: dirname
real, intent(out)              :: seis_snglcomp(nt,3,4)
integer                        :: it
character(len=4)               :: ind_recchar

call define_io_appendix(ind_recchar,ind_rec)!
!write(6,*)'reading seismogram ',trim(dirname)//'/MZZ/Data/recfile_'//ind_recchar//'_disp.dat'
open(unit=60,file=trim(dirname)//'/MZZ/Data/recfile_'//ind_recchar//'_disp.dat')
open(unit=61,file=trim(dirname)//'/MXX_P_MYY/Data/recfile_'//ind_recchar//'_disp.dat')
open(unit=62,file=trim(dirname)//'/MXZ_MYZ/Data/recfile_'//ind_recchar//'_disp.dat')
open(unit=63,file=trim(dirname)//'/MXY_MXX_M_MYY/Data/recfile_'//ind_recchar//'_disp.dat')
do it=1,nt
   read(60,*)seis_snglcomp(it,1,1),seis_snglcomp(it,3,1)
   read(61,*)seis_snglcomp(it,1,2),seis_snglcomp(it,3,2)
   read(62,*)seis_snglcomp(it,1,3),seis_snglcomp(it,2,3),seis_snglcomp(it,3,3)
   read(63,*)seis_snglcomp(it,1,4),seis_snglcomp(it,2,4),seis_snglcomp(it,3,4)
enddo
close(60); close(61); close(62); close(63)  
end subroutine read_seis
!--------------------------------------------------------------------

!--------------------------------------------------------------------
subroutine rotate_src_rec(th_src,ph_src,nrec,th_rec,ph_rec,th_rec_rot,ph_rec_rot,rot_mat)

implicit none
integer, intent(in) :: nrec
real, intent(in)    :: th_src,ph_src,th_rec(nrec),ph_rec(nrec)
real, intent(out)   :: th_rec_rot(nrec),ph_rec_rot(nrec),rot_mat(3,3)
real                :: x_vec(3),x_vec_rot(3),r,trans_rot_mat(3,3)
integer             :: ircv
real                :: cosph_rec_rot, costh_rec_rot
real, parameter     :: smallval = 1.e-10
real, parameter     :: pi = 3.14159265

! This is the rotation matrix of Nissen-Meyer, Dahlen, Fournier, GJI 2007.
rot_mat(1,1)=cos(th_src)*cos(ph_src)
rot_mat(2,2)=cos(ph_src)
rot_mat(3,3)=cos(th_src)
rot_mat(2,1)=cos(th_src)*sin(ph_src)
rot_mat(3,1)=-sin(th_src)
rot_mat(3,2)=0.
rot_mat(1,2)=-sin(ph_src)
rot_mat(1,3)=sin(th_src)*cos(ph_src)
rot_mat(2,3)=sin(th_src)*sin(ph_src)
trans_rot_mat=transpose(rot_mat)

do ircv=1,nrec
   x_vec(1)=sin(th_rec(ircv))*cos(ph_rec(ircv))
   x_vec(2)=sin(th_rec(ircv))*sin(ph_rec(ircv))
   x_vec(3)=cos(th_rec(ircv))
   x_vec_rot=matmul(trans_rot_mat,x_vec)
   r = sqrt(x_vec_rot(1)**2 + x_vec_rot(2)**2 + x_vec_rot(3)**2)
   costh_rec_rot = x_vec_rot(3) / (r +smallval)
   if (costh_rec_rot.lt.-1.) costh_rec_rot = -1.
   if (costh_rec_rot.gt.1.)  costh_rec_rot = 1.
   th_rec_rot(ircv) = acos(costh_rec_rot)

   cosph_rec_rot = x_vec_rot(1) / (r * sin(th_rec_rot(ircv)) + smallval)
   if (cosph_rec_rot.lt.-1.) cosph_rec_rot = -1.
   if (cosph_rec_rot.gt.1.)  cosph_rec_rot = 1.

   if (x_vec_rot(2) >= 0.) then
      ph_rec_rot(ircv) = acos(cosph_rec_rot)
   else
      ph_rec_rot(ircv) = 2.*pi - acos(cosph_rec_rot)
   end if   
enddo

end subroutine rotate_src_rec
!--------------------------------------------------------------------

!--------------------------------------------------------------------
subroutine compute_radiation_prefactor(Mij,phi,mij_prefact)

implicit none
real, intent(in)  :: Mij(6)
real, intent(in)  :: phi
real, intent(out) :: mij_prefact(3,4)

mij_prefact(:,1) = Mij(1)
mij_prefact(2,1) = 0.
mij_prefact(:,2) = 2.*(Mij(2)+Mij(3))
mij_prefact(2,2) = 0.
mij_prefact(:,3) =  Mij(4)*cos(phi)+Mij(5)*sin(phi)
mij_prefact(2,3) = -Mij(4)*sin(phi)+Mij(5)*cos(phi)
mij_prefact(:,4) = (Mij(2)-Mij(3))*cos(2.*phi)+2.*Mij(6)*sin(2.*phi) 
mij_prefact(2,4) = (Mij(3)-Mij(2))*sin(2.*phi)+2.*Mij(6)*cos(2.*phi)

end subroutine compute_radiation_prefactor
!------------------------------------------------------------------------

!------------------------------------------------------------------------
subroutine sum_individual_wavefields(field_sum,field_in,n,mij_prefact)

implicit none
integer, intent(in)                :: n
real, dimension(n,3,4), intent(in) :: field_in
real, dimension(3,4), intent(in)   :: mij_prefact
real, dimension(n,3), intent(out)  :: field_sum
integer  :: i

field_sum = 0.
do i=1,4
   field_sum(:,1) = field_sum(:,1) + mij_prefact(1,i)*field_in(:,1,i)
   field_sum(:,2) = field_sum(:,2) + mij_prefact(2,i)*field_in(:,2,i)
   field_sum(:,3) = field_sum(:,3) + mij_prefact(3,i)*field_in(:,3,i)
enddo

end subroutine sum_individual_wavefields
!--------------------------------------------------------------------

!--------------------------------------------------------------------
subroutine rotate_receiver_comp(rec_comp_sys,srccolat,srclon,th_rot,ph_rot,th_orig,ph_orig,nt,rot,seis)

implicit none
character(len=3), intent(in) :: rec_comp_sys
real, intent(in)             :: th_rot,ph_rot ! coordinates in the rotated (src at pole) system
real, intent(in)             :: th_orig,ph_orig ! coordinates in the unrotated (actual src) system
real, intent(in)             :: srccolat,srclon,rot(3,3) ! orginal source coordinates
integer, intent(in)          :: nt
real, intent(inout)          :: seis(nt,3)
real                         :: seis_tmp(nt,3),seisvec(3)
integer                      :: i
real, parameter              :: pi = 3.14159265

!write(6,*)'ROTATIONS'
!write(6,*)th_orig*180./pi,ph_orig*180./pi
!write(6,*)th_rot*180./pi,ph_rot*180./pi

! Need to consider *local* spherical geometry in the first place,
! THEN rotate the source-receiver frame to the north pole in the solver.
! E.g., consider the difference between source-projected and spherical coordinates for a source 
! away from the north pole: they are not the same, but in the framework below would 
! be identified as the same.

! Source projected frame: transform to spherical without any further rotations
if (rec_comp_sys=='src') then  
   seis_tmp(:,1) = cos(th_rot) * seis(:,1) - sin(th_rot) * seis(:,3)
   seis_tmp(:,2) = seis(:,2)
   seis_tmp(:,3) = sin(th_rot) * seis(:,1) + cos(th_rot) * seis(:,3)

! Rotate from rotated u_sphiz to rotated u_xyz (both in reference, source-at-pole system) 
else 
   seis_tmp(:,1) = cos(ph_rot) * seis(:,1) - sin(ph_rot) * seis(:,2) 
   seis_tmp(:,2) = sin(ph_rot) * seis(:,1) + cos(ph_rot) * seis(:,2)
   seis_tmp(:,3) = seis(:,3)

   ! Rotate to the original (i.e. real src-rec coordinate-based) u_xyz
   if (srccolat>1.E-10 .or. srclon>1.E-10) then 
      do i=1,nt
         seisvec = seis_tmp(i,:)
         seis_tmp(i,:) = matmul(rot,seisvec)
      enddo
   endif
endif 

! Rotate to coordinate system of choice
if (rec_comp_sys=='enz') then
   seis(:,1) = - cos(th_orig) * cos(ph_orig) * seis_tmp(:,1) &
             & - cos(th_orig) * sin(ph_orig) * seis_tmp(:,2) &
             & + sin(th_orig) * seis_tmp(:,3)
   seis(:,2) = - sin(ph_orig) * seis_tmp(:,1) &
             & + cos(ph_orig) * seis_tmp(:,2)
   seis(:,3) =   sin(th_orig) * cos(ph_orig) * seis_tmp(:,1) & 
             & + sin(th_orig) * sin(ph_orig) * seis_tmp(:,2) &
             & + cos(th_orig) * seis_tmp(:,3)
   

elseif (rec_comp_sys=='sph') then 
   seis(:,1) =   cos(th_orig) * cos(ph_orig) * seis_tmp(:,1) &
             & + cos(th_orig) * sin(ph_orig) * seis_tmp(:,2) &
             & - sin(th_orig) * seis_tmp(:,3)
   seis(:,2) = - sin(ph_orig) * seis_tmp(:,1) & 
             & + cos(ph_orig) * seis_tmp(:,2)
   seis(:,3) =   sin(th_orig) * cos(ph_orig) * seis_tmp(:,1) & 
             & + sin(th_orig) * sin(ph_orig) * seis_tmp(:,2) &
             & + cos(th_orig) * seis_tmp(:,3)

elseif (rec_comp_sys=='cyl') then 
   seis(:,1) =   cos(ph_orig) * seis_tmp(:,1) + sin(ph_orig) * seis_tmp(:,2)
   seis(:,2) = - sin(ph_orig) * seis_tmp(:,1) + cos(ph_orig) * seis_tmp(:,2)
   seis(:,3) =   seis_tmp(:,3)

elseif (rec_comp_sys=='xyz') then
   seis = seis_tmp

elseif (rec_comp_sys=='src') then
   seis = seis_tmp ! taken from above

else
   write(6,*)'unknown component system',rec_comp_sys
   stop
endif

end subroutine rotate_receiver_comp
!--------------------------------------------------------------------

!-----------------------------------------------------------------------------
subroutine convolve_with_stf(t_0,dt,nt,stf,seis,seis_fil)          

implicit none
integer, intent(in)            :: nt
real, intent(in)               :: t_0,dt
real                           :: source,sqrt_pi_inv
integer                        :: i,j,N_j,irec,lffile,ishift
real, intent(in)               :: seis(nt,3)
real, intent(out)              :: seis_fil(nt,3)
real                           :: temp_expo,alpha
character(len=7), intent(in)   :: stf
character(len=4)               :: appidur,appirec
real, parameter                :: decay=3.5d0
real, parameter                :: pi = 3.14159265
real                           :: shiftfact

  shiftfact =  1.5*t_0
!  write(6,*)
!  write(6,*)'Convolving with period=',t_0
!  write(6,*)'convolve:',stf,maxval(seis)
!  write(6,*)'shift factor:',shiftfact
  N_j=int(3.*t_0/dt)
  call define_io_appendix(appidur,int(t_0))
  alpha=decay/t_0
  sqrt_pi_inv=alpha*dt/sqrt(pi)/pi
  do i=1,nt
    seis_fil(i,:)=0.
    do j=1,N_j
       temp_expo=alpha*(real(j)*dt-shiftfact)
       source = exp(-temp_expo**2 )*sqrt_pi_inv 
       if(i > j .and. i-j <= nt) seis_fil(i,:) = seis_fil(i,:)+seis(i-j,:)*source
    enddo
  enddo
  seis_fil=seis_fil*pi

end subroutine convolve_with_stf
!--------------------------------------------------------------------------


!-----------------------------------------------------------------------------
subroutine save_google_earth_kml(srccolat1,srclon1,srcdepth,Mij,per,rcvcolat,rcvlon,interp_method, &
                                appidur,rcvcolat_rot,rcvlon_rot,reccomp,num_rec_glob,receiver_name)

implicit none

integer, intent(in) :: num_rec_glob
real, intent(in) :: srccolat1,srclon1,srcdepth
real, intent(in) :: rcvcolat(1:num_rec_glob),rcvlon(1:num_rec_glob)
real, intent(in) :: rcvcolat_rot(1:num_rec_glob),rcvlon_rot(1:num_rec_glob)
character(len=200), intent(in) :: interp_method
real, intent(in) :: Mij(6),per
character(len=30), intent(in) :: receiver_name(1:num_rec_glob)
character(len=1),dimension(3), intent(in) :: reccomp
real :: slon,slat,rlon(1:num_rec_glob),rlat(1:num_rec_glob)
real :: rlon_rot(1:num_rec_glob),rlat_rot(1:num_rec_glob)
integer :: i
character(len=4) :: app,appidur
character(len=2) :: comp(3)
character(len=200) :: fname2
real, parameter :: pi = 3.14159265

write(6,*)'writing google earth kml file for plotting earthquake and receiver locations/seismograms...'
write(6,*)'Check it out: googleearth_src_rec_seis.kml'

slat=90.-srccolat1*180./pi
slon=srclon1*180./pi
if (slon>180.) slon=slon-360.

rlat=90.-rcvcolat*180./pi
rlon=rcvlon*180./pi
do i=1,num_rec_glob
   if (rlon(i)>180.) rlon(i)=rlon(i)-360.
enddo

rlat_rot=90.-rcvcolat_rot*180./pi
rlon_rot=rcvlon_rot*180./pi
do i=1,num_rec_glob
   if (rlon_rot(i)>180.) rlon_rot(i)=rlon_rot(i)-360.
enddo
open(unit=88,file='googleearth_src_rec_seis.kml')

write(88,14)'<?xml version="1.0" encoding="UTF-8"?> '
write(88,15)'<kml xmlns="http://earth.google.com/kml/2.0"> '
write(88,16)'<Document> '

write(88,*)
write(88,*)'  <name> earthquake-receiver configuration</name>'
write(88,*)'    <LookAt>'
write(88,12)'     <longitude>',slon,'</longitude><latitude>',slat,'</latitude>'
write(88,*)'     <range>7000000</range><tilt>0</tilt><heading>0</heading>'
write(88,*)'    </LookAt>'
write(88,*)
write(88,*)'......'
write(88,*)'  <Placemark>'
write(88,*)'     <Style id="earthquake">'
write(88,*)'       <IconStyle>'
 write(88,*)'       <scale>5</scale>'
write(88,*)'         <Icon>'
write(88,*)' <href>http://maps.google.com/mapfiles/kml/shapes/earthquake.png</href>'
write(88,*)'             </Icon>'
write(88,*)'           </IconStyle>'
write(88,*)'                  <LabelStyle>'
write(88,*)'                      <scale>5</scale>'
 write(88,*)'                 </LabelStyle>'
write(88,*)'        </Style>'
write(88,*)'    <name> earthquake</name>'
write(88,*) ' <description> Event details:'
write(88,20) ' colat,lon [deg]:',srccolat1*180./pi,srclon1*180./pi
write(88,21)' source depth [km]',srcdepth
write(88,23)'Mrr=',Mij(1)
write(88,23)'Mtt=',Mij(2)
write(88,23)'Mpp=',Mij(3)
write(88,23)'Mtr=',Mij(4)
write(88,23)'Mpr=',Mij(5)
write(88,23)'Mtp=',Mij(6)
write(88,21)'source period [s]:',per
write(88,*)'</description>'
write(88,13)'   <Point><coordinates>',slon,',',slat,'</coordinates></Point>'
write(88,*)'   </Placemark>'

do i=1,num_rec_glob
   write(88,*)
   write(88,*) ' <Placemark>'
   write(88,*) '     <Style id="cam">'
   write(88,*) '       <IconStyle>'
 write(88,*)'       <scale>2</scale>'
   write(88,*) '         <Icon>'
      write(88,*) '<href>http://maps.google.com/mapfiles/kml/shapes/camera.png</href>'
   write(88,*) '         </Icon>'
   write(88,*) '       </IconStyle>'
write(88,*)'                  <LabelStyle>'
write(88,*)'                      <scale>2</scale>'
write(88,*)'                 </LabelStyle>'
   write(88,*) '     </Style>'

   fname2=trim(receiver_name(i))//'_'//trim(interp_method)//'_filtperiod'//appidur//&
        '_['//reccomp(1)//reccomp(2)//reccomp(3)//'].ascii'

   write(88,17) ' <name> ',trim(fname2),'  # ',i,'</name>'
   call define_io_appendix(app,i)
   write(88,119) ' <description> rotated ',trim(fname2)
   write(88,20) ' colat,lon [deg]:',rcvcolat_rot(i)*180./pi,rcvlon_rot(i)*180./pi

   write(88,*) '  </description>'
   write(88,13) '   <Point><coordinates>',rlon_rot(i),',',rlat_rot(i),'</coordinates></Point>'
   write(88,*) ' </Placemark>'

   write(88,*) ' <Placemark>'
   write(88,*) '     <Style id="cam">'
   write(88,*) '       <IconStyle>'
 write(88,*)'       <scale>2</scale>'
   write(88,*) '         <Icon>'
      write(88,*)' <href>http://maps.google.com/mapfiles/kml/pushpin/ylw-pushpin.png</href>'
   write(88,*) '         </Icon>'
   write(88,*) '       </IconStyle>'
write(88,*)'                  <LabelStyle>'
write(88,*)'                      <scale>2</scale>'
write(88,*)'                 </LabelStyle>'
   write(88,*) '     </Style>'
   write(88,17) ' <name> ',trim(fname2),'  # ',i,'</name>'
   call define_io_appendix(app,i)
   write(88,119) ' <description> original ',trim(fname2)
   write(88,20) ' colat,lon [deg]:',rcvcolat(i)*180./pi,rcvlon(i)*180./pi
   write(88,*) '  </description>'
   write(88,13) '   <Point><coordinates>',rlon(i),',',rlat(i),'</coordinates></Point>'
   write(88,*) ' </Placemark>'

enddo

write(88,*)'......'
write(88,*)
write(88,*)'</Document>'
write(88,*)'</kml>'

close(88)

12 format(a16,f14.2,a23,f14.2,a12)
13 format(a23,f14.2,a1,f14.2,a23)
14 format(a39)
15 format(a46)
16 format(a11)
17 format(a7,a50,a10,i4,a7)
18 format(a36,a4,a14)
19 format(a24,a30)
119 format(a24,a50)
20 format(A18,f14.2,f14.2)
21 format(A18,f14.2)
23 format(a5,1pe14.2)

end subroutine save_google_earth_kml
!=============================================================================

!-----------------------------------------------------------------------------
subroutine define_io_appendix(app,iproc)

implicit none
integer :: iproc
character(len=4) :: app
character(len=1) :: milp,cenp,dizp,unip

  milp = char(48+    iproc/1000)
  cenp = char(48+mod(iproc/100,10))
  dizp = char(48+mod(iproc/10,10))
  unip = char(48+mod(iproc,10))
  app = milp//cenp//dizp//unip

end subroutine define_io_appendix
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
subroutine depths(srcdepth,model,chosen_dir)

real, intent(in)                :: srcdepth
character(len=34), intent(in)   :: model
character(len=200), intent(out) :: chosen_dir
integer, parameter              :: ndepths_prem = 1
integer, parameter              :: ndepths_iasp91 = 1
integer, parameter              :: ndepths_ak135 = 1
integer, parameter              :: ndepths = max(ndepths_prem,ndepths_ak135,ndepths_ak135)
character(len=200)              :: rundir(ndepths), rundirtemp
real, dimension(ndepths)        :: depth
integer                         :: dirind,ndepths_chosen

  select case(model)

  case ('prem') !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
    ndepths_chosen = ndepths_prem
    depth(1) = 100.; rundir(1) = '../PREM_40S_MIJ100KM_DIRAC_2000S_DB/ '
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

  case ('iasp91') !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
    ndepths_chosen = ndepths_iasp91
    depth(1) = 100.; rundir(1) = '../PREM_40S_MIJ100KM_DIRAC_2000S_DB/ '
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

  case ('ak135') !vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
    ndepths_chosen = ndepths_ak135
    depth(1) = 100.; rundir(1) = '../PREM_40S_MIJ100KM_DIRAC_2000S_DB/ '
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

  case default 
     write(6,*)'no simulation for model ',trim(model),' available'; stop
  end select

  write(6,*) 'Rundir?'
  read(*,*) rundirtemp
  dirind = minloc(abs(srcdepth-depth(1:ndepths_chosen)),1)
  chosen_dir = '../'//trim(rundirtemp)//'/'
  write(6,*)'Background model: ',trim(model)
  write(6,*)'Desired, offered depth [km]:',srcdepth,depth(dirind)
  write(6,*)'Directory name: ',trim(chosen_dir)

end subroutine depths
!--------------------------------------------------------------------------
