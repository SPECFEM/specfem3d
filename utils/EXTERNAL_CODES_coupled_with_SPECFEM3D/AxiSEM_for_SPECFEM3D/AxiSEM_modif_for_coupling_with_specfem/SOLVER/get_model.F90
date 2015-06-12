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

!========================
module get_model
!========================

  use global_parameters
  use data_mesh
  use data_io
  use data_proc

  use background_models
  use utlity

  implicit none

  public :: read_model
  private
  contains

!-----------------------------------------------------------------------------
!> First define array ieldom that specifically appoints a respective domain 
!! between discontinuities for each element (to avoid issues very close to 
!! discontinuities).
!! Fill each GLL(GLJ) grid point with rho,lambda,mu on domain basis for given 
!! background model as defined in background_models.f90 (NOTE: that routine 
!! is shared between solver and mesher, do not edit for one without the other!).
!! Compute global min/max values.
!!
!! Runs various checks and outputs related to the background model:
!!
!! 1) Discontinuities (given by discont(ndisc) from the mesher):
!!    Count points to check on conformity with coarsening and total number of 
!!    lateral GLL points, check sphericity and lateral elastic homogeneity, 
!!    output velocities just above/below
!! 2) Elastic background model:
!!    check whether all GLL points assumed reasonable values, whether all 
!!    lateral neighbor points share the same elastic properties and 
!!    compute maximal variation across an element (issue warning if large)
!! 3) Resolution test:
!!    Define radial trigonometric functions with wavelengths corresponding 
!!    to the source period, or half, or double of it.
!!    Compute accuracy of same for each element, using max/min of vp/vs 
!!    elementally and globally.
!!    This test gives a somewhat reasonable handle on how well the desired 
!!    source frequencies are resolved at worst, globally and locally.
!!    (of course only spatially: time discretization is usually worse such that 
!!     these values should not be taken absolutely, but rather as a relative 
!!     heuristic feel for the grid-model melange).
!!    Also outputs various forms of this test to files (called resolutionsine).
!! 4) Numerical resolution: 
!!    computes the characteristic lead time (grid spacing)/velocity which 
!!    is the determining factor for the numerical performance, stability 
!!    and mesh efficiency. Outputs various descendants such as time step and 
!!    source period as a function of radius, and separately writes min/max 
!!    characteristic lead times for P- and S-wave velocities of all elements 
!!    per processor, and for critical subdomains such as the central cube 
!!    and coarsening layers. 
!!    On-the-fly verification of respective radial averages:
!!       xmgrace timestep_rad.dat or xmgrace period_rad.dat
!-----------------------------------------------------------------------------
subroutine read_model(rho, lambda, mu, xi_ani, phi_ani, eta_ani, &
                          fa_ani_theta, fa_ani_phi, Q_mu, Q_kappa)

  use commun, ONLY : barrier
  use lateral_heterogeneities
  use data_source, ONLY : rot_src
  use data_mesh, only: npol, nelem, nel_solid, ielsolid
  use nc_routines, only: nc_dump_elastic_parameters

  real(kind=dp), dimension(0:npol,0:npol,nelem), intent(out) :: rho
  real(kind=dp), dimension(0:npol,0:npol,nelem), intent(out) :: lambda, mu
  real(kind=dp), dimension(0:npol,0:npol,nelem), intent(out) :: xi_ani, phi_ani, eta_ani
  real(kind=dp), dimension(0:npol,0:npol,nelem), intent(out) :: fa_ani_theta, fa_ani_phi

  real(kind=realkind), dimension(nel_solid), intent(out), optional :: Q_mu, Q_kappa

  real(kind=dp)    :: s,z,r,theta,r1
  real(kind=dp)    :: vphtmp, vpvtmp, vshtmp, vsvtmp
  integer :: iel,ipol,jpol,iidom,ieldom(nelem),domcount(ndisc),iel_count
  logical :: foundit
  character(len=100) :: modelstring

  if (make_homo ) then 
       write(6,*)'  '
       write(6,*)'ERROR: homogeneous AND anisotropic model does not make '
       write(6,*)'       sense, check input file'
       write(6,*)'  '
       stop
  endif

  ! Set elastic parameters to crazy values to later check if all have been filled
  rho(0:npol,0:npol,1:nelem) = -1.E30
  lambda(0:npol,0:npol,1:nelem) = -1.E30
  mu(0:npol,0:npol,1:nelem) = -1.E30
  xi_ani(0:npol,0:npol,1:nelem) = -1.E30
  phi_ani(0:npol,0:npol,1:nelem) = -1.E30
  eta_ani(0:npol,0:npol,1:nelem) = -1.E30

  ! check for each element which domain it belongs to
  if (diagfiles) open(unit=65,file=infopath(1:lfinfo)//'/elems_bkgrdmodel_domain.dat'&
                                   //appmynum)
  
  if (do_mesh_tests) open(60000+mynum,file='Data/model_r_th_rho_vp_vs_'&
                                          //appmynum//'.dat')
  
  do iel=1,nelem
      foundit=.false.
      r1 = rcoord(int(npol/2),int(npol/2),iel) 
      do iidom=1,ndisc-1
          if (r1<discont(iidom) .and. r1> discont(iidom+1)) then
              ieldom(iel) = iidom
              foundit = .true.
              if (diagfiles) write(65,10) iel, r1, iidom, discont(iidom), discont(iidom+1)
          endif
      enddo
      if (r1 < discont(ndisc)) then
          ieldom(iel) = ndisc
          foundit = .true.
          if (diagfiles) write(65,10) iel, r1, ndisc, discont(ndisc)
      endif
      if (.not. foundit) then 
          write(6,*)'Processor', mynum, ' reports problem:'
          write(6,*)'havent found domain for element', iel
          write(6,*)'...of radius', r1
          write(6,*)'Disconts are located at:'
          write(6,11) discont(1:ndisc)
          stop
      endif
  enddo
  if (diagfiles) close(65)

10 format(i9,1pe15.7,i3,2(1pe11.3))
11 format(e11.3)
  if (do_mesh_tests) then
      if (lpr .and. verbose > 1) write(6,*)'    checking discontinuity discretization...' 
      call check_mesh_discontinuities(ieldom,domcount)
  endif

  if (lpr .and. verbose > 1) write(6,*)'   filling mesh with elastic properties...'   

  modelstring = bkgrdmodel

  iel_count = 0
  !========================
  do iel=1, nelem
  !========================

     iidom = ieldom(iel)
     do ipol=0,npol      !ipol
        do jpol=0,npol  !jpol
           call compute_coordinates(s,z,r,theta,iel,ipol,jpol)

           vphtmp = velocity(r,'vph',iidom,modelstring,lfbkgrdmodel)
           vpvtmp = velocity(r,'vpv',iidom,modelstring,lfbkgrdmodel)
           vshtmp = velocity(r,'vsh',iidom,modelstring,lfbkgrdmodel)
           vsvtmp = velocity(r,'vsv',iidom,modelstring,lfbkgrdmodel)
           eta_ani(ipol,jpol,iel) = velocity(r,'eta',iidom,modelstring,lfbkgrdmodel)
           rho(ipol,jpol,iel) = velocity(r,'rho',iidom,modelstring,lfbkgrdmodel)
           
           lambda(ipol,jpol,iel) = rho(ipol,jpol,iel) * (vphtmp**2 - two*vshtmp**2)
           mu(ipol,jpol,iel) = rho(ipol,jpol,iel) * vshtmp**2

           if (vsvtmp > (smallval_sngl * vphtmp)) then
              xi_ani(ipol,jpol,iel) = vshtmp**2 / vsvtmp**2
           else
              xi_ani(ipol,jpol,iel) = one
           endif
           
           phi_ani(ipol,jpol,iel) = vpvtmp**2 / vphtmp**2
            
           ! radial anisotropy (otherwise lateral heterogeneity!)
           fa_ani_theta(ipol,jpol,iel) = thetacoord(ipol, jpol, iel)
           fa_ani_phi(ipol,jpol,iel) = 0.


           ! @TODO generalize tests for anisotropic parameters:
           ! test round-off errors for elastic parameters & velocities
           if ( .not. dblreldiff_small(vphtmp,dsqrt( (lambda(ipol,jpol,iel) + &
                two*mu(ipol,jpol,iel) )/rho(ipol,jpol,iel) ) ) ) then 
              write(6,*)
              write(6,*)procstrg,&
                   'PROBLEM: Elastic params and p-vel conversion erroneous!'
              write(6,*)procstrg,'r,theta,vphtmp:',&
                   r/1000.,theta*180./pi,vphtmp/1000.
              write(6,*)procstrg,'rho,lambda,mu:',rho(ipol,jpol,iel)/1000., & 
                   lambda(ipol,jpol,iel)/1000.,mu(ipol,jpol,iel)/1000.
              stop
           endif
           
           if ( .not. dblreldiff_small(vshtmp,dsqrt( mu(ipol,jpol,iel)/ &
                rho(ipol,jpol,iel) ))) then 
              write(6,*)
              write(6,*)procstrg,&
                   'PROBLEM: Elastic params and s-vel conversion erroneous!'
              write(6,*)procstrg,'r,theta,vshtmp:',&
                   r/1000.,theta*180./pi,vshtmp/1000.
              write(6,*)procstrg,'rho,lambda,mu:',rho(ipol,jpol,iel)/1000., & 
                   lambda(ipol,jpol,iel)/1000.,mu(ipol,jpol,iel)/1000.
              stop
           endif
        enddo
     enddo

     ! @TODO generalize tests for anisotropic parameters:
     ! write out for later snaps
     if (do_mesh_tests) then
        do ipol=ibeg, iend
           do jpol=ibeg, iend
              call compute_coordinates(s,z,r,theta,iel,ipol,jpol)
              write(60000+mynum,14)r,theta,rho(ipol,jpol,iel),&
                 sqrt( (lambda(ipol,jpol,iel)+2.*mu(ipol,jpol,iel))/rho(ipol,jpol,iel)), &
                 sqrt( mu(ipol,jpol,iel)/rho(ipol,jpol,iel))
           enddo
        enddo
     endif

  !========================
  enddo ! nelem
  !========================

  ! Fill up Q arrays with values from the backgroundmodel
  if (anel_true) then
      if (lpr .and. verbose > 1) print *, '...filling up Q arrays'
      do iel=1, nel_solid
          iidom = ieldom(ielsolid(iel))
          call compute_coordinates(s, z, r, theta, ielsolid(iel), npol/2 - 1, npol/2 - 1)
          Q_mu(iel) = velocity(r, 'Qmu', iidom, modelstring, lfbkgrdmodel)
          Q_kappa(iel) = velocity(r, 'Qka', iidom, modelstring, lfbkgrdmodel)
      enddo
  endif

  if (lpr .and. verbose > 1) write(6,*) 'done with big mesh loop to define model'
  if (do_mesh_tests) close(60000+mynum)

14 format(5(1pe13.4))

  if (add_hetero) call compute_heterogeneities(rho, lambda, mu, xi_ani, phi_ani, &
                                                 eta_ani, fa_ani_theta, fa_ani_phi, ieldom)

  if (diagfiles) then
      ! plot final velocity model in vtk
      if (lpr .and. verbose > 1) write(6,*) 'plotting vtks for the model properties....'

      if (anel_true) then
         call plot_model_vtk(rho, lambda, mu, xi_ani, phi_ani, eta_ani, fa_ani_theta, &
                             fa_ani_phi, Q_mu, Q_kappa)
      else
         call plot_model_vtk(rho, lambda, mu, xi_ani, phi_ani, eta_ani, fa_ani_theta, &
                                 fa_ani_phi)
      endif
  end if

  if (use_netcdf) then
      if (anel_true) then
          call nc_dump_elastic_parameters(rho, lambda, mu, xi_ani, phi_ani, eta_ani, &
                                          fa_ani_theta, fa_ani_phi, Q_mu, Q_kappa)
      else
          call nc_dump_elastic_parameters(rho, lambda, mu, xi_ani, phi_ani, eta_ani, &
                                          fa_ani_theta, fa_ani_phi)
      endif
  end if


  ! Some tests....
  if (do_mesh_tests) then
     call barrier

     if (.not. add_hetero) then
         if (lpr) write(6,*) &
            '    checking elastic properties on discontinuities...' 
         call check_elastic_discontinuities(ieldom,domcount,lambda,mu,rho)
     else
         if (lpr) write(6,*) &
            '    NOT checking elastic properties on discontinuities since we added hetero...'   
     endif

     if (lpr) write(6,*) '    checking the background model discretization...' 
     call check_background_model(lambda,mu,rho)

     if (lpr) write(6,*) '    testing the mesh/background model resolution...' 
     call test_mesh_model_resolution(lambda,mu,rho)
  endif 

  ! MvD: Do we need the following for each anisotropic velocity can we leave it as is? Do we need it at all? 
  !  - 2.2.12: just outputted to standard out in parameters.f90

  ! compute min/max velocities in whole domain
  vpmin    = minval(dsqrt((lambda+two*mu)/rho))
  vpminloc = minloc(dsqrt((lambda+two*mu)/rho))
  vpmax    = maxval(dsqrt((lambda+two*mu)/rho))
  vpmaxloc = maxloc(dsqrt((lambda+two*mu)/rho))
  
  vsmin    = minval(dsqrt(mu/rho))
  vsminloc = minloc(dsqrt(mu/rho))
  vsmax    = maxval(dsqrt(mu/rho))
  vsmaxloc = maxloc(dsqrt(mu/rho))
    
  ! since minloc/maxloc start counting at 1...
  vpminloc(1) = vpminloc(1)-1; vpminloc(2) = vpminloc(2) - 1;
  vsminloc(1) = vsminloc(1)-1; vsminloc(2) = vsminloc(2) - 1;
  vpmaxloc(1) = vpmaxloc(1)-1; vpmaxloc(2) = vpmaxloc(2) - 1;
  vsmaxloc(1) = vsmaxloc(1)-1; vsmaxloc(2) = vsmaxloc(2) - 1;
  
  call compute_coordinates(s, z, vpminr, theta, vpminloc(3),&
                           vpminloc(1), vpminloc(2))
  call compute_coordinates(s, z, vsminr, theta, vsminloc(3),&
                           vsminloc(1), vsminloc(2))
  call compute_coordinates(s, z, vpmaxr, theta, vpmaxloc(3),&
                           vpmaxloc(1), vpmaxloc(2))
  call compute_coordinates(s, z, vsmaxr, theta, vsmaxloc(3),&
                           vsmaxloc(1), vsmaxloc(2))
end subroutine read_model
!=============================================================================

!-----------------------------------------------------------------------------
!> file-based, step-wise model in terms of domains separated by disconts.
!! format:
!! ndisc
!! r vp vs rho
!subroutine arbitr_sub_solar_arr(s,z,v_p,v_s,rho,bkgrdmodel2)
!
!  use data_mesh
!  real(kind=dp)   , intent(in) :: s(0:npol,0:npol,1:nelem),z(0:npol,0:npol,1:nelem)
!  character(len=100), intent(in) :: bkgrdmodel2
!  real(kind=dp)   , dimension(:,:,:), intent(out) :: rho(0:npol,0:npol,1:nelem)
!  real(kind=dp)   , dimension(:,:,:), intent(out) :: v_s(0:npol,0:npol,1:nelem)
!  real(kind=dp)   , dimension(:,:,:), intent(out) :: v_p(0:npol,0:npol,1:nelem)
!  real(kind=dp)   , allocatable, dimension(:) :: disconttmp,rhotmp,vstmp,vptmp
!  integer :: ndisctmp,i,ind(2),ipol,jpol,iel
!  logical :: bkgrdmodelfile_exists
!  real(kind=dp)    :: w(2),wsum,r0
!
!  ! Does the file bkgrdmodel".bm" exist?
!  !@TODO: Change to new name convention scheme. Should start in the MESHER.
!  inquire(file=bkgrdmodel2(1:index(bkgrdmodel2,' ')-1)//'.bm', &
!          exist=bkgrdmodelfile_exists)
!  if (bkgrdmodelfile_exists) then
!      open(unit=77,file=bkgrdmodel2(1:index(bkgrdmodel2,' ')-1)//'.bm')
!      read(77,*)ndisctmp
!      allocate(disconttmp(1:ndisctmp))
!      allocate(vptmp(1:ndisctmp),vstmp(1:ndisctmp),rhotmp(1:ndisctmp))
!      do i=1, ndisctmp
!          read(77,*)disconttmp(i),rhotmp(i),vptmp(i),vstmp(i)
!      enddo
!      close(77)
!      do iel=1,nelem
!          do jpol=0,npol
!              do ipol=0,npol
!                  r0 = dsqrt(s(ipol,jpol,iel)**2 +z(ipol,jpol,iel)**2 )
!                  call interp_vel(r0,disconttmp(1:ndisctmp),ndisctmp,ind,w,wsum)
!                  rho(ipol,jpol,iel)=sum(w*rhotmp(ind))*wsum
!                  v_p(ipol,jpol,iel)=(w(1)*vptmp(ind(1))+w(2)*vptmp(ind(2)))*wsum
!                  v_s(ipol,jpol,iel)=sum(w*vstmp(ind))*wsum
!              enddo
!          enddo
!      enddo
!      deallocate(disconttmp,vstmp,vptmp,rhotmp)
!  else 
!      write(6,*)'Background model file', &
!                trim(bkgrdmodel2)//'.bm','does not exist!!!'
!      stop
!  endif
!
!end subroutine arbitr_sub_solar_arr
!!=============================================================================
!
!!-----------------------------------------------------------------------------
!!> Calculate interpolation parameters w to interpolate velocity at radius r0
!!! from a model defined at positions r(1:n)
!subroutine interp_vel(r0,r,n,ind,w,wsum)
!
!  integer, intent(in)           :: n      !< number of supporting points
!  real(kind=dp)   , intent(in)  :: r(1:n) !< supporting points in depth
!  real(kind=dp)   , intent(in)  :: r0     !< Target depth
!  integer, intent(out)          :: ind(2) !< Indizes of supporting points 
!                                          !! between which r0 is found
!  real(kind=dp)   , intent(out) :: w(2),wsum !< Weighting factors
!  integer                       :: i,p
!  real(kind=dp)                 :: dr1,dr2
!
!  p = 1
!
!  i = minloc(dabs(r-r0),1)
!
!  if (r0>0.d0) then
!     if ((r(i)-r0)/r0> 1.d-8) then ! closest discont. at larger radius
!        ind(1)=i
!        ind(2)=i+1
!        dr1=r(ind(1))-r0
!        dr2=r0-r(ind(2))
!     elseif ((r0-r(i))/r0> 1.d-8) then  ! closest discont. at smaller radius
!        if (r0>maxval(r)) then ! for round-off errors where mesh is above surface
!           ind(1)=i
!           ind(2)=i
!           dr1=1.d0
!           dr2=1.d0
!        else
!           ind(1)=i-1
!           ind(2)=i
!           dr1=r(ind(1))-r0
!           dr2=r0-r(ind(2))
!        endif
!     elseif (dabs((r(i)-r0)/r0)< 1.d-8) then ! closest discont identical
!        ind(1)=i
!        ind(2)=i
!        dr1=1.d0
!        dr2=1.d0
!     else
!        write(6,*)'problem with round-off errors in interpolating......'
!        write(6,*)'r0,r(i),i',r0,r(i),abs((r(i)-r0)/r0),i
!        stop
!     endif
!  else !r0=0
!     if (r(i)==0.d0) then ! center of the sun
!        ind(1)=i
!        ind(2)=i
!        dr1=1.d0
!        dr2=1.d0
!     else
!        ind(1)=i
!        ind(2)=i+1
!        dr1=r(ind(1))-r0
!        dr2=r0-r(ind(2))        
!     endif
!  endif
!
!  ! inverse distance weighting
!  w(1) = (dr1)**(-p)
!  w(2) = (dr2)**(-p)
!  wsum = 1.d0 / sum(w)
!
!end subroutine interp_vel
!=============================================================================

!-----------------------------------------------------------------------------
subroutine check_mesh_discontinuities(ieldom,domcount)

  use data_mesh, only: ndisc, npol, nelem

  integer, intent(in)  :: ieldom(:)
  integer, intent(out) :: domcount(ndisc)
  integer              :: iel,ipol,jpol,iidom,eldomcount
  real(kind=dp)        :: s,z,r,theta

  domcount(:) = 0
  ! Count points on discontinuities and check whether right domains are assigned
  do iel=1,nelem
      iidom=ieldom(iel)
      do jpol=0,npol
          eldomcount=-1
          do ipol=0,npol
              call compute_coordinates(s,z,r,theta,iel,ipol,jpol)
              if (r>zero) then 
                  if ( (r-discont(iidom))/r > smallval_dble ) then 
                      write(6,*)procstrg,'PROBLEM with domains and discontinuities!'
                      write(6,*)procstrg,'radius > associated discont.:',&
                           iidom,r,discont(iidom)
                      stop
                  endif
              endif

              if ( dblreldiff_small(r,discont(iidom)) ) then
                 domcount(iidom)=domcount(iidom)+1
                 eldomcount=eldomcount+1
              endif
          enddo
         
          ! Check number of GLL points on discontinuities
          if (dblreldiff_small(r,discont(iidom)) .and. (eldomcount /= npol)) then
              write(6,*) procstrg,&
                         'PROBLEM: not every GLL along xi is on the discont!'
              write(6,*) procstrg, 'iel,idom,r,theta:', iel, ieldom(iel), &
                         r/1000.d0, theta*180.d0/pi
              write(6,*) procstrg, 'ipol count, npol:', eldomcount, npol
              stop
          endif
         
      enddo
  enddo
  
  ! make sure on each discontinuity we have a fraction of # surface points
  ! (in agreement with conformal coarsening/size doubling)

  do iidom=2, ndisc
     if (mod(domcount(1),domcount(iidom)) /= 0 ) then
        write(6,*) ' ' 
        write(6,*) procstrg,&
                   'PR0BLEM: # points on discontinuity not fraction of surface '
        write(6,*) procstrg, 'Domain, discontinuity:', iidom, discont(iidom)
        write(6,*) procstrg, 'number of points at discont, surface:', &
                   domcount(iidom), domcount(1)
        stop
     endif
  enddo
  
  ! same as above, but make sure each level is fraction of its above neighbor
  do iidom=2,ndisc
     if (domcount(iidom-1) /= domcount(iidom)) then
        if (mod(domcount(iidom-1),2*domcount(iidom)) /= 0 ) then
           write(6,*) ' '
           write(6,*) procstrg, 'PR0BLEM: # points discont. not even fraction', &
                      ' of discont. above'
           write(6,*) procstrg, 'Domain, discontinuity:', iidom, discont(iidom)
           write(6,*) procstrg, '#points at discont above,here:', &
                      domcount(iidom-1), domcount(iidom)
           stop
        endif
     endif
  enddo

end subroutine check_mesh_discontinuities
!=============================================================================

!-----------------------------------------------------------------------------
subroutine check_elastic_discontinuities(ieldom,domcount,lambda,mu,rho)

  use data_mesh, only : ndisc, nelem

  integer, intent(in)          :: ieldom(nelem),domcount(ndisc)
  real(kind=dp)   , intent(in) :: rho(0:npol,0:npol,nelem)
  real(kind=dp)   , intent(in) :: lambda(0:npol,0:npol,nelem)
  real(kind=dp)   , intent(in) :: mu(0:npol,0:npol,nelem)
  integer                      :: iel,jpol,iidom
  real(kind=dp)                :: minvpdomabove(ndisc),maxvpdomabove(ndisc)
  real(kind=dp)                :: minvsdomabove(ndisc),maxvsdomabove(ndisc)
  real(kind=dp)                :: minrodomabove(ndisc),maxrodomabove(ndisc)
  real(kind=dp)                :: minvpdombelow(ndisc),maxvpdombelow(ndisc)
  real(kind=dp)                :: minvsdombelow(ndisc),maxvsdombelow(ndisc)
  real(kind=dp)                :: minrodombelow(ndisc),maxrodombelow(ndisc)
  real(kind=dp)                :: s,z,r,theta

  ! count points on discontinuities and check whether right domains are assigned

  open(unit=6464, file=infopath(1:lfinfo)//&
       '/velocities_below_discont.dat'//appmynum)
  open(unit=6465, file=infopath(1:lfinfo)//&
       '/velocities_above_discont.dat'//appmynum)
  
  minvpdombelow = 1.d6; maxvpdombelow = 0.d0; minvpdomabove = 1.d6; maxvpdomabove = 0.d0
  minvsdombelow = 1.d6; maxvsdombelow = 0.d0; minvsdomabove = 1.d6; maxvsdomabove = 0.d0
  minrodombelow = 1.d6; maxrodombelow = 0.d0; minrodomabove = 1.d6; maxrodomabove = 0.d0
  
  do iel=1,nelem
     iidom=ieldom(iel)
     do jpol=0,npol
       
        ! write out velocities for each element that shares a discontinuity
        call compute_coordinates(s,z,r,theta,iel,int(npol/2),jpol)
        if ( dblreldiff_small(r,discont(iidom)) ) then 
         
           if ( minvpdombelow(iidom) > sqrt((lambda(int(npol/2),jpol,iel) &
                 + two*mu(int(npol/2),jpol,iel))/rho(int(npol/2),jpol,iel)) ) &
              minvpdombelow(iidom) = sqrt((lambda(int(npol/2),jpol,iel) &
                 + two*mu(int(npol/2),jpol,iel))/rho(int(npol/2),jpol,iel))
           
           if ( maxvpdombelow(iidom)< sqrt((lambda(int(npol/2),jpol,iel) &
                 + two*mu(int(npol/2),jpol,iel))/rho(int(npol/2),jpol,iel)) ) & 
              maxvpdombelow(iidom)=sqrt((lambda(int(npol/2),jpol,iel) &
                 + two*mu(int(npol/2),jpol,iel))/rho(int(npol/2),jpol,iel)) 
           
           if ( minvsdombelow(iidom) > & 
                 sqrt(mu(int(npol/2),jpol,iel)/rho(int(npol/2),jpol,iel)) ) & 
              minvsdombelow(iidom)=&
              sqrt(mu(int(npol/2),jpol,iel)/rho(int(npol/2),jpol,iel))
           
           if ( maxvsdombelow(iidom) < &
                 sqrt(mu(int(npol/2),jpol,iel)/rho(int(npol/2),jpol,iel)) ) & 
              maxvsdombelow(iidom) = &
                 sqrt(mu(int(npol/2),jpol,iel)/rho(int(npol/2),jpol,iel))
           
           if ( minrodombelow(iidom) > rho(int(npol/2),jpol,iel) ) & 
              minrodombelow(iidom) = rho(int(npol/2),jpol,iel)
           
           if ( maxrodombelow(iidom) < rho(int(npol/2),jpol,iel) ) & 
                maxrodombelow(iidom) = rho(int(npol/2),jpol,iel)
         
           ! r,theta,vp,vs,rho for one point of each element on the discontinuity
           write(6464,14) r/1000.d0, theta*180.d0/pi,                          &
                          sqrt((lambda(npol/2,jpol,iel)                        &
                          + 2 * mu(npol/2,jpol,iel)) / rho(npol/2,jpol,iel)),  &
                          sqrt(mu(npol/2,jpol,iel) / rho(npol/2,jpol,iel)),    &
                          rho(npol/2,jpol,iel)
        endif
        
        if (iidom < ndisc) then 
           if ( dblreldiff_small(r,discont(iidom+1)) ) then 
          
              if ( minvpdomabove(iidom+1) > sqrt((lambda(int(npol/2),jpol,iel) &
                    + two*mu(int(npol/2),jpol,iel))/rho(int(npol/2),jpol,iel)) ) &
                 minvpdomabove(iidom+1) = sqrt((lambda(int(npol/2),jpol,iel) &
                    + two*mu(int(npol/2),jpol,iel))/rho(int(npol/2),jpol,iel))
              
              if ( maxvpdomabove(iidom+1) < sqrt((lambda(int(npol/2),jpol,iel) &
                    + two*mu(int(npol/2),jpol,iel))/rho(int(npol/2),jpol,iel)) ) & 
                 maxvpdomabove(iidom+1) = sqrt((lambda(int(npol/2),jpol,iel) &
                    + two*mu(int(npol/2),jpol,iel))/rho(int(npol/2),jpol,iel)) 
              
              if ( minvsdomabove(iidom+1) > & 
                    sqrt(mu(int(npol/2),jpol,iel)/rho(int(npol/2),jpol,iel)) ) & 
                 minvsdomabove(iidom+1) = &
                    sqrt(mu(int(npol/2),jpol,iel)/rho(int(npol/2),jpol,iel))
              
              if ( maxvsdomabove(iidom+1) < &
                    sqrt(mu(int(npol/2),jpol,iel)/rho(int(npol/2),jpol,iel)) ) &
                 maxvsdomabove(iidom+1) = &
                    sqrt(mu(int(npol/2),jpol,iel)/rho(int(npol/2),jpol,iel))
              
              if ( minrodomabove(iidom+1) > rho(int(npol/2),jpol,iel) ) & 
                 minrodomabove(iidom+1) = rho(int(npol/2),jpol,iel)
              
              if ( maxrodomabove(iidom+1)<rho(int(npol/2),jpol,iel) ) & 
                   maxrodomabove(iidom+1)=rho(int(npol/2),jpol,iel)
              
              ! r,theta,vp,vs,rho for one point of each element on the discontinuity
              write(6464,14) r/1000.d0, theta*180.d0/pi,                          &
                             sqrt((lambda(npol/2,jpol,iel)                        &
                             + 2 * mu(npol/2,jpol,iel)) / rho(npol/2,jpol,iel)),  &
                             sqrt(mu(npol/2,jpol,iel) / rho(npol/2,jpol,iel)),    &
                             rho(npol/2,jpol,iel)
           endif
        endif
     enddo !jpol
  enddo ! ielem
  close(6464); close(6465)
14  format(1pe14.5,1pe14.5,3(1pe14.6))
  
  ! Affix "above" values on the surface to corresponding "below" values
  minvpdomabove(1) = minvpdombelow(1)
  minvsdomabove(1) = minvsdombelow(1)
  minrodomabove(1) = minrodombelow(1)

  maxvpdomabove(1) = maxvpdombelow(1)
  maxvsdomabove(1) = maxvsdombelow(1)
  maxrodomabove(1) = maxrodombelow(1)  

  ! check whether the right number of points resides on each discontinuity
  ! and compute the velocities above/below the discontinuities

  if (verbose > 1) then
     write(69,*)' '
     write(69,*)'# points and elastic info around discontinuities'
     write(69,*)
     write(69,133)'','r disc. [km]','# pts','vp [km/s]','vs [km/s]','rho [g/cm^3]'
     write(69,13) 'Below', discont(1)/1.d3, domcount(1), minvpdombelow(1)/1.d3,&
                   minvsdombelow(1)/1.d3, minrodombelow(1)/1.d3
  endif
133 format(a7,a13,a7,3(a13))
 
 
  if (.not. dblreldiff_small(minvpdombelow(1),maxvpdombelow(1))) then
     write(6,*)
     write(6,*) procstrg,&
                'PROBLEM: Vp at the surface is not spherically symmetric!'
     write(6,*) procstrg, 'Min/max vp along surface:',&
                minvpdombelow(1), maxvpdombelow(1)
     stop
  endif
 
  if (.not. dblreldiff_small(minvsdombelow(1),maxvsdombelow(1))) then 
     write(6,*)
     write(6,*) procstrg,&
               'PROBLEM: Vs at the surface is not spherically symmetric!'
     write(6,*) procstrg, 'Min/max vs along surface:',&
                minvsdombelow(1), maxvsdombelow(1)
     stop
  endif
 
  if (.not. dblreldiff_small(minrodombelow(1),maxrodombelow(1))) then 
     write(6,*)
     write(6,*) procstrg,&
                'PROBLEM: density at the surface is not spherically symmetric!'
     write(6,*) procstrg,&
                'Min/max rho along surface:',minrodombelow(1),maxrodombelow(1)
     stop
  endif
 
  do iidom=2,ndisc
     if (verbose > 1) write(69,13) 'Above', discont(iidom)/1.d3, domcount(iidom), &
                  minvpdomabove(iidom)/1.d3, minvsdomabove(iidom)/1.d3, &
                  minrodomabove(iidom)/1.d3  
     
     if (.not. dblreldiff_small(minvpdomabove(iidom),maxvpdomabove(iidom)))then
        write(6,*)
        write(6,*)procstrg,'PROBLEM: Vp not spherically symm. above discont',&
             discont(iidom)
        write(6,*)procstrg,'Min/max vp along discont:',minvpdomabove(iidom), &
             maxvpdomabove(iidom)
        stop
     endif
     
     if (.not. dblreldiff_small(minvsdomabove(iidom),maxvsdomabove(iidom)))then
        write(6,*)
        write(6,*)procstrg,'PROBLEM: Vs not spherically symm. above discont',&
             discont(iidom)
        write(6,*)procstrg,'Min/max vs along discont:',minvsdomabove(iidom), &
             maxvsdomabove(iidom)
        stop 
     endif
     
     if (.not. dblreldiff_small(minrodomabove(iidom),maxrodomabove(iidom)))then
        write(6,*)
        write(6,*)procstrg,&
             'PROBLEM: density not spherically symm. above discont', &
             discont(iidom)
        write(6,*)procstrg,'Min/max rho along discont:',minrodomabove(iidom), &
             maxrodomabove(iidom)
        stop
     endif
     
     if (verbose > 1) write(69,13)  'Below', discont(iidom)/1.d3, domcount(iidom), &
                   minvpdombelow(iidom)/1.d3, minvsdombelow(iidom)/1.d3, &
                   minrodombelow(iidom)/1.d3
     
     if (.not. dblreldiff_small(minvpdombelow(iidom),maxvpdombelow(iidom)))then
        write(6,*)
        write(6,*)procstrg,'PROBLEM: Vp not spherically symm. below discont', &
             discont(iidom)
        write(6,*)procstrg,'Min/max vp along discont:',minvpdombelow(iidom), &
             maxvpdombelow(iidom)
        stop
     endif
     
     if (.not. dblreldiff_small(minvsdombelow(iidom),maxvsdombelow(iidom)))then
        write(6,*)
        write(6,*)procstrg,'PROBLEM: Vs not spherically symm. below discont', &
             discont(iidom)
        write(6,*)procstrg,'Min/max vs along discont:',minvsdombelow(iidom), &
             maxvsdombelow(iidom)
        stop
     endif
     
     if (.not.dblreldiff_small(minrodombelow(iidom),maxrodombelow(iidom))) then
        write(6,*)
        write(6,*)procstrg,&
             'PROBLEM: density not spherically symmetric below discont', &
             discont(iidom)
        write(6,*)procstrg,'Min/max rho along discont:',minrodombelow(iidom), &
             maxrodombelow(iidom)
        stop
     endif
  enddo !idom=2:ndisc
13 format(a7,1pe12.4,i7,3(1pe13.3))

end subroutine check_elastic_discontinuities
!=============================================================================

!-----------------------------------------------------------------------------
!> Do sanity checks on background model:
!! - no invalid (negative) values on lambda, mu, rho
!! - no lateral heterogeneities, unless it is a lathet model
!! - no insane variations of model parameters within one cell.
subroutine check_background_model(lambda,mu,rho)

  use data_mesh, ONLY : eltype, coarsing
  use data_mesh, only : npol, nelem

  real(kind=dp)   , intent(in)  :: rho(0:npol,0:npol,nelem) 
  real(kind=dp)   , intent(in)  :: lambda(0:npol,0:npol,nelem)
  real(kind=dp)   , intent(in)  :: mu(0:npol,0:npol,nelem)
  integer                       :: iel,ipol,jpol,i
  integer, dimension (2)        :: rhominloc,lamminloc,muminloc
  integer, dimension (2)        :: rhomaxloc,lammaxloc,mumaxloc
  real(kind=dp)   , allocatable :: maxdiffrho(:)
  real(kind=dp)   , allocatable :: maxdifflam(:),maxdiffmu(:)

  allocate(maxdiffrho(nelem),maxdifflam(nelem),maxdiffmu(nelem))

  ! check whether all grid points assumed non-crazy elastic properties
  if (minval(rho) < zero) then 
     write(6,*)
     write(6,*) procstrg, 'PROBLEM: Not every grid point has a valid density!'
     vsminloc = minloc(rho)
     write(6,*) procstrg, 'rho=',minval(rho)
     write(6,*) procstrg, 'r,theta:',&
                rcoord(vsminloc(1)-1,vsminloc(2)-1,vsminloc(3))/1000.,&
                thetacoord(vsminloc(1)-1,vsminloc(2)-1,vsminloc(3))*180./pi
     stop
  endif
  
  if (minval(mu) < zero) then 
     write(6,*)
     write(6,*) procstrg,&
          'PROBLEM: Not every grid point has a valid mu (Lame) parameter!'
     vsminloc = minloc(mu)
     write(6,*) procstrg, 'mu=',minval(mu)
     write(6,*) procstrg, 'r,theta:',&
                rcoord(vsminloc(1)-1,vsminloc(2)-1,vsminloc(3))/1000.,&
                thetacoord(vsminloc(1)-1,vsminloc(2)-1,vsminloc(3))*180./pi
     stop
  endif
  
  if (minval(lambda) < zero) then 
     write(6,*)
     write(6,*) procstrg,&
               'PROBLEM: Not every grid point has a valid lambda (Lame) parameter'
     vsminloc = minloc(lambda)
     write(6,*) procstrg, 'lambda=', minval(lambda)
     write(6,*) procstrg, 'r,theta:', &
                rcoord(vsminloc(1)-1,vsminloc(2)-1,vsminloc(3))/1000.,&
                thetacoord(vsminloc(1)-1,vsminloc(2)-1,vsminloc(3))*180./pi
     stop
  endif

  ! Test lateral homogeneity: Does each point have same elastic properties 
  ! as its neighbor in the xi-direction?
  do iel = 1, nelem
     if (.not. add_hetero) then 
        if ( eltype(iel)=='curved' .and. .not. coarsing(iel) ) then
           do jpol = 0, npol
              do ipol = 1, npol-1
      
                 if (.not. dblreldiff_small(rho(ipol,jpol,iel), &
                     rho(ipol-1,jpol,iel)) .or. &
                     .not. dblreldiff_small(rho(ipol,jpol,iel), &
                     rho(ipol+1,jpol,iel))) then
                    write(6,*)
                    write(6,*) procstrg,'PROBLEM: lateral density inhomogeneity!'
                    write(6,*) procstrg,' at r[km], theta[deg]:',&
                               rcoord(ipol,jpol,iel)/1000., &
                               thetacoord(ipol,jpol,iel)*180/pi
                    write(6,*) procstrg,'Lateral ipol,theta,density:'
                    do i=0,npol
                        write(6,*) procstrg, i, thetacoord(i,jpol,iel)*180./pi,&
                                   rho(i,jpol,iel)
                    enddo
                    stop
                 endif
                 
                 if (.not. dblreldiff_small(lambda(ipol,jpol,iel), &
                     lambda(ipol-1,jpol,iel)) .or. &
                     .not. dblreldiff_small(lambda(ipol,jpol,iel), &
                     lambda(ipol+1,jpol,iel)) ) then 
                    write(6,*)
                    write(6,*) procstrg, 'PROBLEM: lateral elastic inhomogeneity!'
                    write(6,*) procstrg, 'at r[km], theta[deg]:',&
                               rcoord(ipol,jpol,iel)/1000., &
                               thetacoord(ipol,jpol,iel)*180/pi
                    write(6,*) procstrg, 'Lateral ipol,theta,lambda:'
                    do i=0,npol
                        write(6,*) procstrg, i, thetacoord(i,jpol,iel)*180./pi,&
                                   lambda(i,jpol,iel)
                    enddo
                    stop
                 endif
                 
                 if (.not. dblreldiff_small(mu(ipol,jpol,iel), &
                     mu(ipol-1,jpol,iel)) .or. &
                     .not. dblreldiff_small(mu(ipol,jpol,iel), &
                     mu(ipol+1,jpol,iel)) ) then 
                    write(6,*)
                    write(6,*) procstrg,'PROBLEM: lateral elastic inhomogeneity!'
                    write(6,*) procstrg,'at r[km], theta[deg]:',&
                               rcoord(ipol,jpol,iel)/1000., &
                               thetacoord(ipol,jpol,iel)*180/pi
                    write(6,*)'Lateral ipol,theta,mu:'
                    do i=0,npol
                        write(6,*) procstrg, i, thetacoord(i,jpol,iel)*180./pi,&
                                   mu(i,jpol,iel)
                    enddo
                    stop
                 endif
      
              enddo !ipol
           enddo !jpol
        endif ! only curved elements
     endif ! add_hetero
          
  
  
     ! Compute maximal elemental variation in medium properties
     rhominloc = minloc(rho(0:npol,0:npol,iel))
     lamminloc = minloc(lambda(0:npol,0:npol,iel))
     muminloc  = minloc(mu(0:npol,0:npol,iel))
     rhomaxloc = maxloc(rho(0:npol,0:npol,iel))
     lammaxloc = maxloc(lambda(0:npol,0:npol,iel))
     mumaxloc  = maxloc(mu(0:npol,0:npol,iel))
     
     maxdiffrho(iel) = dbleabsreldiff(rho(rhomaxloc(1)-1, rhomaxloc(2)-1, iel), &
                                      rho(rhominloc(1)-1, rhominloc(2)-1, iel))
     maxdifflam(iel) = dbleabsreldiff(lambda(lammaxloc(1)-1, lammaxloc(2)-1, iel), &
                                      lambda(lamminloc(1)-1, lamminloc(2)-1, iel))
     maxdiffmu(iel)  = dbleabsreldiff(mu(mumaxloc(1)-1, mumaxloc(2)-1, iel), & 
                                     mu(muminloc(1)-1, muminloc(2)-1, iel))
  enddo
  
  ! Maximal elemental variation in medium properties
  open(unit=6464,file=infopath(1:lfinfo)//&
       '/model_variation_across_elem.dat'//appmynum)
  open(unit=4,file=infopath(1:lfinfo)//&
       '/model_var_morethan_12percent.dat'//appmynum)
  do iel = 1,nelem
     if (maxdiffrho(iel) >= 0.12) then
        write(4,*)
        write(4,*)'WARNING: Density varies by >= 12 percent across elem: r=',&
             rcoord(npol/2,npol/2,iel)/1000.
     endif
     
     if (maxdifflam(iel) >= 0.12) then 
        write(4,*)
        write(4,*)'WARNING: Lambda varies by >= 12 percent across elem: r=',&
             rcoord(npol/2,npol/2,iel)/1000.
     endif
     
     if (maxdiffmu(iel) >= 0.12) then 
        write(4,*)
        write(4,*)'WARNING: Mu varies by >= 12 percent across elem: r=',&
             rcoord(npol/2,npol/2,iel)/1000.
     endif
     
     write(6464,64) rcoord(npol/2,npol/2,iel) / 1000., &
                    thetacoord(npol/2,npol/2,iel) * 180 / pi, &
                    maxdiffrho(iel), maxdifflam(iel), maxdiffmu(iel)
  enddo
  close(4)
  close(6464)
  
64 format(1pe14.4,1pe13.3,3(1pe15.5))
    
  if (verbose > 1) then
     write(69,*)
     write(69,*) 'Maximal variation in density across one element:', &
                 maxval(maxdiffrho)
     write(69,*) 'at r [km], theta[deg]:', &
                 rcoord(npol/2,npol/2,maxloc(maxdiffrho,1))/1000., &
                 thetacoord(npol/2,npol/2,maxloc(maxdiffrho,1))*180./pi
     write(69,*)
     write(69,*) 'Maximal variation in lambda across one element:',&
                 maxval(maxdifflam)
     write(69,*) 'at r [km], theta[deg]:', &
                 rcoord(npol/2,npol/2,maxloc(maxdifflam,1))/1000., &
                 thetacoord(npol/2,npol/2,maxloc(maxdifflam,1))*180./pi
     write(69,*)
     write(69,*) 'Maximal variation in mu across one element:',maxval(maxdiffmu)
     write(69,*) 'at r [km], theta[deg]:', &
                 rcoord(npol/2,npol/2,maxloc(maxdiffmu,1))/1000., &
                 thetacoord(npol/2,npol/2,maxloc(maxdiffmu,1))*180./pi
  endif

  deallocate(maxdiffrho, maxdifflam, maxdiffmu)

end subroutine check_background_model
!=============================================================================

!-----------------------------------------------------------------------------
!> Resolution test: Insert radial sine function with wavenumbers according to 
!! actual wavelengths, compute numerical and analytical integrals for various 
!! setups by changing the wavenumber by source period and seismic velocities.
subroutine test_mesh_model_resolution(lambda,mu,rho)

  use def_grid,  ONLY : massmatrix,massmatrix_dble
  use data_time, ONLY : period
  use data_mesh, ONLY : eltype, coarsing, north, axis
  use commun, ONLY : psum_dble,broadcast_int,broadcast_dble


  real(kind=dp)   , intent(in)  :: rho(0:npol,0:npol,nelem)
  real(kind=dp)   , intent(in)  :: lambda(0:npol,0:npol,nelem)
  real(kind=dp)   , intent(in)  :: mu(0:npol,0:npol,nelem)
  real(kind=dp)   , allocatable :: fp(:,:,:),fs(:,:,:)
  real(kind=dp)   , allocatable :: mass(:,:,:)
  real(kind=dp)                 :: intp_num,ints_num
  real(kind=dp)                 :: intp_ana,ints_ana,arg
  integer                       :: iel,ipol,jpol,i,irad,iirad,iidom
  real(kind=dp)                 :: tmpradii(naxel,2)
  real(kind=dp)                 :: vptmp,vstmp
  real(kind=dp)   , allocatable :: radii(:,:),vel(:,:)

  allocate(fp(0:npol,0:npol,nelem),fs(0:npol,0:npol,nelem))
  allocate(mass(0:npol,0:npol,nelem))

  do iidom = 1,2
     if (iidom == 1) then
        vptmp=minval(dsqrt( (lambda+2.d0*mu )/rho ))
        vstmp = vptmp/dsqrt(3.d0)
        if (verbose > 1) write(69,'(/,a,/,a)') &
                ':::::::::Resolution test for globally MINIMAL vp, vs::::::::::::::::::', &
                '     sin(2 pi/(T0*minv) r) integrated over 3-D volume'
     elseif (iidom ==2) then
        vptmp=maxval(dsqrt( (lambda+2.d0*mu )/rho ))
        vstmp = vptmp/dsqrt(3.d0)
        if (verbose > 1) write(69,'(/,a,/,a)') &
            ':::::::::Resolution test for globally MAXIMAL vp, vs::::::::::::::::::', &
            '     sin(2 pi/(T0*maxv) r) integrated over 3-D volume'
     endif

     do i = 1,4 ! 0.5*period, period, 1.5*period, 2*period
        do iel = 1, nelem
           do jpol = 0, npol
              do ipol = 0,npol
                 fp(ipol,jpol,iel) = dsin(two*pi*rcoord(ipol,jpol,iel)/ & 
                                    (dble(i)/2.d0*period*vptmp))
                 if (vstmp > zero) then 
                    fs(ipol,jpol,iel) = dsin(two*pi*rcoord(ipol,jpol,iel)/ & 
                                        (dble(i)/2.d0*period*vstmp))
                 else
                    fs(ipol,jpol,iel) = fp(ipol,jpol,iel)
                 endif
              enddo
           enddo
        enddo

        if (iidom ==3) fp = one

        !Numerical integration:
        call massmatrix_dble(mass,nelem,'total')
        intp_num = zero; ints_num = zero
        do iel = 1, nelem
           do ipol = 0, npol
              do jpol = 0, npol
                  intp_num = intp_num + fp(ipol,jpol,iel)*mass(ipol,jpol,iel)
                  ints_num = ints_num + fs(ipol,jpol,iel)*mass(ipol,jpol,iel)
              end do
           end do
        end do
        intp_num = 2.d0 * pi * intp_num
        ints_num = 2.d0 * pi * ints_num

        intp_num = psum_dble(intp_num)
        ints_num = psum_dble(ints_num)

        !Analytical integration: 
        !\int sin(ar) r^2 dr = 2r/a^2 sin(ar) - (r^2/a - 2/a^3) cos(ar)
        !  where a = 2 pi/(T0 v_{s,p})
        !Here: analytically for whole earth if v_{s,p} = const

        !Constant velocities:
        arg = two*pi/(real(i)/2.d0*period * vptmp)
        intp_ana = four*pi*( two*router/arg**2 * dsin(arg*router) - &
             (router**2/arg-two/arg**3)*dcos(arg*router) - two/arg**3 )
        arg = two*pi/(real(i)/2.*period * vstmp)
        ints_ana = four*pi*( two*router/arg**2 * dsin(arg*router) - &
             (router**2/arg-two/arg**3)*dcos(arg*router) - two/arg**3 )

        if (verbose > 1) then
           write(69,*) 'T0=period*', real(i)/2.
           write(69,1234) ' Num., ana., diff vp:', intp_num, intp_ana, &
                          dbleabsreldiff(intp_num,intp_ana)
           write(69,1234) ' Num., ana., diff vs:', ints_num, ints_ana, &
                          dbleabsreldiff(ints_num,ints_ana)
        endif
1234    format(a22,2(1pe18.8),1pe13.3)

     enddo ! 4 periods
  enddo ! min/max velocities

  ! Same but for constant function, i.e. the pure volume of the sphere
  if (verbose > 1) then
     write(69,*)
     write(69,*) ':::::::::Resolution test for constant function::::::::::::::::::::::::'
     write(69,*) '             3-D spherical volume itself'
  endif

  call massmatrix_dble(mass,nelem,'total')
  intp_num = zero
  do iel = 1, nelem
     do ipol = 0, npol
        do jpol = 0, npol
           intp_num = intp_num + mass(ipol,jpol,iel)
        end do
     end do
  end do
  intp_num = 2.d0 * pi * intp_num
  intp_num = psum_dble(intp_num)

  intp_ana = four / three * pi * router**3

  if (verbose > 1) write(69,1234) ' Num., ana., diff vp:',intp_num,intp_ana, &
                                  dbleabsreldiff(intp_num,intp_ana)

  ! Piecewise constant velocities over each element,
  ! only for spheroidal element shapes and leaving out coarsening levels

  irad = 0
  if (mynum==0) then
     do iel = 1, nelem
        if (eltype(iel)=='curved' .and. .not. coarsing(iel) ) then
           ! radii needed for analytical integration
           if (axis(iel) .and. north(iel)) then
              irad = irad + 1
              tmpradii(irad,1) = rcoord(0,npol,iel)
              tmpradii(irad,2) = rcoord(0,0,iel)
           endif
        endif !eltype curved
     enddo
  endif
  call broadcast_int(irad,0) 
  allocate(radii(irad,2),vel(irad,2))
  if (mynum==0) radii(1:irad,1:2) = tmpradii(1:irad,1:2)
  do iirad=1,irad
     call broadcast_dble(radii(iirad,1),0)
     call broadcast_dble(radii(iirad,2),0)
  enddo

  open(unit=779,file=infopath(1:lfinfo)//'/resolutionsine.dat'//appmynum)

  do iidom = 1,2
     if (iidom == 1) then
        open(unit=779,file=infopath(1:lfinfo)// &
             '/resolutionsine_elminvp.dat'//appmynum)      
        if (verbose > 1) write(69,'(/,a,/,a)')&
            ':::::::::Resolution test for elementally MINIMAL vp, vs:::::::::::::::', &
            '     sin(2 pi/(T0*minv) r) integrated over 3-D volume'
     elseif (iidom ==2) then
        open(unit=779,file=infopath(1:lfinfo)//'/resolutionsine_elmaxvp.dat' &
             //appmynum)  
        if (verbose > 1) write(69,'(/,a,/,a)') &
            ':::::::::Resolution test for elementally MAXIMAL vp, vs:::::::::::::::', &
            '     sin(2 pi/(T0*maxv) r) integrated over 3-D volume'
     endif

     fp = zero
     fs = zero
     iirad = 0
     do iel = 1, nelem
        if (eltype(iel)=='curved' .and. .not. coarsing(iel) ) then

           if (iidom == 1) then
              vptmp=minval(dsqrt((lambda(:,:,iel)+&
                                 two*mu(:,:,iel))/rho(:,:,iel)))
              vstmp= minval( dsqrt( mu(:,:,iel) / rho(:,:,iel) ) )
              if (vstmp==zero) vstmp = vptmp  
           elseif (iidom ==2) then
              vptmp=maxval(dsqrt((lambda(:,:,iel)+&
                                  two*mu(:,:,iel))/rho(:,:,iel)))
              vstmp= maxval( dsqrt( mu(:,:,iel) / rho(:,:,iel) )) 
              if (vstmp==zero) vstmp = vptmp  
           endif
           do jpol = 0, npol
              do ipol = 0,npol
                 fp(ipol,jpol,iel)= dsin(two*pi*rcoord(ipol,jpol,iel)/ &
                      (period*vptmp))
                 fs(ipol,jpol,iel)= dsin(two*pi*rcoord(ipol,jpol,iel)/ &
                      (period*vstmp))
              enddo
              if (axis(iel) .and. north(iel)) then
                 write(779,*) rcoord(0,jpol,iel),fp(0,jpol,iel),fs(0,jpol,iel)
              endif
           enddo
           ! radii needed for analytical integration
           if (axis(iel) .and. north(iel) .and. mynum==0) then
              iirad = iirad + 1
              vel(iirad,1) = vptmp
              vel(iirad,2) = vstmp
           endif

        endif !eltype curved
     enddo

     do iirad=1,irad
        call broadcast_dble(vel(iirad,1),0)
        call broadcast_dble(vel(iirad,2),0)
     enddo

     close(779)

     ! Numerical integration:
     call massmatrix_dble(mass,nelem,'total')
     intp_num = zero; ints_num = zero
     do iel = 1, nelem
        do ipol = 0, npol
           do jpol = 0, npol
              intp_num = intp_num + fp(ipol,jpol,iel)*mass(ipol,jpol,iel)
              ints_num = ints_num + fs(ipol,jpol,iel)*mass(ipol,jpol,iel)
           end do
        end do
     end do
     intp_num = 2.d0 * pi * intp_num
     ints_num = 2.d0 * pi * ints_num

     intp_num = psum_dble(intp_num)
     ints_num = psum_dble(ints_num)

     ! analytical integration
     intp_ana = zero
     ints_ana = zero; 
     do i=1,irad
        arg = two*pi/(period * vel(i,1))
        intp_ana = intp_ana + &
             four*pi*(two*radii(i,1)/arg**2*dsin(arg*radii(i,1)) - &
             (radii(i,1)**2/arg-two/arg**3)*dcos(arg*radii(i,1)) - &
             two*radii(i,2)/arg**2*dsin(arg*radii(i,2)) + &
             (radii(i,2)**2/arg-two/arg**3)*dcos(arg*radii(i,2)) )
        arg = two*pi/(period * vel(i,2))
        ints_ana = ints_ana + &
             four*pi*(two*radii(i,1)/arg**2*dsin(arg*radii(i,1)) - &
             (radii(i,1)**2/arg-two/arg**3)*dcos(arg*radii(i,1)) - &
             two*radii(i,2)/arg**2*dsin(arg*radii(i,2)) + &
             (radii(i,2)**2/arg-two/arg**3)*dcos(arg*radii(i,2)) )
     enddo
    
     if (verbose > 1) then
        write(69,1234)' Num., ana., diff vp:',intp_num,intp_ana, &
             dbleabsreldiff(intp_num,intp_ana)
        write(69,1234)' Num., ana., diff vs:',ints_num,ints_ana, &
             dbleabsreldiff(ints_num,ints_ana)
     endif

  enddo ! min/max velocities

  deallocate(fp,fs)
  deallocate(mass)

end subroutine test_mesh_model_resolution
!=============================================================================

!-----------------------------------------------------------------------------
!> Write scalar values into a binary VTK files
subroutine write_VTK_bin_scal(x,y,z,u1,elems,filename)
  implicit none
  integer                                   :: i,t
  !> Number of points
  integer, intent(in)                       :: elems
  !> Point coordinates
  real(sp), dimension(1:elems*4), intent(in) :: x,y,z 
  !> Values at coordinates x,y,z
  real(sp), dimension(1:elems*4), intent(in) :: u1 
  !> Filename of produced VTK file
  character (len=200), intent(in)           :: filename
  integer, dimension(1:elems*5)             :: cell
  integer, dimension(1:elems)               :: cell_type
  character (len=50)                        :: ss
  !points structure

  do i=5, elems*5, 5
      cell(i-4) = 4
  enddo

  t = 0

  do i=5,elems*5,5
      t = t + 4
      cell(i-3) = t - 4
      cell(i-2) = t - 3
      cell(i-1) = t - 2
      cell(i) = t - 1
  enddo

  cell_type = 9

#if defined(_CRAYFTN)
   open(110, file=trim(filename)//'.vtk', access='stream', status='replace', &
        form='unformatted')
#else
   open(110, file=trim(filename)//'.vtk', access='stream', status='replace', &
        form='unformatted', convert='big_endian')
#endif

  write(110) '# vtk DataFile Version 4.0'//char(10)
  write(110) 'mittico'//char(10)
  write(110) 'BINARY'//char(10)
  write(110) 'DATASET UNSTRUCTURED_GRID'//char(10)
  write(ss,fmt='(A6,I10,A5)') 'POINTS',elems*4,'float'
  write(110) ss//char(10)

  !points
  write(110) (x(i),y(i),z(i),i=1,elems*4)
  write(110) char(10)

  !cell topology
  write(ss,fmt='(A5,2I10)') 'CELLS',elems,elems*5
  write(110) char(10)//ss//char(10)
  write(110) cell
  write(110) char(10)

  !cell type
  write(ss,fmt='(A10,2I10)') 'CELL_TYPES',elems
  write(110) char(10)//ss//char(10)
  write(110) cell_type
  write(110) char(10)

  !data
  write(ss,fmt='(A10,I10)') 'POINT_DATA',elems*4
  write(110) char(10)//ss//char(10)
  write(110) 'SCALARS data float 1'//char(10)
  write(110) 'LOOKUP_TABLE default'//char(10) !color table?
  write(110) u1

  close(110)
  if (verbose > 1) write(6,*)'...saved ',trim(filename)//'.vtk'

end subroutine write_VTK_bin_scal
!=============================================================================

!-----------------------------------------------------------------------------
!> Plot the background model into vtk files. 
!! Difference to the ones generated in the MESHER: these include lateral
!! heterogeneities. Also, they are written per processor, hence easier to 
!! handle for high frequencies/many cpus
subroutine plot_model_vtk(rho, lambda, mu, xi_ani, phi_ani, eta_ani, &
                          fa_ani_theta, fa_ani_phi, Q_mu, Q_kappa)

  use data_mesh

  real(kind=dp), dimension(0:,0:,:), intent(in) :: rho 
  real(kind=dp), dimension(0:,0:,:), intent(in) :: lambda, mu
  real(kind=dp), dimension(0:,0:,:), intent(in), optional :: &
                                xi_ani, phi_ani, eta_ani
  real(kind=dp), dimension(0:,0:,:), intent(in), optional :: &
                                          fa_ani_theta, fa_ani_phi

  real(kind=realkind), dimension(nel_solid), intent(in), optional :: Q_mu, Q_kappa

  real, dimension(:), allocatable :: vp1,vs1,rho1
  real, dimension(:), allocatable :: vpv1,vsv1,eta1
  real, dimension(:), allocatable :: xi1,phi1
  real, dimension(:), allocatable :: fa_ani_theta1, fa_ani_phi1
  real, dimension(:), allocatable :: Q_mu1, Q_kappa1
  character(len=200) :: fname
  integer :: npts_vtk,ct,iel
  real, allocatable ::  x(:),y(:),z0(:)
  logical :: plot_ani, plot_anel

  npts_vtk = nelem * 4

  allocate(vp1(npts_vtk),vs1(npts_vtk),rho1(npts_vtk))
  if (present(xi_ani) .and. present(phi_ani) .and.  present(eta_ani) &
     .and. present(fa_ani_theta) .and. present(fa_ani_phi)) then
     allocate(vpv1(npts_vtk), vsv1(npts_vtk), eta1(npts_vtk))
     allocate(xi1(npts_vtk), phi1(npts_vtk))
     allocate(fa_ani_theta1(npts_vtk), fa_ani_phi1(npts_vtk))
     plot_ani = .true.
  else
     plot_ani = .false.
  endif

  if (present(Q_mu) .and. present(Q_kappa)) then
     allocate(Q_mu1(npts_vtk), Q_kappa1(npts_vtk))
     ! initializing to zero, as for the outer core has no attenuation
     Q_mu1 = 0.
     Q_kappa1 = 0.
     plot_anel = .true.
  else
     plot_anel = .false.
  endif

  allocate(x(npts_vtk), y(npts_vtk), z0(npts_vtk))

  z0 = 0.d0
  ct = 0
  do iel=1, nelem

     x(ct+1) = scoord(0,0,iel)
     x(ct+2) = scoord(npol,0,iel)
     x(ct+3) = scoord(npol,npol,iel)
     x(ct+4) = scoord(0,npol,iel)
     y(ct+1) = zcoord(0,0,iel)
     y(ct+2) = zcoord(npol,0,iel)
     y(ct+3) = zcoord(npol,npol,iel)
     y(ct+4) = zcoord(0,npol,iel)
 
 
     rho1(ct+1) = rho(0,0,iel)
     vp1(ct+1) = sqrt( (lambda(0,0,iel)+2.*mu(0,0,iel) ) / rho1(ct+1)  )
     vs1(ct+1) = sqrt( mu(0,0,iel)  / rho1(ct+1)  )
 
     if (plot_ani) then
        vpv1(ct+1) = sqrt(phi_ani(0,0,iel)) * vp1(ct+1)
        vsv1(ct+1) = vs1(ct+1) / sqrt(xi_ani(0,0,iel))
        xi1(ct+1) = xi_ani(0,0,iel)
        phi1(ct+1) = phi_ani(0,0,iel)
        eta1(ct+1) = eta_ani(0,0,iel)
        fa_ani_theta1(ct+1) = fa_ani_theta(0,0,iel)
        fa_ani_phi1(ct+1) = fa_ani_phi(0,0,iel)
     endif
 
     rho1(ct+2) = rho(npol,0,iel)
     vp1(ct+2) = sqrt( (lambda(npol,0,iel)+2.*mu(npol,0,iel) ) / rho1(ct+2)  )
     vs1(ct+2) = sqrt( mu(npol,0,iel)  / rho1(ct+2)  )
 
     if (plot_ani) then
        vpv1(ct+2) = sqrt(phi_ani(npol,0,iel)) * vp1(ct+2)
        vsv1(ct+2) = vs1(ct+2) / sqrt(xi_ani(npol,0,iel))
        xi1(ct+2) = xi_ani(npol,0,iel)
        phi1(ct+2) = phi_ani(npol,0,iel)
        eta1(ct+2) = eta_ani(npol,0,iel)
        fa_ani_theta1(ct+2) = fa_ani_theta(npol,0,iel)
        fa_ani_phi1(ct+2) = fa_ani_phi(npol,0,iel)
     endif
 
     rho1(ct+3) = rho(npol,npol,iel)
     vp1(ct+3) = sqrt( (lambda(npol,npol,iel)+2.*mu(npol,npol,iel) ) / rho1(ct+3)  )
     vs1(ct+3) = sqrt( mu(npol,npol,iel)  / rho1(ct+3)  )
 
     if (plot_ani) then
        vpv1(ct+3) = sqrt(phi_ani(npol,npol,iel)) * vp1(ct+3)
        vsv1(ct+3) = vs1(ct+3) / sqrt(xi_ani(npol,npol,iel))
        xi1(ct+3) = xi_ani(npol,npol,iel)
        phi1(ct+3) = phi_ani(npol,npol,iel)
        eta1(ct+3) = eta_ani(npol,npol,iel)
        fa_ani_theta1(ct+3) = fa_ani_theta(npol,npol,iel)
        fa_ani_phi1(ct+3) = fa_ani_phi(npol,npol,iel)
     endif
 
     rho1(ct+4) = rho(0,npol,iel)
     vp1(ct+4) = sqrt( (lambda(0,npol,iel)+2.*mu(0,npol,iel) ) / rho1(ct+4)  )
     vs1(ct+4) = sqrt( mu(0,npol,iel)  / rho1(ct+4)  )
 
     if (plot_ani) then
        vpv1(ct+4) = sqrt(phi_ani(0,npol,iel)) * vp1(ct+4)
        vsv1(ct+4) = vs1(ct+4) / sqrt(xi_ani(0,npol,iel))
        xi1(ct+4) = xi_ani(0,npol,iel)
        phi1(ct+4) = phi_ani(0,npol,iel)
        eta1(ct+4) = eta_ani(0,npol,iel)
        fa_ani_theta1(ct+4) = fa_ani_theta(0,npol,iel)
        fa_ani_phi1(ct+4) = fa_ani_phi(0,npol,iel)
     endif
 
     ct = ct + 4
  enddo
     
  if (plot_anel) then
     do iel=1, nel_solid
        ct = ielsolid(iel) * 4 - 4
        Q_mu1(ct+1) = Q_mu(iel)
        Q_mu1(ct+2) = Q_mu(iel)
        Q_mu1(ct+3) = Q_mu(iel)
        Q_mu1(ct+4) = Q_mu(iel)
        
        Q_kappa1(ct+1) = Q_kappa(iel)
        Q_kappa1(ct+2) = Q_kappa(iel)
        Q_kappa1(ct+3) = Q_kappa(iel)
        Q_kappa1(ct+4) = Q_kappa(iel)
     enddo
  endif


  if (plot_ani) then
     fname=trim(infopath(1:lfinfo)//'/model_vph_'//appmynum)
     call write_VTK_bin_scal(x,y,z0,vp1,npts_vtk/4,fname)
 
     fname=trim(infopath(1:lfinfo)//'/model_vsh_'//appmynum)
     call write_VTK_bin_scal(x,y,z0,vs1,npts_vtk/4,fname)
 
     fname=trim(infopath(1:lfinfo)//'/model_rho_'//appmynum)
     call write_VTK_bin_scal(x,y,z0,rho1,npts_vtk/4,fname)
 
     fname=trim(infopath(1:lfinfo)//'/model_vpv_'//appmynum)
     call write_VTK_bin_scal(x,y,z0,vpv1,npts_vtk/4,fname)
     
     fname=trim(infopath(1:lfinfo)//'/model_vsv_'//appmynum)
     call write_VTK_bin_scal(x,y,z0,vsv1,npts_vtk/4,fname)
     
     fname=trim(infopath(1:lfinfo)//'/model_eta_'//appmynum)
     call write_VTK_bin_scal(x,y,z0,eta1,npts_vtk/4,fname)
     
     fname=trim(infopath(1:lfinfo)//'/model_xi_'//appmynum)
     call write_VTK_bin_scal(x,y,z0,xi1,npts_vtk/4,fname)
     
     fname=trim(infopath(1:lfinfo)//'/model_phi_'//appmynum)
     call write_VTK_bin_scal(x,y,z0,phi1,npts_vtk/4,fname)
 
     fname=trim(infopath(1:lfinfo)//'/model_fa_theta_'//appmynum)
     call write_VTK_bin_scal(x,y,z0,fa_ani_theta1,npts_vtk/4,fname)
 
     fname=trim(infopath(1:lfinfo)//'/model_fa_phi_'//appmynum)
     call write_VTK_bin_scal(x,y,z0,fa_ani_phi1,npts_vtk/4,fname)
 
     deallocate(vp1,vs1,rho1, vsv1, vpv1, eta1, xi1, phi1, fa_ani_theta1, fa_ani_phi1)
  else
     fname=trim(infopath(1:lfinfo)//'/model_vp_'//appmynum)
     call write_VTK_bin_scal(x,y,z0,vp1,npts_vtk/4,fname)
 
     fname=trim(infopath(1:lfinfo)//'/model_vs_'//appmynum)
     call write_VTK_bin_scal(x,y,z0,vs1,npts_vtk/4,fname)
 
     fname=trim(infopath(1:lfinfo)//'/model_rho_'//appmynum)
     call write_VTK_bin_scal(x,y,z0,rho1,npts_vtk/4,fname)
     deallocate(vp1,vs1,rho1)
  endif

  if (plot_anel) then
     fname=trim(infopath(1:lfinfo)//'/model_Qmu_'//appmynum)
     call write_VTK_bin_scal(x,y,z0,Q_mu1,npts_vtk/4,fname)
 
     fname=trim(infopath(1:lfinfo)//'/model_Qkappa_'//appmynum)
     call write_VTK_bin_scal(x,y,z0,Q_kappa1,npts_vtk/4,fname)
     deallocate(Q_mu1, Q_kappa1)
  endif

end subroutine plot_model_vtk
!=============================================================================

!========================
end module get_model
!========================
