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
module model_discontinuities

  use data_bkgrdmodel

  implicit none

  public :: define_discont
  private

  contains

!-----------------------------------------------------------------------------------------
subroutine define_discont
   use data_diag, only: dump_1dmodel

! wrapper routine to call different model types for number of layers ndisc,
! discontinuity radii discont(1:ndisc), and corresponding velocities
! vp(1:ndisc,1:2) and vs(1:ndisc,1:2), where e.g. vs(2,1) is the shear
! velocity at the top of the 2nd domain (just below discont(2), and
! vs(2,2) is the shear velocity at the bottom of the 2nd domain discont(3).

  select case(bkgrdmodel)
     case('ak135')
        write(*,*)'Reading AK135 discontinuities...'
        call ak135_discont
     case('ak135f')
        write(*,*)'Reading AK135f discontinuities...'
        call ak135f_discont
     case('prem_iso')
        write(*,*)'Reading PREM discontinuities...'
        call prem_discont
     case('prem_ani')
        write(*,*)'Reading PREM_ANI discontinuities...'
        call prem_ani_discont
     case('prem_iso_solid')
        write(*,*)'Reading PREM_SOLID discontinuities...'
        call prem_solid_discont
     case('prem_iso_light')
        write(*,*)'Reading PREM_LIGHT discontinuities...'
        call prem_light_discont
     case('prem_ani_light')
        write(*,*)'Reading PREM_LIGHT_ANI discontinuities...'
        call prem_light_ani_discont
     case('prem_iso_onecrust')
        write(*,*)'Reading PREM_ONECRUST discontinuities...'
        call prem_onecrust_discont
     case('prem_ani_onecrust')
        write(*,*)'Reading PREM_ONECRUST_ANI discontinuities...'
        call prem_onecrust_ani_discont
     case('prem_iso_solid_light')
        write(*,*)'Reading PREM_SOLID_LIGHT discontinuities...'
        call prem_solid_light_discont
     case('iasp91')
        write(*,*)'Reading IASP91 discontinuities...'
        call iasp91_discont
     case('prem_iso_custom')
        write(*,*)'Reading PREM discontinuities...'
        call prem_discont
     case('external')
        write(*,*)'Reading step-wise model from file: ', trim(fnam_ext_model)
        call arbitrmodel_discont
     case default
        write(*,*) 'Unknown model' ,bkgrdmodel
        stop
  end select

  if (dump_1dmodel) then
     print *, 'Writing out the current model to Diags/1dmodel.bm'
     call write_1Dmodel(discont)
  endif

end subroutine define_discont
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine ak135f_discont
! Montagner and Kennett 1996

    ndisc = 12

    allocate(discont(ndisc),vp(ndisc,2),vs(ndisc,2))

    ! 0 to 10
    discont(1) = 6371.000000
    vp(1,1) = 5.800000
    vs(1,1) = 3.200000
    vp(1,2) = 5.800000
    vs(1,2) = 3.900000

    ! 10 to 18
    discont(2) = 6361.000000
    vp(2,1) = 6.800000
    vs(2,1) = 3.900000
    vp(2,2) = 6.800000
    vs(2,2) = 3.900000

    ! 18 to 80
    discont(3) = 6353.000000
    vp(3,1) = 8.035500
    vs(3,1) = 4.483900
    vp(3,2) = 8.040000
    vs(3,2) = 4.480000

    ! 80 to 120
    discont(4) = 6291.000000
    vp(4,1) = 8.045000
    vs(4,1) = 4.490000
    vp(4,2) = 8.050500
    vs(4,2) = 4.500000

    ! 120 to 210
    discont(5) = 6251.000000
    vp(5,1) = 8.050500
    vs(5,1) = 4.500000
    vp(5,2) = 8.300700
    vs(5,2) = 4.518400

    ! 210 to 410
    discont(6) = 6161.000000
    vp(6,1) = 8.300700
    vs(6,1) = 4.518400
    vp(6,2) = 9.030200
    vs(6,2) = 4.870200

    ! 410 to 660
    discont(7) = 5961.000000
    vp(7,1) = 9.360100
    vs(7,1) = 5.080600
    vp(7,2) = 10.200000
    vs(7,2) = 5.610400

    ! 660 to 760
    discont(8) = 5711.000000
    vp(8,1) = 10.790900
    vs(8,1) = 5.960700
    vp(8,2) = 11.055300
    vs(8,2) = 6.210000

    ! 760 to 2740
    discont(9) = 5611.000000
    vp(9,1) = 11.055300
    vs(9,1) = 6.210000
    vp(9,2) = 13.649800
    vs(9,2) = 7.248500

    ! 2740 to 2891.5
    discont(10) = 3631.000000
    vp(10,1) = 13.649800
    vs(10,1) = 7.248500
    vp(10,2) = 13.660100
    vs(10,2) = 7.281700

    ! 2891.5 to 5153.5
    discont(11) = 3479.500000
    vp(11,1) = 8.000000
    vs(11,1) = 0.000000
    vp(11,2) = 10.289000
    vs(11,2) = 0.000000

    ! 5153.5 to center
    discont(12) = 1217.500000
    vp(12,1) = 11.042700
    vs(12,1) = 3.504300
    vp(12,2) = 11.262200
    vs(12,2) = 3.667800

    ! numbering relates to regions within, i.e. counting numbers as in discont
    ! for regions above the respective discontinuities

    vp = vp * 1000.
    vs = vs * 1000.
    discont = discont * 1000.

end subroutine ak135f_discont
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine ak135_discont
! Kennett et al., 1995: Constrains on seismic velocities in the Earth

  ndisc = 9

  allocate(discont(ndisc),vp(ndisc,2),vs(ndisc,2))

  ! upper crust
  discont(1) = 6371.
  vp(1,:) = 5.8
  vs(1,:) = 3.46

  ! lower crust
  discont(2) = 6351.
  vp(2,:) = 6.5
  vs(2,:) = 3.85

  ! moho -> 210
  discont(3) = 6336.
  vp(3,1) = 8.04
  vs(3,1) = 4.48
  vp(3,2) = 8.3
  vs(3,2) = 4.518

  ! 210 -> 410
  discont(4) = 6161.
  vp(4,1) = 8.3
  vs(4,1) = 4.523
  vp(4,2) = 9.03
  vs(4,2) = 4.870

  ! 410 -> 660
  discont(5) = 5961.
  vp(5,1) = 9.36
  vs(5,1) = 5.08
  vp(5,2) = 10.2
  vs(5,2) = 5.61

  ! 660 -> D''
  discont(6) = 5711.
  vp(6,1) = 10.79
  vs(6,1) = 5.96
  vp(6,2) = 13.649
  vs(6,2) = 7.249

  ! D'' -> outer core
  discont(7) = 3631.
  vp(7,1) = 13.649
  vs(7,1) = 7.249
  vp(7,2) = 13.66
  vs(7,2) = 7.281

  ! outer core -> inner core
  discont(8) = 3479.5
  vp(8,1) = 8.0
  vs(8,1) = 0.0
  vp(8,2) = 10.289
  vs(8,2) = 0.0

  ! inner core
  discont(9) = 1217.5
  vp(9,1) = 11.043
  vs(9,1) = 3.504
  vp(9,2) = 11.262
  vs(9,2) = 3.668

  ! numbering relates to regions within, i.e. counting numbers as in discont
  ! for regions above the respective discontinuities

  vp = vp * 1000.
  vs = vs * 1000.
  discont = discont * 1000.

end subroutine ak135_discont
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine prem_discont
! PREM discontinuities to be honored by the mesh
! each index represents the layer *below* its corresponding discontinuity

  ndisc = 11

  allocate(discont(ndisc),vp(ndisc,2),vs(ndisc,2))

  ! UPPER CRUST
  discont(1) = 6371.
  vp(1,:) = 5.8
  vs(1,:) = 3.2

  ! LOWER CRUST
  discont(2) = 6356.
  vp(2,:) = 6.8
  vs(2,:) = 3.9

  ! MOHO --> 220
  discont(3) = 6346.6
  vp(3,1) = 8.1106
  vs(3,1) = 4.4910
  vp(3,2) = 7.9897
  vs(3,2) = 4.4189

  ! TRANSITION ZONE: 220 --> 400
  discont(4) = 6151.
  vp(4,1) = 8.55895
  vs(4,1) = 4.64390
  vp(4,2) = 8.90524
  vs(4,2) = 4.76990

  ! 400 --> 600
  discont(5) = 5971.
  vp(5,1) = 9.133917
  vs(5,1) = 4.932487
  vp(5,2) = 10.15783
  vs(5,2) = 5.515931

  ! 600 --> 670
  discont(6) = 5771.
  vp(6,1) = 10.15776
  vs(6,1) = 5.516017
  vp(6,2) = 10.26617
  vs(6,2) = 5.570211

  ! 670 --> 770
  discont(7) = 5701.
  vp(7,1) = 10.7513
  vs(7,1) = 5.9451
  vp(7,2) = 11.0656
  vs(7,2) = 6.2405

  !LOWER MANTLE: 770 --> TOP D"
  discont(8) = 5600.
  vp(8,1) = 11.0656
  vs(8,1) = 6.2404
  vp(8,2) = 13.6804
  vs(8,2) = 7.2659

  ! D" LAYER
  discont(9) = 3630.
  vp(9,1) = 13.6805
  vs(9,1) = 7.2660
  vp(9,2) = 13.7166
  vs(9,2) = 7.2647

  ! FLUID OUTER CORE: CMB --> ICB
  discont(10) = 3480.
  vp(10,1) = 8.0650
  vs(10,1) = 0.0
  vp(10,2) = 10.3557
  vs(10,2) = 0.0

  ! SOLID INNER CORE: ICB --> CENTER
  discont(11) = 1221.5
  vp(11,1) = 11.0283
  vs(11,1) = 3.5043
  vp(11,2) = 11.2622
  vs(11,2) = 3.6678

  ! numbering relates to regions within, i.e. counting numbers as in discont
  ! for regions above the respective discontinuities

  vp = vp * 1000.
  vs = vs * 1000.
  discont = discont * 1000.

end subroutine prem_discont
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine prem_ani_discont

! PREM discontinuities to be honored by the mesh
! each index represents the layer *below* its corresponding discontinuity

! for the anisotropic case: vp = max(vpv, vph), vs=min(vsv, vsh)

  ndisc = 11

  allocate(discont(ndisc),vp(ndisc,2),vs(ndisc,2))

  ! UPPER CRUST
  discont(1) = 6371.
  vp(1,:) = 5.8
  vs(1,:) = 3.2

  ! LOWER CRUST
  discont(2) = 6356.
  vp(2,:) = 6.8
  vs(2,:) = 3.9

  ! MOHO --> 220
  discont(3) = 6346.6
  vp(3,1) = 8.19032
  vs(3,1) = 4.39602
  vp(3,2) = 8.04856
  vs(3,2) = 4.43626

  ! TRANSITION ZONE: 220 --> 400
  discont(4) = 6151.
  vp(4,1) = 8.55895
  vs(4,1) = 4.64390
  vp(4,2) = 8.90524
  vs(4,2) = 4.76990

  ! 400 --> 600
  discont(5) = 5971.
  vp(5,1) = 9.133917
  vs(5,1) = 4.932487
  vp(5,2) = 10.15783
  vs(5,2) = 5.515931

  ! 600 --> 670
  discont(6) = 5771.
  vp(6,1) = 10.15776
  vs(6,1) = 5.516017
  vp(6,2) = 10.26617
  vs(6,2) = 5.570211

  ! 670 --> 770
  discont(7) = 5701.
  vp(7,1) = 10.7513
  vs(7,1) = 5.9451
  vp(7,2) = 11.0656
  vs(7,2) = 6.2405

  !LOWER MANTLE: 770 --> TOP D"
  discont(8) = 5600.
  vp(8,1) = 11.0656
  vs(8,1) = 6.2404
  vp(8,2) = 13.6804
  vs(8,2) = 7.2659

  ! D" LAYER
  discont(9) = 3630.
  vp(9,1) = 13.6805
  vs(9,1) = 7.2660
  vp(9,2) = 13.7166
  vs(9,2) = 7.2647

  ! FLUID OUTER CORE: CMB --> ICB
  discont(10) = 3480.
  vp(10,1) = 8.0650
  vs(10,1) = 0.0
  vp(10,2) = 10.3557
  vs(10,2) = 0.0

  ! SOLID INNER CORE: ICB --> CENTER
  discont(11) = 1221.5
  vp(11,1) = 11.0283
  vs(11,1) = 3.5043
  vp(11,2) = 11.2622
  vs(11,2) = 3.6678

  ! numbering relates to regions within, i.e. counting numbers as in discont
  ! for regions above the respective discontinuities

  vp = vp * 1000.
  vs = vs * 1000.
  discont = discont * 1000.

end subroutine prem_ani_discont
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine prem_solid_discont

! PREM discontinuities to be honored by the mesh
! each index represents the layer *below* its corresponding discontinuity

  ndisc = 11

  allocate(discont(ndisc),vp(ndisc,2),vs(ndisc,2))

  ! UPPER CRUST

  discont(1) = 6371.
  vp(1,:) = 5.8
  vs(1,:) = 3.2

  ! LOWER CRUST
  discont(2) = 6356.
  vp(2,:) = 6.8
  vs(2,:) = 3.9

  ! MOHO --> 220
  discont(3) = 6346.6
  vp(3,1) = 8.1106
  vs(3,1) = 4.4910
  vp(3,2) = 7.9897
  vs(3,2) = 4.4189

  ! TRANSITION ZONE: 220 --> 400
  discont(4) = 6151.
  vp(4,1) = 8.55895
  vs(4,1) = 4.64390
  vp(4,2) = 8.90524
  vs(4,2) = 4.76990

  ! 400 --> 600
  discont(5) = 5971.
  vp(5,1) = 9.133917
  vs(5,1) = 4.932487
  vp(5,2) = 10.15783
  vs(5,2) = 5.515931

  ! 600 --> 670
  discont(6) = 5771.
  vp(6,1) = 10.15776
  vs(6,1) = 5.516017
  vp(6,2) = 10.26617
  vs(6,2) = 5.570211

  ! 670 --> 770
  discont(7) = 5701.
  vp(7,1) = 10.7513
  vs(7,1) = 5.9451
  vp(7,2) = 11.0656
  vs(7,2) = 6.2405

  !LOWER MANTLE: 770 --> TOP D"
  discont(8) = 5600.
  vp(8,1) = 11.0656
  vs(8,1) = 6.2404
  vp(8,2) = 13.6804
  vs(8,2) = 7.2659

  ! D" LAYER
  discont(9) = 3630.
  vp(9,1) = 13.6805
  vs(9,1) = 7.2660
  vp(9,2) = 13.7166
  vs(9,2) = 7.2647

  ! FLUID OUTER CORE: CMB --> ICB
  discont(10) = 3480.
  vp(10,1) = 8.0650
  vs(10,1) = vp(10,1)/sqrt(3.)
  vp(10,2) = 10.3557
  vs(10,2) = vp(10,2)/sqrt(3.)

  ! SOLID INNER CORE: ICB --> CENTER
  discont(11) = 1221.5
  vp(11,1) = 11.0283
  vs(11,1) = 3.5043
  vp(11,2) = 11.2622
  vs(11,2) = 3.6678

  ! numbering relates to regions within, i.e. counting numbers as in discont
  ! for regions above the respective discontinuities

  vp = vp * 1000.
  vs = vs * 1000.
  discont = discont * 1000.

end subroutine prem_solid_discont
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine prem_onecrust_discont

! PREM discontinuities to be honored by the mesh
! but with lower crust extended to the surface
! each index represents the layer *below* its corresponding discontinuity

  ndisc = 10

  allocate(discont(ndisc),vp(ndisc,2),vs(ndisc,2))

  ! ONE CRUST
  discont(1) = 6371.
  vp(1,:) = 5.8
  vs(1,:) = 3.2

  ! MOHO --> 220
  discont(2) = 6346.6
  vp(2,1) = 8.1106
  vs(2,1) = 4.4910
  vp(2,2) = 7.9897
  vs(2,2) = 4.4189

  ! TRANSITION ZONE: 220 --> 400
  discont(3) = 6151.
  vp(3,1) = 8.55895
  vs(3,1) = 4.64390
  vp(3,2) = 8.90524
  vs(3,2) = 4.76990

  ! 400 --> 600
  discont(4) = 5971.
  vp(4,1) = 9.133917
  vs(4,1) = 4.932487
  vp(4,2) = 10.15783
  vs(4,2) = 5.515931

  ! 600 --> 670
  discont(5) = 5771.
  vp(5,1) = 10.15776
  vs(5,1) = 5.516017
  vp(5,2) = 10.26617
  vs(5,2) = 5.570211

  ! 670 --> 770
  discont(6) = 5701.
  vp(6,1) = 10.7513
  vs(6,1) = 5.9451
  vp(6,2) = 11.0656
  vs(6,2) = 6.2405

  !LOWER MANTLE: 770 --> TOP D"
  discont(7) = 5600.
  vp(7,1) = 11.0656
  vs(7,1) = 6.2404
  vp(7,2) = 13.6804
  vs(7,2) = 7.2659

  ! D" LAYER
  discont(8) = 3630.
  vp(8,1) = 13.6805
  vs(8,1) = 7.2660
  vp(8,2) = 13.7166
  vs(8,2) = 7.2647

  ! FLUID OUTER CORE: CMB --> ICB
  discont(9) = 3480.
  vp(9,1) = 8.0650
  vs(9,1) = 0.0
  vp(9,2) = 10.3557
  vs(9,2) = 0.0

  ! SOLID INNER CORE: ICB --> CENTER
  discont(10) = 1221.5
  vp(10,1) = 11.0283
  vs(10,1) = 3.5043
  vp(10,2) = 11.2622
  vs(10,2) = 3.6678

  ! numbering relates to regions within, i.e. counting numbers as in discont
  ! for regions above the respective discontinuities

  vp = vp * 1000.
  vs = vs * 1000.
  discont = discont * 1000.

end subroutine prem_onecrust_discont
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine prem_onecrust_ani_discont

! PREM discontinuities to be honored by the mesh
! each index represents the layer *below* its corresponding discontinuity

! for the anisotropic case: vp = max(vpv, vph), vs=min(vsv, vsh)

  ndisc = 10

  allocate(discont(ndisc),vp(ndisc,2),vs(ndisc,2))

  ! UPPER CRUST
  discont(1) = 6371.
  vp(1,:) = 5.8
  vs(1,:) = 3.2

  ! MOHO --> 220
  discont(2) = 6346.6
  vp(2,1) = 8.19032
  vs(2,1) = 4.39602
  vp(2,2) = 8.04856
  vs(2,2) = 4.43626

  ! TRANSITION ZONE: 220 --> 400
  discont(3) = 6151.
  vp(3,1) = 8.55895
  vs(3,1) = 4.64390
  vp(3,2) = 8.90524
  vs(3,2) = 4.76990

  ! 400 --> 600
  discont(4) = 5971.
  vp(4,1) = 9.133917
  vs(4,1) = 4.932487
  vp(4,2) = 10.15783
  vs(4,2) = 5.515931

  ! 600 --> 670
  discont(5) = 5771.
  vp(5,1) = 10.15776
  vs(5,1) = 5.516017
  vp(5,2) = 10.26617
  vs(5,2) = 5.570211

  ! 670 --> 770
  discont(6) = 5701.
  vp(6,1) = 10.7513
  vs(6,1) = 5.9451
  vp(6,2) = 11.0656
  vs(6,2) = 6.2405

  !LOWER MANTLE: 770 --> TOP D"
  discont(7) = 5600.
  vp(7,1) = 11.0656
  vs(7,1) = 6.2404
  vp(7,2) = 13.6804
  vs(7,2) = 7.2659

  ! D" LAYER
  discont(8) = 3630.
  vp(8,1) = 13.6805
  vs(8,1) = 7.2660
  vp(8,2) = 13.7166
  vs(8,2) = 7.2647

  ! FLUID OUTER CORE: CMB --> ICB
  discont(9) = 3480.
  vp(9,1) = 8.0650
  vs(9,1) = 0.0
  vp(9,2) = 10.3557
  vs(9,2) = 0.0

  ! SOLID INNER CORE: ICB --> CENTER
  discont(10) = 1221.5
  vp(10,1) = 11.0283
  vs(10,1) = 3.5043
  vp(10,2) = 11.2622
  vs(10,2) = 3.6678

  ! numbering relates to regions within, i.e. counting numbers as in discont
  ! for regions above the respective discontinuities

  vp = vp * 1000.
  vs = vs * 1000.
  discont = discont * 1000.

end subroutine prem_onecrust_ani_discont
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine prem_light_discont

! PREM LIGHT discontinuities to be honored by the mesh:
! isotropic PREM without the crust, extending the upper mantle to the surface
! each index represents the layer *below* its corresponding discontinuity

  ndisc = 9

  allocate(discont(ndisc),vp(ndisc,2),vs(ndisc,2))

  ! 0 --> 220
  discont(1) = 6371.
  vp(1,1) = 8.1257
  vs(1,1) = 4.5000
  vp(1,2) = 7.9897
  vs(1,2) = 4.4189

  ! TRANSITION ZONE: 220 --> 400
  discont(2) = 6151.
  vp(2,1) = 8.55895
  vs(2,1) = 4.64390
  vp(2,2) = 8.90524
  vs(2,2) = 4.76990

  ! 400 --> 600
  discont(3) = 5971.
  vp(3,1) = 9.133917
  vs(3,1) = 4.932487
  vp(3,2) = 10.15783
  vs(3,2) = 5.515931

  ! 600 --> 670
  discont(4) = 5771.
  vp(4,1) = 10.15776
  vs(4,1) = 5.516017
  vp(4,2) = 10.26617
  vs(4,2) = 5.570211

  ! 670 --> 770
  discont(5) = 5701.
  vp(5,1) = 10.7513
  vs(5,1) = 5.9451
  vp(5,2) = 11.0656
  vs(5,2) = 6.2405

  !LOWER MANTLE: 770 --> TOP D"
  discont(6) = 5600.
  vp(6,1) = 11.0656
  vs(6,1) = 6.2404
  vp(6,2) = 13.6804
  vs(6,2) = 7.2659

  ! D" LAYER
  discont(7) = 3630.
  vp(6,1) = 13.6805
  vs(7,1) = 7.2660
  vp(7,2) = 13.7166
  vs(7,2) = 7.2647

  ! FLUID OUTER CORE: CMB --> ICB
  discont(8) = 3480.
  vp(8,1) = 8.0650
  vs(8,1) = 0.0
  vp(8,2) = 10.3557
  vs(8,2) = 0.0

  ! SOLID INNER CORE: ICB --> CENTER
  discont(9) = 1221.5
  vp(9,1) = 11.0283
  vs(9,1) = 3.5043
  vp(9,2) = 11.2622
  vs(9,2) = 3.6678

  ! numbering relates to regions within, i.e. counting numbers as in discont
  ! for regions above the respective discontinuities

  vp = vp * 1000.
  vs = vs * 1000.
  discont = discont * 1000.

end subroutine prem_light_discont
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine prem_light_ani_discont

! PREM discontinuities to be honored by the mesh
! each index represents the layer *below* its corresponding discontinuity

! for the anisotropic case: vp = max(vpv, vph), vs=min(vsv, vsh)

  ndisc = 9

  allocate(discont(ndisc),vp(ndisc,2),vs(ndisc,2))

  ! MOHO --> 220
  discont(1) = 6371.
  vp(1,1) = 8.20800
  vs(1,1) = 4.39040
  vp(1,2) = 8.04856
  vs(1,2) = 4.43626

  ! TRANSITION ZONE: 220 --> 400
  discont(2) = 6151.
  vp(2,1) = 8.55895
  vs(2,1) = 4.64390
  vp(2,2) = 8.90524
  vs(2,2) = 4.76990

  ! 400 --> 600
  discont(3) = 5971.
  vp(3,1) = 9.133917
  vs(3,1) = 4.932487
  vp(3,2) = 10.15783
  vs(3,2) = 5.515931

  ! 600 --> 670
  discont(4) = 5771.
  vp(4,1) = 10.15776
  vs(4,1) = 5.516017
  vp(4,2) = 10.26617
  vs(4,2) = 5.570211

  ! 670 --> 770
  discont(5) = 5701.
  vp(5,1) = 10.7513
  vs(5,1) = 5.9451
  vp(5,2) = 11.0656
  vs(5,2) = 6.2405

  !LOWER MANTLE: 770 --> TOP D"
  discont(6) = 5600.
  vp(6,1) = 11.0656
  vs(6,1) = 6.2404
  vp(6,2) = 13.6804
  vs(6,2) = 7.2659

  ! D" LAYER
  discont(7) = 3630.
  vp(7,1) = 13.6805
  vs(7,1) = 7.2660
  vp(7,2) = 13.7166
  vs(7,2) = 7.2647

  ! FLUID OUTER CORE: CMB --> ICB
  discont(8) = 3480.
  vp(8,1) = 8.0650
  vs(8,1) = 0.0
  vp(8,2) = 10.3557
  vs(8,2) = 0.0

  ! SOLID INNER CORE: ICB --> CENTER
  discont(9) = 1221.5
  vp(9,1) = 11.0283
  vs(9,1) = 3.5043
  vp(9,2) = 11.2622
  vs(9,2) = 3.6678

  ! numbering relates to regions within, i.e. counting numbers as in discont
  ! for regions above the respective discontinuities

  vp = vp * 1000.
  vs = vs * 1000.
  discont = discont * 1000.

end subroutine prem_light_ani_discont
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine prem_solid_light_discont
! PREM LIGHT discontinuities to be honored by the mesh:
! isotropic PREM without the crust, extending the upper mantle to the surface
! No fluid outer core, but instead vs=vp/sqrt(3)
! each index represents the layer *below* its corresponding discontinuity

  ndisc = 9

  allocate(discont(ndisc),vp(ndisc,2),vs(ndisc,2))

  ! 0 --> 220
  discont(1) = 6371.
  vp(1,1) = 8.1257
  vs(1,1) = 4.5000
  vp(1,2) = 7.9897
  vs(1,2) = 4.4189

  ! TRANSITION ZONE: 220 --> 400
  discont(2) = 6151.
  vp(2,1) = 8.55895
  vs(2,1) = 4.64390
  vp(2,2) = 8.90524
  vs(2,2) = 4.76990

  ! 400 --> 600
  discont(3) = 5971.
  vp(3,1) = 9.133917
  vs(3,1) = 4.932487
  vp(3,2) = 10.15783
  vs(3,2) = 5.515931

  ! 600 --> 670
  discont(4) = 5771.
  vp(4,1) = 10.15776
  vs(4,1) = 5.516017
  vp(4,2) = 10.26617
  vs(4,2) = 5.570211

  ! 670 --> 770
  discont(5) = 5701.
  vp(5,1) = 10.7513
  vs(5,1) = 5.9451
  vp(5,2) = 11.0656
  vs(5,2) = 6.2405

  !LOWER MANTLE: 770 --> TOP D"
  discont(6) = 5600.
  vp(6,1) = 11.0656
  vs(6,1) = 6.2404
  vp(6,2) = 13.6804
  vs(6,2) = 7.2659

  ! D" LAYER
  discont(7) = 3630.
  vp(6,1) = 13.6805
  vs(7,1) = 7.2660
  vp(7,2) = 13.7166
  vs(7,2) = 7.2647

  ! FLUID OUTER CORE: CMB --> ICB
  discont(8) = 3480.
  vp(8,1) = 8.0650
  vs(8,1) = vp(8,1)/sqrt(3.)
  vp(8,2) = 10.3557
  vs(8,2) = vp(8,2)/sqrt(3.)

  ! SOLID INNER CORE: ICB --> CENTER
  discont(9) = 1221.5
  vp(9,1) = 11.0283
  vs(9,1) = 3.5043
  vp(9,2) = 11.2622
  vs(9,2) = 3.6678

  ! numbering relates to regions within, i.e. counting numbers as in discont
  ! for regions above the respective discontinuities

  vp = vp * 1000.
  vs = vs * 1000.
  discont = discont * 1000.

end subroutine prem_solid_light_discont
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine iasp91_discont
! IASP91 discontinuities to be honored by the mesh
! each index represents the layer *below* its corresponding discontinuity

  ! Tabulated all media parameters just above and below discontinuities
  !!$    radius         v_p              v_s            rho
  !!$   6371.000       5.800000       3.360000       2.720000
  !!$   6351.000       6.500000       3.750000       2.920000
  !!$   6336.010       6.500000       3.750000       2.920000
  !!$   6335.990       8.039999       4.470003       3.319806
  !!$   6291.010       8.045291       4.485878       3.347059
  !!$   6290.990       8.045293       4.485885       3.347071
  !!$   6251.010       8.049996       4.499995       3.371294
  !!$   6250.990       8.050032       4.500002       3.370357
  !!$   6161.010       8.299976       4.517998       3.360578
  !!$   6160.990       8.300035       4.522027       3.429809
  !!$   5961.010       9.029963       4.869992       3.549229
  !!$   5960.990       9.360033       5.070020       3.931578
  !!$   5711.010       10.19997       5.599978       3.989790
  !!$   5710.990       10.79003       5.950025       4.374491
  !!$   5611.010       11.05577       6.209453       4.436460
  !!$   5610.990       11.05585       6.209506       4.436472
  !!$   3631.010       13.65642       7.264518       5.490973
  !!$   3630.990       13.65640       7.264504       5.490983
  !!$   3482.010       13.69080       7.301500       5.565448
  !!$   3481.990       8.008780      0.0000000E+00   9.900254
  !!$   1217.010       10.25781      0.0000000E+00   12.16863
  !!$   1216.990       11.09145       3.438566       12.76601
  !!$  9.9999998E-03   11.24094       3.564540       13.08850

  ndisc = 11

  allocate(discont(ndisc),vp(ndisc,2),vs(ndisc,2))

  ! UPPER CRUST
  discont(1) = 6371.
  vp(1,:) = 5.8
  vs(1,:) = 3.36

  ! LOWER CRUST
  discont(2) = 6351.
  vp(2,:) = 6.5
  vs(2,:) = 3.75

  ! MOHO --> 120
  discont(3) = 6336.
  vp(3,1) = 8.04
  vs(3,1) = 4.47
  vp(3,2) = 8.049996
  vs(3,2) = 4.499995

  ! 120 --> 220
  discont(4) = 6251.
  vp(4,1) = 8.050032
  vs(4,1) = 4.500002
  vp(4,2) = 8.299976
  vs(4,2) = 4.517998

  ! 220 --> 400
  discont(5) = 6161.
  vp(5,1) = 8.300035
  vs(5,1) = 4.522027
  vp(5,2) = 9.029963
  vs(5,2) = 4.869992

  ! 400 --> 660
  discont(6) = 5961.
  vp(6,1) = 9.360033
  vs(6,1) = 5.070020
  vp(6,2) = 10.19997
  vs(6,2) = 5.599978

  ! 660 --> 771
  discont(7) = 5711.
  vp(7,1) = 10.79003
  vs(7,1) = 5.950025
  vp(7,2) = 11.05577
  vs(7,2) = 6.209453

  ! 771 --> D"
  discont(8) = 5611.
  vp(8,1) = 11.05585
  vs(8,1) = 6.209506
  vp(8,2) = 13.65642
  vs(8,2) = 7.264518

  ! D" --> CMB
  discont(9) = 3631.
  vp(9,1) = 13.65640
  vs(9,1) = 7.264504
  vp(9,2) = 13.69080
  vs(9,2) = 7.301500

  ! CMB --> ICB
  discont(10) = 3482.
  vp(10,1) = 8.008780
  vs(10,1) = 0.0
  vp(10,2) = 10.25781
  vs(10,2) = 0.0

  ! ICB --> CENTER
  discont(11) = 1217.
  vp(11,1) = 11.09145
  vs(11,1) = 3.438566
  vp(11,2) = 11.24094
  vs(11,2) = 3.564540

  ! numbering relates to regions within, i.e. counting numbers as in discont
  ! for regions above the respective discontinuities

  vp = vp * 1000.
  vs = vs * 1000.
  discont = discont * 1000.

end subroutine iasp91_discont
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine arbitrmodel_discont

  use background_models, only: get_ext_disc
  use data_grid, only: router
! discontinuities (read in from a file) to be honored by the mesh

  integer :: idom

  call get_ext_disc(fnam_ext_model, ndisc, discont, vp, vs, rho)

  router = discont(1)

  do idom = 1, ndisc
     write(1001,*) discont(idom), discont(1) - discont(idom), vp(idom,:), vs(idom,:)
  enddo

end subroutine arbitrmodel_discont
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine write_1Dmodel(discontinuities)
   ! Write out the current model in a .bm file, which can be reused by the mesher.
   use global_parameters, only: smallval_dble
   use background_models, only: velocity, model_is_ani, model_is_anelastic
   use data_diag, only: diagpath

   real(kind=dp), intent(in) :: discontinuities(:)
   integer, parameter        :: maxlayers = 10000
   real(kind=dp), dimension(0:maxlayers)  :: vpv, vsv, rho, radius
   real(kind=dp), dimension(0:maxlayers)  :: qka, qmu, vph, vsh, eta
   real(kind=dp)             :: vp_tmp, vs_tmp
   integer                   :: ndom, irad, idom, ilayer, nlayer, step, nic, noc
   integer, allocatable      :: disc_layer(:)
   character(len=256)        :: fnam, fmtstring
   character(len=8)          :: mydate
   character(len=10)         :: mytime

   ndom = size(discontinuities)
   allocate(disc_layer(ndom))

   ilayer = 0

   ! Domains from first to second-last
   do idom = 1, ndom-1
      ! Adapt layer spacing according to overall size of domain. Should be at least 5
      ! layers in each domain
      if (discontinuities(idom) - discontinuities(idom+1) < 25000.) then
         step = -1000
      else if (discontinuities(idom) - discontinuities(idom+1) < 250000.) then
         step = -5000
      else
         step = -10000
      endif
      fmtstring = "(' Domain:', I3, ', width: ', F12.1, ', step:', I7)"
      print fmtstring, idom, discontinuities(idom) - discontinuities(idom+1), step
      ! Layers within the domain
      do irad = nint(discontinuities(idom)), nint(discontinuities(idom+1)), step
         vp_tmp = velocity(real(irad, kind=dp), 'v_p', idom, bkgrdmodel, lfbkgrdmodel)
         vs_tmp = velocity(real(irad, kind=dp), 'v_s', idom, bkgrdmodel, lfbkgrdmodel)
         if (vp_tmp /= vpv(ilayer) .or. vs_tmp /= vsv(ilayer)) then
            ilayer = ilayer + 1
            radius(ilayer) = real(irad, kind=dp)
            vpv(ilayer)   = vp_tmp
            vsv(ilayer)   = vs_tmp
            rho(ilayer)   = velocity(real(irad, kind=dp), 'rho', idom, bkgrdmodel, lfbkgrdmodel)
            if (model_is_anelastic(bkgrdmodel)) then
               qka(ilayer)   = velocity(real(irad, kind=dp), 'Qka', idom, bkgrdmodel, lfbkgrdmodel)
               qmu(ilayer)   = velocity(real(irad, kind=dp), 'Qmu', idom, bkgrdmodel, lfbkgrdmodel)
               vph(ilayer)   = velocity(real(irad, kind=dp), 'vph', idom, bkgrdmodel, lfbkgrdmodel)
               vsh(ilayer)   = velocity(real(irad, kind=dp), 'vsh', idom, bkgrdmodel, lfbkgrdmodel)
               eta(ilayer)   = velocity(real(irad, kind=dp), 'eta', idom, bkgrdmodel, lfbkgrdmodel)
            endif

         endif
         !write(2000,*) real(irad, kind=dp), rho, vp, vs

      enddo

      ! Layer at the bottom of the domain
      if ((radius(ilayer)-discontinuities(idom+1)) > smallval_dble*radius(ilayer)) then
         ilayer = ilayer + 1
         radius(ilayer) = discontinuities(idom+1)
         vpv(ilayer)   = velocity(discontinuities(idom+1), 'v_p', idom, bkgrdmodel, lfbkgrdmodel)
         vsv(ilayer)   = velocity(discontinuities(idom+1), 'v_s', idom, bkgrdmodel, lfbkgrdmodel)
         rho(ilayer)   = velocity(discontinuities(idom+1), 'rho', idom, bkgrdmodel, lfbkgrdmodel)
         if (model_is_anelastic(bkgrdmodel)) then
             qka(ilayer)   = velocity(discontinuities(idom+1), 'Qka', idom, bkgrdmodel, lfbkgrdmodel)
             qmu(ilayer)   = velocity(discontinuities(idom+1), 'Qmu', idom, bkgrdmodel, lfbkgrdmodel)
             vph(ilayer)   = velocity(discontinuities(idom+1), 'vph', idom, bkgrdmodel, lfbkgrdmodel)
             vsh(ilayer)   = velocity(discontinuities(idom+1), 'vsh', idom, bkgrdmodel, lfbkgrdmodel)
             eta(ilayer)   = velocity(discontinuities(idom+1), 'eta', idom, bkgrdmodel, lfbkgrdmodel)
         endif
      endif

      disc_layer(idom) = ilayer+1
   enddo !idom = 1, ndom-1

   ! Layers within the last domain
   if (discontinuities(ndom) < 25000.) then
      step = -1000
   else if (discontinuities(ndom) < 250000.) then
      step = -5000
   else
      step = -10000
   endif

   fmtstring = "(' Domain:', I3, ', width: ', F12.1, ', step:', I7)"
   print fmtstring, idom, discontinuities(ndom), step

   do irad = nint(discontinuities(ndom)), 0, step
      vp_tmp = velocity(real(irad, kind=dp), 'v_p', idom, bkgrdmodel, lfbkgrdmodel)
      vs_tmp = velocity(real(irad, kind=dp), 'v_s', idom, bkgrdmodel, lfbkgrdmodel)
      if (vp_tmp /= vpv(ilayer) .or. vs_tmp /= vsv(ilayer)) then
         ilayer = ilayer + 1
         radius(ilayer) = real(irad, kind=dp)
         vpv(ilayer)    = vp_tmp
         vsv(ilayer)    = vs_tmp
         rho(ilayer)   = velocity(real(irad, kind=dp), 'rho', idom, bkgrdmodel, lfbkgrdmodel)
         if (model_is_anelastic(bkgrdmodel)) then
            qka(ilayer)   = velocity(real(irad, kind=dp), 'Qka', idom, bkgrdmodel, lfbkgrdmodel)
            qmu(ilayer)   = velocity(real(irad, kind=dp), 'Qmu', idom, bkgrdmodel, lfbkgrdmodel)
            vph(ilayer)   = velocity(real(irad, kind=dp), 'vph', idom, bkgrdmodel, lfbkgrdmodel)
            vsh(ilayer)   = velocity(real(irad, kind=dp), 'vsh', idom, bkgrdmodel, lfbkgrdmodel)
            eta(ilayer)   = velocity(real(irad, kind=dp), 'eta', idom, bkgrdmodel, lfbkgrdmodel)
         endif
      endif
   enddo

   ! Layer at the bottom of the model
   if (radius(ilayer) > smallval_dble) then
      ilayer = ilayer + 1
      radius(ilayer) = 0.0d0
      vpv(ilayer)    = velocity(0.0d0, 'v_p', idom, bkgrdmodel, lfbkgrdmodel)
      vsv(ilayer)    = velocity(0.0d0, 'v_s', idom, bkgrdmodel, lfbkgrdmodel)
      rho(ilayer)   = velocity(0.0d0, 'rho', idom, bkgrdmodel, lfbkgrdmodel)
      if (model_is_anelastic(bkgrdmodel)) then
          qka(ilayer)   = velocity(0.0d0, 'Qka', idom, bkgrdmodel, lfbkgrdmodel)
          qmu(ilayer)   = velocity(0.0d0, 'Qmu', idom, bkgrdmodel, lfbkgrdmodel)
          vph(ilayer)   = velocity(0.0d0, 'vph', idom, bkgrdmodel, lfbkgrdmodel)
          vsh(ilayer)   = velocity(0.0d0, 'vsh', idom, bkgrdmodel, lfbkgrdmodel)
          eta(ilayer)   = velocity(0.0d0, 'eta', idom, bkgrdmodel, lfbkgrdmodel)
      endif
   endif

   print *
   nlayer = ilayer

   ! Write input file for AxiSEM
   fnam = trim(diagpath)//'/1dmodel_axisem.bm'

   call date_and_time(mydate,mytime)

   open(2000, file=fnam, action='write')
11 format('# Input file for AXISEM created from external model on ', &
            A2,'/',A2,'/',A4,', at ',A2,'h ',A2,'min')
   write(2000,11) mydate(7:8), mydate(5:6), mydate(1:4), mytime(1:2), mytime(3:4)
   write(2000,'("ANELASTIC    ", L4)') model_is_anelastic(bkgrdmodel)
   write(2000,'("ANISOTROPIC  ", L4)') model_is_ani(bkgrdmodel)
   write(2000,'("UNITS        m")')
   if (model_is_anelastic(bkgrdmodel)) then
      if (model_is_ani(bkgrdmodel)) then !ANI=true, ANE=true
         write(2000,'(A11, 9(A9))') 'COLUMNS    ', 'radius', 'rho', 'vpv', 'vsv', 'qka', &
                                               'qmu', 'vph', 'vsh', 'eta'
         idom = 1
         do ilayer = 1, nlayer
            if (ilayer == disc_layer(idom)) then
                write(2000, '("#          Discontinuity ", I3, ", depth: ", F10.2, " km")') idom, &
                             (radius(1)-radius(ilayer))*0.001
                idom = idom + 1
            endif
            write(2000, '("           ", f9.0, 3f9.2, 2f11.1, 2f9.2, f9.5)') &
                        radius(ilayer), rho(ilayer), vpv(ilayer), vsv(ilayer), &
                        qka(ilayer), qmu(ilayer), vph(ilayer), vsh(ilayer), eta(ilayer)
         enddo

      else !ANI=false, ANE=true
         write(2000,'(A11, 6(A9))') 'COLUMNS    ', 'radius', 'rho', 'vpv', 'vsv', 'qka', &
                                               'qmu'
         idom = 1
         do ilayer = 1, nlayer
            if (ilayer == disc_layer(idom)) then
                write(2000, '("#          Discontinuity ", I3, ", depth: ", F10.2, " km")') idom, &
                             (radius(1)-radius(ilayer))*0.001
                idom = idom + 1
            endif
            write(2000, '("           ", f9.0, 3f9.2, 2f11.1)') &
                        radius(ilayer), rho(ilayer), vpv(ilayer), vsv(ilayer), &
                        qka(ilayer), qmu(ilayer)
         enddo
      endif

   else
      if (model_is_ani(bkgrdmodel)) then !ANI=true, ANE=false
         write(2000,'(A11, 7(A9))') 'COLUMNS    ', 'radius', 'rho', 'vpv', 'vsv', 'vph', &
                                    'vsh', 'eta'
         idom = 1
         do ilayer = 1, nlayer
            if (ilayer == disc_layer(idom)) then
                write(2000, '("#          Discontinuity ", I3, ", depth: ", F10.2, " km")') idom, &
                             (radius(1)-radius(ilayer))*0.001
                idom = idom + 1
            endif
            write(2000, '("           ", f9.0, 3f9.2, 2f9.2, f9.5)') &
                        radius(ilayer), rho(ilayer), vpv(ilayer), vsv(ilayer), &
                        vph(ilayer), vsh(ilayer), eta(ilayer)
         enddo
      else !ANI=false, ANE=false
         write(2000,'(A11, 6(A9))') 'COLUMNS    ', 'radius', 'rho', 'vpv', 'vsv'
         idom = 1
         do ilayer = 1, nlayer
            if (ilayer == disc_layer(idom)) then
                write(2000, '("#          Discontinuity ", I3, ", depth: ", F10.2, " km")') idom, &
                             (radius(1)-radius(ilayer))*0.001
                idom = idom + 1
            endif
            write(2000, '("           ", f9.0, 3f9.2, 2f9.1)') &
                        radius(ilayer), rho(ilayer), vpv(ilayer), vsv(ilayer)
         enddo
      endif
   endif

   close(2000)

   ! Write input file for YSPEC
   fnam = trim(diagpath)//'/1dmodel_yspec.bm'
   noc = count(radius(1:nlayer) <= discontinuities(ndom-1) .and. radius(1:nlayer) > discontinuities(ndom))
   nic = count(radius(1:nlayer) <= discontinuities(ndom)) - 1
   open(2000, file=fnam, action='write')
   write (2000,*) 'AXISEM model for YSPEC: ', bkgrdmodel(1:lfbkgrdmodel)
   if (model_is_ani(bkgrdmodel)) then
       write (2000,'(i2,f4.1)') 1, 1.
   else
       write (2000,'(i2,f4.1)') 0, 1.
   endif
   write(2000,*) nlayer, nic, noc
   do ilayer = nlayer, 1, -1
      if (model_is_anelastic(bkgrdmodel)) then
          write(2000, '(f8.0, 3f9.2, 2f9.1, 2f9.2, f9.5)') &
                      radius(ilayer), rho(ilayer), vpv(ilayer), vsv(ilayer), &
                      qka(ilayer), qmu(ilayer), vph(ilayer), vsh(ilayer), eta(ilayer)
      else
          write(2000, '(f8.0, 3f9.2, 2f9.1, 2f9.2, f9.5)') &
                      radius(ilayer), rho(ilayer), vpv(ilayer), vsv(ilayer)
      endif
   enddo
   close(2000)

end subroutine write_1Dmodel
!-----------------------------------------------------------------------------------------

end module model_discontinuities
!=========================================================================================
