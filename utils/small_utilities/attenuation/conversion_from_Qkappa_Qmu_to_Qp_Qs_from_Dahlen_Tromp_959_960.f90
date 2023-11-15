
  program conversion

! Dimitri Komatitsch, CNRS Marseille, France, October 2015 and March 2017

! see formulas 9.59 and 9.60 in the book of Dahlen and Tromp, 1998
! (in that book, P is called alpha and S is called beta).
! The formulas in Dahlen and Tromp, 1998 are for the 3D case,
! in the 2D plane strain case the 4/3 coefficient must be changed to 1,
! see more details in file doc/Qkappa_Qmu_versus_Qp_Qs_relationship_in_2D_plane_strain.pdf

  implicit none

! coefficient for the 3D case is 4/3, for the 2D plane strain case it is 1,
! see more details in file doc/Qkappa_Qmu_versus_Qp_Qs_relationship_in_2D_plane_strain.pdf
  double precision, parameter :: coefficient = 4.d0/3.d0   !! 1.d0

  double precision :: Qkappa,Qmu,Qp,Qs,inverse_of_Qp,cp,cs

!!! this is for the Carcione et al. 1988 example

! enter your Qkappa and Qmu here
  Qkappa = 40.d0
  Qmu = 20.d0

! enter the cp and cs velocities of the medium here, at the frequency at which you want this conversion to be performed
  cp = 3000.d0
  cs = 2000.d0

!!! this is for the Carcione 1993 example

! enter your Qkappa and Qmu here
! Qkappa = 20.d0
! Qmu = 10.d0

! enter the cp and cs velocities of the medium here, at the frequency at which you want this conversion to be performed
! cp = 3249.d0
! cs = 2235.d0

! Qs is the same as Qmu
  Qs = Qmu

! for Qp the formula is more complex
  inverse_of_Qp = (1.d0 - coefficient*(cs**2)/(cp**2))/Qkappa + coefficient*(cs**2)/(cp**2)/Qmu
  Qp = 1.d0/inverse_of_Qp

! print the result
  print *,'Qp = ',Qp
  print *,'Qs = ',Qs

  end program conversion

