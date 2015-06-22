
Demo codes written by Dimitri Komatitsch, CNRS, France, around 2010.
--------------------------------------------------------------------

Useful to learn how the spectral-element method works, and how to write or modify a code to implement it.
Also useful to test new ideas by modifying these simple codes to run some tests.

Here is how to compile and use them:

1/ cd mesher_for_serial

2/ edit the Makefile if needed to change the compiler or compiler options

3/ make clean ; make all

4/ ./xmeshfem3D

5/ cd ..

6/ edit make_all_Fortran.csh and make_all_C.csh if needed to change the compiler or compiler options

7/ ./make_all_Fortran.csh

8/ ./make_all_C.csh

9/ ./xspecfem3D_F90_normal

10/ ./xspecfem3D_F90_faster

11/ ./xspecfem3D_C

and that's it.

