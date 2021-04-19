!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
!               ---------------------------------------
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
!
! United States and French Government Sponsorship Acknowledged.

  program xgenerate_databases

  use adios_manager_mod
  use generate_databases_par

  implicit none

!=============================================================================!
!                                                                             !
!  generate_databases produces a spectral element grid                        !
!  for a local or regional model.                                             !
!  The mesher uses the UTM projection                                         !
!                                                                             !
!=============================================================================!
!
! If you use this code for your own research, please cite at least one article
! written by the developers of the package, for instance:
!
! @ARTICLE{TrKoLi08,
! author = {Jeroen Tromp and Dimitri Komatitsch and Qinya Liu},
! title = {Spectral-Element and Adjoint Methods in Seismology},
! journal = {communications in Computational Physics},
! year = {2008},
! volume = {3},
! pages = {1-32},
! number = {1}}
!
! or
!
! @ARTICLE{LiPoKoTr04,
! author = {Qinya Liu and Jascha Polet and Dimitri Komatitsch and Jeroen Tromp},
! title = {Spectral-element moment tensor inversions for earthquakes in {S}outhern {C}alifornia},
! journal={Bull. Seismol. Soc. Am.},
! year = {2004},
! volume = {94},
! pages = {1748-1761},
! number = {5},
! doi = {10.1785/012004038}}
!
! @INCOLLECTION{ChKoViCaVaFe07,
! author = {Emmanuel Chaljub and Dimitri Komatitsch and Jean-Pierre Vilotte and
! Yann Capdeville and Bernard Valette and Gaetano Festa},
! title = {Spectral Element Analysis in Seismology},
! booktitle = {Advances in Wave Propagation in Heterogeneous Media},
! publisher = {Elsevier - Academic Press},
! year = {2007},
! editor = {Ru-Shan Wu and Val\'erie Maupin},
! volume = {48},
! series = {Advances in Geophysics},
! pages = {365-419}}
!
! @ARTICLE{KoVi98,
! author={D. Komatitsch and J. P. Vilotte},
! title={The spectral-element method: an efficient tool to simulate the seismic response of 2{D} and 3{D} geological structures},
! journal={Bull. Seismol. Soc. Am.},
! year=1998,
! volume=88,
! number=2,
! pages={368-392}}
!
! @ARTICLE{KoTr99,
! author={D. Komatitsch and J. Tromp},
! year=1999,
! title={Introduction to the spectral-element method for 3-{D} seismic wave propagation},
! journal={Geophys. J. Int.},
! volume=139,
! number=3,
! pages={806-822},
! doi={10.1046/j.1365-246x.1999.00967.x}}
!
! @ARTICLE{KoLiTrSuStSh04,
! author={Dimitri Komatitsch and Qinya Liu and Jeroen Tromp and Peter S\"{u}ss
!   and Christiane Stidham and John H. Shaw},
! year=2004,
! title={Simulations of Ground Motion in the {L}os {A}ngeles {B}asin
!   based upon the Spectral-Element Method},
! journal={Bull. Seism. Soc. Am.},
! volume=94,
! number=1,
! pages={187-206}}
!
! and/or another article from http://web.univ-pau.fr/~dkomati1/publications.html
!
!
! If you use the kernel capabilities of the code, please cite at least one article
! written by the developers of the package, for instance:
!
! @ARTICLE{TrKoLi08,
! author = {Jeroen Tromp and Dimitri Komatitsch and Qinya Liu},
! title = {Spectral-Element and Adjoint Methods in Seismology},
! journal = {communications in Computational Physics},
! year = {2008},
! volume = {3},
! pages = {1-32},
! number = {1}}
!
! or
!
! @ARTICLE{LiTr06,
! author={Qinya Liu and Jeroen Tromp},
! title={Finite-frequency kernels based on adjoint methods},
! journal={Bull. Seismol. Soc. Am.},
! year=2006,
! volume=96,
! number=6,
! pages={2383-2397},
! doi={10.1785/0120060041}}
!
!
! Reference frame - convention:
! ----------------------------
!
! The code uses the following convention for the reference frame:
!
!  - X axis is East
!  - Y axis is North
!  - Z axis is up
!
! Note that this convention is different from both the Aki-Richards convention
! and the Harvard CMT convention.
!
! Let us recall that the Aki-Richards convention is:
!
!  - X axis is North
!  - Y axis is East
!  - Z axis is down
!
! and that the Harvard CMT convention is:
!
!  - X axis is South
!  - Y axis is East
!  - Z axis is up
!
! To report bugs or suggest improvements to the code, please use our online
! bug tracking system at http://www.geodynamics.org/roundup .
!
! Evolution of the code:
! ---------------------
!
! MPI v. 2.0, November 2010:
! Dimitri Komatitsch, Nicolas Le Goff, Roland Martin and Pieyre Le Loher, University of Pau, France,
! Daniel Peter, Jeroen Tromp and the Princeton group of developers, Princeton University, USA,
! and Emanuele Casarotti, INGV Roma, Italy:
!  support for CUBIT meshes decomposed by SCOTCH;
!  much faster solver using Michel Deville's inlined matrix products.
!
! MPI v. 1.4 Dimitri Komatitsch, University of Pau, Qinya Liu and others, Caltech, September 2006:
!  better adjoint and kernel calculations, faster and better I/Os
!  on very large systems, many small improvements and bug fixes
!
! MPI v. 1.3 Dimitri Komatitsch, University of Pau, and Qinya Liu, Caltech, July 2005:
!  serial version, regular mesh, adjoint and kernel calculations, ParaView support
!
! MPI v. 1.2 Min Chen and Dimitri Komatitsch, Caltech, July 2004:
!  full anisotropy, volume movie
!
! MPI v. 1.1 Dimitri Komatitsch, Caltech, October 2002: Zhu's Moho map, scaling
!  of Vs with depth, Hauksson's regional model, attenuation, oceans, movies
!
! MPI v. 1.0 Dimitri Komatitsch, Caltech, May 2002: first MPI version based on global code


  include 'version.fh'

  ! local parameters
  ! timing
  double precision, external :: wtime

  ! MPI initialization
  call init_mpi()

  ! sizeprocs returns number of processes started (should be equal to NPROC).
  ! myrank is the rank of each process, between 0 and NPROC-1.
  ! as usual in MPI, process 0 is in charge of coordinating everything
  ! and also takes care of the main output
  call world_size(sizeprocs)
  call world_rank(myrank)

  ! open main output file, only written to by process 0
  if (myrank == 0 .and. IMAIN /= ISTANDARD_OUTPUT) &
    open(unit=IMAIN,file=trim(OUTPUT_FILES)//'/output_generate_databases.txt',status='unknown')

  ! get MPI starting time
  time_start = wtime()

  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '*****************************************'
    write(IMAIN,*) '*** Specfem3D MPI database generation ***'
    write(IMAIN,*) '*****************************************'
    write(IMAIN,*)
    write(IMAIN,*) 'Running Git package version of the code: ', git_package_version
    write(IMAIN,*) 'which is Git ', git_commit_version
    write(IMAIN,*) 'dating ', git_date_version
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! read the parameter file
  call read_parameters()

  ! reads topography and bathymetry file
  call read_topography()

  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '************************************'
    write(IMAIN,*) 'reading partition files in the model'
    write(IMAIN,*) '************************************'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! Initialize ADIOS I/O
  if (ADIOS_ENABLED) then
    call adios_setup()
  endif

  ! reads Databases files
  if (ADIOS_FOR_DATABASES) then
    call read_partition_files_adios()
  else
    call read_partition_files()
  endif

  ! external mesh creation
  call setup_mesh()

  ! finalize mesher
  call finalize_databases()

  if (ADIOS_ENABLED) then
    call adios_cleanup()
  endif

  ! MPI finish
  call finalize_mpi()

  end program xgenerate_databases

