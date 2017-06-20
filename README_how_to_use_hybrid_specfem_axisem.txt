INSTALLATION of Hybrid Specfem / axisem


  1-- go to specfem3d_Git_devel and run configure script

     ./configure FC=ifort MPIFC=mpif90 CC=icc


  2-- add following line : (just before FC = ...)

     DEBUG_COUPLED_FLAG=-DDEBUG_COUPLED -DUSE_VTK_INSTEAD_OF_MESH

     in Makefile file created by configure script


  3-- Do 2 modifications in

         file ./setup/constants.h line 191 :
         logical, parameter :: USE_SOURCES_RECEIVERS_Z = .true.

         file ./src/meshfem3D/save_databases.F90 line 112
         logical, parameter :: SAVE_MESH_AS_CUBIT = .true.

    **** NOTE THAT configure script wil discart your changes in Makefile and in constants.h when your run it
    **** you need to do all modification after each run of configure script (configure script ned to be run once only).


  4-- Now compile the code for coupled mode :

      in specfem3d_Git_devel directory

           make all (or make all -j 8 )


    **** NOTE when doing any modifications in files add_to*.*90 you must compile after doing a make clean

         make clean
   make all (or make all -J 8)

    ****** you can use 'make realclean' to clean also scotch directory.



  5-- compile tools for coupling Axisem/Specfem

     go to ./DSM_FOR_SPECFEM3D/EXTERNAL_CODES_coupled_with_SPECFEM3D/AxiSEM_for_SPECFEM3D/UTILS_COUPLING_SpecFEM
     and type

            make all

     **** NOTE the compilers are defined in config.h file here you must apdate to yours (by default it's intel)


  6-- installing the matlab gui

       launch matlab in DSM_FOR_SPECFEM3D/inverse_problem/matlab_gui directory and run

         INSTALL_MATLAB_GUI_FOR_SPECFEM


  7-- install m_map tools for displaying maps in gui

       go to /DSM_FOR_SPECFEM3D/inverse_problem/matlab_gui/m_map/m_map
       and set :

       PATHNAME=' **** your path to DSM_FOR_SPECFEM **** /inverse_problem/matlab_gui/m_map/m_map/';

    *** NOTE this is the procedure needed by m_map


  7-- using matlab gui : launch matlab and add the path to gui tools :

          Hybrid_Axisem_Specfem


