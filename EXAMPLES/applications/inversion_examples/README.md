# Inversion examples

Folder for examples of full-waveform inversion (FWI) approaches.


The examples in this folder will use very simple workflows, based mostly on the fortran tools provided in this SPECFEM3D package.


For more serious approaches, please consider additional external packages built to run FWI with SPECFEM solvers.
A recommended landing page which hosts the latest inversion tools is [adjTomo](https://github.com/adjtomo).

Training material can be [found here](https://specfem.org/training).

Further links to inversion frameworks:
  - [Pyatoa](https://github.com/adjtomo/pyatoa)
  - [Seisflows](https://github.com/adjtomo/seisflows)

and somewhat older frameworks which at some point worked with SPECFEM3D:
  - [Lasif](http://lasif.net)
  - [ASKI](http://www.gmg.ruhr-uni-bochum.de/geophysik/seismology/aski.html)
  - [specfem_FWI_workflow](https://github.com/alanschiemenz/specfem_FWI_workflow)


More specific tools:
* for adjoint inversions:
  - [pytomo3d](https://github.com/computational-seismology/pytomo3d)

* time-window selection:
  - [pyflex](https://github.com/adjtomo/pyflex)
  - FLEXWIN: original version from [FLEXWIN](https://github.com/geodynamics/flexwin)
             provided as submodule in utils/ADJOINT_TOMOGRAPHY_TOOLS/

* adjoint source creation:
  - [pyadjoint](https://github.com/adjtomo/pyadjoint)
  - measure_adj: provided as submodule in utils/ADJOINT_TOMOGRAPHY_TOOLS/

* moment-tensor inversions:
  - [pycmt3d](https://github.com/wjlei1990/pycmt3d)
    provided in src/inverse_problem_for_source/
  - CMT3D: original version from [GRD_CMT3D](https://github.com/UTCompSeismo/GRD_CMT3D) 
           


Please consider contributing your own examples.

