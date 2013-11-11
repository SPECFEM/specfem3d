GEOCUBIT

GeoCubit is a Python library wrapping around the Cubit Python Interface.

It aims to facilitate the meshing process in some common problems in seismic wave propagation. In particular, it is focused on the meshing requests of Specfem3D and it is helpful for such tedious tasks as:

• Creation of geophysical surfaces and volumes (ex. topography).

• Mesh of layered volumes with hexahedral.

• Creation of an anisotropic mesh suitable for cases where some alluvial basins (or slow velocity zones) are present.

• It can be used as a serial or parallel process. The parallel meshing capabilities are fun- damental for large geophysical problems (ex. mesh of Southern California using SRTM topography).

GeoCubit can be used inside the graphical interface of Cubit (i.e. as a Python object in the script tab) or as unix command.

see Cubit here [cubit.sandia.gov]
