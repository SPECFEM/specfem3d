--------------------
PaToH - version 3.2
--------------------

downloadable from:
https://www.cc.gatech.edu/~umit/software.html
(accessed Jan, 2019)


Setup:

1. unzip PaToH package for your OS in this folder:
   > tar -xvzf patoh-***-x86_64.tar

2. link folder with a symbolic link:
   > cd ~/SPECFEM3D/external_libs
   > ln -s patoh-3.2/build/***-x86_64/ patoh

3. edit and uncomment the following Makefile entries in root directory SPECFEM3D/ to add PATOH support:
   ..
   PATOH_INC = $(FC_DEFINE)USE_PATOH -I./external_libs/patoh/
   PATOH_LIBS = -L./external_libs/patoh/ -lpatoh -lconfig++ -lstdc++
   ..
   PART_FLAGS += $(PATOH_INC)
   PART_LIBS  += $(PATOH_LIBS)

4. re-compile decompose_mesh tool:
   > make xdecompose_mesh



requirements:

  - libconfig++
  Make sure to have libconfig.h++ installed (http://www.hyperrealm.com/libconfig/) in your system.


