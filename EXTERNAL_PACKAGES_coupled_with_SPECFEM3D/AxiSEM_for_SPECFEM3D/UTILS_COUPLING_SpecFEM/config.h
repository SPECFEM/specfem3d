CC=icc
FC=mpif90
CCFLAGS = -O3
# for debugging: change -O3 -check nobounds to      -check all -debug -g -O0 -fp-stack-check -traceback -ftrapuv

# for full vectorization on CURIE, change -xHost to -xAVX (since the frontend node has Nehalem processors, not Sandy Bridge)

FFLAGS = -O3 -check nobounds -xAVX -ftz -assume buffered_io -assume byterecl -vec-report0 -implicitnone -warn truncated_source -warn argument_checking -warn declarations -warn alignments -warn ignore_loc -warn usage -mcmodel=large -shared-intel

#FFLAGS = -check all -debug -g -O0 -fp-stack-check -traceback -ftrapuv -implicitnone -warn truncated_source -warn argument_checking -warn declarations -warn alignments -warn ignore_loc -warn usage -mcmodel=medium -shared-intel


################################################
#   Intel ifort
################################################

# it is crucial to use -xAVX here, so that big loops that contain double precision complex numbers are vectorized;
# this means that only Intel Sandy Bridge processors can vectorize these loops, not Intel Nehalem processors,
# since they only have SSE4.2 vector instructions but not AVX (Advanced Vector Extensions).
# option "-assume buffered_io" is important especially on
# parallel file systems like SFS 3.2 / Lustre 1.8. If omitted
# I/O throughput lingers at 2.5 MB/s, with it it can increase to ~44 MB/s
# However it does not make much of a difference on NFS mounted volumes or with SFS 3.1.1 / Lustre 1.6.7.1
#FC = mpif90
#FFLAGS = -O3 -check nobounds -xAVX -ftz -assume buffered_io -assume byterecl -vec-report3 -implicitnone -warn truncated_source -warn argument_checking -warn declarations -warn alignments -warn ignore_loc -warn usage -mcmodel=medium -shared-intel

# useful for debugging:
#  change   -O3 -check nobounds     to      -check all -debug -g -O0 -fp-stack-check -traceback -ftrapuv

# change    -vec-report0      to      -vec-report3     to get a vectorization report

################################################
#   GNU gfortran
################################################

#FC = mpif90
#FFLAGS = -std=gnu -fimplicit-none -frange-check -O2 -pedantic -pedantic-errors -Waliasing -Wampersand -Wline-truncation -Wsurprising -Wunderflow -fbounds-check

################################################
#   IBM Blue Gene
################################################

# at IDRIS (France) maybe change -qarch=auto to -qarch=450d
#
# you will probably need to add " module load bgq-xl " or similar to your .bash_profile to load the compilers
#
# It could also help to put this in your .bash_profile: export XLFRTEOPTS=aggressive_array_io=yes:buffering=enable
#
# On some (but not all) IBM machines one might need to add -qsave otherwise the IBM compiler allocates the
# arrays in the stack and the code crashes if the stack size is too
# small (which is sometimes the case, but less often these days on large machines)
#
# to debug with IBM xlf, one can add this: -g -O0 -C -qddim -qfullpath -qflttrap=overflow:zerodivide:invalid:enable -qfloat=nans -qinitauto=7FBFFFFF
#
# options -qreport -qsource -qlist create a *.lst file containing detailed information about vectorization
#
#FC = mpixlf95_r
#FFLAGS = -O4 -qnostrict -qhot -qsimd=auto -qassert=contig -g -Q -qarch=auto -q64 -qfree=f90 -qsuffix=f=f90 -qsuppress=1500-036 -qreport -qsource -qlist

