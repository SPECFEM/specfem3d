SPECFEM3D input data files to reproduce the pressure and seismic waveforms presented in Figures 4 and 8 of the paper

Bishop, J. W., Fee, D., Modrak, R., Tape, C., & Kim, K. (2022). Spectral element modeling of acoustic to seismic coupling over topography. Journal of Geophysical Research: Solid Earth, 127, e2021JB023142.
https://doi.org/10.1029/2021JB023142

Figure 4 depicts an acoustic wave (from a source in the acoustic medium) coupling into a planar elastic halfspace. Figure 8 shows an acoustic wave (from a source in the acoustic medium) coupling into an elastic halfspace with complex, realistic terrain. These examples were initially created for infrasound propagation, so the acoustic media properties are representative of a homogeneous atmosphere.

Note: The STATIONS files in the DATA folders only contain the reference receiver locations. The STATIONS_LARGE files contain the station information for all the receivers in Figures 4 and 8, respectively.

Note: This example is designed to generate results close to the aforementioned paper, but the steps here might not be exactly the same as that were followed originally.

Step-by-step tutorial:

1. Check that all software is available (or that modules are loaded):

    - openmpi: > which mpirun

    - optional:
        - gnuplot: > which gnuplot
        - xmgrace: > which xmgrace

2. Configure package:

   - From the SPECFEM3D root directory SPECFEM3D/
     configure the package, e.g. using intel's ifort compiler:
     > cd SPECFEM3D
     >
     > ./configure FC=ifort MPIFC=mpif90

     If successful, this will generate the files in setup/:
     config.h, constants.h, precision.h, among others

   - Copy run scripts for your particular cluster from utils/Cluster/ into SPECFEM3D/, e.g., go_mesher_slurm.bash, go_generate_databases_slurm.bash, and go_solver_slurm.bash

     > cp utils/Cluster/slurm/*.bash .

     Note: you may need to adjust the commands for your particular cluster (e.g., -q or -l)
     Note: in this example you do not need go_decomposer_slurm.bash

3. Copy input files from examples directory into SPECFEM3D/DATA/

     Note: this may not be needed if the files are already symbolically linked.
     Note: This line is for the Figure 4 example. Replace `DATA_FIG4` with `DATA_FIG8` to run the example with terrain.

     > cp -r EXAMPLES/seismoacoustic_bishop2022/DATA_FIG4/ DATA

4. Make mesh internally:

   - compile internal mesher in directory SPECFEM3D/:
     > make xmeshfem3D

   - modify the run mesher script to use 72 cores, e.g.
     > #SBATCH --ntasks=4

     to:

     > #SBATCH --ntasks=72

   - run mesher:
     > sbatch go_mesher_slurm.bash

     Note: this script will need to be tailored to your cluster, e.g.,
     > bsub < go_mesher_lsf.bash
     >
     > qsub go_mesher_pbs.bash

     This command creates mesh partitions "proc000***_Database" in directory DATABASES_default/

     Note 1: the program xmeshfem3D is a parallel program (runs on $NPROC$ cores).

5. Generate databases:

   - compile generate_databases in directory SPECFEM3D/:
     > make xgenerate_databases

   - modify the run generate databases script to use 72 cores

   - submit job script:
     > sbatch go_generate_databases_slurm.bash

     Note: this script will need to be tailored to your cluster, e.g.,
     > bsub < go_generate_databases_lsf.bash
     >
     > qsub go_generate_databases_pbs.bash

     This command will create binary mesh files, e.g. "proc000***_external_mesh.bin"
     in directory DATABASES_default/


6. Run simulation:

   - compile specfem3D in directory SPECFEM3D/:
     > make xspecfem3D

   - modify the run solver script to use 72 cores

   - submit job script:
     > sbatch go_solver_slurm.bash

     Note: this script will need to be tailored to your cluster, e.g.,
     > bsub < go_solver_lsf.bash
     >
     > qsub go_solver_pbs.bash

     Note: This simulation runs on 72 processes and will take approximately 2 hours and 56 minutes (halfspace example) or 3 hours 29 minutes (topography example).
           You can track the progress with the timestamp files
           generated in OUTPUT_FILES/

   - when the job is complete, you should have 6 velocity seismogram files and 2 pressure seismogram files in the directory OUTPUT_FILES. There will also be 330 timestamp****** files (halfspace example) or 67 timestamp****** files (topography example).
     > ls OUTPUT_FILES/\*semv | wc -l  
     > ls OUTPUT_FILES/\*semp | wc -l  
     > ls OUTPUT_FILES/timestamp\* | wc -l  

7. Check with the four reference seismograms in SPECFEM3D/EXAMPLES/seismoacoustic_bishop2022/REF_SEIS_FIG[48]/

     Note: For the halfspace simulation (Figure 4; directory REF_SEIS_FIG4), compare receivers
     1. JP.STN83P.FXP.semp
     2. JP.STN83S.FXX.semv
     3. JP.STN83S.FXY.semv
     4. JP.STN83S.FXZ.semv

     Note: For the simulation with topography (Figure 8; directory REF_SEIS_FIG8), compare the receivers
     1. JP.STN939.FXP.semp
     2. JP.STN939S.FXX.semv
     3. JP.STN939S.FXY.semv
     4. JP.STN939S.FXZ.semv

   - Example: from SPECFEM3D/, quick viewing using xmgrace (if available):

     > xmgrace EXAMPLES/seismoacoustic_bishop2022/REF_SEIS_FIG4/JP.STN83S.FXZ.semv &  
     > xmgrace OUTPUT_FILES/JP.STN83S.FXZ.semv &