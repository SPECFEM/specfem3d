#!/bin/bash
bsub -W 60 -n 15 -o ./%J.out -a mpich_gm -q normal -J parallel_batch.py mpirun.lsf `which mpipython.exe` /usr/bin/parallel_batch.py run_file_001
