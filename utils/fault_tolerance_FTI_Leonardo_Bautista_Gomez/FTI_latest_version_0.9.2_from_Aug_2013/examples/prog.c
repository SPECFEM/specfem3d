/*
 * =====================================================================================
 *
 * Dummy program to show how to use FTI library
 *
 * =====================================================================================
 */

#include "fti.h"
#include <stdio.h>
#include <stdlib.h>


// Dummy computation function with collective communications
float compute(float *grid, unsigned long gridSize, int ind) {
    unsigned long j, start = rand()%(gridSize/2);
    float res;

    for (j = start; j < gridSize/4; j++) {
        res = ((rand() % 10) + 1)/10;
        grid[j] = grid[j-1] + ((rand()%3)-1)*res;
    }

    MPI_Allreduce(grid+ind, &res, 1, MPI_FLOAT, MPI_SUM, FTI_COMM_WORLD);
    //usleep(100000);
    return res;
}


int main(int argc, char **argv)
{

    // MPI initialization
    MPI_Init(&argc, &argv);

    // FTI initialization
    FTI_Init(argv[3], MPI_COMM_WORLD);

    int rank, nbProc, steps, i;
    unsigned long gridSize, j;
    float *grid, ranDif, res = 1;
    double t;

    // Rank and size of FTI communicator
    MPI_Comm_rank(FTI_COMM_WORLD, &rank);
    MPI_Comm_size(FTI_COMM_WORLD, &nbProc);
    srand(time(NULL)*rank);

    // Get parameters
    gridSize = atoi(argv[1])*262150; // 262150 x 4 bytes = 1MB
    steps = atoi(argv[2]);
    grid = malloc(gridSize * sizeof(float));
    float memsize = (float) gridSize*sizeof(float);
    memsize = memsize /(1024*1024);
    if (rank == 0) printf("Gridsize : %ldK and memory size : %fMB \n", gridSize/1000, memsize);

    // Initialize grid
    if (rank == 0) printf("Initializing vector...\n");
    grid[0] = 300.0;
    for (j = 1; j < gridSize; j++) {
        ranDif = ((rand() % 10) + 1)/10;
        grid[j] = grid[j-1] + ((rand()%3)-1)*ranDif;
    }

    // Defining which variables to protect
    FTI_Protect(0, &i, 1*sizeof(int));
    FTI_Protect(1, &res, 1*sizeof(float));
    FTI_Protect(2, grid, gridSize*sizeof(float));

    // Main loop
    t = MPI_Wtime();
    for (i = 0; i < steps; i++) {

        // FTI snapshot
        FTI_Snapshot();

        // Computation
        res = compute(grid, gridSize, rank);

        // Printing time steps
        if (rank == 0 && i%10 == 0) printf("Step : %d and res : %f \n", i, res);

    }
    // Final result
    if (rank == 0) printf("Final result = %f computed in %f seconds\n", res, MPI_Wtime() - t);

    // Free memory
    free(grid);

    FTI_Finalize();
    MPI_Finalize();

    return 0;
}
