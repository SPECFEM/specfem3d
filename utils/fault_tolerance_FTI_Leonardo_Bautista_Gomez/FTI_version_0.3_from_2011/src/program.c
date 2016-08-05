/*
 * =====================================================================================
 *
 *       Filename:  program.c
 *
 *    Description:  Just for testing FTI library
 *
 *        Version:  1.0
 *        Created:  12/13/2010 04:54:13 PM JST
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Leonardo Bautista (leobago@matsulab.is.titech.ac.jp),
 *        Company:  Tokyo Tech
 *
 * =====================================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include "fti.h"


// Read local checkpoint function
void myRestore(int gridSize, int *i, float *res, float* grid, char* filename) {

    FILE *fd;
    char fn[30];
    int num;
    sprintf(fn,"%s/%s",FTI_Cdir,filename);
    fd = fopen(fn, "r");
    fread(i,sizeof(int),1,fd);
    fread(res,sizeof(float),1,fd);
    num = fread(grid,sizeof(float),gridSize,fd);
    if (num != gridSize) {
        printf("Error reading the ckpt file: %d reads\n", num);
    }
    fclose(fd);
}


// Checkpoint to local file funtion
void myCheckpoint(int gridSize, int i, float res, float* grid, char* filename) {

    FILE *fd;
    char fn[30];
    int num;
    sprintf(fn,"%s/%s",FTI_Cdir,filename);
    fd = fopen(fn, "w");
    fwrite(&i,sizeof(int),1,fd);
    fwrite(&res,sizeof(float),1,fd);
    num = fwrite(grid,sizeof(float),gridSize,fd);
    if (num != gridSize) {
        printf("Error writting the ckpt file: %d writes\n", num);
    }
    fclose(fd);
}

// Computation function
float compute(float *grid,  int gridSize, int ind) {

    int i;
    float res, red = ind;
    for(i = 0; i < gridSize; i++) {
        grid[i] = (grid[i] * (ind + 1)) / ind;
        red = (red + grid[i]);
        if (red > gridSize)
            red = red - gridSize;
    }
    red = red / (float)gridSize;
    MPI_Allreduce(&red, &res, 1, MPI_FLOAT, MPI_SUM, FTI_COMM_WORLD);
    return res;
}


int main(int argc, char **argv)
{

    // MPI initialization
    MPI_Init(&argc, &argv);

    // FTI initialization
    FTI_Init("config.fti");

    int rank, nbProc, gridSize, steps, offset, i, ex = 6000;
    float *grid, res = 1;
    char filename[25];
    double t1, t2;

    // Rank and size of FTI communicator
    MPI_Comm_rank(FTI_COMM_WORLD, &rank);
    MPI_Comm_size(FTI_COMM_WORLD, &nbProc);

    // Parameters
    gridSize = atoi(argv[1])/nbProc;
    steps = atoi(argv[2]);
    offset = gridSize * rank;
    grid = malloc(gridSize * sizeof(float));
    //printf("Rank: %d, nbProc: %d, gridZise: %d, steps: %d, offset: %d\n",rank,nbProc,gridSize,steps,offset);

    // Setup grid
    for (i = 0; i < gridSize; i++) {
        grid[i] = i * 0.1;
    }

    // Main loop
    t1 = MPI_Wtime();
    t2 = MPI_Wtime();
    for (i = 0; i < steps; i++) {

        // Recovering in case of failure
        if (FTI_Fail) {
            myRestore(gridSize, &i, &res, grid, FTI_File);
            FTI_Restarted();
        }

        // Computation
        res = compute(grid, gridSize, res+1) / (float)nbProc;

        // Printing time steps
        if (rank == 0) {
            printf("Step : %d and res : %f \n", i, res);
        }

        // Every 30 steps checkpoint
        if ((i+1)%30 == 0) {
            myCheckpoint(gridSize, i+1, res, grid, FTI_File);
            FTI_Checkpointed();
        }
    }

    // Final result
    if (rank == 0) {
        printf("Final result = %f computed in %f seconds\n", res, MPI_Wtime() - t2);
    }

    // Free memory
    free(grid);

    FTI_Finalize();
    MPI_Finalize();

    return 0;
}
