/*
 * =====================================================================================
 *
 *       Filename:  tools.c
 *
 *    Description:  Initialization and configuracion functions for FTI library
 *
 *        Version:  1.0
 *        Created:  09/28/2010 03:58:30 PM JST
 *       Revision:  none
 *       Compiler:  mpicc
 *
 *         Author:  Leonardo BAUTISTA GOMEZ (leobago@matsulab.is.titech.com),
 *        Company:  Tokyo Institue of Technology
 *
 * =====================================================================================
 */


// Including FTI header
#include "fti.h"


// This function print an error on stderr and abort the execution
void FTI_Prerror(char *s)
{
    if (s != NULL) fprintf(stderr, "\n%s\n\n", s);
    MPI_Finalize();
    exit(1);
}


// This function creates the encoding matrix
void FTI_CreateMatrix() {
    int i,j;
    FTI_Mtrx = talloc(int, FTI_Grsz*FTI_Grsz);
    // Creating matrix
    for (i = 0; i < FTI_Grsz; i++) {
        for (j = 0; j < FTI_Grsz; j++) {
            FTI_Mtrx[i*FTI_Grsz+j] = galois_single_divide(1, i ^ (FTI_Grsz + j), FTI_Wdsz);
        }
    }
}


// This function reads the configuration file
void FTI_ReadConf(char *filename) {
    // Getting configuration
    Py_Initialize();
    FILE*        exp_file;
    struct stat fileStatus;
    PyObject     *obj, *main_module, * global_dict, * expression;

    // Open and execute the Python file
    exp_file = fopen("src/conf.py", "r");
    PyRun_SimpleFile(exp_file, "src/conf.py");

    // Get a reference to the main module and global dictionary
    main_module = PyImport_AddModule("__main__");
    global_dict = PyModule_GetDict(main_module);

    // Extract a reference to the function "getConf" from the global dictionary
    expression = PyDict_GetItemString(global_dict, "getConf");

    // Make a call to the function referenced by "expression"
    if(stat("config.fti", &fileStatus) != 0)
        FTI_Prerror("Error with stat on the FTI configuration file");

    obj = PyObject_CallFunction(expression, "s", filename);

    MPI_Comm_rank(MPI_COMM_WORLD, &FTI_Rank);
    MPI_Comm_size(MPI_COMM_WORLD, &FTI_Nbpr);
    FTI_Grsz = (int) PyInt_AsLong(PyList_GetItem(obj, 0));
    FTI_Bksz = (int) PyInt_AsLong(PyList_GetItem(obj, 1));
    FTI_Wdsz = (int) PyInt_AsLong(PyList_GetItem(obj, 2));
    FTI_Ndsz = (int) PyInt_AsLong(PyList_GetItem(obj, 3));
    FTI_Mtag = (int) PyInt_AsLong(PyList_GetItem(obj, 4));
    FTI_Fail = (int) PyInt_AsLong(PyList_GetItem(obj, 5));
    FTI_Nbhe = (int) PyInt_AsLong(PyList_GetItem(obj, 6));
    FTI_Cdir = PyString_AsString(PyList_GetItem(obj, 7));
    FTI_Mdir = PyString_AsString(PyList_GetItem(obj, 8));
    FTI_Nbnd = FTI_Nbpr / FTI_Ndsz;
    FTI_Nbgr = FTI_Ndsz - FTI_Nbhe;
    FTI_Ckpt = -9991;
    FTI_Endw = -9992;
    FTI_Idck = 0;
    sprintf(FTI_File,"Ckpt%d-Rank%d.fti", FTI_Idck, FTI_Rank);

    if (FTI_Grsz <= 2)
        FTI_Prerror("Wrong value for FTI_Grsz (k and m). It must be bigger than 2");
    if (FTI_Wdsz != 8 && FTI_Wdsz != 16)
        FTI_Prerror("Wrong value for m. It must be equal to 8 or 16");
    if (2*FTI_Grsz > (1 << FTI_Wdsz))
        FTI_Prerror("k + m must be <= 2 ^ w");
    if (FTI_Bksz % sizeof(long) != 0)
        FTI_Prerror("FTI_Bksz must be multiple of sizeof(long)");

    Py_Finalize();
}


// This function get the metadata in order to
// recover the lost data after a failure
void FTI_GetMeta(int *fs, int *mfs, int group) {
    // Getting configuration
    Py_Initialize();
    char mfn[FTI_BUFS], *cfn;
    FILE*        exp_file;
    struct stat fileStatus;
    PyObject     *obj, *main_module, * global_dict, * expression;

    // Open and execute the Python file
    exp_file = fopen("src/conf.py", "r");
    PyRun_SimpleFile(exp_file, "src/conf.py");

    // Get a reference to the main module and global dictionary
    main_module = PyImport_AddModule("__main__");
    global_dict = PyModule_GetDict(main_module);

    // Extract a reference to the function "getMeta" from the global dictionary
    expression = PyDict_GetItemString(global_dict, "getMeta");

    // Make a call to the function referenced by "expression"
    sprintf(mfn,"%s/sector%d-group%d.fti",FTI_Mdir, FTI_Idsc, group);
    if(stat(mfn, &fileStatus) != 0)
        FTI_Prerror("Error with stat on the FTI metadata file");

    cfn = talloc(char, FTI_BUFS);
    obj = PyObject_CallFunction(expression, "si", mfn, FTI_Idgr);
    cfn = PyString_AsString(PyList_GetItem(obj, 0));
    *fs = (int) PyInt_AsLong(PyList_GetItem(obj, 1));
    *mfs = (int) PyInt_AsLong(PyList_GetItem(obj, 2));
    strncpy(FTI_File,cfn,FTI_BUFS);
    free(cfn);

    Py_Finalize();
}


// This function writes the metadata to use in case of restart
void FTI_CreateMetadata(int *fs, int *mfs, int group) {

    int i;
    char mfn[FTI_BUFS], cfn[FTI_BUFS], *fnl;
    FILE* exp_file;
    struct stat fileStatus;

    // Getting size of files
    fnl = talloc(char, FTI_Grsz*FTI_BUFS);
    sprintf(cfn,"%s/%s",FTI_Cdir, FTI_File);
    if(stat(cfn, &fileStatus) == 0)
        fs[FTI_Idgr] = (int) fileStatus.st_size;
    else
        FTI_Prerror("Error with stat on the checkpoint file");

    // Gather all the file sizes
    sprintf(fnl+(FTI_Idgr*FTI_BUFS),"%s",FTI_File);
    MPI_Allgather(fs+FTI_Idgr, 1, MPI_INT, fs, 1, MPI_INT, FTI_Cenc);
    MPI_Allgather(fnl+(FTI_Idgr*FTI_BUFS), FTI_BUFS, MPI_CHAR, fnl, FTI_BUFS, MPI_CHAR, FTI_Cenc);

    // Write the file sizes in metadata file and look for the max
    for(i = 0; i < FTI_Grsz; i++) {
        if (fs[i] > *mfs)
            *mfs = fs[i];
    }

    if (FTI_Idgr == 0) {
        sprintf(mfn, "%s/sector%d-group%d.fti", FTI_Mdir, FTI_Idsc, group);

        // Getting configuration
        Py_Initialize();
        PyObject *obj, *main_module, *global_dict, *expression;

        // Open and execute the Python file
        exp_file = fopen("src/conf.py", "r");
        PyRun_SimpleFile(exp_file, "src/conf.py");

        // Get a reference to the main module and global dictionary
        main_module = PyImport_AddModule("__main__");
        global_dict = PyModule_GetDict(main_module);

        // Extract a reference to the function "rmMeta" from the global dictionary
        expression = PyDict_GetItemString(global_dict, "rmMeta");

        // Calling the function with one string parameter
        obj = PyObject_CallFunction(expression, "s", mfn);

        // Extract a reference to the function "setMeta" from the global dictionary
        expression = PyDict_GetItemString(global_dict, "setMeta");

        for (i = 0; i < FTI_Grsz; i++) {
            strncpy(cfn,fnl+(i*FTI_BUFS),FTI_BUFS);
            obj = PyObject_CallFunction(expression, "sisii", mfn, i, cfn, fs[i], *mfs);
        }

        Py_Finalize();
    }
    free(fnl);

}


// This functions informs FTI that the job has succesfully
// restarted and prepare the next checkpoint
void FTI_Restarted() {
    FTI_Idck++;
    sprintf(FTI_File,"Ckpt%d-Rank%d.fti", FTI_Idck, FTI_Rank);
    FTI_Fail = 0;
}


// This functions informs the head of the node that a checkpoint
// has been taken and that the encoding work can start
void FTI_Checkpointed() {
    MPI_Send(&FTI_Ckpt, 1, MPI_INT, FTI_Idhe, FTI_Mtag, MPI_COMM_WORLD);
    FTI_Idck++;
    sprintf(FTI_File,"Ckpt%d-Rank%d.fti", FTI_Idck, FTI_Rank);
}


// This functions is to inform the head of the node that
// the application has finished in order to close correctly
void FTI_Finalize() {

    if (!FTI_Head) {
        MPI_Send(&FTI_Endw, 1, MPI_INT, FTI_Idhe, FTI_Mtag, MPI_COMM_WORLD);
        //free(FTI_Endw);
        free(FTI_Body);
        free(FTI_Mtrx);
        free(FTI_Cdir);
        free(FTI_Mdir);
        free(FTI_File);
    }
}


// This function initialize the FTI context and
// prepare the heads to wait for checkpoints
void FTI_Init(char *configFile) {
    MPI_Status status;

    FTI_Cdir = talloc(char, FTI_BUFS);
    FTI_Mdir = talloc(char, FTI_BUFS);
    FTI_File = talloc(char, FTI_BUFS);

    // Reading the configuration
    FTI_ReadConf(configFile);

    // Building the topology and groups
    FTI_Topology();

    // Initializing matrix
    FTI_CreateMatrix();

    if (FTI_Head) {

        if (FTI_Fail) {
            int maxFs, fs, buf, i;
            for(i = 0; i < FTI_Nbgr; i++) {
                FTI_GetMeta(&fs, &maxFs, i);
                FTI_Decode(FTI_File, fs, maxFs);
                sscanf(FTI_File,"Ckpt%d-Rank%d.fti", &FTI_Idck, &buf);
                MPI_Send(FTI_File, FTI_BUFS, MPI_CHAR, FTI_Body[i], FTI_Mtag, MPI_COMM_WORLD);
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }

        int end = FTI_Nbgr;
        while (1) {
            int i, buf;
            for(i = 0; i < FTI_Nbgr; i++) {
                MPI_Recv(&buf, 1, MPI_INT, FTI_Body[i], FTI_Mtag, MPI_COMM_WORLD, &status);
                if (buf == FTI_Endw) {
                    end--;
                } else {
                    if (buf == FTI_Ckpt) {
                        int maxFs, *fs, src;
                        sprintf(FTI_File,"Ckpt%d-Rank%d.fti", FTI_Idck, FTI_Body[i]);
                        fs = talloc(int, FTI_Grsz);
                        // Writting metadata file and getting file size and max file size
                        FTI_CreateMetadata(fs, &maxFs, i);
                        // Starting the encoding
                        FTI_Encode(FTI_File, fs[FTI_Idgr], maxFs);
                        // Freeing memory
                        free(fs);
                    }
                }
            }
            if (end == 0) {
                break;
            }
        }

        // Freeing memory
        //free(FTI_Endw);
        free(FTI_Body);
        free(FTI_Mtrx);
        free(FTI_Cdir);
        free(FTI_Mdir);

        // Closing MPI
        MPI_Finalize();

        // Exit from program
        exit(0);
    } else {
        if (FTI_Fail) {
            int buf;
            MPI_Recv(FTI_File, FTI_BUFS, MPI_CHAR, FTI_Idhe, FTI_Mtag, MPI_COMM_WORLD, &status);
            sscanf(FTI_File,"Ckpt%d-Rank%d.fti", &FTI_Idck, &buf);
            if (buf != FTI_Rank)
                FTI_Prerror("Error with the recovered filename\n");
            MPI_Barrier(MPI_COMM_WORLD);
        }
    }
}

// Function developped by James Planck
// Copied here under GPL terms
int FTI_InvertMatrix(int *mat, int *inv, int rows, int w)
{
  int cols, i, j, k, x, rs2;
  int row_start, tmp, inverse;

  cols = rows;

  k = 0;
  for (i = 0; i < rows; i++) {
    for (j = 0; j < cols; j++) {
      inv[k] = (i == j) ? 1 : 0;
      k++;
    }
  }

  /* First -- convert into upper triangular  */
  for (i = 0; i < cols; i++) {
    row_start = cols*i;

    /* Swap rows if we ave a zero i,i element.  If we can't swap, then the
       matrix was not invertible  */

    if (mat[row_start+i] == 0) {
      for (j = i+1; j < rows && mat[cols*j+i] == 0; j++) ;
      if (j == rows) return -1;
      rs2 = j*cols;
      for (k = 0; k < cols; k++) {
        tmp = mat[row_start+k];
        mat[row_start+k] = mat[rs2+k];
        mat[rs2+k] = tmp;
        tmp = inv[row_start+k];
        inv[row_start+k] = inv[rs2+k];
        inv[rs2+k] = tmp;
      }
    }

    /* Multiply the row by 1/element i,i  */
    tmp = mat[row_start+i];
    if (tmp != 1) {
      inverse = galois_single_divide(1, tmp, w);
      for (j = 0; j < cols; j++) {
        mat[row_start+j] = galois_single_multiply(mat[row_start+j], inverse, w);
        inv[row_start+j] = galois_single_multiply(inv[row_start+j], inverse, w);
      }
    }

    /* Now for each j>i, add A_ji*Ai to Aj  */
    k = row_start+i;
    for (j = i+1; j != cols; j++) {
      k += cols;
      if (mat[k] != 0) {
        if (mat[k] == 1) {
          rs2 = cols*j;
          for (x = 0; x < cols; x++) {
            mat[rs2+x] ^= mat[row_start+x];
            inv[rs2+x] ^= inv[row_start+x];
          }
        } else {
          tmp = mat[k];
          rs2 = cols*j;
          for (x = 0; x < cols; x++) {
            mat[rs2+x] ^= galois_single_multiply(tmp, mat[row_start+x], w);
            inv[rs2+x] ^= galois_single_multiply(tmp, inv[row_start+x], w);
          }
        }
      }
    }
  }

  /* Now the matrix is upper triangular.  Start at the top and multiply down  */
  for (i = rows-1; i >= 0; i--) {
    row_start = i*cols;
    for (j = 0; j < i; j++) {
      rs2 = j*cols;
      if (mat[rs2+i] != 0) {
        tmp = mat[rs2+i];
        mat[rs2+i] = 0;
        for (k = 0; k < cols; k++) {
          inv[rs2+k] ^= galois_single_multiply(tmp, inv[row_start+k], w);
        }
      }
    }
  }
  return 0;
}

