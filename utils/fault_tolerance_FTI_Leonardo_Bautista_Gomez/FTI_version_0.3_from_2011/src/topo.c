/*
 * =====================================================================================
 *
 *       Filename:  topo.c
 *
 *    Description:  MPI topology funtions for FTI library
 *
 *        Version:  1.0
 *        Created:  01/14/2011 09:41:11 PM JST
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Leonardo BAUTISTA GOMEZ (leobago@matsulab.is.titech.ac.jp),
 *        Company:  Tokyo Institute of Technology
 *
 * =====================================================================================
 */


#include "fti.h"


// This function write the node topology in a file
// so it can be restored after failure
void FTI_CreateTopo(char *nameList) {

    int i;
    char hname[FTI_BUFS], mfn[FTI_BUFS];
    FILE* exp_file;
    sprintf(mfn,"%s/Topology.fti", FTI_Mdir);

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

    // Extract a reference to the function "setTopo" from the global dictionary
    expression = PyDict_GetItemString(global_dict, "setTopo");

    for (i = 0; i < FTI_Nbnd; i++) {
        strncpy(hname,nameList+(i*FTI_BUFS),FTI_BUFS);
        obj = PyObject_CallFunction(expression, "sis", mfn, i, hname);
    }

    Py_Finalize();
}


// This function reorder the nodes accordingly
// with the previous topology before failure
void FTI_ReorderNodes(int *nodeList, char *nameList) {

    int i, j, *nl, *old, *new;
    char *nmls, mfn[FTI_BUFS], *tmp;
    FILE* exp_file;
    sprintf(mfn,"%s/Topology.fti", FTI_Mdir);
    nl = talloc(int, FTI_Nbpr);
    nmls = talloc(char, FTI_Nbnd*FTI_BUFS);
    tmp = talloc(char, FTI_BUFS);

    // Getting configuration
    Py_Initialize();
    PyObject *obj, *main_module, *global_dict, *expression;

    // Open and execute the Python file
    exp_file = fopen("src/conf.py", "r");
    PyRun_SimpleFile(exp_file, "src/conf.py");

    // Get a reference to the main module and global dictionary
    main_module = PyImport_AddModule("__main__");
    global_dict = PyModule_GetDict(main_module);

    // Extract a reference to the function "setTopo" from the global dictionary
    expression = PyDict_GetItemString(global_dict, "getTopo");

    // Initializing transformation vectors
    old = talloc(int, FTI_Nbnd);
    new = talloc(int, FTI_Nbnd);
    for (i = 0; i < FTI_Nbnd; i++) {
        old[i] = -1;
        new[i] = -1;
    }

    // Get the old order of nodes
    for (i = 0; i < FTI_Nbnd; i++) {
        obj = PyObject_CallFunction(expression, "si", mfn, i);
        tmp = PyString_AsString(PyList_GetItem(obj, 0));
        strncpy(nmls+(i*FTI_BUFS),tmp,FTI_BUFS);
        for (j = 0; j < FTI_Nbnd; j++) {
            if (strncmp(tmp,nameList+(j*FTI_BUFS),FTI_BUFS) == 0) {
                old[j] = i;
                new[i] = j;
                j = FTI_Nbnd;
            }
        }
    }

    // Introducing missing nodes
    j = 0;
    for (i = 0; i < FTI_Nbnd; i++) {
        if (new[i] == -1) {
            while(old[j] != -1) {
                j++;
            }
            old[j] = i;
            new[i] = j;
            j++;
        }
    }

    // Copying nodeList in nl
    for (i = 0; i < FTI_Nbpr; i++) {
        nl[i] = nodeList[i];
    }

    // Creating the new nodeList with the old order
    for (i = 0; i < FTI_Nbnd; i++) {
        for (j = 0; j < FTI_Ndsz; j++) {
            nodeList[(i*FTI_Ndsz)+j] = nl[(new[i]*FTI_Ndsz)+j];
        }
    }

    Py_Finalize();

}


// This function creates the MPI topology
// either for first launch jobs or restarted jobs
void FTI_Topology() {

    MPI_Status status;
    MPI_Group newGroup, origGroup;
    int nn, pos, found, c1, c2, p, i, j;

    // Memory allocation
    char *lhn = talloc(char, FTI_BUFS*FTI_Nbpr);
    char *nameList = talloc(char, FTI_Nbnd * FTI_BUFS);
    int *headProcList = talloc(int, (FTI_Nbnd*FTI_Nbhe));
    int *userProcList = talloc(int, FTI_Nbpr-(FTI_Nbnd*FTI_Nbhe));
    int *nodeList = talloc(int, FTI_Nbnd * FTI_Ndsz);
    int *group = talloc(int, FTI_Grsz);

    char hname[FTI_BUFS];
    gethostname(lhn+(FTI_Rank*FTI_BUFS),FTI_BUFS);

    // Distributing host names
    MPI_Allgather(lhn+(FTI_Rank*FTI_BUFS), FTI_BUFS, MPI_CHAR, lhn, FTI_BUFS, MPI_CHAR, MPI_COMM_WORLD);

    // Preparing the node list
    for (i = 0; i < FTI_Nbpr; i++) {
        nodeList[i] = -1;
    }

    nn = 0;
    // Creating the node list
    for (i = 0; i < FTI_Nbpr; i++) {

        strncpy(hname,lhn+(i*FTI_BUFS),FTI_BUFS);
        pos = 0;
        found = 0;

        // Search the name in the list
        while ((pos < nn) && (found == 0)) {
            // If we find it, we add the process id in the group list
            if (strncmp(&(nameList[pos*FTI_BUFS]),hname,FTI_BUFS) == 0) {
                found = 1;
                p = pos*FTI_Ndsz;
                while (p < pos*FTI_Ndsz + FTI_Ndsz) {
                    if (nodeList[p] == -1) {
                        nodeList[p] = i;
                        p = pos*FTI_Ndsz + FTI_Ndsz;
                    } else {
                        p++;
                    }
                }
            } else {
                // Else we compare with the next name in the list
                pos++;
            }
        }

        // If we didnt find it, we add it to the end of the list
        if (found == 0) {
            strncpy(&(nameList[pos*FTI_BUFS]),hname,FTI_BUFS);
            nodeList[pos*FTI_Ndsz] = i;
            nn++;
        }
    }

    // Reorder the nodes in case of recovering
    if (FTI_Fail) {
        FTI_ReorderNodes(nodeList, nameList);
    }

    c1 = 0;
    c2 = 0;
    pos = 0;
    nn = 0;
    // Creating the head and user lists
    for (i = 0; i < FTI_Nbpr; i++) {
        if (i % FTI_Ndsz == 0) {
            if (FTI_Rank == nodeList[i]) {
                FTI_Idnd = nn;
                FTI_Idhe = nodeList[i];
                FTI_Head = 1;
            }
            headProcList[c1] = nodeList[i];
            c1++;
            nn++;
        } else {
            if (FTI_Rank == nodeList[i]) {
                FTI_Idnd = nn-1;
                FTI_Idhe = nodeList[(nn-1)*FTI_Ndsz];
                FTI_Head = 0;
            }
            userProcList[c2] = nodeList[i];
            c2++;
        }
    }

    if (FTI_Rank == 0) {
        FTI_CreateTopo(nameList);
    }

    // Creating the old group
    MPI_Comm_group(MPI_COMM_WORLD, &origGroup);

    if (FTI_Head) {
        // Receiving the list of processes and creating the group
        MPI_Group_incl(origGroup, FTI_Nbnd*FTI_Nbhe, headProcList, &newGroup);

        // The head get the body of the node
        FTI_Body = talloc(int, FTI_Ndsz-FTI_Nbhe);
        for (j = FTI_Nbhe; j < FTI_Ndsz; j++) {
            int src = nodeList[(FTI_Idnd*FTI_Ndsz)+j];
            MPI_Recv(&found, 1, MPI_INT, src, FTI_Mtag, MPI_COMM_WORLD, &status);
            if (found == src) {
                FTI_Body[j-FTI_Nbhe] = src;
            }
        }

        // Creating the group and communicator
        MPI_Comm_create(MPI_COMM_WORLD, newGroup, &FTI_COMM_WORLD);

        MPI_Comm_rank(FTI_COMM_WORLD, &c1);
        c2 = (c1/FTI_Grsz)*FTI_Grsz;
        for (i = 0; i < FTI_Grsz; i++) {
            group[i] = c2 + i;
        }

        MPI_Comm_group(FTI_COMM_WORLD, &origGroup);
        MPI_Group_incl(origGroup, FTI_Grsz, group, &newGroup);
        MPI_Comm_create(FTI_COMM_WORLD, newGroup, &FTI_Cenc);
        MPI_Group_rank (newGroup, &FTI_Idgr);
        FTI_Idsc = c1 / FTI_Grsz;

    } else {
        MPI_Group_incl(origGroup, FTI_Nbpr-(FTI_Nbnd*FTI_Nbhe), userProcList, &newGroup);
        // Communication with the head
        MPI_Send(&(FTI_Rank), 1, MPI_INT, FTI_Idhe, FTI_Mtag, MPI_COMM_WORLD);
        // Creating the group and communicator
        MPI_Comm_create(MPI_COMM_WORLD, newGroup, &FTI_COMM_WORLD);
    }

    // Free memory
    free(headProcList);
    free(userProcList);
    free(nameList);
    free(nodeList);
    free(group);
    free(lhn);

    MPI_Group_free(&origGroup);
    MPI_Group_free(&newGroup);

}



