/*
 * =====================================================================================
 *
 *       Filename:  fti.h
 *
 *    Description:  Header file of FTI library
 *
 *        Version:  1.0
 *        Created:  09/13/2010 06:13:27 PM JST
 *       Revision:  none
 *       Compiler:  mipcc
 *
 *         Author:  Leonardo BAUTISTA GOMEZ (leobago@matsulab.is.titech.ac.jp),
 *        Company:  Tokyo Institue of Technology
 *
 * =====================================================================================
 */

#ifndef  _FTI_H
#define  _FTI_H

#include <sys/stat.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#include <python2.4/Python.h>
#include "mpi.h"

#include "galois.h"

#define talloc(type, num) (type *) malloc(sizeof(type)*(num))
#define FTI_BUFS    50

extern MPI_Comm FTI_COMM_WORLD;
extern MPI_Comm FTI_Cenc;
extern int FTI_Rank;
extern int FTI_Nbpr;
extern int FTI_Nbnd;
extern int FTI_Grsz;
extern int FTI_Wdsz;
extern int FTI_Ndsz;
extern int FTI_Bksz;
extern int FTI_Idgr;
extern int FTI_Idsc;
extern int FTI_Idck;
extern int FTI_Mtag;
extern int FTI_Idnd;
extern int FTI_Idhe;
extern int FTI_Head;
extern int FTI_Nbhe;
extern int FTI_Nbgr;
extern int FTI_Fail;
extern int FTI_Endw;
extern int FTI_Ckpt;
extern int* FTI_Body;
extern int* FTI_Mtrx;
extern char* FTI_Cdir;
extern char* FTI_Mdir;
extern char* FTI_File;

MPI_Comm FTI_COMM_WORLD;
MPI_Comm FTI_Cenc;
int FTI_Rank;
int FTI_Nbpr;
int FTI_Nbnd;
int FTI_Grsz;
int FTI_Wdsz;
int FTI_Ndsz;
int FTI_Bksz;
int FTI_Idgr;
int FTI_Idsc;
int FTI_Idck;
int FTI_Mtag;
int FTI_Idnd;
int FTI_Idhe;
int FTI_Head;
int FTI_Nbhe;
int FTI_Nbgr;
int FTI_Fail;
int FTI_Endw;
int FTI_Ckpt;
int* FTI_Body;
int* FTI_Mtrx;
char* FTI_Cdir;
char* FTI_Mdir;
char* FTI_File;

// Function in file tools.c
void FTI_Prerror(char *s);
void FTI_CreateMatrix();
void FTI_ReadConf(char *filename);
void FTI_GetMeta(int *fs, int *mfs, int group);
void FTI_CreateMetadata(int *fs, int *mfs, int group);
void FTI_Restarted();
void FTI_Checkpointed();
void FTI_Finalize();
void FTI_Init(char *configFile);
int FTI_InvertMatrix(int *mat, int *inv, int rows, int w);

// Function in file topo.c
void FTI_CreateTopo(char *nameList);
void FTI_ReorderNodes(int *nodeList, char *nameList);
void FTI_Topology();

// Functions in file enc.c
void FTI_Encode(char *filename, int fs, int maxFs);

// Functions in file dec.c
void FTI_Decode(char *filename, int fs, int maxFs);

#endif   /* ----- #ifndef _FTI_H  ----- */

