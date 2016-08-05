/**
 *  @file   fti.h
 *  @author Leonardo A. Bautista Gomez (leobago@gmail.com)
 *  @date   July, 2013
 *  @brief  Header file for the FTI library.
 */

#ifndef  _FTI_H
#define  _FTI_H

#include <sys/stat.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <time.h>
#include <errno.h>

#include "mpi.h"

#include "iniparser.h"
#include "galois.h"
#include "jerasure.h"


/*---------------------------------------------------------------------------
                                  Defines
---------------------------------------------------------------------------*/


/** Malloc macro.                                                          */
#define talloc(type, num) (type *) malloc(sizeof(type)*(num))

/** Standard size of buffer and mas node size.                             */
#define FTI_BUFS    256
/** Word size used during RS encoding.                                     */
#define FTI_WORD    16
/** Token returned if a FTI function succeeds.                             */
#define FTI_SCES    0
/** Token returned if a FTI function fails.                                */
#define FTI_NSCS    -1

/** Verbosity level to print only errors.                                  */
#define FTI_EROR    3
/** Verbosity level to print main information.                             */
#define FTI_INFO    2
/** Verbosity level to print debug messages.                               */
#define FTI_DBUG    1

/** Token for checkpoint Baseline.                                         */
#define FTI_BASE    990
/** Token for checkpoint Level 1.                                          */
#define FTI_CKTW    991
/** Token for checkpoint Level 2.                                          */
#define FTI_XORW    992
/** Token for checkpoint Level 3.                                          */
#define FTI_RSEW    993
/** Token for checkpoint Level 4.                                          */
#define FTI_PFSW    994
/** Token for end of the execution.                                        */
#define FTI_ENDW    995


/*---------------------------------------------------------------------------
                                  New types
---------------------------------------------------------------------------*/


/*-------------------------------------------------------------------------*/
/** @typedef    FTIT_dataset
    @brief      Dataset metadata.

    This type stores the metadata related with a dataset.
*/
/*-------------------------------------------------------------------------*/
typedef struct FTIT_dataset {           /** Dataset metadata.              */
    int             id;                 /** ID to search/update dataset.   */
    void            *ptr;               /** Pointer to the dataset.        */
    long            size;               /** Size in bytes of the dataset.  */
} FTIT_dataset;

/*-------------------------------------------------------------------------*/
/** @typedef    FTIT_execution
    @brief      Execution metadata

    This type stores all the dynamic metadata related to the current execution
*/
/*-------------------------------------------------------------------------*/
typedef struct FTIT_execution {         /** Execution metadata.            */
    char            id[FTI_BUFS];       /** Execution ID.                  */
    char            ckptFile[FTI_BUFS]; /** Checkpoint file name.          */
    int             ckpt;               /** Checkpoint flag.               */
    int             reco;               /** Recovery flag.                 */
    int             ckptLvel;           /** Checkpoint level.              */
    int             ckptIntv;           /** Ckpt. interval in minutes.     */
    int             lastCkptLvel;       /** Last checkpoint level.         */
    int             wasLastOffline;     /** TRUE if last ckpt. offline.    */
    double          iterTime;           /** Current wall time.             */
    double          lastIterTime;       /** Time spent in the last iter.   */
    double          meanIterTime;       /** Mean iteration time.           */
    double          globMeanIter;       /** Global mean iteration time.    */
    double          totalIterTime;      /** Total main loop time spent.    */
    unsigned int    syncIter;           /** To check mean iter. time.      */
    unsigned int    ckptCnt;            /** Checkpoint loop counter.       */
    unsigned int    ckptID;             /** Checkpoint ID.                 */
    unsigned int    ckptNext;           /** Iteration for next checkpoint. */
    unsigned int    ckptLast;           /** Iteration for last checkpoint. */
    unsigned int    nbVar;              /** Number of protected variables. */
    MPI_Comm        globalComm;         /** Global communicator.           */
    MPI_Comm        groupComm;          /** Group communicator.            */
} FTIT_execution;

/*-------------------------------------------------------------------------*/
/** @typedef    FTIT_configuration
    @brief      Configuration metadata.

    This type stores the general configuration metadata.
 */
/*-------------------------------------------------------------------------*/
typedef struct FTIT_configuration {     /** Configuration metadata.        */
    char            cfgFile[FTI_BUFS];  /** Configuration file name.       */
    int             saveLastCkpt;       /** TRUE to save last checkpoint.  */
    int             verbosity;          /** Verbosity level.               */
    int             blockSize;          /** Communication block size.      */
    int             tag;                /** Tag for MPI messages in FTI.   */
    int             test;               /** TRUE if local test.            */
    int             l3WordSize;         /** RS encoding word size.         */
    char            localDir[FTI_BUFS]; /** Local directory.               */
    char            glbalDir[FTI_BUFS]; /** Global directory.              */
    char            metadDir[FTI_BUFS]; /** Metadata directory.            */
    char            lTmpDir[FTI_BUFS];  /** Local temporary directory.     */
    char            gTmpDir[FTI_BUFS];  /** Global temporary directory.    */
    char            mTmpDir[FTI_BUFS];  /** Metadata temporary directory.  */
} FTIT_configuration;

/*-------------------------------------------------------------------------*/
/** @typedef    FTIT_topology
    @brief      Topology metadata.

    This type stores the topology metadata.
 */
/*-------------------------------------------------------------------------*/
typedef struct FTIT_topology {          /** Topology metadata.             */
    int             nbProc;             /** Total global number of proc.   */
    int             nbNodes;            /** Total global number of nodes.  */
    int             myRank;             /** My rank on the global comm.    */
    int             splitRank;          /** My rank on the FTI comm.       */
    int             nodeSize;           /** Total number of pro. per node. */
    int             nbHeads;            /** Number of FTI proc. per node.  */
    int             nbApprocs;          /** Number of app. proc. per node. */
    int             groupSize;          /** Group size for L2 and L3.      */
    int             sectorID;           /** Sector ID in the system.       */
    int             nodeID;             /** Node ID in the system.         */
    int             groupID;            /** Group ID in the node.          */
    int             amIaHead;           /** TRUE if FTI process.           */
    int             headRank;           /** Rank of the head in this node. */
    int             groupRank;          /** My rank in the group comm      */
    int             right;              /** Proc. on the right of the ring.*/
    int             left;               /** Proc. on the left of the ring. */
    int             body[FTI_BUFS];     /** List of app. proc. in the node.*/
} FTIT_topology;

/*-------------------------------------------------------------------------*/
/** @typedef    FTIT_checkpoint
    @brief      Checkpoint metadata.

    This type stores all the checkpoint metadata.
 */
/*-------------------------------------------------------------------------*/
typedef struct FTIT_checkpoint {        /** Checkpoint metadata.           */
    char            dir[FTI_BUFS];      /** Checkpoint directory.          */
    char            metaDir[FTI_BUFS];  /** Metadata directory.            */
    int             isInline;           /** TRUE if work is inline.        */
    int             ckptIntv;           /** Checkpoint interval.           */
    int             lastCkpt;           /** ID of last checkpoint.         */
} FTIT_checkpoint;


/*---------------------------------------------------------------------------
                                  Global variables
---------------------------------------------------------------------------*/


/** MPI communicator that splits the global one into app and FTI appart.   */
MPI_Comm            FTI_COMM_WORLD;
/** Topology of the system.                                                */
FTIT_topology       FTI_Topo;
/** Dynamic information for this execution.                                */
FTIT_execution      FTI_Exec;
/** Checkpoint information for all levels of checkpoint.                   */
FTIT_checkpoint     FTI_Ckpt[5];
/** General configuration information used by FTI.                         */
FTIT_configuration  FTI_Conf;
/** Array of datasets and all their internal information.                  */
FTIT_dataset        FTI_Data[FTI_BUFS];


/*---------------------------------------------------------------------------
                            FTI public functions
---------------------------------------------------------------------------*/


int FTI_Init(char *configFile, MPI_Comm globalComm);
int FTI_Protect(int id, void *ptr, long size);
int FTI_Snapshot();
int FTI_Finalize();


#endif   /* ----- #ifndef _FTI_H  ----- */

