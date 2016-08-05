/**
 *  @file   api.c
 *  @author Leonardo A. Bautista Gomez (leobago@gmail.com)
 *  @date   July, 2013
 *  @brief  API functions for the FTI library.
 */


#include "fti.h"


/*-------------------------------------------------------------------------*/
/**
    @brief      Initializes FTI.
    @param      configFile      FTI configuration file.
    @param      globalComm      Main MPI communicator of the application.
    @return     integer         FTI_SCES if successful.

    This function initialize the FTI context and prepare the heads to wait
    for checkpoints. FTI processes should never get out of this function. In
    case of restart, checkpoint files should be recovered and in place at the
    end of this function.

 **/
/*-------------------------------------------------------------------------*/
int FTI_Init(char *configFile, MPI_Comm globalComm) {
    FTI_Exec.globalComm = globalComm;
    MPI_Comm_rank(FTI_Exec.globalComm, &FTI_Topo.myRank);
    MPI_Comm_size(FTI_Exec.globalComm, &FTI_Topo.nbProc);
    snprintf(FTI_Conf.cfgFile, FTI_BUFS, "%s", configFile);
    FTI_Conf.verbosity = 1;
    FTI_COMM_WORLD = globalComm; // Temporary assignement before initialization
    if (FTI_ReadConf() == FTI_SCES) FTI_Print("Configuration loaded.", FTI_DBUG);
    if (FTI_TestConfig() == FTI_SCES) FTI_Print("Configuration tested", FTI_DBUG);
    if (FTI_Topo.myRank == 0) if (FTI_UpdateConf(1) == FTI_SCES) FTI_Print("Configuration updated.", FTI_DBUG);
    if (FTI_CreateDirs() == FTI_SCES) FTI_Print("Checkpoint and metadata directories created.", FTI_DBUG);
    if (FTI_Topology() == FTI_SCES) FTI_Print("Topology built.", FTI_DBUG);
    if (FTI_Topo.amIaHead) { // If I am a FTI dedicated process...
        if (FTI_Exec.reco) if (FTI_RecoverFiles() == FTI_SCES) FTI_Print("Ckpt. files recovered.", FTI_DBUG);
        if (FTI_Listen() == FTI_SCES) FTI_Print("Head stopped listening for notifications.", FTI_DBUG);
        FTI_Finalize();
    } else { // If I am an application process...
        if (FTI_Exec.reco) if (FTI_RecoverFiles() == FTI_SCES) FTI_Print("Ckpt. files recovered.", FTI_DBUG);
    }
    return FTI_SCES;
}


/*-------------------------------------------------------------------------*/
/**
    @brief      It stores a pointer to a variable that needs to be protected.
    @param      id              ID for searches and update.
    @param      ptr             Pointer to the data structure.
    @param      size            Size of the data structure in bytes.
    @return     integer         FTI_SCES if successful.

    This function stores a pointer to a data structure and its size. It only
    takes into accounts calls for this function that are done before the
    first FTI_Snapshot call. Any later call to this function will be ignored.
    This list of structures is the data that will be stored during a
    checkpoint and loaded during a recovery.

 **/
/*-------------------------------------------------------------------------*/
int FTI_Protect(int id, void *ptr, long size) {
    if (FTI_Exec.ckptCnt == 0) {
        if (FTI_Exec.nbVar > FTI_BUFS) FTI_Print("Too many variables registered.", FTI_EROR);
        FTI_Data[FTI_Exec.nbVar].id = id;
        FTI_Data[FTI_Exec.nbVar].ptr = ptr;
        FTI_Data[FTI_Exec.nbVar].size = size;
        FTI_Exec.nbVar = FTI_Exec.nbVar + 1;
    }
    return FTI_SCES;
}


/*-------------------------------------------------------------------------*/
/**
    @brief      Takes an FTI snapshot or recover the data if it is a restart.
    @return     integer         FTI_SCES if successful.

    This function loads the checkpoint data from the checkpoint file in case
    of restart. Otherwise, it checks if the current iteration requires
    checkpointing, if it does it checks which checkpoint level, write the
    data in the files and it communicates with the head of the node to inform
    that a checkpoint has been taken. Checkpoint ID and counters are updated.

 **/
/*-------------------------------------------------------------------------*/
int FTI_Snapshot() {
    char fn[FTI_BUFS], str[FTI_BUFS];
    MPI_Status status;
    int buff, tres, res, i, nbProcs;
    FILE *fd;
    if (FTI_Exec.reco) { // If this is a recovery
        sprintf(fn,"%s/%s" ,FTI_Ckpt[FTI_Exec.ckptLvel].dir, FTI_Exec.ckptFile);
        sprintf(str, "Trying to load FTI checkpoint file (%s)...", fn); FTI_Print(str, FTI_DBUG);
        if (access(fn, F_OK) != 0) FTI_Print("FTI checkpoint file is NOT accesible.", FTI_EROR);
        else FTI_Print("FTI checkpoint file is accesible.", FTI_DBUG);
        fd = fopen(fn, "rb"); // Open, read and close checkpoint file
        if (fd == NULL) FTI_Print("Could not open FTI checkpoint file.", FTI_EROR);
        for(i = 0; i < FTI_Exec.nbVar; i++) fread(FTI_Data[i].ptr, 1, FTI_Data[i].size, fd);
        if (fclose(fd) == 0) FTI_Print("Checkpoint data successfully loaded.", FTI_DBUG);
        FTI_Exec.reco = 0;
        FTI_Exec.ckptID++;
        sprintf(FTI_Exec.ckptFile,"Ckpt%d-Rank%d.fti", FTI_Exec.ckptID, FTI_Topo.myRank);
        FTI_Exec.iterTime = MPI_Wtime();
    } else{ // If it is a checkpoint test
        double last = FTI_Exec.iterTime;
        FTI_Exec.iterTime = MPI_Wtime();
        if (FTI_Exec.ckptCnt > 0) {
            FTI_Exec.lastIterTime = FTI_Exec.iterTime - last;
            FTI_Exec.totalIterTime = FTI_Exec.totalIterTime + FTI_Exec.lastIterTime;
            FTI_Exec.meanIterTime = FTI_Exec.totalIterTime / FTI_Exec.ckptCnt;
            if (FTI_Exec.ckptCnt % FTI_Exec.syncIter == 0) {
                sprintf(str, "Local mean iteration time: %f seconds.", FTI_Exec.meanIterTime); FTI_Print(str, FTI_INFO);
                MPI_Allreduce(&FTI_Exec.meanIterTime, &FTI_Exec.globMeanIter, 1, MPI_DOUBLE, MPI_SUM, FTI_COMM_WORLD);
                MPI_Comm_size(FTI_COMM_WORLD, &nbProcs);
                FTI_Exec.globMeanIter = FTI_Exec.globMeanIter/nbProcs;
                sprintf(str, "Global mean iteration time: %f seconds.", FTI_Exec.globMeanIter); FTI_Print(str, FTI_INFO);
                FTI_Ckpt[1].ckptIntv = (FTI_Exec.ckptIntv*60)/FTI_Exec.globMeanIter;
                res = FTI_Exec.ckptLast + FTI_Ckpt[1].ckptIntv;
                if (res >= FTI_Exec.ckptCnt) FTI_Exec.ckptNext = res;
                FTI_Exec.ckptLast + FTI_Ckpt[1].ckptIntv;
                sprintf(str, "Checkpoint interval is %d iterations.", FTI_Ckpt[1].ckptIntv); FTI_Print(str, FTI_INFO);
                if (FTI_Exec.syncIter < (FTI_Ckpt[1].ckptIntv/2)) FTI_Exec.syncIter = FTI_Exec.syncIter * 2;
            }
        }
        if (FTI_Exec.ckptNext == FTI_Exec.ckptCnt) { // If it is time to checkpoint
            if (FTI_Exec.wasLastOffline == 1) {
                MPI_Recv(&buff, 1, MPI_INT, FTI_Topo.headRank, FTI_Conf.tag, FTI_Exec.globalComm, &status);
                if (buff == FTI_SCES) FTI_UpdateCkptInfo(1);
            }

            FTI_Exec.ckptLast = FTI_Exec.ckptNext;
            FTI_Exec.ckptNext = FTI_Exec.ckptCnt + FTI_Ckpt[1].ckptIntv;
            FTI_Exec.ckptLvel = 1; // Set checkpoint level ...
            if (FTI_Exec.ckptID % FTI_Ckpt[2].ckptIntv == 0) FTI_Exec.ckptLvel = 2;
            if (FTI_Exec.ckptID % FTI_Ckpt[3].ckptIntv == 0) FTI_Exec.ckptLvel = 3;
            if (FTI_Exec.ckptID % FTI_Ckpt[4].ckptIntv == 0) FTI_Exec.ckptLvel = 4;
            sprintf(str, "Ordering Ckpt. ID %d (L%d) at iteration %d.",
                    FTI_Exec.ckptID, FTI_Exec.ckptLvel, FTI_Exec.ckptCnt);
            FTI_Print(str, FTI_INFO);
            sprintf(str, "Next Ckpt. will be at iteration %d.", FTI_Exec.ckptNext); FTI_Print(str, FTI_INFO);
            if (FTI_Ckpt[4].isInline && FTI_Exec.ckptLvel == 4) sprintf(fn,"%s/%s",FTI_Conf.gTmpDir, FTI_Exec.ckptFile);
            else sprintf(fn,"%s/%s",FTI_Conf.lTmpDir, FTI_Exec.ckptFile);

            if (access(FTI_Conf.lTmpDir, F_OK) != 0) mkdir(FTI_Conf.lTmpDir, 0777);
            if (access(FTI_Conf.gTmpDir, F_OK) != 0) mkdir(FTI_Conf.gTmpDir, 0777);
            fd = fopen(fn, "wb"); // Open, write, flush, close and stat checkpoint file
            if (fd != NULL) FTI_Print("FTI checkpoint file successfully opened.", FTI_DBUG);
            else return FTI_NSCS;
            for(i = 0; i < FTI_Exec.nbVar; i++) fwrite(FTI_Data[i].ptr, 1, FTI_Data[i].size, fd);
            if (fflush(fd) == 0) FTI_Print("FTI checkpoint file flushed.", FTI_DBUG);
            else return FTI_NSCS;
            if (fclose(fd) == 0) FTI_Print("FTI checkpoint file closed.", FTI_DBUG);
            else return FTI_NSCS;

            res = FTI_CreateMetadata();
            MPI_Allreduce(&res, &tres, 1, MPI_INT, MPI_SUM, FTI_COMM_WORLD);
            if (tres == FTI_SCES) {
                FTI_Print("Metadata created.", FTI_DBUG);
                if (    (!FTI_Ckpt[2].isInline && FTI_Exec.ckptLvel == 2) ||
                        (!FTI_Ckpt[3].isInline && FTI_Exec.ckptLvel == 3) ||
                        (!FTI_Ckpt[4].isInline && FTI_Exec.ckptLvel == 4)  ) { // If offline work...
                    buff = FTI_BASE + FTI_Exec.ckptLvel;
                    MPI_Send(&buff, 1, MPI_INT, FTI_Topo.headRank, FTI_Conf.tag, FTI_Exec.globalComm);
                } else {
                    FTI_CkptNotified();
                }
            } else {
                FTI_Print("Error while creating metadata. Discarding checkpoint...", FTI_DBUG);
                FTI_Clean(0, FTI_Topo.groupID, FTI_Topo.myRank);
            }

        }
    }
    FTI_Exec.ckptCnt++; // Increment checkpoint loop counter
    return FTI_SCES;
}


/*-------------------------------------------------------------------------*/
/**
    @brief      It returns the current status of the recovery flag.
    @return     integer         FTI_Exec.reco

    This function returns the current status of the recovery flag.

 **/
/*-------------------------------------------------------------------------*/
int FTI_Status() {
    return FTI_Exec.reco;
}


/*-------------------------------------------------------------------------*/
/**
    @brief      It closes FTI properly on the application processes.
    @return     integer         FTI_SCES if successful.

    This function notify the FTI processes that the execution is over, frees
    some data structures and it closes. If this function is not called on the
    application processes the FTI processes will never finish (deadlock).

 **/
/*-------------------------------------------------------------------------*/
int FTI_Finalize() {
    if (!FTI_Topo.amIaHead) {
        int buff = FTI_ENDW;
        MPI_Status status;
        if (FTI_Exec.wasLastOffline == 1)
            MPI_Recv(&buff, 1, MPI_INT, FTI_Topo.headRank, FTI_Conf.tag, FTI_Exec.globalComm, &status);
        buff = FTI_ENDW;
        if (FTI_Topo.nbHeads == 1) MPI_Send(&buff, 1, MPI_INT, FTI_Topo.headRank, FTI_Conf.tag, FTI_Exec.globalComm);
        if (FTI_Conf.saveLastCkpt) {
            if (FTI_Exec.lastCkptLvel != 4) {
                if (FTI_Flush(FTI_Topo.groupID, 0) == 0) FTI_Print("Last checkpoint saved in the PFS.", FTI_DBUG);
            } else {
                rename(FTI_Ckpt[4].dir, FTI_Conf.gTmpDir); // To avoid erasing during cleaning
            }
            if (FTI_Topo.splitRank == 0) if (FTI_UpdateConf(2) == 0) FTI_Print("Configuration updated to 2.",FTI_DBUG);
            buff = 6;
        } else {
            if (FTI_Topo.splitRank == 0) if (FTI_UpdateConf(0) == 0) FTI_Print("Configuration updated to 0.",FTI_DBUG);
            buff = 5;
        }
        MPI_Barrier(FTI_Exec.globalComm); FTI_Print("Final barrier passed.", FTI_DBUG);
        if (FTI_Clean(buff, FTI_Topo.groupID, FTI_Topo.myRank) == FTI_SCES) FTI_Print("Final clean done.", FTI_DBUG);
        if (FTI_Conf.saveLastCkpt) rename(FTI_Conf.gTmpDir, FTI_Ckpt[4].dir);
        FTI_Print("FTI finalized.", FTI_DBUG);
    } else {
        MPI_Barrier(FTI_Exec.globalComm); FTI_Print("Final barrier passed.", FTI_DBUG);
        FTI_Print("FTI finalized.", FTI_DBUG);
        MPI_Finalize();
        exit(0);
    }
    return FTI_SCES;
}


