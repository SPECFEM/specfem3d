/**
 *  @file   recover.c
 *  @author Leonardo A. Bautista Gomez (leobago@gmail.com)
 *  @date   July, 2013
 *  @brief  Recovery functions for the FTI library.
 */


#include "fti.h"


/*-------------------------------------------------------------------------*/
/**
    @brief      Decides wich action take depending on the restart level.
    @return     integer         FTI_SCES if successful.

    This function launchs the required action dependeing on the recovery
    level. The recovery level is detected from the checkpoint ID of the
    last checkpoint taken.

 **/
/*-------------------------------------------------------------------------*/
int FTI_RecoverFiles() {
    int     f, i, fs, maxFs, tres, id, maxid, start = 1;
    char    str[FTI_BUFS];
    if (FTI_Topo.nbHeads == 1) f = 1; else f = 0;
    while (start < 5) {
        maxid = -1; tres = -1;
        for (i = start; i < 5; i ++) {
            if (FTI_GetMeta(&fs, &maxFs, f, i) == FTI_SCES) {
                sscanf(FTI_Exec.ckptFile,"Ckpt%d-Rank%d.fti", &id, &fs);
                if (id > maxid) {
                    maxid = id;
                    FTI_Exec.ckptID = id;
                    FTI_Exec.ckptLvel = i;
                }
                FTI_Ckpt[i].lastCkpt = id;
            }
        }
        if (maxid == -1) break;
        FTI_Exec.lastCkptLvel = FTI_Exec.ckptLvel;
        if (FTI_Exec.reco == 2) FTI_Exec.ckptLvel = 4;
        if (FTI_Exec.reco == 2) FTI_Ckpt[4].lastCkpt = FTI_Exec.ckptID;
        if (FTI_Exec.ckptLvel == 4) {
            FTI_GetMeta(&fs, &maxFs, f, 1);
            sscanf(FTI_Exec.ckptFile,"Ckpt%d-Rank%d.fti", &id, &fs);
            FTI_Ckpt[1].lastCkpt = id;
            FTI_Clean(1, FTI_Topo.groupID, FTI_Topo.myRank);
            FTI_Ckpt[1].lastCkpt = FTI_Exec.ckptID;
        }
        sprintf(str, "Trying recovery with Ckpt. %d at level %d.", FTI_Exec.ckptID, FTI_Exec.ckptLvel);
        FTI_Print(str, FTI_DBUG);
        if (!FTI_Topo.amIaHead) {
            if (FTI_Exec.ckptLvel == 4) fs = FTI_RecoverL4(FTI_Topo.groupID);
            if (FTI_Exec.ckptLvel == 3) fs = FTI_RecoverL3(FTI_Topo.groupID);
            if (FTI_Exec.ckptLvel == 2) fs = FTI_RecoverL2(FTI_Topo.groupID);
            if (FTI_Exec.ckptLvel == 1) fs = FTI_RecoverL1(FTI_Topo.groupID);
            MPI_Allreduce(&fs, &tres, 1, MPI_INT, MPI_SUM, FTI_COMM_WORLD);
            if (tres ==  FTI_SCES) { FTI_Print("Recovery function finished successfully.", FTI_DBUG); break; }
            else { FTI_Print("Recovery did NOT finish correctly.", FTI_DBUG); start = FTI_Exec.ckptLvel+1; }
        }
    }
    if (tres == FTI_NSCS) FTI_Print("Impossible to recover from this failure.", FTI_EROR);
    MPI_Barrier(FTI_Exec.globalComm); sleep(1); // Global barrier and sleep for clearer output
    return FTI_SCES;
}


/*-------------------------------------------------------------------------*/
/**
    @brief      Check if a file exist and that its size is 'correct'.
    @param      fn              The ckpt. file name to check.
    @param      fs              The ckpt. file size tocheck.
    @return     integer         FTI_SCES if successful.

    This function checks whether a file exist or not and if its size is
    the expected one.

 **/
/*-------------------------------------------------------------------------*/
FTI_CheckFile(char *fn, unsigned long fs) {
    struct stat fileStatus;
    if (access(fn, F_OK) == 0) {
        if (stat(fn, &fileStatus) == 0) {
            if (fileStatus.st_size == fs) {
                return FTI_SCES;
            } else {
                return FTI_NSCS;
            }
        } else {
            return FTI_NSCS;
        }
    } else {
        return FTI_NSCS;
    }
}


/*-------------------------------------------------------------------------*/
/**
    @brief      Detects all the erasures for a particular level.
    @param      fs              The ckpt. file size for this process.
    @param      maxFs           The max. ckpt. file size in the group.
    @param      group           The group ID.
    @param      erased          The array of erasures to fill.
    @param      level           The ckpt. level to check for erasures.
    @return     integer         FTI_SCES if successful.

    This function detects all the erasures for L1, L2 and L3. It return the
    results in the erased array. The search for erasures is done at the
    three levels independently on the current recovery level.

 **/
/*-------------------------------------------------------------------------*/
int FTI_CheckErasures(unsigned long *fs, unsigned long *maxFs, int group, int *erased, int level) {
    int         buf, j;
    char        fn[FTI_BUFS];
    for (j = 0; j < FTI_Topo.groupSize*4; j++) erased[j] = 1;                   // Initialize erasures table
    if (FTI_GetMeta(fs, maxFs, group, level) == FTI_SCES) FTI_Print("Metadata obtained.", FTI_DBUG);
    else { FTI_Print("Error getting metadata.", FTI_DBUG); return FTI_NSCS; }
    switch(level) {
        case 1: {
                    sprintf(fn, "%s/%s", FTI_Ckpt[1].dir, FTI_Exec.ckptFile);
                    if (FTI_CheckFile(fn, *fs) == FTI_SCES) erased[FTI_Topo.groupRank] = 0;
                    buf = erased[FTI_Topo.groupRank];
                    MPI_Allgather(&buf, 1, MPI_INT, erased, 1, MPI_INT, FTI_Exec.groupComm);
                    break;
                }
        case 2: {
                    sprintf(fn, "%s/%s", FTI_Ckpt[2].dir, FTI_Exec.ckptFile);
                    if (FTI_CheckFile(fn, *fs) == FTI_SCES) erased[FTI_Topo.groupRank] = 0;
                    buf = erased[FTI_Topo.groupRank];
                    MPI_Allgather(&buf, 1, MPI_INT, erased, 1, MPI_INT, FTI_Exec.groupComm);
                    sscanf(FTI_Exec.ckptFile,"Ckpt%d-Rank%d.fti", &FTI_Exec.ckptID, &buf);
                    sprintf(fn,"%s/Ckpt%d-Pcof%d.fti", FTI_Ckpt[2].dir, FTI_Exec.ckptID, buf);
                    if (FTI_CheckFile(fn, *fs) == FTI_SCES) erased[FTI_Topo.groupSize+FTI_Topo.groupRank] = 0;
                    buf = erased[FTI_Topo.groupSize+FTI_Topo.groupRank];
                    MPI_Allgather(&buf, 1, MPI_INT, erased+FTI_Topo.groupSize, 1, MPI_INT, FTI_Exec.groupComm);
                    break;
                }
        case 3: {
                    sprintf(fn, "%s/%s", FTI_Ckpt[3].dir, FTI_Exec.ckptFile);
                    if (FTI_CheckFile(fn, *fs) == FTI_SCES) erased[FTI_Topo.groupRank] = 0;
                    buf = erased[FTI_Topo.groupRank];
                    MPI_Allgather(&buf, 1, MPI_INT, erased, 1, MPI_INT, FTI_Exec.groupComm);
                    sscanf(FTI_Exec.ckptFile,"Ckpt%d-Rank%d.fti", &FTI_Exec.ckptID, &buf);
                    sprintf(fn,"%s/Ckpt%d-RSed%d.fti", FTI_Ckpt[3].dir, FTI_Exec.ckptID, buf);
                    if (FTI_CheckFile(fn, *fs) == FTI_SCES) erased[FTI_Topo.groupSize+FTI_Topo.groupRank] = 0;
                    buf = erased[FTI_Topo.groupRank+FTI_Topo.groupSize];
                    MPI_Allgather(&buf, 1, MPI_INT, erased+FTI_Topo.groupSize, 1, MPI_INT, FTI_Exec.groupComm);
                    break;
                }
        case 4: {
                    sprintf(fn, "%s/%s", FTI_Ckpt[4].dir, FTI_Exec.ckptFile);
                    if (FTI_CheckFile(fn, *fs) == FTI_SCES) erased[FTI_Topo.groupRank] = 0;
                    buf = erased[FTI_Topo.groupRank];
                    MPI_Allgather(&buf, 1, MPI_INT, erased, 1, MPI_INT, FTI_Exec.groupComm);
                    break;
                }
    }
    return FTI_SCES;
}


/*-------------------------------------------------------------------------*/
/**
    @brief      Checks that all L1 ckpt. files are present.
    @param      group           The group ID.
    @return     integer         FTI_SCES if successful.

    This function detects all the erasures for L1. If there is at least one,
    L1 is not considered as recoverable.

 **/
/*-------------------------------------------------------------------------*/
int FTI_RecoverL1(int group) {
    int         erased[FTI_BUFS], buf, i, j; // FTI_BUFS > 32*3
    unsigned long fs, maxFs;
    if (FTI_CheckErasures(&fs, &maxFs, group, erased, 1) != FTI_SCES)
        { FTI_Print("Error checking erasures.", FTI_DBUG); return FTI_NSCS; }
    buf = 0; for(j = 0; j < FTI_Topo.groupSize; j++) if(erased[j]) buf++; // Counting erasures
    if (buf > 0) { FTI_Print("Checkpoint files missing at L1.", FTI_DBUG); return FTI_NSCS; }
    return FTI_SCES;
}


/*-------------------------------------------------------------------------*/
/**
    @brief      Recover L2 ckpt. files using the partner copy.
    @param      group           The group ID.
    @return     integer         FTI_SCES if successful.

    This function tries to recover the L2 ckpt. files missing using the
    partner copy. If a ckpt. file and its copy are both missing, then we
    consider this checkpoint unavailable.

 **/
/*-------------------------------------------------------------------------*/
int FTI_RecoverL2(int group) {
    int         erased[FTI_BUFS], gs, buf, j, src, dest;
    char        str[FTI_BUFS], lfn[FTI_BUFS], pfn[FTI_BUFS], jfn[FTI_BUFS], qfn[FTI_BUFS];
    char        *blBuf1, *blBuf2, *blBuf3, *blBuf4;
    unsigned long ps, fs, maxFs, pos = 0;
    FILE        *lfd, *pfd, *jfd, *qfd;
    MPI_Request reqSend1, reqRecv1, reqSend2, reqRecv2;
    MPI_Status  status;
    blBuf1 = talloc(char, FTI_Conf.blockSize);
    blBuf2 = talloc(char, FTI_Conf.blockSize);
    blBuf3 = talloc(char, FTI_Conf.blockSize);
    blBuf4 = talloc(char, FTI_Conf.blockSize);
    gs = FTI_Topo.groupSize;
    src = FTI_Topo.left;
    dest = FTI_Topo.right;
    if (access(FTI_Ckpt[2].dir, F_OK) != 0) mkdir(FTI_Ckpt[2].dir, 0777);
    if ( FTI_CheckErasures(&fs, &maxFs, group, erased, 2) != FTI_SCES) // Checking erasures
        { FTI_Print("Error checking erasures.", FTI_DBUG); return FTI_NSCS; }
    buf = -1; for(j = 0; j < gs; j++) if(erased[j] && erased[((j+1)%gs)+gs]) buf=j; // Counting erasures
    sprintf(str, "A checkpoint file and its partner copy (ID in group : %d) have been lost", buf);
    if (buf > -1) { FTI_Print(str, FTI_DBUG); return FTI_NSCS; }
    buf = 0; for(j = 0; j < gs*2; j++) if(erased[j]) buf++; // Counting erasures
    if (buf > 0) {
        ps = (maxFs/FTI_Conf.blockSize)*FTI_Conf.blockSize; pos = 0; // For the logic
        if (ps < maxFs) ps = ps + FTI_Conf.blockSize; // Calculating padding size
        sprintf(str,"File size: %ld, max. file size : %ld and padding size : %ld.", fs, maxFs, ps);
        FTI_Print(str, FTI_DBUG);
        if (erased[FTI_Topo.groupRank]) { // Open checkpoint file to recover
            sprintf(lfn,"%s/%s", FTI_Ckpt[2].dir, FTI_Exec.ckptFile);
            sscanf(FTI_Exec.ckptFile,"Ckpt%d-Rank%d.fti", &FTI_Exec.ckptID, &buf);
            sprintf(jfn,"%s/Ckpt%d-Pcof%d.fti", FTI_Ckpt[2].dir, FTI_Exec.ckptID, buf);
            sprintf(str,"Opening checkpoint file (%s) to recover (L2).", lfn); FTI_Print(str, FTI_DBUG);
            sprintf(str,"Opening partner ckpt. file (%s) to recover (L2).", jfn); FTI_Print(str, FTI_DBUG);
            lfd = fopen(lfn, "wb"); jfd = fopen(jfn, "wb");
            if (lfd == NULL) { FTI_Print("R2 cannot open the checkpoint file.", FTI_DBUG); return FTI_NSCS; }
            if (jfd == NULL) { FTI_Print("R2 cannot open the partner ckpt. file.", FTI_DBUG); return FTI_NSCS; }
        }
        if (erased[src] && !erased[gs+FTI_Topo.groupRank]) { // Truncate and open partner file to transfer
            sscanf(FTI_Exec.ckptFile,"Ckpt%d-Rank%d.fti", &FTI_Exec.ckptID, &buf);
            sprintf(pfn,"%s/Ckpt%d-Pcof%d.fti", FTI_Ckpt[2].dir, FTI_Exec.ckptID, buf);
            sprintf(str,"Opening partner ckpt. file (%s) to transfer (L2).", pfn); FTI_Print(str, FTI_DBUG);
            if (truncate(pfn,ps) == -1) { FTI_Print("R2 cannot truncate the partner ckpt. file.", FTI_DBUG); return FTI_NSCS; }
            pfd = fopen(pfn, "rb");
            if (pfd == NULL) { FTI_Print("R2 cannot open partner ckpt. file.", FTI_DBUG); return FTI_NSCS; }
        }
        if (erased[dest] && !erased[gs+FTI_Topo.groupRank]) { // Truncate and open partner file to transfer
            sprintf(qfn,"%s/%s", FTI_Ckpt[2].dir, FTI_Exec.ckptFile);
            sprintf(str,"Opening ckpt. file (%s) to transfer (L2).", qfn); FTI_Print(str, FTI_DBUG);
            if (truncate(qfn,ps) == -1) { FTI_Print("R2 cannot truncate the ckpt. file.", FTI_DBUG); return FTI_NSCS; }
            qfd = fopen(qfn, "rb");
            if (qfd == NULL) { FTI_Print("R2 cannot open ckpt. file.", FTI_DBUG); return FTI_NSCS; }
        }
        while(pos < ps) { // Checkpoint files exchange
            if (erased[src] && !erased[gs+FTI_Topo.groupRank]) {
                fread(blBuf1, sizeof(char), FTI_Conf.blockSize, pfd);
                MPI_Isend(blBuf1, FTI_Conf.blockSize, MPI_CHAR, src, FTI_Conf.tag, FTI_Exec.groupComm, &reqSend1);
            }
            if (erased[dest] && !erased[gs+FTI_Topo.groupRank]) {
                fread(blBuf3, sizeof(char), FTI_Conf.blockSize, qfd);
                MPI_Isend(blBuf3, FTI_Conf.blockSize, MPI_CHAR, dest, FTI_Conf.tag, FTI_Exec.groupComm, &reqSend2);
            }
            if (erased[FTI_Topo.groupRank]) {
                MPI_Irecv(blBuf2, FTI_Conf.blockSize, MPI_CHAR, dest, FTI_Conf.tag, FTI_Exec.groupComm, &reqRecv1);
                MPI_Irecv(blBuf4, FTI_Conf.blockSize, MPI_CHAR, src, FTI_Conf.tag, FTI_Exec.groupComm, &reqRecv2);
            }
            if (erased[src] && !erased[gs+FTI_Topo.groupRank]) MPI_Wait(&reqSend1, &status);
            if (erased[dest] && !erased[gs+FTI_Topo.groupRank]) MPI_Wait(&reqSend2, &status);
            if (erased[FTI_Topo.groupRank]) {
                MPI_Wait(&reqRecv1, &status);
                MPI_Wait(&reqRecv2, &status);
                if (fwrite(blBuf2, sizeof(char), FTI_Conf.blockSize, lfd) != FTI_Conf.blockSize)
                    { FTI_Print("Errors writting the data in the R2 checkpoint file.", FTI_DBUG); return FTI_NSCS; }
                if (fwrite(blBuf4, sizeof(char), FTI_Conf.blockSize, jfd) != FTI_Conf.blockSize)
                    { FTI_Print("Errors writting the data in the R2 partner ckpt. file.", FTI_DBUG); return FTI_NSCS; }
            }
            pos = pos + FTI_Conf.blockSize;
        }
        if (erased[FTI_Topo.groupRank]) { // Close files
            if (fclose(lfd) != 0) { FTI_Print("R2 cannot close the checkpoint file.", FTI_DBUG); return FTI_NSCS; }
            if (truncate(lfn,fs) == -1) { FTI_Print("R2 cannot re-truncate the checkpoint file.", FTI_DBUG); return FTI_NSCS; }
            if (fclose(jfd) != 0) { FTI_Print("R2 cannot close the partner ckpt. file.", FTI_DBUG); return FTI_NSCS; }
            if (truncate(jfn,fs) == -1) { FTI_Print("R2 cannot re-truncate the partner ckpt. file.", FTI_DBUG); return FTI_NSCS; }
        }
        if (erased[src] && !erased[gs+FTI_Topo.groupRank]) {
            if (fclose(pfd) != 0) { FTI_Print("R2 cannot close the partner ckpt. file", FTI_DBUG); return FTI_NSCS; }
            if (truncate(pfn,fs) == -1) { FTI_Print("R2 cannot re-truncate the partner ckpt. file.", FTI_DBUG); return FTI_NSCS; }
        }
        if (erased[dest] && !erased[gs+FTI_Topo.groupRank]) {
            if (fclose(qfd) != 0) { FTI_Print("R2 cannot close the ckpt. file", FTI_DBUG); return FTI_NSCS; }
            if (truncate(qfn,fs) == -1) { FTI_Print("R2 cannot re-truncate the ckpt. file.", FTI_DBUG); return FTI_NSCS; }
        }
    }
    free(blBuf1); free(blBuf2); // Free memory
    return FTI_SCES;
}


/*-------------------------------------------------------------------------*/
/**
    @brief      Recover L3 ckpt. files ordering the RS decoding algorithm.
    @param      group           The group ID.
    @return     integer         FTI_SCES if successful.

    This function tries to recover the L3 ckpt. files missing using the
    RS decoding. If to many files are missing in the group, then we
    consider this checkpoint unavailable.

 **/
/*-------------------------------------------------------------------------*/
int FTI_RecoverL3(int group) {
    int         erased[FTI_BUFS], gs, j, l = 0;
    unsigned long fs, maxFs;
    char        str[FTI_BUFS];
    gs = FTI_Topo.groupSize;
    if (access(FTI_Ckpt[3].dir, F_OK) != 0) mkdir(FTI_Ckpt[3].dir, 0777);
    if ( FTI_CheckErasures(&fs, &maxFs, group, erased, 3) != FTI_SCES) // Checking erasures
        { FTI_Print("Error checking erasures.", FTI_DBUG); return FTI_NSCS; }
    l = 0; for(j = 0; j < gs; j++) { if(erased[j]) l++; if(erased[j+gs]) l++; } // Counting erasures
    if (l > gs) { FTI_Print("Too many erasures at L3.", FTI_DBUG); return FTI_NSCS; }
    if (l > 0) {
        sprintf(str, "There are %d encoded/checkpoint files missing in this group.", l); FTI_Print(str, FTI_DBUG);
        if (FTI_Decode(fs, maxFs, erased) == FTI_NSCS)
            { FTI_Print("RS-decoding could not regenerate the missing data."); return FTI_NSCS; }
    } // Reed-Solomon decoding
    return FTI_SCES;
}


/*-------------------------------------------------------------------------*/
/**
    @brief      Recover L4 ckpt. files from the PFS.
    @param      group           The group ID.
    @return     integer         FTI_SCES if successful.

    This function tries to recover the ckpt. files using the L4 ckpt. files
    stored in the PFS. If at least one ckpt. file is missing in the PFS, we
    consider this checkpoint unavailable.

 **/
/*-------------------------------------------------------------------------*/
int FTI_RecoverL4(int group) {
    unsigned long maxFs, fs, ps, pos = 0;
    int         j, l, gs, erased[FTI_BUFS];
    char        gfn[FTI_BUFS], lfn[FTI_BUFS], *blBuf1;
    FILE        *gfd, *lfd;
    blBuf1 = talloc(char, FTI_Conf.blockSize); // Allocate memory
    gs = FTI_Topo.groupSize;
    if (access(FTI_Ckpt[1].dir, F_OK) != 0) mkdir(FTI_Ckpt[1].dir, 0777);
    if ( FTI_CheckErasures(&fs, &maxFs, group, erased, 4) != FTI_SCES) // Checking erasures
        { FTI_Print("Error checking erasures.", FTI_DBUG); return FTI_NSCS; }
    l = 0; for(j = 0; j < gs; j++) { if(erased[j]) l++; } // Counting erasures
    if (l > 0) { FTI_Print("Checkpoint file missing at L4.", FTI_DBUG); return FTI_NSCS; }
    ps = (fs/FTI_Conf.blockSize)*FTI_Conf.blockSize; pos = 0; // For the logic
    if (ps < fs) ps = ps + FTI_Conf.blockSize; // Calculating padding size
    sprintf(gfn,"%s/%s", FTI_Ckpt[4].dir, FTI_Exec.ckptFile); // Open and resize files
    sprintf(lfn,"%s/%s", FTI_Ckpt[1].dir, FTI_Exec.ckptFile);
    if (access(gfn, R_OK) != 0) { FTI_Print("R4 cannot read the checkpoint file in the PFS.", FTI_DBUG); return FTI_NSCS; }
    if (truncate(gfn,ps) == -1) { FTI_Print("R4 cannot truncate the ckpt. file in the PFS.", FTI_DBUG); return FTI_NSCS; }
    gfd = fopen(gfn, "rb"); lfd = fopen(lfn, "wb");
    if (gfd == NULL) { FTI_Print("R4 cannot open the ckpt. file in the PFS.", FTI_DBUG); return FTI_NSCS; }
    if (lfd == NULL) { FTI_Print("R4 cannot open the local ckpt. file.", FTI_DBUG); return FTI_NSCS; }
    while(pos < ps) { // Checkpoint files transfer from PFS
        fread(blBuf1, sizeof(char), FTI_Conf.blockSize, gfd);
        fwrite(blBuf1, sizeof(char), FTI_Conf.blockSize, lfd);
        pos = pos + FTI_Conf.blockSize;
    }
    fclose(gfd); fclose(lfd); // Close files
    if (truncate(gfn,fs) == -1) { FTI_Print("R4 cannot re-truncate the checkpoint file in the PFS.", FTI_DBUG); return FTI_NSCS; }
    if (truncate(lfn,fs) == -1) { FTI_Print("R4 cannot re-truncate the local checkpoint file.",  FTI_DBUG); return FTI_NSCS; }
    free(blBuf1);
    return FTI_SCES;
}


/*-------------------------------------------------------------------------*/
/**
    @brief      Recover a set of ckpt. files using RS decoding.
    @return     integer         FTI_SCES if successful.

    This function tries to recover the L3 ckpt. files missing using the
    RS decoding.

 **/
/*-------------------------------------------------------------------------*/
int FTI_Decode(int fs, int maxFs, int *erased) {
    int *matrix, *decMatrix, *dm_ids, *tmpmat, i, j, k, m, ps, bs, pos = 0;
    char **coding, **data, *dataTmp, fn[FTI_BUFS], efn[FTI_BUFS], str[FTI_BUFS];
    FILE *fd, *efd;
    bs = FTI_Conf.blockSize; k = FTI_Topo.groupSize; m = k;
    ps = ((maxFs/FTI_Conf.blockSize))*FTI_Conf.blockSize;
    if (ps < maxFs) ps = ps + FTI_Conf.blockSize; // Calculating padding size
    if (access(FTI_Ckpt[3].dir, F_OK) != 0) mkdir(FTI_Ckpt[3].dir, 0777);
    sscanf(FTI_Exec.ckptFile,"Ckpt%d-Rank%d.fti", &FTI_Exec.ckptID, &i);
    sprintf(fn,"%s/%s",FTI_Ckpt[3].dir, FTI_Exec.ckptFile);
    sprintf(efn,"%s/Ckpt%d-RSed%d.fti", FTI_Ckpt[3].dir, FTI_Exec.ckptID, i);
    data = talloc(char *, k); coding = talloc(char *, m); dataTmp = talloc(char, FTI_Conf.blockSize*k);
    dm_ids = talloc(int, k); decMatrix = talloc(int, k*k); tmpmat = talloc(int, k*k); matrix =  talloc(int, k*k);
    if (FTI_CreateMatrix(matrix) == FTI_SCES) FTI_Print("Matrix created.", FTI_DBUG);
    for (i = 0; i < m; i++) {
        coding[i] = talloc(char, FTI_Conf.blockSize);
        data[i] = talloc(char, FTI_Conf.blockSize);
    }
    j = 0; for (i = 0; j < k; i++) { if (erased[i] == 0) {dm_ids[j] = i; j++;} }
    for (i = 0; i < k; i++) { // Building the matrix
        if (dm_ids[i] < k) {
            for (j = 0; j < k; j++) tmpmat[i*k+j] = 0;
            tmpmat[i*k+dm_ids[i]] = 1;
        } else for (j = 0; j < k; j++) { tmpmat[i*k+j] = matrix[(dm_ids[i]-k)*k+j]; }
    } // Inversing the matrix
    if (jerasure_invert_matrix(tmpmat, decMatrix, k, FTI_Conf.l3WordSize) < 0)
        { FTI_Print("Error inversing matrix", FTI_DBUG); return FTI_NSCS; }
    if(erased[FTI_Topo.groupRank] == 0) { // Resize and open files
        if (truncate(fn,ps) == -1) { FTI_Print("Error with truncate on checkpoint file", FTI_DBUG); return FTI_NSCS; }
        fd = fopen(fn, "rb"); efd = fopen(efn, "rb");
    } else { fd = fopen(fn, "wb"); efd = fopen(efn, "wb"); }
    if (fd == NULL) { FTI_Print("R3 cannot open checkpoint file.", FTI_DBUG); return FTI_NSCS; }
    if (efd == NULL) { FTI_Print("R3 cannot open encoded ckpt. file.", FTI_DBUG); return FTI_NSCS; }
    while(pos < ps) { // Main loop, block by block
        if(erased[FTI_Topo.groupRank] == 0) { // Reading the data
            fread(data[FTI_Topo.groupRank]+0, sizeof(char), bs, fd);
            fread(coding[FTI_Topo.groupRank]+0, sizeof(char), bs, efd);
        } else { bzero(data[FTI_Topo.groupRank], bs); bzero(coding[FTI_Topo.groupRank], bs); } // Erasure found
        MPI_Allgather(data[FTI_Topo.groupRank]+0, bs, MPI_CHAR, dataTmp, bs, MPI_CHAR, FTI_Exec.groupComm);
        for (i = 0; i < k; i++) memcpy(data[i]+0, &(dataTmp[i*bs]), sizeof(char)*bs);
        MPI_Allgather(coding[FTI_Topo.groupRank]+0, bs, MPI_CHAR, dataTmp, bs, MPI_CHAR, FTI_Exec.groupComm);
        for (i = 0; i < k; i++) memcpy(coding[i]+0, &(dataTmp[i*bs]), sizeof(char)*bs);
        if (erased[FTI_Topo.groupRank]) // Decoding the lost data work
            jerasure_matrix_dotprod(k, FTI_Conf.l3WordSize, decMatrix+(FTI_Topo.groupRank*k), dm_ids, FTI_Topo.groupRank, data, coding, bs);
        MPI_Allgather(data[FTI_Topo.groupRank]+0, bs, MPI_CHAR, dataTmp, bs, MPI_CHAR, FTI_Exec.groupComm);
        for (i = 0; i < k; i++) memcpy(data[i]+0, &(dataTmp[i*bs]), sizeof(char)*bs);
        if (erased[FTI_Topo.groupRank + k]) // Finally, re-encode any erased encoded checkpoint file
            jerasure_matrix_dotprod(k, FTI_Conf.l3WordSize, matrix+(FTI_Topo.groupRank*k), NULL, FTI_Topo.groupRank+k, data, coding, bs);
        if (erased[FTI_Topo.groupRank]) fwrite(data[FTI_Topo.groupRank]+0, sizeof(char), bs, fd);
        if (erased[FTI_Topo.groupRank + k]) fwrite(coding[FTI_Topo.groupRank]+0, sizeof(char), bs, efd);
        pos = pos + bs;
    }
    fclose(fd); fclose(efd); // Closing files
    if (truncate(fn,fs) == -1) { FTI_Print("R3 cannot re-truncate checkpoint file.", FTI_DBUG); return FTI_NSCS; }
    if (truncate(efn,fs) == -1) { FTI_Print("R3 cannot re-truncate encoded ckpt. file.", FTI_DBUG); return FTI_NSCS; }
    free(tmpmat); free(dm_ids); free(decMatrix); free(matrix); free(data); free(dataTmp); free(coding);
    return FTI_SCES;
}

