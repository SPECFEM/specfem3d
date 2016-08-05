/**
 *  @file   checkpoint.c
 *  @author Leonardo A. Bautista Gomez (leobago@gmail.com)
 *  @date   July, 2013
 *  @brief  Checkpointing functions for the FTI library.
 */


#include "fti.h"


/*-------------------------------------------------------------------------*/
/**
    @brief      Creates a distribution matrix for the RS encoding.
    @param      matrix          Matrix to be filled.
    @return     integer         FTI_SCES if successfull.

    This function creates a distribution matrix for the RS encoding and
    fill the values of the one passed in parameter. It uses the FTI
    configuration and the jerasure library.

 **/
/*-------------------------------------------------------------------------*/
int FTI_CreateMatrix(int *matrix) {
    int i,j;
    for (i = 0; i < FTI_Topo.groupSize; i++) {
        for (j = 0; j < FTI_Topo.groupSize; j++) {
            matrix[i*FTI_Topo.groupSize+j] = galois_single_divide(1, i ^ (FTI_Topo.groupSize + j), FTI_Conf.l3WordSize);
        }
    }
    return FTI_SCES;
}


/*-------------------------------------------------------------------------*/
/**
    @brief      It gets the metadata to recover the data after a failure.
    @param      fs              Pointer to fill the checkpoint file size.
    @param      mfs             Pointer to fill the maximum file size.
    @param      group           The group in the node.
    @param      level           The level of the ckpt or 0 if tmp.
    @return     integer         FTI_SCES if successfull.

    This function read the metadata file created during checkpointing and
    recover the checkpoint file name, file size and the size of the largest
    file in the group (for padding if ncessary during decoding).

 **/
/*-------------------------------------------------------------------------*/
int FTI_GetMeta(unsigned long *fs, unsigned long *mfs, int group, int level) {
    dictionary *ini;
    char mfn[FTI_BUFS], str[FTI_BUFS], *cfn;
    if(level == 0) sprintf(mfn,"%s/sector%d-group%d.fti",FTI_Conf.mTmpDir, FTI_Topo.sectorID, group);
    else sprintf(mfn,"%s/sector%d-group%d.fti",FTI_Ckpt[level].metaDir, FTI_Topo.sectorID, group);
    sprintf(str, "Getting FTI metadata file (%s)...",mfn);
    FTI_Print(str, FTI_DBUG);
    if (access(mfn, R_OK) == 0) FTI_Print("FTI metadata file accessible.", FTI_DBUG);
    else return FTI_NSCS;
    ini = iniparser_load(mfn);
    if (ini != NULL) FTI_Print("Iniparser parsed the FTI metadata file.", FTI_DBUG);
    else return FTI_NSCS;
    sprintf(str, "%d:Ckpt_file_name", FTI_Topo.groupRank);
    cfn = iniparser_getstring(ini, str, NULL);
    snprintf(FTI_Exec.ckptFile, FTI_BUFS, "%s", cfn);
    sprintf(str, "%d:Ckpt_file_size", FTI_Topo.groupRank);
    *fs = (int) iniparser_getint(ini, str, -1);
    sprintf(str, "%d:Ckpt_file_maxs", FTI_Topo.groupRank);
    *mfs = (int) iniparser_getint(ini, str, -1);
    iniparser_freedict(ini);
    return FTI_SCES;
}


/*-------------------------------------------------------------------------*/
/**
    @brief      It writes the metadata to recover the data after a failure.
    @return     integer         FTI_SCES if successfull.

    This function gathers information about the checkpoint files in the
    group (name and sizes), and creates the metadata file used to recover in
    case of failure.

 **/
/*-------------------------------------------------------------------------*/
int FTI_CreateMetadata() {
    char str[FTI_BUFS], buf[FTI_BUFS], *fnl = talloc(char, FTI_Topo.groupSize*FTI_BUFS);
    unsigned long fs[FTI_BUFS], mfs;
    struct stat fileStatus;
    unsigned long tmpo;
    dictionary *ini;
    int i, start;
    if (FTI_Ckpt[4].isInline && FTI_Exec.ckptLvel == 4) sprintf(buf,"%s/%s",FTI_Conf.gTmpDir, FTI_Exec.ckptFile);
    else sprintf(buf,"%s/%s",FTI_Conf.lTmpDir, FTI_Exec.ckptFile); // Getting size of files
    if(stat(buf, &fileStatus) == 0) fs[FTI_Topo.groupRank] = (unsigned long) fileStatus.st_size;
    else { FTI_Print("Error with stat on the checkpoint file.", FTI_DBUG); return FTI_NSCS; }
    sprintf(str, "Checkpoint file size : %ld bytes.", fs[FTI_Topo.groupRank]); FTI_Print(str, FTI_DBUG);
    sprintf(fnl+(FTI_Topo.groupRank*FTI_BUFS),"%s",FTI_Exec.ckptFile);
    tmpo = fs[FTI_Topo.groupRank]; // Gather all the file sizes
    MPI_Allgather(&tmpo, 1, MPI_UNSIGNED_LONG, fs, 1, MPI_UNSIGNED_LONG, FTI_Exec.groupComm);
    strncpy(str,fnl+(FTI_Topo.groupRank*FTI_BUFS),FTI_BUFS); // Gather all the file names
    MPI_Allgather(str, FTI_BUFS, MPI_CHAR, fnl, FTI_BUFS, MPI_CHAR, FTI_Exec.groupComm);
    mfs = 0; for(i = 0; i < FTI_Topo.groupSize; i++) if (fs[i] > mfs) mfs = fs[i]; // Search max. size
    sprintf(str, "Max. file size %ld.", mfs); FTI_Print(str, FTI_DBUG);
    if (FTI_Topo.groupRank == 0) { // Only one process in the group create the metadata
        snprintf(buf, FTI_BUFS, "%s/Topology.fti",FTI_Conf.metadDir);
        sprintf(str, "Temporary load of topology file (%s)...", buf);
        FTI_Print(str, FTI_DBUG);
        ini = iniparser_load(buf); // To bypass iniparser bug while empty dictionary
        if (ini != NULL) FTI_Print("Temporary topology file parsed", FTI_DBUG);
        else { FTI_Print("Temporary topology file could NOT be parsed", FTI_DBUG); return FTI_NSCS; }
        //if (FTI_Topo.nbHeads == 1) start = 0; else start = 0;
        for (i = 0; i < FTI_Topo.groupSize; i++) { // Add metadata to dictionary
            strncpy(buf,fnl+(i*FTI_BUFS),FTI_BUFS);
            sprintf(str,"%d", i);
            iniparser_set(ini, str, NULL);
            sprintf(str,"%d:Ckpt_file_name", i);
            iniparser_set(ini, str, buf);
            sprintf(str,"%d:Ckpt_file_size", i);
            sprintf(buf,"%ld", fs[i]);
            iniparser_set(ini, str, buf);
            sprintf(str,"%d:Ckpt_file_maxs", i);
            sprintf(buf,"%ld", mfs);
            iniparser_set(ini, str, buf);
        }
        iniparser_unset(ini, "topology"); // Remove topology section
        if (access(FTI_Conf.mTmpDir, F_OK) != 0) mkdir(FTI_Conf.mTmpDir, 0777);
        sprintf(buf, "%s/sector%d-group%d.fti", FTI_Conf.mTmpDir, FTI_Topo.sectorID, FTI_Topo.groupID);
        remove(buf);
        sprintf(str, "Creating metadata file (%s)...", buf);
        FTI_Print(str, FTI_DBUG);
        FILE *fd = fopen(buf, "w");
        if (fd != NULL) FTI_Print("Metadata file opened.", FTI_DBUG);
        else return FTI_NSCS;
        iniparser_dump_ini(ini, fd); // Write metadata
        if (fflush(fd) == 0) FTI_Print("Metadata file flushed.", FTI_DBUG);
        else return FTI_NSCS;
        if (fclose(fd) == 0) FTI_Print("Metadata file closed.", FTI_DBUG);
        else return FTI_NSCS;
        iniparser_freedict(ini);
    }
    free(fnl);
    return FTI_SCES;
}


/*-------------------------------------------------------------------------*/
/**
    @brief      Updates the checkpoint information.
    @param      offline         TRUE if the last ckpt. was offline.
    @return     integer         FTI_SCES if successful.

    This function updates the checkpoint information for the next ckpt.

 **/
/*-------------------------------------------------------------------------*/
int FTI_UpdateCkptInfo(int offline) {
    FTI_Exec.wasLastOffline = offline;
    FTI_Ckpt[FTI_Exec.ckptLvel].lastCkpt = FTI_Exec.ckptID;
    FTI_Exec.lastCkptLvel = FTI_Exec.ckptLvel;
    FTI_Exec.ckptID++;
    sprintf(FTI_Exec.ckptFile,"Ckpt%d-Rank%d.fti", FTI_Exec.ckptID, FTI_Topo.myRank);
    return FTI_SCES;
}


/*-------------------------------------------------------------------------*/
/**
    @brief      Decides wich action start depending on the ckpt. level.
    @return     integer         FTI_SCES if successful.

    This function launchs the required action dependeing on the ckpt. level.
    It does thast for each group (application process in the node).

 **/
/*-------------------------------------------------------------------------*/
int FTI_CkptNotified() {
    char str[FTI_BUFS];
    int i, tres, res = 0;
    unsigned long maxFs, fs;
    double timer = MPI_Wtime(); // For timing checkpoint post-processing
    if (FTI_Topo.amIaHead) { // If I am the Head
        for(i = 0; i < FTI_Topo.nbApprocs; i++) { // I post-process every ckpt. in the node
            switch(FTI_Exec.ckptLvel) { // Action required depending on the ckpt. level
                case 4 : if (!FTI_Ckpt[4].isInline) res += FTI_Flush(i+1, 0); break;
                case 3 : if (!FTI_Ckpt[3].isInline) res += FTI_RSencode(i+1); break;
                case 2 : if (!FTI_Ckpt[2].isInline) res += FTI_Partner(i+1); break;
                case 1 : res += FTI_SCES; break;
            }
        }
        if (res == FTI_SCES) FTI_Print("Ckpt. post-processing finished successfully.", FTI_DBUG);
        else FTI_Print("Ckpt. post-processing did NOT finish successfully.", FTI_DBUG);
        MPI_Allreduce(&res, &tres, 1, MPI_INT, MPI_SUM, FTI_COMM_WORLD);
        for(i = 0; i < FTI_Topo.nbApprocs; i++) {
            if (tres == FTI_SCES) FTI_Clean(FTI_Exec.ckptLvel, i+1, FTI_Topo.body[i]);
            else FTI_Clean(0, i+1, FTI_Topo.body[i]);
        }
        FTI_Print("Barrier after cleaning previous checkpoint.", FTI_DBUG); MPI_Barrier(FTI_COMM_WORLD);
        sprintf(str, "Ckpt. post-processing at L%d took %f seconds.", FTI_Exec.ckptLvel, MPI_Wtime()-timer);
        if (tres == FTI_SCES) FTI_UpdateCkptInfo(1);
    } else { // If I am not Head I only post-process my ckpt.
        switch(FTI_Exec.ckptLvel) { // Action required depending on the ckpt. level
            case 4 : if (FTI_Ckpt[4].isInline) res = FTI_SCES; break;
            case 3 : if (FTI_Ckpt[3].isInline) res = FTI_RSencode(FTI_Topo.groupID); break;
            case 2 : if (FTI_Ckpt[2].isInline) res = FTI_Partner(FTI_Topo.groupID); break;
            case 1 : res = FTI_SCES; break;
        }
        if (res == FTI_SCES) FTI_Print("Ckpt. post-processing finished successfully.", FTI_DBUG);
        else FTI_Print("Ckpt. post-processing did NOT finish successfully.", FTI_DBUG);
        MPI_Allreduce(&res, &tres, 1, MPI_INT, MPI_SUM, FTI_COMM_WORLD);
        if (tres == FTI_SCES) FTI_Clean(FTI_Exec.ckptLvel, FTI_Topo.groupID, FTI_Topo.myRank);
        else FTI_Clean(0, FTI_Topo.groupID, FTI_Topo.myRank);
        FTI_Print("Barrier after cleaning previous checkpoint.", FTI_DBUG); MPI_Barrier(FTI_COMM_WORLD);
        sprintf(str, "Ckpt. post-processing at L%d took %f seconds.", FTI_Exec.ckptLvel, MPI_Wtime()-timer);
        if (tres == FTI_SCES) FTI_UpdateCkptInfo(0);
    }
    if (res == FTI_SCES) FTI_Print(str, FTI_DBUG);
    switch(FTI_Exec.ckptLvel) { // Action required depending on the ckpt. level
        case 4 : rename(FTI_Conf.gTmpDir, FTI_Ckpt[4].dir); rename(FTI_Conf.mTmpDir, FTI_Ckpt[4].metaDir); break;
        case 3 : rename(FTI_Conf.lTmpDir, FTI_Ckpt[3].dir); rename(FTI_Conf.mTmpDir, FTI_Ckpt[3].metaDir); break;
        case 2 : rename(FTI_Conf.lTmpDir, FTI_Ckpt[2].dir); rename(FTI_Conf.mTmpDir, FTI_Ckpt[2].metaDir); break;
        case 1 : rename(FTI_Conf.lTmpDir, FTI_Ckpt[1].dir); rename(FTI_Conf.mTmpDir, FTI_Ckpt[1].metaDir); break;
    }
    if (FTI_Topo.amIaHead)  // If I am the Head
        for(i = 0; i < FTI_Topo.nbApprocs; i++) // Send msg. to avoid checkpoint collision
            MPI_Send(&tres, 1, MPI_INT, FTI_Topo.body[i], FTI_Conf.tag, FTI_Exec.globalComm);
    return FTI_SCES;
}


/*-------------------------------------------------------------------------*/
/**
    @brief      It copies ckpt. files in to the partner node.
    @param      group           The group ID.
    @return     integer         FTI_SCES if successful.

    This function copies the checkpoint files into the pertner node. It
    follows a ring, where the ring size is the group size given in the FTI
    configuration file.

 **/
/*-------------------------------------------------------------------------*/
int FTI_Partner(int group) {
    char        *blBuf1, *blBuf2, lfn[FTI_BUFS], pfn[FTI_BUFS], str[FTI_BUFS];
    unsigned long maxFs, fs, ps, pos = 0; // FTI_BUFS > 32
    int         dest, src;
    MPI_Request reqSend, reqRecv;
    FILE        *lfd, *pfd;
    MPI_Status  status;
    blBuf1 = talloc(char, FTI_Conf.blockSize);
    blBuf2 = talloc(char, FTI_Conf.blockSize);
    if (FTI_GetMeta(&fs, &maxFs, group, 0) == FTI_SCES) FTI_Print("Metadata obtained.", FTI_DBUG);
    else { FTI_Print("Error getting metadata.", FTI_DBUG); return FTI_NSCS; }
    ps = (maxFs/FTI_Conf.blockSize)*FTI_Conf.blockSize;
    if (ps < maxFs) ps = ps + FTI_Conf.blockSize;
    sprintf(str, "Max. file size %ld and padding size %ld.", maxFs, ps); FTI_Print(str, FTI_DBUG);
    sscanf(FTI_Exec.ckptFile,"Ckpt%d-Rank%d.fti", &FTI_Exec.ckptID, &src);
    sprintf(lfn,"%s/%s",FTI_Conf.lTmpDir, FTI_Exec.ckptFile);
    sprintf(pfn,"%s/Ckpt%d-Pcof%d.fti", FTI_Conf.lTmpDir, FTI_Exec.ckptID, src);
    sprintf(str, "L2 trying to access local ckpt. file (%s).", lfn); FTI_Print(str, FTI_DBUG);
    dest = FTI_Topo.right; // Defining source and destination partners
    src = FTI_Topo.left;
    if (access(lfn, R_OK) != 0) { FTI_Print("L2 cannot access the checkpoint file.", FTI_DBUG); return FTI_NSCS; }
    if (truncate(lfn,ps) == -1) { FTI_Print("L2 cannot truncate the checkpoint file.", FTI_DBUG); return FTI_NSCS; }
    lfd = fopen(lfn, "rb"); pfd = fopen(pfn, "wb");
    if (lfd == NULL) { FTI_Print("L2 cannot open the checkpoint file.", FTI_DBUG); return FTI_NSCS; }
    if (pfd == NULL) { FTI_Print("L2 cannot open the partner ckpt. file.", FTI_DBUG); return FTI_NSCS; }
    while(pos < ps) { // Checkpoint files partner copy
        fread(blBuf1, sizeof(char), FTI_Conf.blockSize, lfd);
        MPI_Isend(blBuf1, FTI_Conf.blockSize, MPI_CHAR, dest, FTI_Conf.tag, FTI_Exec.groupComm, &reqSend);
        MPI_Irecv(blBuf2, FTI_Conf.blockSize, MPI_CHAR, src, FTI_Conf.tag, FTI_Exec.groupComm, &reqRecv);
        MPI_Wait(&reqSend, &status);
        MPI_Wait(&reqRecv, &status);
        fwrite(blBuf2, sizeof(char), FTI_Conf.blockSize, pfd);
        pos = pos + FTI_Conf.blockSize;
    }
    fclose(lfd); fclose(pfd); // Close files
    if (truncate(lfn,fs) == -1) { FTI_Print("L2 cannot re-truncate the checkpoint file.", FTI_DBUG); return FTI_NSCS; }
    if (truncate(pfn,fs) == -1) { FTI_Print("L2 cannot re-truncate the partner ckpt. file.", FTI_DBUG); return FTI_NSCS; }
    return FTI_SCES;
}


/*-------------------------------------------------------------------------*/
/**
    @brief      It performs RS encoding with the ckpt. files in to the group.
    @param      group           The group ID.
    @return     integer         FTI_SCES if successful.

    This function performs the Reed-Solomon encoding for a given group. The
    checkpoint files are padded to the maximum size of the largest checkpoint
    file in the group +- the extra space to be a multiple of block size.

 **/
/*-------------------------------------------------------------------------*/
int FTI_RSencode(int group) {
    char *myData, *data, *coding, lfn[FTI_BUFS], efn[FTI_BUFS], str[FTI_BUFS];
    int *matrix, cnt, i, init, src, offset, dest, matVal, bs = FTI_Conf.blockSize;
    unsigned long maxFs, fs, ps, pos = 0; // FTI_BUFS > 32
    MPI_Request reqSend, reqRecv;
    MPI_Status status;
    FILE *lfd, *efd;
    if (FTI_GetMeta(&fs, &maxFs, group, 0) == FTI_SCES) FTI_Print("Metadata obtained.", FTI_DBUG);
    else { FTI_Print("Error getting metadata.", FTI_DBUG); return FTI_NSCS; }
    ps = ((maxFs/bs))*bs;
    if (ps < maxFs) ps = ps + bs; // Calculating padding size
    myData = talloc(char, bs);
    data = talloc(char, 2*bs);
    coding = talloc(char, bs);
    matrix =  talloc(int, FTI_Topo.groupSize*FTI_Topo.groupSize);
    if (FTI_CreateMatrix(matrix) == FTI_SCES) FTI_Print("Matrix created.", FTI_DBUG);
    sscanf(FTI_Exec.ckptFile,"Ckpt%d-Rank%d.fti", &FTI_Exec.ckptID, &i);
    sprintf(lfn,"%s/%s",FTI_Conf.lTmpDir, FTI_Exec.ckptFile);
    sprintf(efn,"%s/Ckpt%d-RSed%d.fti", FTI_Conf.lTmpDir, FTI_Exec.ckptID, i);
    sprintf(str, "L3 trying to access local ckpt. file (%s).", lfn); FTI_Print(str, FTI_DBUG);
    if (access(lfn, R_OK) != 0) { FTI_Print("L3 cannot access the checkpoint file.", FTI_DBUG); return FTI_NSCS; }
    if (truncate(lfn,ps) == -1) { FTI_Print("L3 cannot truncate the checkpoint file.", FTI_DBUG); return FTI_NSCS; }
    lfd = fopen(lfn, "rb"); efd = fopen(efn, "wb");
    if (lfd == NULL) { FTI_Print("L3 cannot open checkpoint file.", FTI_DBUG); return FTI_NSCS; }
    if (efd == NULL) { FTI_Print("L3 cannot open encoded ckpt. file.", FTI_DBUG); return FTI_NSCS; }
    while(pos < ps) { // For each block
        fread(myData, sizeof(char), bs, lfd); // Reading checkpoint files
        dest = FTI_Topo.groupRank; i = FTI_Topo.groupRank; offset = 0; init = 0; cnt = 0; // For the logic
        while(cnt < FTI_Topo.groupSize) { // For each encoding
            if (cnt == 0) memcpy(&(data[offset*bs]), myData, sizeof(char)*bs);
            else { MPI_Wait(&reqSend, &status); MPI_Wait(&reqRecv, &status); }
            if (cnt != FTI_Topo.groupSize-1) { // At every loop *but* the last one we send the data
                dest = (dest+FTI_Topo.groupSize-1)%FTI_Topo.groupSize;
                src = (i+1)%FTI_Topo.groupSize;
                MPI_Isend(myData, bs, MPI_CHAR, dest, FTI_Conf.tag, FTI_Exec.groupComm, &reqSend);
                MPI_Irecv(&(data[(1-offset)*bs]), bs, MPI_CHAR, src, FTI_Conf.tag, FTI_Exec.groupComm, &reqRecv);
            }
            matVal = matrix[FTI_Topo.groupRank*FTI_Topo.groupSize+i];
            if (matVal == 1) { // First copy or xor any data that does not need to be multiplied by a factor
                if (init == 0) { memcpy(coding, &(data[offset*bs]), bs); init = 1; }
                else galois_region_xor(&(data[offset*bs]), coding, coding, bs);
            }
            if (matVal != 0 && matVal != 1) // Then the data that needs to be multiplied by a factor
                { galois_w16_region_multiply(&(data[offset*bs]), matVal, bs, coding, init); init = 1; }
            i = (i+1)%FTI_Topo.groupSize; offset = 1 - offset; cnt++; // For the logic
        }
        fwrite(coding, sizeof(char), bs, efd); // Writting encoded checkpoints
        pos = pos + bs; // Next block
    }
    fclose(lfd); fclose(efd); // Closing and resizing files
    if (truncate(lfn,fs) == -1) { FTI_Print("L3 cannot re-truncate the checkpoint file.", FTI_DBUG); return FTI_NSCS; }
    if (truncate(efn,fs) == -1) { FTI_Print("L3 cannot re-truncate the encoded ckpt. file.", FTI_DBUG); return FTI_NSCS; }
    free(data); free(matrix); free(coding); free(myData); // Freeing memory
    return FTI_SCES;
}


/*-------------------------------------------------------------------------*/
/**
    @brief      It flushes the local ckpt. files in to the PFS.
    @param      group           The group ID.
    @param      level           The level from which ckpt. files are flushed.
    @return     integer         FTI_SCES if successful.

    This function flushes the local checkpoint files in to the PFS. The files
    are padded to match the block size and then truncated.

 **/
/*-------------------------------------------------------------------------*/
int FTI_Flush(int group, int level) {
    char        lfn[FTI_BUFS], gfn[FTI_BUFS], str[FTI_BUFS], *blBuf1 = talloc(char, FTI_Conf.blockSize);
    unsigned long maxFs, fs, ps, pos = 0; // FTI_BUFS > 32
    FILE        *lfd, *gfd;
    if (FTI_GetMeta(&fs, &maxFs, group, level) == FTI_SCES) FTI_Print("Metadata obtained.", FTI_DBUG);
    else { FTI_Print("Error getting metadata.", FTI_DBUG); return FTI_NSCS; }
    if (access(FTI_Conf.gTmpDir, F_OK) != 0) mkdir(FTI_Conf.gTmpDir, 0777);
    ps = (maxFs/FTI_Conf.blockSize)*FTI_Conf.blockSize;
    if (ps < maxFs) ps = ps + FTI_Conf.blockSize;
    switch(level) {
        case 0: sprintf(lfn,"%s/%s", FTI_Conf.lTmpDir, FTI_Exec.ckptFile); break;
        case 1: sprintf(lfn,"%s/%s", FTI_Ckpt[1].dir, FTI_Exec.ckptFile); break;
        case 2: sprintf(lfn,"%s/%s", FTI_Ckpt[2].dir, FTI_Exec.ckptFile); break;
        case 3: sprintf(lfn,"%s/%s", FTI_Ckpt[3].dir, FTI_Exec.ckptFile); break;
    }
    sprintf(gfn,"%s/%s", FTI_Conf.gTmpDir, FTI_Exec.ckptFile); // Open and resize files
    sprintf(str, "L4 trying to access local ckpt. file (%s).", lfn); FTI_Print(str, FTI_DBUG);
    if (access(lfn, R_OK) != 0) { FTI_Print("L4 cannot access the checkpoint file.", FTI_DBUG); return FTI_NSCS; }
    if (truncate(lfn,ps) == -1) { FTI_Print("L4 cannot truncate the checkpoint file.", FTI_DBUG); return FTI_NSCS; }
    lfd = fopen(lfn, "rb"); gfd = fopen(gfn, "wb");
    if (lfd == NULL) { FTI_Print("L4 cannot open the checkpoint file.", FTI_DBUG); return FTI_NSCS; }
    if (gfd == NULL) { FTI_Print("L4 cannot open ckpt. file in the PFS.", FTI_DBUG); return FTI_NSCS; }
    while(pos < ps) { // Checkpoint files exchange
        fread(blBuf1, sizeof(char), FTI_Conf.blockSize, lfd);
        fwrite(blBuf1, sizeof(char), FTI_Conf.blockSize, gfd);
        pos = pos + FTI_Conf.blockSize;
    }
    fclose(lfd); fclose(gfd); // Close files
    if (truncate(lfn,fs) == -1) { FTI_Print("L4 cannot re-truncate the checkpoint file.", FTI_DBUG); return FTI_NSCS; }
    if (truncate(gfn,fs) == -1) { FTI_Print("L4 cannot re-truncate the ckpt. file on PFS.", FTI_DBUG); return FTI_NSCS; }
    return FTI_SCES;
}


