/**
 *  @file   tools.c
 *  @author Leonardo A. Bautista Gomez (leobago@gmail.com)
 *  @date   July, 2013
 *  @brief  Utility functions for the FTI library.
 */


#include "fti.h"


/*-------------------------------------------------------------------------*/
/**
    @brief      Prints FTI messages.
    @param      s               Message to print.
    @param      priority        Priority of the message to be printed.
    @return     void

    This function prints messages depending on their priority and the
    verbosity level set by the user. DEBUG messages are printed by all
    processes with their rank. INFO messages are printed by one process.
    ERROR messages are printed and then the application is killed.

 **/
/*-------------------------------------------------------------------------*/
void FTI_Print(char *s, int priority) {
    if (priority >= FTI_Conf.verbosity) {
        if (s != NULL) {
            switch(priority) {
                case FTI_EROR:
                    fprintf(stderr, "[FTI Error - %06d] : %s : %s \n", FTI_Topo.myRank, s, strerror(errno));
                    FTI_Clean(5, FTI_Topo.groupID, FTI_Topo.myRank);
                    MPI_Abort(MPI_COMM_WORLD, -1);
                    MPI_Finalize();
                    exit(1);
                case FTI_INFO:
                    if (FTI_Topo.splitRank == 0)
                        fprintf(stdout, "[FTI Infos - %06d] : %s \n", FTI_Topo.myRank, s);
                    break;
                case FTI_DBUG:
                    fprintf(stdout, "[FTI Debug - %06d] : %s \n", FTI_Topo.myRank, s);
                    break;
                default:
                    break;
            }
        }
    }
    fflush(stdout);
}


/*-------------------------------------------------------------------------*/
/**
    @brief      Set the exec. ID and failure parameters in the conf. file.
    @param      restart         Value to set in the conf. file (0 or 1).
    @return     integer         FTI_SCES if successfull.

    This function sets the execution ID and failure parameters in the
    configuration file. This is to avoid forcing the user to change these
    values manually in case of recovery needed. In this way, relaunching the
    execution in the same way as the initial time will make FTI detect that
    it is a restart. It also allows to set the failure parameter back to 0
    at the end of a successful execution.

 **/
/*-------------------------------------------------------------------------*/
int FTI_UpdateConf(int restart) {
    char str[FTI_BUFS];
    dictionary *ini;
    ini = iniparser_load(FTI_Conf.cfgFile); // Load dictionary
    sprintf(str, "Updating configuration file (%s)...", FTI_Conf.cfgFile);
    FTI_Print(str, FTI_DBUG);
    if (ini != NULL) FTI_Print("Iniparser parsed the configuration file.", FTI_DBUG);
    else return FTI_NSCS;
    sprintf(str, "%d", restart);
    iniparser_set(ini, "Restart:failure", str); // Set failure to 'restart'
    iniparser_set(ini, "Restart:exec_id", FTI_Exec.id); // Set the execution ID
    FILE *fd = fopen(FTI_Conf.cfgFile, "w");
    if (fd != NULL) FTI_Print("Configuration file opened.", FTI_DBUG);
    else return FTI_NSCS;
    iniparser_dump_ini(ini, fd); // Write new configuration
    if (fflush(fd) == 0) FTI_Print("Configuration file flushed.", FTI_DBUG);
    else return FTI_NSCS;
    if (fclose(fd) == 0) FTI_Print("Configuration file closed.", FTI_DBUG);
    else return FTI_NSCS;
    iniparser_freedict(ini); // Free dictionary
    return FTI_SCES;
}


/*-------------------------------------------------------------------------*/
/**
    @brief      It creates and broadcast a global execution ID.
    @return     integer         FTI_SCES if successfull.

    This function creates and broadcast an execution ID, so that all ranks
    have the same execution ID.

 **/
/*-------------------------------------------------------------------------*/
int FTI_CreateExecID() {
    time_t tim = time(NULL);
    struct tm *now = localtime(&tim);
    snprintf(FTI_Exec.id, FTI_BUFS, "%d-%02d-%02d_%02d-%02d-%02d",
            now->tm_year+1900, now->tm_mon+1, now->tm_mday, now->tm_hour, now->tm_min, now->tm_sec);
    MPI_Bcast(FTI_Exec.id, FTI_BUFS, MPI_CHAR, 0, FTI_Exec.globalComm);
    return FTI_SCES;
}


/*-------------------------------------------------------------------------*/
/**
    @brief      It reads the configuration given in the configuration file.
    @return     integer         FTI_SCES if successfull.

    This function reads the configuration given in the FTI configuration
    file and sets other required parameters.

 **/
/*-------------------------------------------------------------------------*/
int FTI_ReadConf() {
    // Check access to FTI configuration file and load dictionary
    dictionary *ini;
    char *par, str[FTI_BUFS];
    sprintf(str, "Reading FTI configuration file (%s)...", FTI_Conf.cfgFile);
    FTI_Print(str, FTI_DBUG);
    if (access(FTI_Conf.cfgFile, F_OK) == 0) FTI_Print("FTI configuration file accessible.", FTI_DBUG);
    else FTI_Print("FTI configuration file NOT accessible.", FTI_EROR);
    ini = iniparser_load(FTI_Conf.cfgFile);
    if (ini != NULL) FTI_Print("Iniparser parsed the FTI configuration file.", FTI_DBUG);
    else FTI_Print("Iniparser could NOT parse the FTI configuration file.", FTI_EROR);

    // Setting/reading checkpoint configuration metadata
    par = iniparser_getstring(ini, "Basic:ckpt_dir", NULL);
    snprintf(FTI_Conf.localDir, FTI_BUFS, "%s", par);
    par = iniparser_getstring(ini, "Basic:glbl_dir", NULL);
    snprintf(FTI_Conf.glbalDir, FTI_BUFS, "%s", par);
    FTI_Ckpt[1].ckptIntv = (int) 1;
    FTI_Ckpt[2].ckptIntv = (int) iniparser_getint(ini, "Basic:ckpt_l2", -1);
    FTI_Ckpt[3].ckptIntv = (int) iniparser_getint(ini, "Basic:ckpt_l3", -1);
    FTI_Ckpt[4].ckptIntv = (int) iniparser_getint(ini, "Basic:ckpt_l4", -1);
    FTI_Ckpt[2].isInline = (int) iniparser_getint(ini, "Basic:inline_l2", 0);
    FTI_Ckpt[3].isInline = (int) iniparser_getint(ini, "Basic:inline_l3", 0);
    FTI_Ckpt[4].isInline = (int) iniparser_getint(ini, "Basic:inline_l4", 0);

    // Reading/setting configuration metadata
    par = iniparser_getstring(ini, "Basic:meta_dir", NULL);
    snprintf(FTI_Conf.metadDir, FTI_BUFS, "%s", par);
    FTI_Conf.verbosity = (int) iniparser_getint(ini, "Basic:verbosity", -1);
    FTI_Conf.saveLastCkpt = (int) iniparser_getint(ini, "Basic:keep_last_ckpt", 0);
    FTI_Conf.blockSize = (int) iniparser_getint(ini, "Advanced:block_size", -1) * 1024;
    FTI_Conf.tag = (int) iniparser_getint(ini, "Advanced:mpi_tag", -1);
    FTI_Conf.test = (int) iniparser_getint(ini, "Advanced:local_test", -1);
    FTI_Conf.l3WordSize = FTI_WORD;

    // Reading/setting execution metadata
    FTI_Exec.nbVar = 0;
    FTI_Exec.ckpt = 0;
    FTI_Exec.ckptCnt = 0;
    FTI_Exec.ckptID = 1;
    FTI_Exec.ckptLvel = 0;
    FTI_Exec.ckptIntv = (int) iniparser_getint(ini, "Basic:ckpt_int", -1);
    FTI_Exec.wasLastOffline = 0;
    FTI_Exec.ckptNext = 1;
    FTI_Exec.ckptLast = 0;
    FTI_Exec.syncIter = 1;
    FTI_Exec.lastIterTime = 0;
    FTI_Exec.totalIterTime = 0;
    FTI_Exec.meanIterTime = 0;
    FTI_Exec.reco = (int) iniparser_getint(ini, "restart:failure", 0);
    if (FTI_Exec.reco == 0) {
        FTI_CreateExecID();
        sprintf(str, "The execution ID is: %s", FTI_Exec.id);
        FTI_Print(str, FTI_INFO);
        snprintf(FTI_Exec.ckptFile, FTI_BUFS, "Ckpt%d-Rank%d.fti", FTI_Exec.ckptID, FTI_Topo.myRank);
    } else {
        par = iniparser_getstring(ini, "restart:exec_id", NULL);
        snprintf(FTI_Exec.id, FTI_BUFS, "%s", par);
        sprintf(str, "This is a restart. The execution ID is: %s", FTI_Exec.id);
        FTI_Print(str, FTI_INFO);
    }

    // Reading/setting topology metadata
    FTI_Topo.nbHeads = (int) iniparser_getint(ini, "Basic:head", 0);
    FTI_Topo.groupSize = (int) iniparser_getint(ini, "Basic:group_size", -1);
    FTI_Topo.nodeSize = (int) iniparser_getint(ini, "Basic:node_size", -1);
    FTI_Topo.nbApprocs = FTI_Topo.nodeSize - FTI_Topo.nbHeads;
    FTI_Topo.nbNodes = FTI_Topo.nbProc / FTI_Topo.nodeSize;

    // Synchronize after config reading and free dictionary
    MPI_Barrier(FTI_Exec.globalComm);
    iniparser_freedict(ini);
    return FTI_SCES;
}


/*-------------------------------------------------------------------------*/
/**
    @brief      It tests that the configuration given is correct.
    @return     integer         FTI_SCES if successfull.

    This function tests the FTI configuration to make sure that all
    parameter's values are correct.

 **/
/*-------------------------------------------------------------------------*/
int FTI_TestConfig() {
    // Checking topology
    char str[FTI_BUFS];

    if (FTI_Topo.nbHeads != 0 && FTI_Topo.nbHeads != 1)
        FTI_Print("The number of heads needs to be set to 0 or 1.", FTI_EROR);
    if (FTI_Topo.nbProc % FTI_Topo.nodeSize != 0)
        FTI_Print("Number of ranks is not a multiple of the node size.", FTI_EROR);
    if (FTI_Topo.nbNodes % FTI_Topo.groupSize != 0)
        FTI_Print("The group size is not multiple of the number of nodes.", FTI_EROR);
    if (FTI_Topo.groupSize <= 2)
        FTI_Print("The group size must be bigger than 2", FTI_EROR);
    if (FTI_Topo.groupSize >= 32)
        FTI_Print("The group size must be lower than 32", FTI_EROR);

    // Checking general configuration
    if (FTI_Conf.verbosity > 3 || FTI_Conf.verbosity < 1)
        FTI_Print("Verbosity needs to be set to 1, 2 or 3.", FTI_EROR);
    if (FTI_Conf.blockSize > (2048*1024) || FTI_Conf.blockSize < (1*1024))
        FTI_Print("Block size needs to be set between 1 and 2048.", FTI_EROR);
    if (FTI_Conf.test != 0 && FTI_Conf.test != 1)
        FTI_Print("Local test size needs to be set to 0 or 1.", FTI_EROR);
    if (FTI_Conf.saveLastCkpt != 0 && FTI_Conf.saveLastCkpt != 1)
        FTI_Print("Keep last ckpt. needs to be set to 0 or 1.", FTI_EROR);


    // Checking checkpoint configuration
    if (FTI_Ckpt[1].ckptIntv == 0) FTI_Ckpt[1].ckptIntv = -1;
    if (FTI_Ckpt[2].ckptIntv == 0) FTI_Ckpt[2].ckptIntv = -1;
    if (FTI_Ckpt[3].ckptIntv == 0) FTI_Ckpt[3].ckptIntv = -1;
    if (FTI_Ckpt[4].ckptIntv == 0) FTI_Ckpt[4].ckptIntv = -1;
    if (FTI_Ckpt[2].isInline != 0 && FTI_Ckpt[2].isInline != 1) FTI_Ckpt[2].isInline = 0;
    if (FTI_Ckpt[3].isInline != 0 && FTI_Ckpt[3].isInline != 1) FTI_Ckpt[3].isInline = 0;
    if (FTI_Ckpt[4].isInline != 0 && FTI_Ckpt[4].isInline != 1) FTI_Ckpt[4].isInline = 0;
    if (FTI_Ckpt[2].isInline == 0 && FTI_Topo.nbHeads != 1)
        FTI_Print("If inline_l2 is set to 0 then head should be set to 1.", FTI_EROR);
    if (FTI_Ckpt[3].isInline == 0 && FTI_Topo.nbHeads != 1)
        FTI_Print("If inline_l3 is set to 0 then head should be set to 1.", FTI_EROR);
    if (FTI_Ckpt[4].isInline == 0 && FTI_Topo.nbHeads != 1)
        FTI_Print("If inline_l4 is set to 0 then head should be set to 1.", FTI_EROR);

    // Checking metadata directory
    sprintf(str,"Checking the metadata directory (%s)...", FTI_Conf.metadDir);
    FTI_Print(str, FTI_DBUG);
    if (access(FTI_Conf.metadDir, W_OK) != 0) {
        FTI_Print("The metadata directory does not exist or has no write access.", FTI_DBUG);
        if (mkdir(FTI_Conf.metadDir, 0777) != 0)
            FTI_Print("The metadata directory could NOT be created.", FTI_EROR);
    }

    // Checking local directory
    sprintf(str,"Checking the local directory (%s)...", FTI_Conf.localDir);
    FTI_Print(str, FTI_DBUG);
    if (access(FTI_Conf.localDir, W_OK) != 0) {
        FTI_Print("The local directory does not exist or has no write access.", FTI_DBUG);
        if (mkdir(FTI_Conf.localDir, 0777) != 0)
            FTI_Print("The local directory could NOT be created.", FTI_EROR);
    }

    // Checking global directory
    sprintf(str,"Checking the global directory (%s)...", FTI_Conf.glbalDir);
    FTI_Print(str, FTI_DBUG);
    if (access(FTI_Conf.glbalDir, W_OK) != 0) {
        FTI_Print("The global directory does not exist or has no write access.", FTI_DBUG);
        if (mkdir(FTI_Conf.glbalDir, 0777) != 0)
            FTI_Print("The global directory could NOT be created.", FTI_EROR);
    }

    return FTI_SCES;
}


/*-------------------------------------------------------------------------*/
/**
    @brief      It creates the directories required for current execution.
    @return     integer         FTI_SCES if successfull.

    This function creates the temporary metadata, local and global
    directories required for the current execution.

 **/
/*-------------------------------------------------------------------------*/
int FTI_CreateDirs() {
    // Create metadata timestamp directory
    char fn[FTI_BUFS];
    snprintf(fn, FTI_BUFS, "%s/%s", FTI_Conf.metadDir, FTI_Exec.id);
    if (access(fn, F_OK) != 0) mkdir(fn, 0777);
    snprintf(FTI_Conf.metadDir, FTI_BUFS, "%s", fn);
    snprintf(FTI_Conf.mTmpDir, FTI_BUFS, "%s/tmp", fn);
    snprintf(FTI_Ckpt[1].metaDir, FTI_BUFS, "%s/l1", fn);
    snprintf(FTI_Ckpt[2].metaDir, FTI_BUFS, "%s/l2", fn);
    snprintf(FTI_Ckpt[3].metaDir, FTI_BUFS, "%s/l3", fn);
    snprintf(FTI_Ckpt[4].metaDir, FTI_BUFS, "%s/l4", fn);
    // Create global checkpoint timestamp directory
    snprintf(fn, FTI_BUFS, "%s", FTI_Conf.glbalDir);
    snprintf(FTI_Conf.glbalDir, FTI_BUFS, "%s/%s", fn, FTI_Exec.id);
    if (access(FTI_Conf.glbalDir, F_OK) != 0) mkdir(FTI_Conf.glbalDir, 0777);
    snprintf(FTI_Conf.gTmpDir, FTI_BUFS, "%s/tmp", FTI_Conf.glbalDir);
    snprintf(FTI_Ckpt[4].dir, FTI_BUFS, "%s/l4", FTI_Conf.glbalDir);
    // Create local checkpoint timestamp directory
    if (FTI_Conf.test) { // If local test generate name by topology
        snprintf(fn, FTI_BUFS, "%s/node%d", FTI_Conf.localDir, FTI_Topo.myRank/FTI_Topo.nodeSize);
        if (access(fn, F_OK) != 0) mkdir(fn, 0777);
    } else {
        snprintf(fn, FTI_BUFS, "%s", FTI_Conf.localDir);
    }
    snprintf(FTI_Conf.localDir, FTI_BUFS, "%s/%s", fn, FTI_Exec.id);
    if (access(FTI_Conf.localDir, F_OK) != 0) mkdir(FTI_Conf.localDir, 0777);
    snprintf(FTI_Conf.lTmpDir, FTI_BUFS, "%s/tmp", FTI_Conf.localDir);
    snprintf(FTI_Ckpt[1].dir, FTI_BUFS, "%s/l1", FTI_Conf.localDir);
    snprintf(FTI_Ckpt[2].dir, FTI_BUFS, "%s/l2", FTI_Conf.localDir);
    snprintf(FTI_Ckpt[3].dir, FTI_BUFS, "%s/l3", FTI_Conf.localDir);
    return FTI_SCES;
}


/*-------------------------------------------------------------------------*/
/**
    @brief      It erases the previous checkpoints and their metadata.
    @param      level           Level of cleaning.
    @param      group           Group ID of the cleanning target process.
    @param      rank            Rank of the cleanning target process.
    @return     integer         FTI_SCES if successfull.

    This function erases previous checkpoint depending on the level of the
    current checkpoint. Level 5 means complete clean up. Level 6 means clean
    up local nodes but keep last checkpoint data and metadata in the PFS.

 **/
/*-------------------------------------------------------------------------*/
int FTI_Clean(int level, int group, int rank) {
    int ckptID;
    char buf[FTI_BUFS];
    sprintf(buf, "Clean level %d, group %d and rank %d.", level, group, rank); FTI_Print(buf, FTI_DBUG);
    if (level == 0) {
        snprintf(buf, FTI_BUFS, "%s/sector%d-group%d.fti", FTI_Conf.mTmpDir, FTI_Topo.sectorID, group); remove(buf);
        snprintf(buf, FTI_BUFS, "%s/Ckpt%d-Rank%d.fti", FTI_Conf.gTmpDir, FTI_Exec.ckptID, rank); remove(buf);
        snprintf(buf, FTI_BUFS, "%s/Ckpt%d-Rank%d.fti", FTI_Conf.lTmpDir, FTI_Exec.ckptID, rank); remove(buf);
        snprintf(buf, FTI_BUFS, "%s/Ckpt%d-Pcof%d.fti", FTI_Conf.lTmpDir, FTI_Exec.ckptID, rank); remove(buf);
        snprintf(buf, FTI_BUFS, "%s/Ckpt%d-RSed%d.fti", FTI_Conf.lTmpDir, FTI_Exec.ckptID, rank); remove(buf);
        rmdir(FTI_Conf.mTmpDir);
        rmdir(FTI_Conf.gTmpDir);
        rmdir(FTI_Conf.lTmpDir);
    }
    if (level >= 1) {
        snprintf(buf, FTI_BUFS, "%s/sector%d-group%d.fti", FTI_Ckpt[1].metaDir, FTI_Topo.sectorID, group); remove(buf);
        snprintf(buf, FTI_BUFS, "%s/Ckpt%d-Rank%d.fti",FTI_Ckpt[1].dir, FTI_Ckpt[1].lastCkpt, rank); remove(buf);
        rmdir(FTI_Ckpt[1].metaDir);
        rmdir(FTI_Ckpt[1].dir);
    }
    if (level >= 2) {
        snprintf(buf, FTI_BUFS, "%s/sector%d-group%d.fti", FTI_Ckpt[2].metaDir, FTI_Topo.sectorID, group); remove(buf);
        snprintf(buf, FTI_BUFS, "%s/Ckpt%d-Rank%d.fti",FTI_Ckpt[2].dir, FTI_Ckpt[2].lastCkpt, rank); remove(buf);
        snprintf(buf, FTI_BUFS, "%s/Ckpt%d-Pcof%d.fti",FTI_Ckpt[2].dir, FTI_Ckpt[2].lastCkpt, rank); remove(buf);
        rmdir(FTI_Ckpt[2].metaDir);
        rmdir(FTI_Ckpt[2].dir);
    }
    if (level >= 3) {
        snprintf(buf, FTI_BUFS, "%s/sector%d-group%d.fti", FTI_Ckpt[3].metaDir, FTI_Topo.sectorID, group); remove(buf);
        snprintf(buf, FTI_BUFS, "%s/Ckpt%d-Rank%d.fti",FTI_Ckpt[3].dir, FTI_Ckpt[3].lastCkpt, rank); remove(buf);
        snprintf(buf, FTI_BUFS, "%s/Ckpt%d-RSed%d.fti",FTI_Ckpt[3].dir, FTI_Ckpt[3].lastCkpt, rank); remove(buf);
        rmdir(FTI_Ckpt[3].metaDir);
        rmdir(FTI_Ckpt[3].dir);
    }
    if (level >= 4) {
        snprintf(buf, FTI_BUFS, "%s/sector%d-group%d.fti", FTI_Ckpt[4].metaDir, FTI_Topo.sectorID, group); remove(buf);
        snprintf(buf, FTI_BUFS, "%s/Ckpt%d-Rank%d.fti",FTI_Ckpt[4].dir, FTI_Ckpt[4].lastCkpt, rank); remove(buf);
        rmdir(FTI_Ckpt[4].metaDir);
        rmdir(FTI_Conf.gTmpDir);
        rmdir(FTI_Ckpt[4].dir);
    }
    if (level == 5) { // If it is the very last cleaning
        rmdir(FTI_Conf.localDir);
        rmdir(FTI_Conf.glbalDir);
        snprintf(buf, FTI_BUFS, "%s/Topology.fti",FTI_Conf.metadDir); remove(buf);
        rmdir(FTI_Conf.metadDir);
    }
    if (level == 6) rmdir(FTI_Conf.localDir);
    return FTI_SCES;
}


/*-------------------------------------------------------------------------*/
/**
    @brief      It listens for checkpoint notifications.
    @return     integer         FTI_SCES if successfull.

    This function listens for notifications from the application processes
    and take the required actions after notification. This function is only
    executed by the head of the nodes.

 **/
/*-------------------------------------------------------------------------*/
int FTI_Listen() {
    MPI_Status status;
    char str[FTI_BUFS];
    int i, buf, flags[6];
    while (1) { // Each loop is a checkpoint interval
        for (i = 0; i < 6; i++) flags[i] = 0;
        FTI_Print("Head listening...", FTI_DBUG);
        for(i = 0; i < FTI_Topo.nbApprocs; i++) { // Iterate on the application processes in the node
            MPI_Recv(&buf, 1, MPI_INT, FTI_Topo.body[i], FTI_Conf.tag, FTI_Exec.globalComm, &status);
            sprintf(str, "The head received a %d message", buf); FTI_Print(str, FTI_DBUG);
            flags[buf-FTI_BASE] = flags[buf-FTI_BASE] + 1;
        }
        for (i = 1; i < 6; i++) if (flags[i] == FTI_Topo.nbApprocs) FTI_Exec.ckptLvel = i;
        if (FTI_Exec.ckptLvel == 5) break;
        if (FTI_CkptNotified() == FTI_SCES) FTI_Print("Checkpoint successfully done.", FTI_DBUG);
    }
    return FTI_SCES;
}


