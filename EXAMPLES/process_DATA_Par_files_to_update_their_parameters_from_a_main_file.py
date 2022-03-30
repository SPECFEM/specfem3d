#!/usr/bin/env python
#
# reads in a Par_file as main (template) and updates parameters and comments in all other Par_files
# in current directory to have a consistent set of Par_files
#
from __future__ import print_function

import sys
import os
import collections

#
#----------------------------------------------------------------------------
#

# USER PARAMETERS

# deprecated parameter names which have been renamed
DEPRECATED_RENAMED_PARAMETERS = [ \
  ("ABSORBING_CONDITIONS", "STACEY_ABSORBING_CONDITIONS"), \
  ("ABSORB_INSTEAD_OF_FREE_SURFACE", "STACEY_INSTEAD_OF_FREE_SURFACE"), \
  ("OCEANS", "APPROXIMATE_OCEAN_LOAD"), \
  ("NZ_DOUGLING_1", "NZ_DOUBLING_1"), \
  ("NZ_DOUGLING_2", "NZ_DOUBLING_2"), \
  ("EXTERNAL_SOURCE_FILE","USE_EXTERNAL_SOURCE_FILE"),
  ("WRITE_SEISMOGRAMS_BY_MASTER","WRITE_SEISMOGRAMS_BY_MAIN") \
]

# exclude other possible files with similar name, but with different format
# (more file names can be appended in this list)
EXCLUDE_NAME_LIST = [ \
  "Par_file_faults", \
  "Par_file.yml" \
]

# exclude old unused directories from search directories
# (more exclude directories can be appended in this list)
EXCLUDE_DIR_LIST = [ \
  "unused_routines", \
  "small_SEM_solvers_in_Fortran_and_C_without_MPI_to_learn" \
]

#
#----------------------------------------------------------------------------
#

# more global parameters

# ordered dictionary: ordering is kept, what is filled in first, will be listed first
#
# each dictionary entry will have format (value,comment,appendix), use for example:
# (value,comment,appendix) = main_parameters[name]
#
main_parameters = collections.OrderedDict()

# Mesh_Par_file data lines (NMATERIALS and NREGIONS sections)
mesh_par_file_data_counter = 0

# parameter file type
is_Mesh_Par_file = False

#
#----------------------------------------------------------------------------
#

def read_Par_file_sections(parameters,file,verbose=False):
    """
    reads a Par_file and fills in parameter sections
    """
    global mesh_par_file_data_counter
    global is_Mesh_Par_file

    # user info
    if verbose:
        print("reading in file: ",file)
        print("")

    # checks file string
    if len(file) <= 1:
        print("invalid file name: %s, please check..." % file)
        sys.tracebacklimit=0
        raise Exception('invalid file name: %s' % file)

    # checks if file exists
    if not os.path.isfile(file):
        print("file %s does not exist, please check..." % file)
        sys.tracebacklimit=0
        raise Exception('file does not exist: %s' % file)

    # opens file
    try:
        f = open(file,'r')
    except:
        print("Error opening file ",file)
        sys.tracebacklimit=0
        raise Exception('file does not open: %s' % file)

    # reads in dictionaries
    nsections = 0
    comment = ''
    mesh_par_file_data_counter = 0

    for line in f:
        dataline = line.strip()
        #print("debug: dataline = ",dataline)

        if dataline:
            # line with some data
            if dataline.startswith('#'):
                # comment line
                # for example: # simulation input parameters
                comment += dataline + '\n'
            else:
                # parameter line (eventually with comment appended)
                # for example: SIMULATION_TYPE                 = 1  # or 2 # or 3
                par_string = ""
                app_string = ""

                # separates parameter from appended comment(s)
                index_app = dataline.find('#')

                if index_app == 0:
                    # should have been a comment, notify user
                    print("Error parameter line: ",dataline)
                    print("A line starting with a # sign should be a comment line")
                    sys.tracebacklimit=0
                    raise Exception('Invalid parameter line: %s' % dataline)
                elif index_app > 0:
                    # separates parameter part from rest
                    par_string = dataline[0:index_app-1]
                    par_string = par_string.strip()
                    # appended part
                    app_string = dataline[index_app:]
                    app_string = app_string.strip()
                else:
                    # no appended comment found
                    par_string = dataline

                # first part holds parameter
                tokens = par_string.split('=')

                # checks format
                wrong_format = False
                if len(tokens) != 2:
                    wrong_format = True
                    # check with Mesh_Par_file format for NMATERIALS and NREGIONS
                    if is_Mesh_Par_file:
                        dataitems = par_string.split()
                        # check number of items
                        # for example 7-entries:
                        # #NEX_XI_BEGIN  #NEX_XI_END  #NEX_ETA_BEGIN  #NEX_ETA_END  #NZ_BEGIN #NZ_END  #material_id
                        # 1              128          1               128           1         5        14
                        # for example 8-entries: (default line)
                        # #material_id  #rho  #vp  #vs  #Q_Kappa  #Q_mu  #anisotropy_flag  #domain_id
                        # 1             2120  2325 1360 9999.0    35.0   0                 2
                        # for example 17-entries: (poroelastic line)
                        # #material_id #rho_s #rho_f #phi #tort #kxx #kxy #kxz #kyy #kyz #kzz #kappa_s #kappa_f #kappa_fr #eta #mu_fr #domain_id
                        # 2            2500   1020   0.4  2.0   1d-11 0.0 0.0  1d-11 0.0 1d-11 16.0554d9 2.295d9 10.0d9   0.0  9.63342d9 3
                        #debug
                        #print("debug: dataitems ",dataitems,"length: ",len(dataitems))
                        if len(dataitems) == 7 or len(dataitems) == 8 or len(dataitems) == 17:
                            # this is a data format line
                            wrong_format = False
                if wrong_format:
                    print("Error invalid format on parameter line: \n",dataline)
                    print("\nPlease use a line format: PARAMETER_NAME = value")
                    sys.tracebacklimit=0
                    raise Exception('Invalid parameter line: %s' % dataline)

                # determines format
                is_Par_file_format = True
                if is_Mesh_Par_file:
                    if par_string.find('=') < 0:
                        is_Par_file_format = False
                        index_app = 0

                # Par_file format: PARAMETER_NAME = value #..
                if is_Par_file_format:
                    name = tokens[0].strip()
                    value = tokens[1].strip()
                else:
                    # sets as a Mesh_Par_file data line
                    mesh_par_file_data_counter += 1
                    name = "MESH_PAR_FILE_DATA" + str(mesh_par_file_data_counter)
                    value = par_string

                # appended comment
                if len(app_string) > 0:
                    appendix = app_string
                else:
                    appendix = ''

                # checks that parameter name is unique
                if name in parameters.keys():
                    print("Error parameter name already used: ",name)
                    print("See section: \n")
                    (value,comment,appendix) = parameters[name]
                    print("%s" % comment)
                    print("%s = %s" % (name,value))
                    print("\nPlease verify that parameter names are unique.")
                    sys.tracebacklimit=0
                    raise Exception('Invalid parameter name: %s' % name)

                # removes last newline character from comment lines
                comment = comment.rstrip()

                # stores this as a new section
                parameters.update({name : (value,comment,appendix)})

                # starts a new section after this parameter
                nsections += 1
                comment = ''

                # user info
                if verbose:
                    print("  done section for parameter: ",name)
                    #print("starting new section      : ",nsections)
                    #print("")

        else:
            # empty line
            comment += '\n'

    f.close()

    # checks number of entries
    if len(parameters.keys()) != nsections:
        print("Error number of sections ",nsections," but number of keys is ",len(parameters.keys()))
        sys.tracebacklimit=0
        raise Exception('number of sections read in invalid: %s' % file)

    # user info
    if verbose:
        print("")
        print("  got %d parameters" % nsections)
        print("")
        #print("main file sections")
        #print("")
        #print(parameters)


#
#----------------------------------------------------------------------------
#

def get_maximum_parameter_name_length(parameters,verbose=False):
    """
    determines maximum parameter name length for formatting

    example:
        max_name_length = get_maximum_parameter_name_length(parameters)
    """
    # determines maximum parameter name length
    max_name_length = 0
    for name in parameters.keys():
        #print("parameter name: ",name)
        if len(name) > max_name_length: max_name_length = len(name)

    # user info
    if verbose:
        print("parameter names:")
        print("  maximum name length: ",max_name_length)

    # restrict to a minimum of 32 character until = sign, we will add 1 space when writing out line below
    if max_name_length <= 32-1:
      max_name_length = 32-1
    else:
      # add additional white space
      max_name_length += 1

    # user info
    if verbose:
        print("  using parameter format length: ",max_name_length," + 1")
        print("")

    return max_name_length

#
#----------------------------------------------------------------------------
#

def write_template_file(parameters,tmp_file,verbose=False):
    """
    creates a new template file with parameter sections
    """
    # checks if anything to do
    if len(parameters.keys()) == 0:
        print("Error no parameters to write out into template file:",parameters)
        sys.tracebacklimit=0
        raise Exception('invalid parameters dictionary: %s' % tmp_file)

    # determines maximum parameter name length
    max_name_length = get_maximum_parameter_name_length(parameters,verbose=verbose)

    # opens file
    try:
        f = open(tmp_file,'w')
    except:
        print("Error opening file ",tmp_file)
        sys.tracebacklimit=0
        raise Exception('file does not open: %s' % tmp_file)

    # writes out parameter sections
    for name in parameters.keys():
        # gets associated values
        (value,comment,appendix) = parameters[name]

        # prints out parameter section lines
        # comments
        if comment:
            f.write( "%s\n" % comment )

        # writes out line with values
        if not "MESH_PAR_FILE_DATA" in name:
            # parameter line
            if appendix:
                f.write( "%s = %s   %s\n" % (name.ljust(max_name_length),value,appendix) )
            else:
                f.write( "%s = %s\n" % (name.ljust(max_name_length),value) )
        else:
            # Mesh_Par_file data line
            if appendix:
                f.write( "%s   %s\n" % (value,appendix) )
            else:
                f.write( "%s\n" % value )

    f.write( "\n" )
    f.close()

    # user info
    if verbose:
        print("written temporary file, see: ",tmp_file)
        print("")

#
#----------------------------------------------------------------------------
#

def get_files_in_subdirectories(dir,files,basename):
    """
    recursive function to list all files in subdirectories
    """
    global EXCLUDE_NAME_LIST
    global EXCLUDE_DIR_LIST

    # loops over all files in this directory
    for subdir in os.listdir(dir):
        path = os.path.join(dir, subdir)
        if not os.path.isdir(path):
            filename = os.path.basename(path)
            # adds filename to file list
            if filename.startswith(basename):
                # checks with exclude list of name
                do_add_me = True
                for name in EXCLUDE_NAME_LIST:
                    #print("name: ",name," filename: ",filename)
                    if filename.startswith(name): do_add_me = False
                # adds to list
                if do_add_me:
                    files.append(path)
        else:
          do_search = True
          # checks with list of exclude directories
          for name in EXCLUDE_DIR_LIST:
              if subdir == name: do_search = False
          if do_search:
              get_files_in_subdirectories(path,files,basename)


#
#----------------------------------------------------------------------------
#

def compare_and_replace_file(main_file,temp_file,verbose=False,replace=False):
    """
    notifies user if template is different than main
    e.g. by different indentation or white space)
    """
    # diff files
    command = "diff " + main_file + " " + temp_file
    #print command
    ret = os.system(command)
    if ret > 0:
        # got some differences
        # replaces main with new file format
        print("")
        print("replacing main with new parameter file format:")
        print("")

        # copy over
        if replace:
            command = "cp -v " + temp_file + " " + main_file
            #command = "echo hello"
            #print command
            ret = os.system(command)
            if ret != 0:
                print("Error replacing file with new format",main_file)
                sys.tracebacklimit=0
                raise Exception('file can not be updated: %s' % main_file)
    else:
        if verbose:
            print("  no differences")
            print("  main file format is okay")
            print("")


def check_and_update_Par_file(my_parameters,file):
    """
    updates parameter entries
    """
    global DEPRECATED_RENAMED_PARAMETERS
    global main_parameters
    global is_Mesh_Par_file

    # checks for old, deprecated parameters
    nold_parameters = 0
    #print("  searching deprecated parameters...")
    my_parameters_new = my_parameters.copy()

    for name in my_parameters.keys():
        if not name in main_parameters.keys():
            if (not "MESH_PAR_FILE_DATA" in name) and (not "NZ_DOUBLING" in name):
                print("  deprecated parameter: ",name)
                nold_parameters += 1

                is_found = False
            else:
                # ignore data line
                is_found = True

            # converts old to new parameter name
            is_unique = False
            for old_name,new_name in DEPRECATED_RENAMED_PARAMETERS:
                if name in old_name :
                    # renames old (e.g. enreg_surf_same_vertical1 to record_at_surface_same_vertical1
                    if name == old_name:
                        name2 = new_name
                        is_unique = True
                    else:
                        name2 = name[0:len(old_name)] + name[len(old_name)+1:len(name)]
                    print("    will be converted to ",new_name)
                    (val_orig,comment_orig,appendix_orig) = my_parameters[name]
                    # stores as new section
                    my_parameters_new[name2] = (val_orig,comment_orig,appendix_orig)
                    is_found = True
                    # removes old section
                    del my_parameters_new[old_name]
                    # done replacing
                    if is_unique: break

            # removes old parameter
            if not is_found:
                # removes old section
                del my_parameters_new[name]

    # replace with new one
    if my_parameters_new != my_parameters: my_parameters = my_parameters_new.copy()
    #print(my_parameters.keys())

    # add missing main parameters and replaces comment lines to compare sections
    nmissing_parameters = 0
    #print("  searching missing parameters...")
    for name in main_parameters.keys():
        if not "MESH_PAR_FILE_DATA" in name:
            # checks if missing
            if not name in my_parameters.keys():
                print("  misses parameter: ",name)
                nmissing_parameters += 1
                # adds from main template record
                (val,comment,appendix) = main_parameters[name]
                my_parameters[name] = (val,comment,appendix)

    # updates comments
    nold_comments = 0
    for name in main_parameters.keys():
        if not "MESH_PAR_FILE_DATA" in name:
            # checks we have this parameter
            if not name in my_parameters.keys():
                print("Error comparing main with current file format parameter",name)
                sys.tracebacklimit=0
                raise Exception('parameter list invalid: %s' % file)

            # compares and replaces comments and appendix
            (val_orig,comment_orig,appendix_orig) = my_parameters[name]
            (val,comment,appendix) = main_parameters[name]
            if comment_orig != comment or appendix != appendix_orig:
                nold_comments += 1
                # replace with new comment/appendix and only keep original value
                my_parameters[name] = (val_orig,comment,appendix)

    # respect ordering
    ordered_parameters = collections.OrderedDict()
    iorder_new = 0
    num_material_entries = 0
    num_region_entries = 0
    for name in main_parameters.keys():
        if not is_Mesh_Par_file:
            # regular Par_file must have a one-to-one match
            # checks that name is available
            if not name in my_parameters.keys():
                print("Error ordering with current file format parameter",name)
                sys.tracebacklimit=0
                raise Exception('parameter list invalid: %s' % file)
            # get values
            (val,comment,appendix) = my_parameters[name]

            # put into same order of appearance
            ordered_parameters[name] = (val,comment,appendix)
        else:
            # Mesh_Par_file can have different number of lines for data ranges
            # new parameter file entry
            if (not "MESH_PAR_FILE_DATA" in name) and (not "NZ_DOUBLING" in name):
                # checks that name is available
                if not name in my_parameters.keys():
                    print("Error ordering with current file format parameter",name)
                    sys.tracebacklimit=0
                    raise Exception('parameter list invalid: %s' % file)

                # get values
                (val,comment,appendix) = my_parameters[name]

                # put into same order of appearance
                ordered_parameters[name] = (val,comment,appendix)
                iorder_new += 1

            else:
                if "NZ_DOUBLING" in name:
                    # get doubling layer values
                    for my_key in my_parameters.keys():
                        if "NZ_DOUBLING" in my_key:
                            # add entries
                            if not my_key in ordered_parameters.keys():
                                (val,comment,appendix) = my_parameters[my_key]
                                # put into same order of appearance
                                ordered_parameters[my_key] = (val,comment,appendix)
                                iorder_new += 1

                else:
                    # first NMATERIALS
                    previous = list(ordered_parameters.keys())[iorder_new-1]
                    (val,comment,appendix) = ordered_parameters[previous]
                    # number of data lines for nmaterials
                    if "NMATERIALS" in previous:
                        num_material_entries = int(val)
                        for i in range(0,num_material_entries):
                            my_key = "MESH_PAR_FILE_DATA" + str(i+1)
                            # check availability
                            if not my_key in my_parameters.keys():
                                print("Error ordering key",my_key," with previous",previous)
                                sys.tracebacklimit=0
                                raise Exception('ordering list failed: %s' % file)
                            # add entries
                            if not my_key in ordered_parameters.keys():
                                (val,comment,appendix) = my_parameters[my_key]
                                # put into same order of appearance
                                ordered_parameters[my_key] = (val,comment,appendix)
                                iorder_new += 1
                    # number of data lines for nregions
                    if "NREGIONS" in previous:
                        # check if materials are done
                        if num_material_entries == 0:
                            print("Error ordering with current file format parameter",name," with previous",previous)
                            print("NREGIONS section must come after NMATERIALS section\n")
                            sys.tracebacklimit=0
                            raise Exception('ordering list failed: %s' % file)
                        # adds regions
                        num_region_entries = int(val)
                        for i in range(0,num_region_entries):
                            my_key = "MESH_PAR_FILE_DATA" + str(i + 1 + num_material_entries)
                            # check availability
                            if not my_key in my_parameters.keys():
                                print("Error ordering key",my_key," with previous",previous)
                                sys.tracebacklimit=0
                                raise Exception('ordering list failed: %s' % file)
                            # add entries
                            if not my_key in ordered_parameters.keys():
                                (val,comment,appendix) = my_parameters[my_key]
                                # put into same order of appearance
                                ordered_parameters[my_key] = (val,comment,appendix)
                                iorder_new += 1


    # check
    if is_Mesh_Par_file:
        if num_material_entries == 0 or num_region_entries == 0:
            print("Error ordering Mesh_Par_file ",file)
            sys.tracebacklimit=0
            raise Exception('ordering list failed: %s' % file)


    # check length
    if len(ordered_parameters.keys()) != len(my_parameters.keys()):
        print("Error ordering with parameters got different lengths:")
        print("new length is ",len(ordered_parameters.keys())," instead of ",len(my_parameters.keys()))
        sys.tracebacklimit=0
        raise Exception('parameter list invalid: %s' % file)

    # checks if order changed
    nold_order = 0
    for i in range(0,len(my_parameters.keys())):
        if list(my_parameters.keys())[i] != list(ordered_parameters.keys())[i]: nold_order += 1

    if nold_order > 0:
        print("  needs re-ordering...",nold_order)

    # replace old file if necessary
    if nold_parameters == 0 and nmissing_parameters == 0 and nold_comments == 0 and nold_order == 0:
        # user info
        print("  file is okay and up-to-date")
    else:
        # user info
        print("  updating parameter file...")

        # opens temporary file with main info
        tmp_file = "_____temp09_____"
        write_template_file(ordered_parameters,tmp_file)

        # notifies user if template is different than main
        # (e.g. by different indentation or white space)
        compare_and_replace_file(file,tmp_file,verbose=True,replace=replace)

        # clean up temporary file
        command = "rm -f " + tmp_file
        os.system(command)

    # frees new order
    del ordered_parameters


def check_parameter_file_type(file):
    """
    determines flag for Par_file or Mesh_Par_file
    """
    global is_Mesh_Par_file

    # determines parameter file typ
    # (line formats may differ)
    basename = os.path.basename(file)
    if "Mesh_Par_file" in basename:
        is_Mesh_Par_file = True
    else:
        is_Mesh_Par_file = False


#
#----------------------------------------------------------------------------
#

def update_Par_files(main_file,replace=False):
    """
    uses a main to update other parameter files
    """
    global main_parameters
    global is_Mesh_Par_file

    # user info
    print("")
    print("main file: ",main_file)
    print("")

    # determines file type
    check_parameter_file_type(main_file)

    # reads in parameters
    read_Par_file_sections(main_parameters,main_file,verbose=True)

    # opens temporary file with main info
    tmp_file = "_____temp01_____"
    write_template_file(main_parameters,tmp_file,verbose=True)

    # notifies user if template is different than main
    # (e.g. by different indentation or white space)
    print("checking differences between new format and main:")
    print("  (different formatting or whitespace can lead to differences)")
    print("")
    compare_and_replace_file(main_file,tmp_file,verbose=True,replace=replace)

    # clean up temporary file
    command = "rm -f " + tmp_file
    os.system(command)

    # finds all Par_files
    basename = os.path.basename(main_file)
    current_dir = os.getcwd()

    # user info
    print("")
    print("finding all files with name: ",basename)
    print("in current directory: ",current_dir)

    # gets all files in subdirectories
    files = []
    get_files_in_subdirectories("./",files,basename)

    nfiles = len(files)
    print("")
    print("found ",nfiles," parameter files")
    print("")

    ifile = 0
    for file in files:
        # updates counter
        ifile += 1

        # determines parameter file typ
        check_parameter_file_type(file)

        # user info
        print("")
        print("file ",ifile," out of ",nfiles)
        print("processing file: ",file)
        if is_Mesh_Par_file:
            print("  file type is ","Mesh_Par_file")
        else:
            print("  file type is ","Par_file")

        # read in parameters
        my_parameters = collections.OrderedDict()
        read_Par_file_sections(my_parameters,file)

        # check and update
        check_and_update_Par_file(my_parameters,file)
        # clean up
        del my_parameters

    # user info
    print("")
    print("done")
    os.system("date")
    print("")

#
#----------------------------------------------------------------------------
#

def usage():
    print("usage:")
    print("    ./process_DATA_Par_files_to_update_their_parameters_from_a_main.py Main-Par-file replace")
    print("  with")
    print("    Main-Par_file - Par_file which serves as main template (e.g. DATA/Par_file)")
    print("    replace         = flag to force replacing of file [0==check-only/1==replace]")

#
#----------------------------------------------------------------------------
#

if __name__ == '__main__':
    # gets arguments
    if len(sys.argv) != 3:
        usage()
        sys.exit(1)
    else:
        main_file = sys.argv[1]
        if int(sys.argv[2]) == 1:
            replace = True
        else:
            replace = False

    update_Par_files(main_file,replace=replace)

