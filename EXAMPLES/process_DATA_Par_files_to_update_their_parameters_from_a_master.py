#!/usr/bin/env python
#
# reads in a Par_file as master (template) and updates parameters and comments in all other Par_files 
# in current directory to have a consistent set of Par_files
# 
import sys
import os
import collections

# ordered dictionary: ordering is kept, what is filled in first, will be listed first
#
# each dictionary entry will have format (value,comment,appendix), use for example:
# (value,comment,appendix) = master_parameters[name]
#
master_parameters = collections.OrderedDict()

# deprecated parameter names which have been renamed
deprecated_renamed_parameters = [ \
  ("ABSORBING_CONDITIONS", "STACEY_ABSORBING_CONDITIONS"), \
  ("ABSORB_INSTEAD_OF_FREE_SURFACE", "STACEY_INSTEAD_OF_FREE_SURFACE"), \
  ("OCEANS", "APPROXIMATE_OCEAN_LOAD"), \
  ]


def read_Par_file_sections(parameters,file,verbose=False):
    """
    reads a Par_file and fills in parameter sections
    """
    # user info
    if verbose:
        print "reading in file: ",file
        print ""

    # checks file string
    if len(file) <= 1:
        print "invalid file name: %s, please check..." % file
        sys.tracebacklimit=0
        raise Exception('invalid file name: %s' % file)

    # checks if file exists
    if not os.path.isfile(file):
        print "file %s does not exist, please check..." % file
        sys.tracebacklimit=0
        raise Exception('file does not exist: %s' % file)

    # opens file
    try:
        f = open(file,'r')
    except:
        print "Error opening file ",file
        sys.tracebacklimit=0
        raise Exception('file does not open: %s' % file)

    # reads in dictionaries
    nsections = 0  
    comment = ''

    for line in f:
        dataline = line.strip()
        #print "line: ",dataline

        if dataline:
            # line with some data
            if dataline.startswith('#'):
                # comment line
                # for example: # simulation input parameters
                comment += dataline + '\n'
            else:
                # parameter line (eventually with comment appended)
                # for example: SIMULATION_TYPE                 = 1  # or 2 # or 3

                # separates parameter from appended comment(s)
                index_app = dataline.find('#')

                if index_app == 0:
                    # should have been a comment, notify user
                    print "Error parameter line: ",dataline
                    print "A line starting with a # sign should be a comment line"
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
                if len(tokens) != 2:
                    print "Error invalid format on parameter line: \n",dataline
                    print "\nPlease use a line format: PARAMETER_NAME = value"
                    sys.tracebacklimit=0         
                    raise Exception('Invalid parameter line: %s' % dataline)

                name = tokens[0].strip()
                value = tokens[1].strip()

                # appended comment
                if index_app > 0:
                    appendix = app_string
                else:
                    appendix = ''

                # checks that parameter name is unique
                if name in parameters.keys():
                    print "Error parameter name already used: ",name
                    print "See section: \n"
                    (value,comment,appendix) = parameters[name]
                    print "%s" % comment
                    print "%s = %s" % (name,value)
                    print "\nPlease verify that parameter names are unique."
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
                    print "  done section for parameter: ",name
                    #print "starting new section      : ",nsections
                    #print ""

        else:
            # empty line
            comment += '\n'

    f.close()

    # checks number of entries
    if len(parameters.keys()) != nsections:
        print "Error number of sections ",nsections," but number of keys is ",len(parameters.keys())
        sys.tracebacklimit=0
        raise Exception('number of sections read in invalid: %s' % file)

    # user info
    if verbose:
        print ""
        print "  got %d parameters" % nsections
        print ""
        #print "master file sections"
        #print ""
        #print parameters



def get_maximum_parameter_name_length(parameters,verbose=False):
    """
    determines maximum parameter name length for formatting
    
    example:
        max_name_length = get_maximum_parameter_name_length(parameters)
    """
    # determines maximum parameter name length
    max_name_length = 0
    for name in parameters.keys():
        #print "parameter name: ",name
        if len(name) > max_name_length: max_name_length = len(name)

    # user info
    if verbose:
        print "parameter names:"
        print "  maximum name length: ",max_name_length

    # restrict to a minimum of 32 character until = sign, we will add 1 space when writing out line below
    if max_name_length < 32-1:
      max_name_length = 32-1
    else:
      # add additional white space
      max_name_length += 1

    # user info
    if verbose:
        print"  using parameter format length: ",max_name_length," + 1"
        print ""

    return max_name_length



def write_template_file(parameters,tmp_file,verbose=False):
    """
    creates a new template file with parameter sections
    """
    # checks if anything to do
    if len(parameters.keys()) == 0:
        print "Error no parameters to write out into template file:",parameters
        sys.tracebacklimit=0
        raise Exception('invalid parameters dictionary: %s' % tmp_file)

    # determines maximum parameter name length
    max_name_length = get_maximum_parameter_name_length(parameters,verbose=verbose)

    # opens file
    try:
        f = open(tmp_file,'w')
    except:
        print "Error opening file ",tmp_file
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

        # parameter line
        if appendix:
            f.write( "%s = %s   %s\n" % (name.ljust(max_name_length),value,appendix) )
        else:
            f.write( "%s = %s\n" % (name.ljust(max_name_length),value) )
    f.write( "\n" )
    f.close()

    # user info
    if verbose:
        print "written temporary file, see: ",tmp_file
        print ""

def get_files_in_subdirectories(dir,files,basename,exclude_dir,exclude_name):
    """
    recursive function to list all files in subdirectories
    """
    for subdir in os.listdir(dir):
        path = os.path.join(dir, subdir)
        if not os.path.isdir(path):
            filename = os.path.basename(path)
            # adds filename to file list
            if filename.startswith(basename):
                # checks with exclude list of name
                do_add_me = True
                for name in exclude_name:
                    if filename.startswith(name): do_add_me = False
                # adds to list
                if do_add_me:
                    files.append(path)
        else:
          do_search = True
          # checks with list of exclude directories
          for name in exclude_dir:
              if subdir == name: do_search = False
          if do_search:
              get_files_in_subdirectories(path,files,basename,exclude_dir,exclude_name)



def compare_and_replace_file(master_file,temp_file,verbose=False,replace=False):
    """
    notifies user if template is different than master
    e.g. by different indentation or white space)
    """
    # diff files
    command = "diff " + master_file + " " + temp_file
    #print command
    ret = os.system(command)
    if ret > 0:
        # got some differences
        # replaces master with new file format
        print ""
        print "replacing master with new parameter file format:"
        print ""

        # copy over
        if replace:
            command = "cp -v " + temp_file + " " + master_file
            #command = "echo hello"
            #print command
            ret = os.system(command)
            if ret != 0:
                print "Error replacing file with new format",master_file
                sys.tracebacklimit=0
                raise Exception('file can not be updated: %s' % master_file)
    else:
        if verbose:
            print "  no differences"
            print "  master file format is okay"
            print ""


def update_Par_files(master_file,replace=False):
    """
    uses a master to update other parameter files  
    """
    global master_parameters
    global deprecated_renamed_parameters

    # user info
    print "master file: ",master_file
    print ""

    # reads in parameters
    read_Par_file_sections(master_parameters,master_file,verbose=True)

    # opens temporary file with master info
    tmp_file = "_____temp01_____"
    write_template_file(master_parameters,tmp_file,verbose=True)

    # notifies user if template is different than master
    # (e.g. by different indentation or white space)
    print "checking differences between new format and master:"
    print "  (different formatting or whitespace can lead to differences)"
    print ""
    compare_and_replace_file(master_file,tmp_file,verbose=True,replace=replace)

    # clean up temporary file
    command = "rm -f " + tmp_file
    os.system(command)

    # finds all Par_files
    basename = os.path.basename(master_file)
    current_dir = os.getcwd()

    # user info
    print ""
    print "finding all files with name: ",basename
    print "in current directory: ",current_dir

    # exclude other possible files with similar name, but with different format
    if basename == "Par_file":
        exclude_name = [ "Par_file_faults" ]
    else:
        exclude_name = []

    # exclude old unused directories from search directories
    exclude_dir = [ "unused_routines", "small_SEM_solvers_in_Fortran_and_C_without_MPI_to_learn" ]

    # gets all files in subdirectories
    files = []
    get_files_in_subdirectories("./",files,basename,exclude_dir,exclude_name)

    nfiles = len(files)
    print ""
    print "found ",nfiles," parameter files"
    print ""

    ifile = 0
    for file in files:
        ifile += 1
        # user info
        print ""
        print "file ",ifile," out of ",nfiles
        print "processing file: ",file

        # read in parameters
        my_parameters = collections.OrderedDict()
        read_Par_file_sections(my_parameters,file)

        # checks for old, deprecated parameters
        nold_parameters = 0
        #print "  searching deprecated parameters..."
        for name in my_parameters.keys():
            if not name in master_parameters.keys():
                print "  deprecated parameter: ",name
                nold_parameters += 1

                # converts old to new parameter name
                is_found = False
                for old_name,new_name in deprecated_renamed_parameters:
                    if name == old_name :
                        print "    will be converted to ",new_name
                        (val_orig,comment_orig,appendix_orig) = my_parameters[old_name]
                        # stores as new section
                        my_parameters[new_name] = (val_orig,comment_orig,appendix_orig)
                        is_found = True
                        # removes old section
                        del my_parameters[old_name]

                # removes old parameter
                if not is_found:
                    # removes old section
                    del my_parameters[name]


        # add missing master parameters and replaces comment lines to compare sections
        nmissing_parameters = 0
        #print "  searching missing parameters..."
        for name in master_parameters.keys():
            if not name in my_parameters.keys():
                print "  misses parameter: ",name
                nmissing_parameters += 1
                # adds from master template record
                (val,comment,appendix) = master_parameters[name]
                my_parameters[name] = (val,comment,appendix)

        # updates comments
        nold_comments = 0
        for name in master_parameters.keys():
            # checks we have this parameter
            if not name in my_parameters.keys():
                print "Error comparing master with current file format parameter",name
                sys.tracebacklimit=0
                raise Exception('parameter list invalid: %s' % file)

            # compares and replaces comments and appendix
            (val_orig,comment_orig,appendix_orig) = my_parameters[name]
            (val,comment,appendix) = master_parameters[name]
            if comment_orig != comment or appendix != appendix_orig:
                nold_comments += 1
                # replace with new comment/appendix and only keep original value
                my_parameters[name] = (val_orig,comment,appendix)

        # replace old file if necessary
        if nold_parameters == 0 and nmissing_parameters == 0 and nold_comments == 0:
            # user info
            print "  parameter file is okay and up-to-date"
        else:
            # user info
            print "  updating parameter file..."

            # opens temporary file with master info
            tmp_file = "_____temp09_____"
            write_template_file(my_parameters,tmp_file)

            # notifies user if template is different than master
            # (e.g. by different indentation or white space)
            compare_and_replace_file(file,tmp_file,verbose=True,replace=replace)

            # clean up temporary file
            command = "rm -f " + tmp_file
            os.system(command)

        # clean up
        del my_parameters


    print ""
    print "done"
    os.system("date")
    print ""


def usage():
    print "usage:"
    print "    ./process_DATA_Par_files_to_update_their_parameters_from_a_master.py Master-Par-file replace"
    print "  with"
    print "    Master-Par_file - Par_file which serves as master template (e.g. DATA/Par_file)"
    print "    replace         = flag to force replacing of file [0==check-only/1==replace]"


if __name__ == '__main__':
    # gets arguments
    if len(sys.argv) != 3:
        usage()
        sys.exit(1)
    else:
        master_file = sys.argv[1]
        if int(sys.argv[2]) == 1:
            replace = True
        else:
            replace = False

    update_Par_files(master_file,replace=replace)

