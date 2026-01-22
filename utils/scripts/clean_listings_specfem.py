#!/usr/bin/env python
#
# fortran source code cleaning script
#
from __future__ import print_function

import sys
import os
import filecmp
import shutil
import re
import subprocess
import glob
import fnmatch

############################################################################################

# shows a diff between original and new content (if any)
show_diff = True

# replaces original with new content (if anything changed)
replace_file_content = False

############################################################################################

# Define the list of file extensions to be processed for Fortran formatting
fortran_file_extensions = ['.fh', '.f90', '.F90', '.fh.in']

# Define the list of file extensions to be processed as general formatting
general_file_extensions = ['.bash','.c','.cpp','.csh','.cu','.h','.h.in','.pl','.py','.tex','.txt','.sh','.rb','.md','.log']

# Define the list of directories to be excluded (these are mostly submodules included in the source repositories)
exclude_dirs = ['.git', 'm4', './utils/ADJOINT_TOMOGRAPHY_TOOLS/flexwin', './src/inverse_problem_for_source/pyCMT3D']

# Define the list of files to be excluded
exclude_files = ['*Par_file*']

# list of regex patterns to be replaced
patterns = [
    # suppress trailing white spaces and carriage return
    (r'\s*$', ''),
    # use new syntax of comparison operators, ignoring case in starting pattern (useful in case of mixed case)
    (r'\.le\.', '<='),
    (r'\.ge\.', '>='),
    (r'\.lt\.', '<'),
    (r'\.gt\.', '>'),
    (r'\.ne\.', '/='),
    (r'\.eq\.', '=='),
    # switch to lowercase for comparison operators
    (r'\.and\.', '.and.'),
    (r'\.or\.', '.or.'),
    (r'\.not\.', '.not.'),
    (r'\.eqv\.', '.eqv.'),
    (r'\.neqv\.', '.neqv.'),
    (r'\.true\.', '.true.'),
    (r'\.false\.', '.false.'),
    # switch to Fortran2008 standard
    (r'call\s*getarg\(', 'call get_command_argument('),
    # constant strings
    (r'endsubroutine', 'end subroutine'),
    (r'if\s*\(', 'if ('),
    (r'\)\s*then', ') then'),
    (r'end\s*if', 'endif'),
    (r'end\s*do', 'enddo'),
    (r'else\s*if', 'else if'),
    # force lowercase keywords
    (r'subroutine', 'subroutine'),
    (r'end\s*subroutine', 'end subroutine'),
    (r'function', 'function'),
    (r'end\s*function', 'end function'),
    (r'continue', 'continue'),
    (r'implicit none', 'implicit none'),
    (r'implicit', 'implicit'),
    (r'return', 'return'),
    (r' go\s*to ', ' goto '),
    (r'use\s*::\s*mpi', 'use mpi'),
    (r',\s*only\s*:\s*', ', only: '),
    (r'NOISE_SOURCE_TIME_FUNCTION_TYPE', 'noise_source_time_function_type'),
    # do not move this before the above line in which we change the keyword "function"
    (r'use_ricker_time_function', 'USE_RICKER_TIME_FUNCTION'),
    (r'print_source_time_function', 'PRINT_SOURCE_TIME_FUNCTION'),
    (r'sourceTimeFunction', 'sourceTimeFunction'),
    #(r'external_source_time_function', 'EXTERNAL_SOURCE_TIME_FUNCTION'),
    #(r'external_stf', 'EXTERNAL_SOURCE_TIME_FUNCTION'),
    #(r'EXTERNAL_SOURCE_TIME_FUNCTION_filename', 'external_source_time_function_filename'),
    #(r'read_EXTERNAL_SOURCE_TIME_FUNCTION', 'read_external_source_time_function'),
    (r'USE_MAP_function', 'USE_MAP_FUNCTION'),
    (r'enddo_LOOP_IJK', 'ENDDO_LOOP_IJK'),
    (r'enddo_LOOP_IJ', 'ENDDO_LOOP_IJ'),
    (r'OMP do', 'OMP DO'),
    (r'OMP enddo', 'OMP ENDDO'),
    (r'print\*', 'print *'),
    (r'print\s*\*', 'print *'),
    (r'spectral-elements', 'spectral elements'),
    (r'gaussian', 'Gaussian'),
    (r'hessian', 'Hessian'),
    (r'cartesian', 'Cartesian'),
    # suppress space between parenthesis and .not. (this can happen when testing logical operators)
    (r'\( \.not\. ', '(.not. '),
    (r'\)call', ') call'),
    # enforce upper case
    (r'CUSTOM_REAL', 'CUSTOM_REAL'),
    # do not use null strings, which are not part of the Fortran standard (and the IBM xlf compiler rejects them for instance)
    (r'print\s*\*\s*,\s*\'\'', 'print *'),
    (r'write\s*\(\s*\*\s*,\s*\*\s*\)\s*\'\'', 'print *'),
    (r'write\s*\(\s*IMAIN\s*,\s*\*\s*\)\s*\'\'', 'write(IMAIN,*)'),
    (r'write\s*\(\s*IOUT\s*,\s*\*\s*\)\s*\'\'', 'write(IOUT,*)'),
    (r'print\s*\*\s*,\s*""', 'print *'),
    (r'write\s*\(\s*\*\s*,\s*\*\s*\)\s*""', 'print *'),
    (r'write\s*\(\s*IMAIN\s*,\s*\*\s*\)\s*""', 'write(IMAIN,*)'),
    (r'write\s*\(\s*IOUT\s*,\s*\*\s*\)\s*""', 'write(IOUT,*)'),
    # unit 6 means standard output, replace it with standard output symbol
    (r'write\s*\(\s*6\s*,\s*\*\s*\)', 'write(*,*)'),
    (r'write\s*\(\s*6\s*,', 'write(*,'),
    # force space in , & at end of line
    (r'\s*\,\s*&\s*$', ', &'),
    # always use upper case for GLL when used as a word
    (r' gll ', ' GLL '),
    (r' mpi ', ' MPI '),
    (r' pml ', ' PML '),
    # fix some typos I have found in the different codes, or non-US spelling.
    # also switch to US spelling in order to have the same standard in all files.
    (r'regularisation', 'regularization'),
    (r'optimisation', 'optimization'),
    (r'analitical', 'analytical'),
    #    (r'communIcation', 'communication'),
    (r' in orfer ', ' in order '),
    (r' stepest ', ' steepest '),
    (r' stepest$', ' steepest'),
    (r'aloow', 'allow'),
    (r'neighbour', 'neighbor'),
    (r'vecotr', 'vector'),
    (r'computse', 'compute'),
    (r'indicies', 'indices'),
    (r'accordig', 'according'),
    (r'paralell', 'parallel'),
    (r'debbug', 'debug'),
    # do not suppress the white space here because it would then change "debugging" for instance
    (r'debugg ', 'debug '),
    (r'debugg$', 'debug'),
    (r'familly', 'family'),
    (r'warnning', 'warning'),
    (r'elemement', 'element'),
    (r'cartesion', 'Cartesian'),
    (r'partiton', 'partition'),
    (r'drection', 'direction'),
    (r'seperation', 'separation'),
    (r'inverision', 'inversion'),
    (r'restauration', 'restoration'),
    (r'restaure', 'restore'),
    (r'memmory', 'memory'),
    (r'convolution formation', 'convolution formulation'),
    (r' fortran', ' Fortran'),
    (r'adress', 'address'),
    (r'gFortran', 'gfortran'),
    (r' usefull ', ' useful '),
    (r' usefull$', ' useful'),
    # enforce upper case
    (r'MAX_neighborS', 'MAX_NEIGHBORS'),
]

# list of regex patterns to be replaced only for selected files w/out excluded files
special_patterns = [
    # operators
    (r'\s*<\s*=\s*', ' <= '),
    (r'\s*>\s*=\s*', ' >= '),
    (r'\s*<\s*', ' < '),
    (r'\s*/=\s*', ' /= '),
    # restore operators that may have been split by the above introduction of white spaces
    (r'<\s*=', '<='),
    (r'>\s*=', '>='),
    (r'=\s*=', '=='),
    (r'/\s*=', '/='),
    # also restore bash file pipes that may appear in some print statements that save bash scripts to disk for future processing
    (r'>\s*&', '>&'),
    (r'<\s*&', '<&'),
    # also restore xml-formatting strings '< and >'
    (r'\'\s*<\s*', '\'<'),
    (r'\s*>\s*\'', '>\''),
    # for pointers
    (r'\s*=\s*>\s*(?!$)', ' => '),
]

# patterns for comment/non-comment lines
comment = [ '!' ]
comment_patterns = [
    (r'-\s*>', '->'),
    (r'<\s*-', '<-'),
]
non_comment_patterns = [
    (r'(?<!\')(?<!=)\s*>(?!=)(?!\')(?!&)\s*', ' > '),
    (r'\s*==(?!=)\s*', ' == '),
    (r'(?<!\s)=\.true\.', ' = .true.'),
    (r'(?<!\s)=\.false\.', ' = .false.'),
]

def format_content_fortran(content):
    """
    applies Fortran code formatting
    """
    # line-by-line
    lines = content.split('\n')  # Split the content into lines

    for i, line in enumerate(lines):
        # first letter on line
        line_nospace = line
        line_nospace = line_nospace.replace(" ", "")
        if len(line_nospace) > 0:
            first_letter = line_nospace[0]
        else:
            first_letter = ''
        #print(f"line {i}: first_letter={first_letter} line: {line}")

        ## general patterns
        for pattern, replacement in patterns:
            #print(f"pattern: {pattern}")
            line = re.sub(pattern, replacement, line, flags=re.IGNORECASE)

        ## special patterns formatting (operators,..)
        # check if line has xml format (contains </ or '< or >' patterns)
        xml_patterns = [ r'</', r'\'<', r'>\'' ]
        has_xml_pattern = False
        for pattern in xml_patterns:
            if re.search(pattern, line):
                has_xml_pattern = True
                break
        if not has_xml_pattern:
            ## Replace special patterns
            for pattern, replacement in special_patterns:
                line = re.sub(pattern, replacement, line, flags=re.IGNORECASE)

            ## Replace patterns on non-comment lines
            if first_letter in comment:
                # comment line
                for pattern, replacement in comment_patterns:
                    #print(f"comment pattern: {pattern}")
                    line = re.sub(pattern, replacement, line, flags=re.IGNORECASE)
            else:
                # non-comment line
                for pattern, replacement in non_comment_patterns:
                    #print(f"non-comment pattern: {pattern}")
                    line = re.sub(pattern, replacement, line, flags=re.IGNORECASE)

        ## special formatting
        # "write(IMAIN,*)'my-comment'" -> "write(IMAIN,*) 'my-comment'"
        if re.search(r'\bwrite\s*\(IMAIN,\*\)\'[^\']*\'', line):
            line = re.sub(r'\)(?=\'[^\']*\')', ') ', line, flags=re.IGNORECASE)

        # "write(IMAIN,*)my-parameter" -> "write(*,*) my-parameter"
        if re.search(r'\bwrite\s*\(IMAIN,\*\)\w+', line):
            line = re.sub(r'\)(?=\w)', ') ', line, flags=re.IGNORECASE)

        # "write(*,*)'my-comment'" -> "write(*,*)'my-comment'"
        if re.search(r'\bwrite\s*\(\*,\*\)\'[^\']*\'', line):
            line = re.sub(r'\)(?=\'[^\']*\')', ') ', line, flags=re.IGNORECASE)

        # "write(*,*)my-parameter" -> "write(*,*) my-parameter"
        if re.search(r'\bwrite\s*\(\*,\*\)\w+', line):
            line = re.sub(r'\)(?=\w)', ') ', line, flags=re.IGNORECASE)

        # "if (a==b)something" -> "if (a==b) something"
        if re.search(r'\bif\s*\(\s*(\w+)\s*==\s*(\w+)\s*\)\w+', line):
            line = re.sub(r'\)(?=\w)', ') ', line, flags=re.IGNORECASE)

        # on non-comment lines
        if not first_letter in comment:
            # do i=1,..   -> do i = 1,..
            if re.search(r'\bdo\s+(\w+)\s*=(\d+)\s*\,', line):
                line = re.sub(r'\bdo\s+(\w+)\s*=(\d+)\s*\,', r'do \1 = \2,', line, flags=re.IGNORECASE)

            # do i=ilat,..   -> do i = ilat,..
            if re.search(r'\bdo\s+(\w+)\s*=(\w+)\s*\,', line):
                line = re.sub(r'\bdo\s+(\w+)\s*=(\w+)\s*\,', r'do \1 = \2,', line, flags=re.IGNORECASE)

            # "myvar==0" -> "myvar == 0"
            if re.search(r'(\w+)==(\d+)', line):
                newline = re.sub(r'(\w+)==(\d+)', r'\1 == \2', line, flags=re.IGNORECASE)
                print("  A newline: ",newline)

            # "myvar==something" -> "myvar == something"
            if re.search(r'\b(\w+)==(\w+)', line):
                newline = re.sub(r'\b(\w+)==(\w+)', r'\1 == \2', line, flags=re.IGNORECASE)
                print("  B newline: ",newline)

            # "a=b" -> "a = b"
            exclude_equal_patterns = [
                r'==',
                r'>=',
                r'<=',
                r'\bopen\s*\(',
                r'\bclose\s*\(',
                r'\binquire\s*\(',
                r'\bread\s*\(',
                r'\bwrite\s*\(',
                r'\brandom_seed\s*\(',
                r'\bminloc\s*\(',
                r'\bminval\s*\(',
                r'\bmaxval\s*\(',
                r'\bcheck_status\s*\(',
                r'\bget_command_argument\s*\(',
                r'\bdate_and_time\s*\(',
                r'\bexit_mpi\s*\(',
                r'\bexit_MPI\s*\(',
                r'\blibxsmm',
                r'\bprint\s*\*',
                r'\brecl=',
                r'\bstat=',
                r'\bexitstat=',
                r'\biostat=',
                r'\blen=',
                r'\bkind=',
                r'\bh5',
                r'^[^)]*\)[^)]*$',
                r'^\s*&',
            ]
            has_equal_pattern = False
            for pattern in exclude_equal_patterns:
                if re.search(pattern, line):
                    has_equal_pattern = True
                    break
            if not has_equal_pattern:
                # "myvar=something" -> "myvar = something"  but not "a==b" or lines with "open(unit=.." etc.
                if re.search(r'(\w+)=(\w+)', line):
                    line = re.sub(r'(\w+)=(\w+)', r'\1 = \2', line, flags=re.IGNORECASE)

        # Replace the original line with the modified line
        lines[i] = line

    # Join the modified lines back together
    content_new = '\n'.join(lines)
    return content_new


def format_content_general(content):
    """
    applies cleaning to general (text) files, for example output_solver.txt files in REF_SEIS/ folders
    """
    # line-by-line
    lines = content.split('\n')  # Split the content into lines

    for i, line in enumerate(lines):
        ## general formatting
        # suppress trailing white spaces and carriage return
        line = re.sub(r'\s*$', '', line)
        # Replace the original line with the modified line
        lines[i] = line

    # Join the modified lines back together
    content_new = '\n'.join(lines)
    return content_new


def clean_code_format(file):
    """
    cleans code format
    """
    # Exclude specified files
    if any(fnmatch.fnmatch(file, filename) for filename in exclude_files):
        return

    # Process files only with specified extensions
    # Fortran files
    is_Fortran_file = False
    is_general_file = False
    if any(file.endswith(ext) for ext in fortran_file_extensions):
        is_Fortran_file = True
    elif any(file.endswith(ext) for ext in general_file_extensions):
        is_general_file = True

    # checks if anything to do
    if not is_Fortran_file and not is_general_file:
        return

    print(f'Processing {file}...')

    # Read the file
    with open(file, 'r') as f:
        content = f.read()

    # content
    if is_Fortran_file:
        # fortran code formatting
        content_new = format_content_fortran(content)
    else:
        # general file cleaning
        content_new = format_content_general(content)

    # output
    if show_diff:
        # show all content
        #print("content:")
        #print(content_new)

        # show differences only line-by-line
        # line-by-line
        lines_org = content.split('\n')  # Split the content into lines
        lines_new = content_new.split('\n')  # Split the content into lines

        len_org = len(lines_org)
        len_new = len(lines_new)
        if len(lines_org) != len(lines_new):
            print("Warning: content number of lines differ: original = {} new = {}".format(len_org,len_new))

        length = min(len_org,len_new)
        for i in range(length):
            line_org = lines_org[i]
            line_new = lines_new[i]
            # show if lines are different
            if line_new != line_org:
                print(f"  line {i}: - {line_org}")
                print(f"  line {i}: + {line_new}")

    if replace_file_content:
        # Write the modified content back to the file if anything changed
        if content_new != content:
            # new content is different
            with open(file, 'w') as f:
                f.write(content_new)

def clean_listings(folder_filename):
    """
    loops over (Fortran) code files and updates code formatting
    """
    # determines whether a folder or a specific file was provided as input
    if os.path.isdir(folder_filename):
        # folder
        # Define the path to the source code
        src_path = folder_filename

        # Iterate over all files in the source directory and its subdirectories
        for root, dirs, files in os.walk(src_path):
            # Exclude specified directories
            dirs[:] = [d for d in dirs if d not in exclude_dirs and os.path.join(root,d) not in exclude_dirs]

            for file in files:
                file_path = os.path.join(root, file)
                # clean code formatting
                clean_code_format(file_path)

    elif os.path.isfile(folder_filename):
        # file
        file_path = folder_filename
        # clean code formatting
        clean_code_format(file_path)

    print("")
    print("all done")
    print("")

# reads in arguments
def usage():
    print("Usage: ./clean_listings_specfem.py filename/folder [--diff/--no-diff] [--replace]")
    print("")
    print("  filename/folder    - required input file or folder containing Fortran source code files, e.g., src/")
    print("  --diff/--no-diff   - show or don't show formatting differences (default is to show differences)")
    print("  --replace          - replace file content with new formatting (default off)")
    sys.exit(1)


if __name__ == '__main__':

    # gets arguments
    if len(sys.argv) < 2:
        usage()

    folder_filename = sys.argv[1] # file or folder with source code files

    # reads arguments
    i = 0
    for arg in sys.argv:
        i += 1
        #print("arg: ",arg)
        # get arguments
        if "--diff" in arg:
            show_diff = True
        elif "--no-diff" in arg:
            show_diff = False
        elif "--replace" in arg:
            replace_file_content = True
        elif i > 2:
            print("argument not recognized: ",arg)
            print("")
            usage()

    # main routine
    clean_listings(folder_filename)
