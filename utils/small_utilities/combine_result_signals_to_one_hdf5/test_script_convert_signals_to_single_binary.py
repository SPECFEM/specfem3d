from __future__ import print_function

import subprocess

def main():
    """
    """

    file_of_signal_list = "list_signal_files.txt" # path to a file which is a list of pathes to .semp files
    output_file_name    = "test_out.h5"           # name of output binary file
    total_timesteps     = 5000                    # total number of timesteps in your calculation
    total_num_sigs      = 5                       # total number of signals (.semp files) to be combined and binarized

    # time window for limited extraction of timesteps
    time_window = [0,total_timesteps-1]

    # binarize_signals
    cmd = './combine ' + file_of_signal_list + ' ' + output_file_name + ' ' + str(total_num_sigs) + \
                   ' ' + str(time_window[0]) + ' ' +  str(time_window[1]) + ' ' + str(total_timesteps)
    subprocess.call(cmd, shell=True)


if __name__ == '__main__':
    main()
