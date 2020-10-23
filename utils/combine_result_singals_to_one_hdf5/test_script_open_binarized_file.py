from __future__ import print_function

import h5py

def open_bin(fname, id_signal):

    # read h5
    data = h5py.File(fname, 'r')['time_pressure']
    print(data)

    # time
    t = data[0,:,id_signal]

    # amplitude of first signal (the first .semp file put on the signal list file)
    amp = data[1,:,id_signal]

    print("length of time array: ",len(t), "length of amplitude array: ", len(amp))

def main():

    fname = "test_out.h5" # file name of binarized signals
    id_signal = 0 # the number of a signal to be opened

    open_bin(fname, id_signal)

if __name__ == "__main__":
    main()
