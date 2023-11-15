#!/usr/bin/env python
#
# class with helper functions
#
#########################################
from __future__ import print_function

import sys
import array
import os.path

import numpy as np
import matplotlib.pyplot as plt


class helper(object):
    """
    class for helper functions
    """
    def __init__(self,verbose=False):
        ## initializations

        # sets default custom_real type
        # defines float 'f' or double precision 'd' for binary values
        self.custom_type = 'f'

        # SU traces
        self.sampling_DT = 0.0 # in microseconds 10-6
        self.source_x = 0      # source position
        self.source_z = 0

        # verbosity
        self.verbose = verbose

        # debugging
        self.debug = False

    #
    #------------------------------------------------------------------------------------------
    #

    def __str__(self):
        info = "helper class"

        return info

    #
    #------------------------------------------------------------------------------------------
    #

    def read_marker(self,file):
        """
        reads array marker from fortran binary file
        (each time a fortran write(*)-routine call is executed, a marker at beginning and end is placed)
        """
        binlength = array.array('i')
        binlength.fromfile(file,1)

        if self.verbose:
            print("marker length = ",binlength[0])

        return binlength[0]

    #
    #------------------------------------------------------------------------------------------
    #

    def read_binary_file_integer(self,file):
        """
        reads integer value from file
        """
        # gets array length in bytes
        binlength = self.read_marker(file)

        # integer (4 bytes)
        num_points = int(binlength / 4)

        if self.verbose:
            print("number of integers: ",num_points)

        int_values = array.array('i')
        int_values.fromfile(file,num_points)

        if self.verbose:
            print("integer values = ",int_values)

        binlength_end = self.read_marker(file)
        if binlength_end != binlength:
            print("Error array markers in fortran binary file:",binlength,binlength_end)
            print("start array length = ",binlength," Bytes")
            print("final array length = ",binlength_end," Bytes")
            print("array lengths should match, please check your file")
            raise Exception('array markers invalid')

        # single value
        ival = int_values[0]
        return ival

    #
    #------------------------------------------------------------------------------------------
    #

    def read_binary_file_integer_array(self,file):
        """
        reads integer data array from file
        """

        # gets array length in bytes
        binlength = self.read_marker(file)

        # integer (4 bytes)
        num_points = int(binlength / 4)

        # user output
        if self.verbose:
            print("  array length = ",binlength," Bytes")
            print("  number of points in array = ",num_points)
            print("")

        # reads in array data
        int_values = array.array('i')
        int_values.fromfile(file,num_points)

        # fortran binary file: file ends with array-size again
        # checks markers
        binlength_end = self.read_marker(file)
        if binlength_end != binlength:
            print("Error array markers in fortran binary file:",binlength,binlength_end)
            print("start array length = ",binlength," Bytes")
            print("final array length = ",binlength_end," Bytes")
            print("array lengths should match, please check your file")
            raise Exception('array markers invalid')

        data = list()
        data.append(int_values)

        # returns data from list-output (first item)
        # converts list to numpy array
        intdata_array = np.array(data[0],dtype='i')

        return intdata_array

    #
    #------------------------------------------------------------------------------------------
    #

    def read_binary_file_custom_real_array(self,file):
        """
        reads data array from file
        """

        # gets array length in bytes
        binlength = self.read_marker(file)

        if self.custom_type == 'f':
            # float (4 bytes) for single precision
            num_points = int(binlength / 4)
        else:
            # double precision
            num_points = int(binlength / 8)

        # user output
        if self.verbose:
            print("  array length = ",binlength," Bytes")
            print("  number of points in array = ",num_points)
            print("")

        # reads in array data
        binvalues = array.array(self.custom_type)
        binvalues.fromfile(file,num_points)

        # fortran binary file: file ends with array-size again
        # checks markers
        binlength_end = self.read_marker(file)
        if binlength_end != binlength:
            print("Error array markers in fortran binary file:",binlength,binlength_end)
            print("start array length = ",binlength," Bytes")
            print("final array length = ",binlength_end," Bytes")
            print("array lengths should match, please check your file")
            raise Exception('array markers invalid')

        data = list()
        data.append(binvalues)

        # returns data from list-output (first item)
        # converts list to numpy array
        data_array = np.array(data[0],dtype=self.custom_type)

        return data_array


    #
    #------------------------------------------------------------------------------------------
    #

    def write_binary_file_custom_real_array(self,filename,data,verbose=False):
        """
        writes data array to file
        """
        with open(filename, mode='wb') as f:
            # gets array length in bytes
            # marker
            binlength = array.array('i')
            num_points = data.size
            if self.custom_type == 'f':
                # float (4 bytes) for single precision
                binlength.fromlist([num_points * 4])
            else:
                # double precision
                binlength.fromlist([num_points * 8])

            # user output
            if verbose:
                print("filename: ",filename)
                print("  array length = ",binlength," Bytes")
                print("  number of points in array = ",num_points)

            # writes array data
            binvalues = array.array(self.custom_type)

            data = np.reshape(data, (num_points), order='F') # fortran-like index ordering
            #print("debug: data ",data.tolist())

            binvalues.fromlist(data.tolist()) #.tolist())

            # fortran binary file: file starts with array-size again
            binlength.tofile(f)
            # data
            binvalues.tofile(f)
            # fortran binary file: file ends with array-size again
            binlength.tofile(f)

            if verbose:
                print("  file written")
                print("")


    #
    #------------------------------------------------------------------------------------------
    #


    def read_binary_NSPEC_ibool_file(self,filename,verbose=False):
        # reads proc000000_NSPEC_ibool.bin,.. arrays

        if verbose: print("binary file: ",filename)

        # checks if file exists
        if not os.path.isfile(filename):
            print("file %s does not exist, please check..." % filename)
            sys.exit(1)

        with open(filename,'rb') as f:
            # fortran file format:
            #   write(888) nspec
            #   write(888) ibool

            # nspec
            nspec = self.read_binary_file_integer(f)

            # ibool
            ibool = self.read_binary_file_integer_array(f)

        return nspec,ibool

    #
    #------------------------------------------------------------------------------------------
    #


    def read_binary_ibool_file(self,filename,verbose=False):
        # reads proc000000_ibool.bin,.. arrays

        if verbose: print("binary file: ",filename)

        # checks if file exists
        if not os.path.isfile(filename):
            print("file %s does not exist, please check..." % filename)
            sys.exit(1)

        with open(filename,'rb') as f:
            # fortran file format:
            #   write(888) ibool

            # ibool
            ibool = self.read_binary_file_integer_array(f)

        return ibool

    #
    #------------------------------------------------------------------------------------------
    #


    def read_binary_SEM_file(self,filename,verbose=False):
        # reads proc000000_x.bin,.. arrays

        if verbose: print("binary file: ",filename)

        # checks if file exists
        if not os.path.isfile(filename):
            print("file %s does not exist, please check..." % filename)
            sys.exit(1)

        with open(filename,'rb') as f:
            # fortran file format: binary or gll
            #   write(172) x_save
            data = self.read_binary_file_custom_real_array(f)

        return data

    #
    #------------------------------------------------------------------------------------------
    #


    def plot_model_image_3D(self,xstore,ystore,zstore,ibool,array,name,plot_kernel=False,verbose=False):
        """
        plots model or kernel image cross-section
        """
        # cross-section in X-Z plane through middle of y-dimension
        # define grid.
        xi = np.linspace(xstore.min(), xstore.max(), 150)
        zi = np.linspace(zstore.min(), zstore.max(), 150)
        ymid = 0.5 * (ystore.max() + ystore.min())

        #print("plot: X-Z plane at Y-midpoint = ",ymid)
        yi = np.ones(xi.size)
        yi *= ymid

        # or example: equation for surface along which to interpolate
        #equation = lambda x,z : 660.0    # (1-x) * z
        #yi = equation(xi, zi)

        # need points location on which array is defined: array(NGLLX*NGLLY*NGLLZ*nspec) vs. xstore(nglob)
        num_points = array.size
        #print("plot: number of data points = ",num_points)

        points = np.empty((num_points,3),dtype='float32')
        for i in range(num_points):
            iglob = ibool[i] - 1 # indexing starts at 0
            x = xstore[iglob]
            y = ystore[iglob]
            z = zstore[iglob]
            points[i] = [x,y,z]

        #print("plot: points min/max = ",points.min(),points.max())
        #print("plot: xi",xi)
        #print("plot: yi",yi)
        #print("plot: zi",zi)

        # grid data
        # Linearly interpolate the data(x,y,z) on a grid defined by (xi, yi, zi).
        if 1 == 1:
            # scipy
            from scipy.interpolate import griddata
            # 'nearest','linear','cubic',..
            data = griddata(points, array, (xi[None,:], yi[:,None], zi[:,None]), method='nearest')
            #data = griddata(points, array, (xi[None,:], yi[:,None], zi[:,None]), method='linear')

            #2D: data = griddata((x,y), z, (xi[None, :], yi[:, None]), method='linear')
        else:
            # 2D only
            import matplotlib.tri as tri
            triang = tri.Triangulation(x, y)
            interpolator = tri.LinearTriInterpolator(triang, z)
            Xi, Yi = np.meshgrid(xi, yi)
            data = interpolator(Xi, Yi)

        #print("plot: data ",data)

        #non-interpolated
        #dimx = int(xstore.max())
        #dimz = int(zstore.max())
        #data = np.zeros((dimx+1,dimz+1))
        #for ispec  in range(0, nspec):
        #    for j in range(0, NGLLZ):
        #        for i in range(0, NGLLX):
        #            # GLL point location (given in m: dimension 2640 m x 2640 x x 1440 m)
        #            iglob = ibool[i,j,ispec]
        #            x = xstore[i,j,ispec]
        #            z = zstore[i,j,ispec]
        #            ix = int(x)
        #            iz = int(z)
        #            data[ix,iz] = vpstore[i,j,ispec]

        #extent = [x.min(),x.max(),y.min(),y.max()]

        # X-Z plane plot
        extent = [xstore.min(),xstore.max(),zstore.min(),zstore.max()]

        # limit size
        total_max = abs(data).max()
        #print("plot: data max = ",total_max)

        plt.clf()
        plt.title(name)

        if plot_kernel:
            if total_max < 1.e-8:
                if total_max != 0.0:
                    total_max = 1.0 * 10**(int(np.log10(total_max))-1)  # example: 2.73e-11 limits to 1.e-11
                else:
                    total_max = 0.0
            if self.debug: print("debug plot: color scale max = ",total_max)

            plt.imshow(data, extent=extent, origin='lower', interpolation='none', cmap="seismic",
                       vmax=total_max, vmin=-total_max)  # colormaps: "RdBu_r", "twilight"
        else:
            plt.imshow(data, extent=extent, origin='lower', interpolation='none', cmap="seismic")  # colormaps: "RdBu_r", "twilight"

        plt.colorbar()  # draw colorbar

        # saves as JPEG file
        filename = name + ".jpg"
        plt.savefig(filename)
        if verbose: print("  plotted as ",filename)

        # show plot
        if 1 == 0: plt.show()


    #
    #------------------------------------------------------------------------------------------
    #

    def read_binary_SU_file(self,file):
        """
        reads SU file data
        """
        # debugging
        DEBUG = False

        # data container
        dat = list()

        # loop over all receiver records
        num_receivers = 0

        #### read NSTEP from seismograms
        ###filename_in=trim(procname)//"_dx_SU"
        ###real(kind=4) :: r4head(60)
        ###!integer(kind=4) :: header4(1)
        ###!integer(kind=2) :: header2(2)
        ###open(111,file="../OUTPUT_FILES/"//trim(filename_in),access='direct',recl=240,iostat = ier)
        ###read(111,rec=1,iostat=ier) r4head
        ###close(111)
        ###header4=r4head(29)
        ###NSTEP=header2(2)
        ###header4=r4head(30)
        ###DT=header2(1)*1.0d-6
        ###print *, 'irec=',r4head(1)
        ###print *, 'xs=',r4head(19)
        ###print *, 'zs=',r4head(20)
        ###print *, 'xr=',r4head(21)
        ###print *, 'zr=',r4head(22)
        ###print *, 'NSTEP=',NSTEP
        ###print *, "DT=",DT
        ###  irec=1
        ###  do while(ios == 0)
        ###    ! reads in data
        ###    read(11,rec=irec,iostat=ios) r4head,dat
        ###    if (ios /= 0) cycle
        ###    ! writes out adjoint source
        ###    write(33,rec=irec,iostat=ios) r4head,adj
        ###    if (ios /= 0) cycle
        ###    irec=irec+1
        ###  enddo

        while True:
            # format:
            # integer, dimension(28) :: header1
            # real(kind=4), dimension(30) :: header4
            # integer(kind=2) :: header2(2),header3(2)
            #
            # ioffset = 4 * ((irec-1) * (NSTEP/subsamp_seismos + 60)) + 1
            # write(12,pos=ioffset) header1,header2,header3,header4

            # read offset 1 - marker
            #binlength = array.array('b')
            #binlength.fromfile(file,1)

            #print("debug: marker binlength ",binlength)

            # 4-bytes for single precision
            #num_points = int(binlength[0] / 4)
            # double precision
            #num_points = int(binlength[0]/8)
            #print("number of data points: ",num_points)

            # header1
            header1_values = array.array('i')   # signed int - 4 bytes
            try:
                header1_values.fromfile(file,28)
            except EOFError as e:
                # no more data
                #print("Error reading header1: ",e)
                # exit while-loop
                break

            # increases receiver count
            num_receivers += 1

            if DEBUG: print("debug: header1 ",header1_values)

            #! write SU headers (refer to Seismic Unix for details)
            # header1(1)  =  irec                          ! receiver ID
            # header1(10) = NINT(st_xval(irec)-x_source)  ! offset
            # header1(19) = NINT(x_source)                ! source location xs
            # header1(20) = NINT(z_source)                ! source location zs
            # header1(21) = NINT(st_xval(irec))           ! receiver location xr
            # header1(22) = NINT(st_zval(irec))           ! receiver location zr

            # receiver ID
            irec = header1_values[0]
            ioffset = header1_values[9]
            isource_x = header1_values[18]
            isource_z = header1_values[19]
            irec_x = header1_values[20]
            irec_z = header1_values[21]

            # header2
            header2_values = array.array('h')    # signed short - 2 bytes
            header2_values.fromfile(file,2)

            if DEBUG: print("debug: header2 ",header2_values)

            # header2(1) = 0  ! dummy
            # header2(2) = int(NSTEP/subsamp_seismos, kind=2)
            num_samples = int(header2_values[1])

            if DEBUG: print("  trace number of samples ",num_samples)

            # header3
            header3_values = array.array('h')    # signed short - 2 bytes
            header3_values.fromfile(file,2)

            if DEBUG: print("debug: header3 ",header3_values)

            # header3(1) = NINT(sampling_deltat*1.0d6, kind=2)  ! deltat (unit: 10^{-6} second)
            sampling_rate = header3_values[0]

            # header4
            header4_values = array.array('f')    # real - 4 bytes
            header4_values.fromfile(file,30)

            if DEBUG: print("debug: header4 ",header4_values)

            # stores trace infos
            if num_receivers == 1:
                self.source_x = isource_x
                self.source_z = isource_z
                self.sampling_DT = sampling_rate

            # trace values

            # ioffset = 4 * ((irec-1) * (NSTEP/subsamp_seismos + 60) + 60 + seismo_offset) + 1
            # write(12,pos=ioffset) single_precision_seismo

            # single header read
            # header r4head 60 * 4 bytes
            #header_skip_values = array.array('f')    # real - 4 bytes
            #header_skip_values.fromfile(file,60)
            #print("debug: header_skip ",header_skip_values)

            # read real values
            binvalues = array.array('f')
            # read double precision values
            #binvalues = array.array('d')

            binvalues.fromfile(file,num_samples)
            if DEBUG: print("debug: binvalues ",binvalues)

            # read marker
            #binlength.fromfile(file,1)

            #print("debug: marker binlength ",binlength)

            # store data values
            dat.append(binvalues)
            #debug
            #print("debug: data ",dat)

        # total number of receivers
        if self.verbose:
            print("  number of receivers    : ",num_receivers)
            print("  trace number of samples: ",num_samples)
            print("")

        # returns data from list-output (first item)
        # converts list to numpy array
        data_array = np.array(dat,dtype='f')

        return data_array

    #
    #------------------------------------------------------------------------------------------
    #


    def read_SU_file(self,filename,verbose=False):
        """
        reads SU-format binary file
        """
        # user output
        if verbose: print("SU file: ",filename)

        # checks if file exists
        if not os.path.isfile(filename):
            print("file %s does not exist, please check..." % filename)
            sys.exit(1)

        # open SU file
        with open(filename,'rb') as f:
            # reads SU data
            data_array = self.read_binary_SU_file(f)

        #debug
        #print("debug: num receivers ",len(data_array))
        #print("debug: samples ",len(data_array[0]))
        #print("debug: total size",data_array.size)

        #sx = dat[2][:nstations]
        #sy = dat[2][nstations:2*nstations]
        #sz = dat[2][2*nstations:3*nstations]
        #station_CPU_index = dat[0][:nstations]

        return data_array

    #
    #------------------------------------------------------------------------------------------
    #

    def write_binary_SU_file(self,data_array,file):
        """
        writes out SU-format file
        """
        #
        num_receivers = len(data_array)
        num_samples = len(data_array[0])

        if self.verbose:
            print("  num receivers = ",num_receivers)
            print("  samples       = ",num_samples)

        ###  irec=1
        ###  do while(ios == 0)
        ###    ! reads in data
        ###    read(11,rec=irec,iostat=ios) r4head,dat
        ###    if (ios /= 0) cycle
        ###    ! writes out adjoint source
        ###    write(33,rec=irec,iostat=ios) r4head,adj
        ###    if (ios /= 0) cycle
        ###    irec=irec+1
        ###  enddo

        # loop over all receiver records
        for irec in range(num_receivers):
            # format:
            # integer, dimension(28) :: header1
            # real(kind=4), dimension(30) :: header4
            # integer(kind=2) :: header2(2),header3(2)
            #
            # ioffset = 4 * ((irec-1) * (NSTEP/subsamp_seismos + 60)) + 1
            # write(12,pos=ioffset) header1,header2,header3,header4

            # read offset 1 - marker
            #binlength = array.array('b')
            #binlength.fromfile(file,1)

            #print("debug: marker binlength ",binlength)

            # dummy header record
            # header r4head 60 * 4 bytes
            #header_skip = array.array('f')    # real - 4 bytes
            #header_skip.fromlist([60 * 4])
            #header_skip.tofile(f)

            # header1
            header1_values = array.array('i')   # signed int - 4 bytes
            header1_values.fromlist(np.zeros(28,dtype='i').tolist())
            #! write SU headers (refer to Seismic Unix for details)
            # header1(1)  =  irec                          ! receiver ID
            # header1(10) = NINT(st_xval(irec)-x_source)  ! offset
            # header1(19) = NINT(x_source)                ! source location xs
            # header1(20) = NINT(z_source)                ! source location zs
            # header1(21) = NINT(st_xval(irec))           ! receiver location xr
            # header1(22) = NINT(st_zval(irec))           ! receiver location zr
            # receiver ID
            header1_values[0] = irec
            header1_values[9] = irec # dummy ioffset
            header1_values[18] = self.source_x
            header1_values[19] = self.source_z
            header1_values[20] = irec # dummy irec_x
            header1_values[21] = 0 # dummy irec_z

            # header2
            header2_values = array.array('h')    # signed short - 2 bytes
            header2_values.fromlist(np.zeros(2,dtype='h').tolist())
            # header2(1) = 0  ! dummy
            # header2(2) = int(NSTEP/subsamp_seismos, kind=2)
            header2_values[0] = 0
            header2_values[1] = num_samples

            # header3
            header3_values = array.array('h')    # signed short - 2 bytes
            header3_values.fromlist(np.zeros(2,dtype='h').tolist())
            # header3(1) = NINT(sampling_deltat*1.0d6, kind=2)  ! deltat (unit: 10^{-6} second)
            header3_values[0] = self.sampling_DT

            # header4
            header4_values = array.array('f')    # real - 4 bytes
            header4_values.fromlist(np.zeros(30,dtype='f').tolist())

            # writes out header
            header1_values.tofile(file)
            header2_values.tofile(file)
            header3_values.tofile(file)
            header4_values.tofile(file)

            # writes real values
            binvalues = array.array('f')
            # read double precision values
            #binvalues = array.array('d')
            #data = np.reshape(data_array[irec], order='F') # fortran-like index ordering
            binvalues.fromlist(data_array[irec].tolist()) #.tolist())
            binvalues.tofile(file)


    #
    #------------------------------------------------------------------------------------------
    #


    def write_SU_file(self,data_array,filename,verbose=False):
        """
        writes out file in binary SU file format
        """
        # user output
        if verbose:
            print("output SU file: ",filename)
            print("  total size    = ",data_array.size)

        # open SU file
        with open(filename,'wb') as f:
            # writes out in SU-format
            self.write_binary_SU_file(data_array,f)

        if verbose:
            print("")
            print("  written: ",filename)
            print("")
            print("")

