#!/usr/bin/env python
#
# based on
# NOISE_TOMOGRAPHY.m
from __future__ import print_function

import os,sys
import numpy as np
from numpy import pi,sin,cos,array,arange,zeros,fft

import matplotlib.pyplot as plt

# noise model defined in PetersonNoiseModel.py
from PetersonNoiseModel import PetersonNoiseModel


def NOISE_TOMOGRAPHY(NSTEP=None,dt=None,Tmin=None,Tmax=None,NOISE_MODEL=None,show_figures=True):
    """
    This is a short program used to sample the spectrum of the
    Peterson's noise model, and convert it to time-domain
    --- the source time function we will use in NOISE TOMOGRAPHY simulations.

    #***********************************************************************
    #*******Please read the manual of SPECFEM3D package for guidance********
    #***********************************************************************

    # Usage:
    NOISE_TOMOGRAPHY(NSTEP,dt,Tmin,Tmax,NOISE_MODEL)

    # Example:
    NOISE_TOMOGRAPHY(2999,0.05,10.,20.,'NLNM')

    # /////////////////////////////////////////////////////////////////////////
    # Parameters:
    # NSTEP       --- number of time steps (always odd for NOISE TOMOGRAPHY)
    # dt          --- time interval of specfem3D solver
    # Tmin, Tmax  --- the period range you are working with (in seconds)
    # NOISE_MODEL --- the Peterson's noise model,
    #                 either 'NLNM' or 'NHNM' (with the quote!)
    #                 'NLNM': New Low  Noise Model (in 1993, the model was New)
    #                 'NHNM': New High Noise Model
    #                 or 'FLAT': for a flat noise spectrum (uniform noise within period range)
    # /////////////////////////////////////////////////////////////////////////

    ATTENTION:
    ***NEVER*** try to calculate "NSTEP" and "dt" by yourself!
    They can be found when you compile the SPECFEM3D package,
    with the correct DATA/Par_file
    If DATA/Par_file is not specified for noise tomography simulations,
    you won't get correct values either!

    It is highly recommended that you compile the package with the correct
    DATA/Par_file you will be using


    # Output from the example above:
    # a figure + following message
    # (the path /data2/yangl/3D_NOSIE/ must be different)
    # *************************************************************
    # the source time function has been saved in:
    # /data2/yangl/3D_NOISE/S_squared
    # S_squared should be put into directory:
    # ./NOISE_TOMOGRAPHY/ in the SPECFEM3D package
    """

    ############################################################
    ## USER PARAMETERS
    ## taper option
    taper_type = 1
    taper_length_percentage = 0.1
    taper_length_min = 40
    ############################################################

    ##
    if show_figures:
        fontsize=10

    ## derived parameters
    T = (NSTEP - 1) * dt
    N_mid = int((NSTEP + 1) / 2)
    fmax = 1.0 / 2.0 / dt
    df = 1.0 / T

    # f should have length == NSTEP
    # using +df/2 to make sure that rounding doesn't lead to length problems
    f = np.concatenate([arange(0.0,fmax+df/2,df),arange(-fmax,-df+df/2,df)])

    #debug
    #print("degug: NSTEP,T,fmax,df ",NSTEP,T,fmax,df)
    #print("debug: range ",len(arange(0.0,fmax+df,df)),len(arange(-fmax,-df+df/2,df)))
    #print("debug: f",f,len(f))

    ## checks length
    if T < Tmax:
        print('The simulation length T('+str(T)+') is smaller than the required maximum period Tmax('+str(Tmax)+')')
        return

    ## checks noise model string
    use_flat_noise_spectrum = False
    if NOISE_MODEL == 'NLNM':
        model_info = 'Peterson New Low Noise Model'
    elif NOISE_MODEL == 'NHNM':
        model_info = 'Peterson New High Noise Model'
    elif NOISE_MODEL == 'FLAT':
        # use a flat noise spectrum instead of Peterson model
        use_flat_noise_spectrum=copy(true)
        model_info = 'uniform noise spectrum'
    else:
        print('Error: noise model %s not recognized, use NLNM, NHNM or FLAT for low, high or flat noise model' %NOISE_MODEL)
        sys.exit(1)

    ## user output
    print('NOISE_TOMOGRAPHY input:')
    print('  number of time steps = %i' % NSTEP)
    print('  time step size       = %f s'% dt)
    print('  period range min/max = %f / %f s' % (Tmin,Tmax))
    print('  noise model          = %s - %s' % (NOISE_MODEL,model_info))
    print('')
    print('  total simulation time = %f s' % T)
    print('  Nyquist frequency: %f Hz' % fmax)

    ## initialize the power spectrum of noise
    accel = zeros(len(f))
    veloc = zeros(len(f))
    displ = zeros(len(f))
    filter_array = zeros(len(f))

    # defined by Tmin & Tmax

    ## calculate the power spectrum of noise from Peterson's model (1993)
    print('  calculating noise power spectrum: number of positive frequencies = %i' % N_mid)
    # only calculate for positive frequencies
    # the negative frequencies will be updated later using symmetry
    for l in arange(0,N_mid):
        if use_flat_noise_spectrum:
            # flat noise spectrum
            accel[l] = 1.0
            veloc[l] = 1.0
            displ[l] = 1.0
        else:
            # Peterson noise model
            if f[l] != 0.0:
                period = 1.0 / f[l]
            else:
                period = 1.e6

            accel[l],veloc[l],displ[l] = PetersonNoiseModel(period,NOISE_MODEL)

    # plots figure
    if show_figures:
        print('  plotting figure: Power spectrum')
        fig = plt.figure(1,figsize=(12, 8),dpi=100)
        plt.subplot(2,1,1)
        plt.semilogx(f[arange(0,N_mid)],accel[arange(0,N_mid)],'b')
        #plt.hold('on')
        plt.semilogx(f[arange(0,N_mid)],veloc[arange(0,N_mid)],'g')
        #plt.hold('on')
        plt.semilogx(f[arange(0,N_mid)],displ[arange(0,N_mid)],'r')
        #plt.hold('on')
        plt.legend(['acceleration','velocity','displacement'])
        plt.xlabel('Frequency (Hz)',size=fontsize)
        plt.ylabel('Amplitude (dB)',size=fontsize)
        plt.title('Power Spectrum of Peterson\'s Noise Model in dB',size=fontsize)
        plt.xlim(0.0001,10.0)
        plt.ylim(-250,- 50)

    ## change power spectrum from dB to real physical unit
    for l in arange(0,N_mid):
        accel[l] = 10.0**(accel[l] / 10.0)
        veloc[l] = 10.0**(veloc[l] / 10.0)
        displ[l] = 10.0**(displ[l] / 10.0)

    #debug
    #n = len(displ); tmp = zeros((n,2));
    #tmp[:,0] = arange(1,n+1); tmp[:,1] = displ[:];
    #np.savetxt('tmp_displ.py.txt',tmp)

    ## constrain the power spectrum only within the range [Tmin Tmax]
    print('  filtering power spectrum:\n    period     min/max = %f / %f\n    frequency  min/max = %f / %f' % (Tmin,Tmax,1.0 / Tmax,1.0 / Tmin))
    for l in arange(0,N_mid):
        if abs(f[l]) >= 1.0/Tmax and abs(f[l]) <= 1.0/Tmin:
            if use_flat_noise_spectrum:
                filter_array[l] = 1.0
            else:
                filter_array[l] = sin((f[l] - 1.0/Tmax) / (1.0/Tmin - 1.0/Tmax) * pi)
        else:
            if abs(f[l]) >= 1.0/1.5 / Tmax and abs(f[l]) < 1.0/Tmax:
                filter_array[l] = sin((f[l] - 1.0/1.5/Tmax) / (1.0/Tmax - 1.0/1.5/Tmax) * pi / 2.0)
            else:
                if abs(f[l]) <= 1.5 / Tmin and abs(f[l]) > 1.0/Tmin:
                    filter_array[l] = sin( -(f[l] - 1.0/Tmin) / (1.5/Tmin - 1.0/Tmin) * pi / 2.0) + 1.0

    #debug
    #n = len(filter_array); tmp = zeros((n,2));
    #tmp[:,0] = arange(1,n+1); tmp[:,1] = filter_array[:];
    #np.savetxt('tmp_filter.py.txt',tmp)

    accel = accel * filter_array
    veloc = veloc * filter_array
    displ = displ * filter_array

    # plots figure
    if show_figures:
        print('  plotting figure: Filtered Power spectrum')
        plt.subplot(2,1,2)
        plt.plot(f[arange(0,N_mid)],accel[arange(0,N_mid)] / max(abs(accel)),'b')
        #hold('on')
        plt.plot(f[arange(0,N_mid)],veloc[arange(0,N_mid)] / max(abs(veloc)),'g')
        #hold('on')
        plt.plot(f[arange(0,N_mid)],displ[arange(0,N_mid)] / max(abs(displ)),'r')
        #hold('on')
        plt.xlabel('Frequency (Hz)',size=fontsize)
        plt.ylabel('Amplitude',size=fontsize)
        plt.title('Power Spectrum filtered between ['+str(Tmin)+' '+str(Tmax)+'] s',size=fontsize)
        plt.legend(['accleration, scaled by '+str(max(abs(accel)))+' m^2/s^4/Hz', \
                    'velocity, scaled by '+str(max(abs(veloc)))+' m^2/s^2/Hz', \
                    'displacement, scaled by '+str(max(abs(displ)))+' m^2/Hz'])
        plt.xlim(0.8 / Tmax, 1.2 / Tmin)
        plt.ylim(-0.1,1.5)

    # using only displacement from here on
    accel = None
    veloc = None

    ## update power spectrum in the negative frequencies, using symmetry
    # note the power spectrum is always REAL, instead of COMPLEX
    for l in arange(N_mid,NSTEP):
        #accel(l)=conj(accel(NSTEP-l+2));
        #veloc(l)=conj(veloc(NSTEP-l+2));
        displ[l] = np.conj(displ[NSTEP - l])

    #debug
    #n = len(displ); tmp = zeros((n,2));
    #tmp[:,0] = arange(1,n+1); tmp[:,1] = displ[:];
    #np.savetxt('tmp_displ2.py.txt',tmp)

    ## prepare source time function for ensemble forward source -- S_squared
    print('  preparing source time function S_squared:\n    NSTEP = %i / dt = %f' % (NSTEP,dt))
    # the file S_squared should be put into directory ./NOISE_TOMOGRAPHY/
    # together with other two files: irec_main_noise & nu_main
    S_squared = zeros((NSTEP,2))

    # second column: source time function
    S_squared[:,0] = (arange(0,NSTEP) - N_mid + 1)*dt
    S_squared[:,1] = fft.ifft(displ[:]).real

    #debug
    #np.savetxt('tmp_S.py.txt',S_squared)

    # change the order of the time series
    # instead of having t=[0 T], we need t=[-T/2 T/2];
    temp = np.zeros(len(S_squared[:,1]))
    temp[arange(N_mid-1,NSTEP)] = np.copy(S_squared[arange(0,N_mid),1])
    temp[arange(0,N_mid-1)] = np.copy(S_squared[arange(N_mid,NSTEP),1])

    #debug
    #n = len(temp); tmp = zeros((n,2));
    #tmp[:,0] = arange(1,n+1); tmp[:,1] = temp[:];
    #np.savetxt('tmp_temp1.py.txt',tmp)


    # tapers ends
    if taper_type == 1:
        print('  tapering source time function S_squared')
        taper_length = int(round(taper_length_percentage*NSTEP))
        if taper_length < taper_length_min:
            taper_length = taper_length_min
        print('  using taper (cosine) length: %i' % taper_length)
        if NSTEP <= 2*taper_length:
            print('Error number of time steps = %i must be bigger than taper window size = %i ! Please retry... ' %(NSTEP,2*taper_length + 1))
            sys.exit(1)
        # cosine taper (both branches, value 0 at index 1, value 1 at index length+1, value 0 at 2*length+1
        taper = zeros(2*taper_length + 1)
        for l in arange(0,2*taper_length + 1):
            taper[l] = (1.0 - cos(pi * 2.0 * l / (2*taper_length))) / 2.0

        # sets up window function for multiplication with signal
        # note: this becomes too memory intensive for large number of time steps, using taper directly
        #window=ones(NSTEP);
        #
        ## adds increasing branch
        #window[0:taper_length] = taper[0:taper_length]
        #
        ## adds decreasing branch
        #window[NSTEP-taper_length:NSTEP] = taper[taper_length+1:2*taper_length+1]
        #
        #for l=0:NSTEP-1
        #  temp(l) = temp(l) * window(l);
        #end

        #debug
        #n = len(taper); tmp = zeros((n,2));
        #tmp[:,0] = arange(1,n+1); tmp[:,1] = taper[:];
        #np.savetxt('tmp_taper.py.txt',tmp)

        # applies taper
        for l in arange(0,taper_length):
            # increasing branch
            temp[l] = temp[l] * taper[l]
            k = NSTEP - taper_length + l
            temp[k] = temp[k] * taper[taper_length + 1 + l]

        #debug
        #n = len(temp); tmp = zeros((n,2));
        #tmp[:,0] = arange(1,n+1); tmp[:,1] = temp[:];
        #np.savetxt('tmp_temp2.py.txt',tmp)

        # plots figure
        if show_figures:
            print('  plotting figure: Taper and Taper window')
            fig = plt.figure(2)
            plt.subplot(2,1,1)
            plt.plot(taper)
            #hold('on')
            plt.xlabel('steps',size=fontsize)
            plt.ylabel('Taper',size=fontsize)
            plt.title('Taper',size=fontsize)
            #subplot(2,1,2);
            #plot(window);
            #xlabel('steps','fontsize',fontsize);ylabel('Taper window','fontsize',fontsize);
            #title('Window','fontsize',fontsize);
            #plt.show()

    # stores source time function values
    S_squared[:,1] = np.copy(temp[:])

    # # filter the source time function if needed
    # Wn=[1/Tmax 1/Tmin]/fmax;
    # [B,A] = butter(4,Wn);
    # S_squared[:,1] = filter(B,A,S_squared[:,1]);

    # plots figure
    if show_figures:
        print('  plotting figure: Source time function S_squared')
        plt.subplot(2,2,3)
        plt.plot(S_squared[:,0] / 60,S_squared[:,1],'r')
        plt.xlabel('Time (min)',size=fontsize)
        plt.ylabel('Amplitude',size=fontsize)
        plt.xlim(- T/2.0/60,T/2.0/60)
        plt.title('Source Time Function for Ensemble Forward Source',size=fontsize)
        plt.subplot(2,2,4)
        plt.plot(S_squared[:,0] / 60,S_squared[:,1],'r')
        plt.xlabel('Time (min)',size=fontsize)
        plt.ylabel('Amplitude',size=fontsize)
        plt.xlim(- Tmax*0.1/60,Tmax*0.1/60)
        plt.title('Zoom-in of the Source Time Function',size=fontsize)
        plt.show()
        plt.savefig("power_spectrum.png")

    ## output the source time function
    print('  saving S_squared as ASCII file')
    filename = 'S_squared'
    np.savetxt(filename,S_squared)

    ## user output
    DIR = os.getcwd()
    print('')
    print('*************************************************************')
    print('the source time function has been saved in:')
    print(''.join([DIR,'/S_squared']))
    print('')
    print('S_squared should be put into directory:')
    print('./NOISE_TOMOGRAPHY/ in the SPECFEM3D package')
    print('*************************************************************')
    print('')


def usage():
    print('usage: NOISE_TOMOGRAPHY.py NSTEP dt Tmin Tmax NOISE_MODEL [show_figures=1]')
    print('')
    print('Parameters:')
    print(' NSTEP       --- number of time steps (always odd for NOISE TOMOGRAPHY)')
    print(' dt          --- time interval of specfem3D solver')
    print(' Tmin, Tmax  --- the period range you are working with (in seconds)')
    print(' NOISE_MODEL --- the Peterson\'s noise model,either NLNM or NHNM:')
    print('                   NLNM: New Low  Noise Model (in 1993, the model was New)')
    print('                   NHNM: New High Noise Model')
    print('                 or')
    print('                   FLAT: for a flat noise spectrum (uniform noise within period range)')
    print('for example:')
    print('     ./NOISE_TOMOGRAPHY.py 1999 0.008 1.0 15.0 NLNM')
    print('')

if __name__ == '__main__':
    # gets arguments
    if len(sys.argv) < 6:
        usage()
        sys.exit(1)
    else:
        NSTEP = int(sys.argv[1])
        DT = float(sys.argv[2])
        T_min = float(sys.argv[3])
        T_max = float(sys.argv[4])
        model = sys.argv[5]
        if len(sys.argv) == 7:
            show_figures = int(sys.argv[6])
        else:
            show_figures = 0

    NOISE_TOMOGRAPHY(NSTEP,DT,T_min,T_max,model,show_figures)

