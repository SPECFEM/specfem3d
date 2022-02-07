#!/usr/bin/env python
#
# perturbs model with a Gaussian perturbation
#
#########################################
from __future__ import print_function

import sys
import numpy as np
import array

from helper_functions import helper

##########################################################

## default setup
NGLLX = 5
NGLLY = 5
NGLLZ = 5

IN_DATA_FILES = './OUTPUT_FILES/DATABASES_MPI/'  #'./DATA/'
KERNEL_FILES  = './KERNELS/'

# plotting model figures
do_plot_figures = False

##########################################################

def get_step_size(NPROC,SIM_TYPE,percent,param):
    """
    reads all kernel files and determines maximum absolute kernel value
    and step size for model update
    """
    # helper functions
    hf = helper()

    # initializes
    step_size = 0.0

    # global maximum value
    max_kernel_val_all = 0.0

    print("step size:")

    # loop over processes
    for myrank in range(0, NPROC):
        # processors name
        prname = IN_DATA_FILES + "proc{:06d}_".format(myrank)

        # K_rhop
        if SIM_TYPE == 1:
            # acoustic domain kernel
            filename = prname + "rhop_acoustic_kernel.bin"
        else:
            # elastic domain kernel
            filename = prname + "rhop_kernel.bin"
        print("  reading kernel {}".format(filename))
        krhop = hf.read_binary_SEM_file(filename)

        # coupled acoustic-elastic
        if SIM_TYPE == 3:
            filename = prname + "rhop_acoustic_kernel.bin"
            print("  reading kernel {}".format(filename))
            ktmp = hf.read_binary_SEM_file(filename)
            krhop += ktmp.copy()

        # K_alpha
        if SIM_TYPE == 1:
            # acoustic domain kernel
            filename = prname + "alpha_acoustic_kernel.bin"  # 2D version: "c_acoustic_kernel.bin"
        else:
            # elastic domain kernel
            filename = prname + "alpha_kernel.bin"
        print("  reading kernel {}".format(filename))
        kalpha = hf.read_binary_SEM_file(filename)

        # coupled acoustic-elastic
        if SIM_TYPE == 3:
            filename = prname + "alpha_acoustic_kernel.bin"
            print("  reading kernel {}".format(filename))
            ktmp = hf.read_binary_SEM_file(filename)
            kalpha += ktmp.copy()

        # K_beta
        if SIM_TYPE == 1:
            # acoustic domain, zero shear for acoustic simulation
            kbeta = np.zeros((kalpha.size))
        else:
            # elastic domain kernel
            filename = prname + "beta_kernel.bin"
            print("  reading kernel {}".format(filename))
            kbeta = hf.read_binary_SEM_file(filename)

        # note: model update is using a gradient in the negative direction to decrease misfits
        #       this is a steepest descent update, using a step size of getting a maximum 1 percent update
        krhop_max = np.abs(krhop).max()
        kalpha_max = np.abs(kalpha).max()
        kbeta_max = np.abs(kbeta).max()

        print("    absolute values:")
        print("    kernel rhop : max = ",krhop_max)
        print("    kernel alpha: max = ",kalpha_max)
        print("    kernel beta : max = ",kbeta_max)
        print("")

        # determines maximum value
        if param == "rho":
            max_val = krhop_max
        elif param == "vp":
            max_val = kalpha_max
        elif param == "vs":
            max_val = kbeta_max
        else:
            max_val = max(krhop_max,kalpha_max,kbeta_max)

        # total maximum over all slices
        if max_val > max_kernel_val_all: max_kernel_val_all = max_val

    # gradient step length factor
    step_factor = percent

    # step size
    if max_kernel_val_all != 0.0:
        step_size = step_factor / max_kernel_val_all
    else:
        print("Error: kernels are all zero, please check...")
        sys.exit(1)

    # user output
    print("  kernel parameterization considered = ",param)
    print("  kernel maximum value               = ",max_kernel_val_all)
    print("  maximum gradient step dln(m)       = ",step_factor)
    print("")
    print("  resulting step size                = ",step_size)
    print("")

    return step_size
#
#------------------------------------------------------------------------------------------
#

def model_update(NPROC,SIM_TYPE,percent,param):
    """
    update model using kernels
    """
    global NGLLX,NGLLY,NGLLZ
    global do_plot_figures

    # helper functions
    hf = helper()

    # user output
    print("model update:")
    print("  NPROC        = ",NPROC)
    print("  SIM_TYPE     = ",SIM_TYPE)
    print("  percent      = ",percent)
    print("  parameter : ",param)
    print("")
    print("  model   input directory: ",IN_DATA_FILES)
    print("  kernels input directory: ",KERNEL_FILES)
    print("")

    # checks parameter
    if not (param == "all" or param == "rho" or param == "vp" or param == "vs"):
        print("Invalid input parameter {}, please choose rho, vp, vs (or all)".format(param))
        sys.exit(1)

    # determines kernel amplitudes and step size for model update
    step_size = get_step_size(NPROC,SIM_TYPE,percent,param)

    # loop over processes
    for myrank in range(0, NPROC):
        print("")
        print("slice {}:".format(myrank))
        print("  reading in model files for rank {}...".format(myrank))
        print("")

        # processors name
        prname = IN_DATA_FILES + "proc{:06d}_".format(myrank)

        # for format of stored arrays:
        #   see routine save_arrays_solver_files() in src/generate_databases
        #
        # we require ibool(NGLLX,NGLLY,NGLLZ,nspec) since xstore(nglob)/ystore(nglob)/zstore(nglobe)
        # are stored on global level, not local level

        # mesh index
        filename = prname + "ibool.bin"
        ibool = hf.read_binary_ibool_file(filename=filename)

        # determines nspec
        total_size = ibool.size
        if total_size % (NGLLX * NGLLY * NGLLZ) == 0:
            nspec = int(total_size / (NGLLX * NGLLY * NGLLZ) )
            print("  ibool: number of elements nspec = {}".format(nspec))
            print("         NGLLX/NGLLY/NGLLZ        = {}/{}/{}".format(NGLLX,NGLLY,NGLLZ))
            print("")
        else:
            print("Error reading ibool: total size {} doesn't match NGLLX/NGLLY/NGLLZ == {}/{}/{}, please check...".format(total_size,NGLLX,NGLLY,NGLLZ))
            sys.exit(1)

        # not needed to reshape array, but would be possible
        #ibool = np.reshape(ibool, (NGLLX, NGLLY, NGLLZ, nspec), order='F') # fortran-like index ordering

        # coordinates
        # global point arrays
        #  xstore_unique(nglob)
        filename = prname + "x.bin"
        xstore = hf.read_binary_SEM_file(filename)
        filename = prname + "y.bin"
        ystore = hf.read_binary_SEM_file(filename)
        filename = prname + "z.bin"
        zstore = hf.read_binary_SEM_file(filename)

        # user output
        print("  model dimension: x min/max = ",xstore.min(),xstore.max())
        print("                   y min/max = ",ystore.min(),ystore.max())
        print("                   z min/max = ",zstore.min(),zstore.max())
        print("")

        # model files
        print("  initial model:")
        print("    directory: ",IN_DATA_FILES)

        # model rho,vp,vs
        filename = prname + "rho.bin"
        rhostore = hf.read_binary_SEM_file(filename)
        filename = prname + "vp.bin"
        vpstore = hf.read_binary_SEM_file(filename)
        filename = prname + "vs.bin"
        vsstore = hf.read_binary_SEM_file(filename)

        print("")
        print("    rho min/max = ",rhostore.min(),rhostore.max())
        print("    vp  min/max = ",vpstore.min(),vpstore.max())
        print("    vs  min/max = ",vsstore.min(),vsstore.max())

        # plot images
        if do_plot_figures:
            hf.plot_model_image_3D(xstore,ystore,zstore,ibool,rhostore,prname + "rho",verbose=True)
            hf.plot_model_image_3D(xstore,ystore,zstore,ibool,vpstore,prname + "vp",verbose=True)
            hf.plot_model_image_3D(xstore,ystore,zstore,ibool,vsstore,prname + "vs",verbose=True)

        # kernels
        print("")
        print("  kernels:")
        print("    input directory: ",KERNEL_FILES)

        prname = KERNEL_FILES + "proc{:06d}_".format(myrank)

        # K_rhop
        if SIM_TYPE == 1:
            # acoustic domain kernel
            filename = prname + "rhop_acoustic_kernel.bin"
        else:
            # elastic domain kernel
            filename = prname + "rhop_kernel.bin"
        krhop = hf.read_binary_SEM_file(filename)

        # coupled acoustic-elastic
        if SIM_TYPE == 3:
            filename = prname + "rhop_acoustic_kernel.bin"
            ktmp = hf.read_binary_SEM_file(filename)
            krhop += ktmp.copy()

        # plotting
        if do_plot_figures:
            hf.plot_model_image_3D(xstore,ystore,zstore,ibool,krhop,filename,plot_kernel=True,verbose=True)

        # K_alpha
        if SIM_TYPE == 1:
            # acoustic domain kernel
            filename = prname + "alpha_acoustic_kernel.bin"  # 2D version: "c_acoustic_kernel.bin"
        else:
            # elastic domain kernel
            filename = prname + "alpha_kernel.bin"
        kalpha = hf.read_binary_SEM_file(filename)

        # coupled acoustic-elastic
        if SIM_TYPE == 3:
            filename = prname + "alpha_acoustic_kernel.bin"
            ktmp = hf.read_binary_SEM_file(filename)
            kalpha += ktmp.copy()

        # plotting
        if do_plot_figures:
            hf.plot_model_image_3D(xstore,ystore,zstore,ibool,kalpha,filename,plot_kernel=True,verbose=True)

        # K_beta
        if SIM_TYPE == 1:
            # acoustic domain, zero shear for acoustic simulation
            kbeta = np.zeros((kalpha.size))
            print("    acoustic simulation: setting zero shear kernel")
        else:
            # elastic domain kernel
            filename = prname + "beta_kernel.bin"
            kbeta = hf.read_binary_SEM_file(filename)
            # plotting
            if do_plot_figures:
                hf.plot_model_image_3D(xstore,ystore,zstore,ibool,kbeta,filename,plot_kernel=True,verbose=True)

        print("")
        print("    kernel rhop : min/max = ",krhop.min(),krhop.max())
        print("    kernel alpha: min/max = ",kalpha.min(),kalpha.max())
        print("    kernel beta : min/max = ",kbeta.min(),kbeta.max())
        print("")
        print("    norm of rhop kernel  = ",np.sum(krhop * krhop))
        print("    norm of alpha kernel = ",np.sum(kalpha * kalpha))
        print("    norm of beta kernel  = ",np.sum(kbeta * kbeta))

        # model update
        print("")
        print("  model update:")
        print("    parameter             = ",param)
        print("    step size             = ",step_size)
        print("")

        # note: model update is using a gradient in the negative direction to decrease misfits
        #       this is a steepest descent update, using a step size of getting a maximum 1 percent update

        # gradients in negative direction
        if param == "rho":
            drho   = - step_size * krhop
            dalpha = 0.0
            dbeta  = 0.0
        elif param == "vp":
            dalpha = - step_size * kalpha
            drho   = 0.0
            dbeta  = 0.0
        elif param == "vs":
            dbeta  = - step_size * kbeta
            drho   = 0.0
            dalpha = 0.0
        else:
            drho   = - step_size * krhop
            dalpha = - step_size * kalpha
            dbeta  = - step_size * kbeta

        # update
        # kernels are for relative perturbations dln(m)
        #
        # dln(m) = G -> m_new = m * exp( dln(m) ) = m * exp( G )
        rhostore = rhostore * np.exp(drho)
        vpstore  = vpstore * np.exp(dalpha)
        vsstore  = vsstore * np.exp(dbeta)
        # or
        # to first order: dln(m) = dm/m = G -> dm = G * m
        #                                      m_new = m + dm = m + G * m = m (1 + G)
        #
        #rhostore = rhostore * (1.0 + drho)
        #vpstore  = vpstore * (1.0 + dalpha)
        #vsstore  = vsstore * (1.0 + dbeta)

        print("  updated model:")
        print("    rho min/max = ",rhostore.min(),rhostore.max())
        print("    vp  min/max = ",vpstore.min(),vpstore.max())
        print("    vs  min/max = ",vsstore.min(),vsstore.max())
        print("")

        # stores perturbed model files
        print("  storing new files...")
        prname = IN_DATA_FILES + "proc{:06d}_".format(myrank)

        # rho
        filename = prname + "rho_new.bin"
        hf.write_binary_file_custom_real_array(filename,rhostore)
        # vp
        filename = prname + "vp_new.bin"
        hf.write_binary_file_custom_real_array(filename,vpstore)
        # vs
        filename = prname + "vs_new.bin"
        hf.write_binary_file_custom_real_array(filename,vsstore)

        print("  updated model files written to:")
        print("    ",prname + "rho_new.bin")
        print("    ",prname + "vp_new.bin")
        print("    ",prname + "vs_new.bin")
        print("")

        # plot images
        if do_plot_figures:
            hf.plot_model_image_3D(xstore,ystore,zstore,ibool,rhostore,prname + "rho_new",verbose=True)
            hf.plot_model_image_3D(xstore,ystore,zstore,ibool,vpstore,prname + "vp_new",verbose=True)
            hf.plot_model_image_3D(xstore,ystore,zstore,ibool,vsstore,prname + "vs_new",verbose=True)


#
#------------------------------------------------------------------------------------------
#

def usage():
    print("Usage: ./model_update.py NPROC SIM_TYPE percent [param]")
    print("")
    print("    NPROC    - number of parallel processes")
    print("    SIM_TYPE - benchmark simulation type (1 == acoustic / 2 == elastic / 3 == coupled acoustic-elastic)")
    print("    percent  - perturbation strength")
    print("    param    - (optional) model parameter, e.g., rho, vp or vs. [default is all]")
    sys.exit(1)

#
#------------------------------------------------------------------------------------------
#

if __name__ == '__main__':
    # gets arguments
    if len(sys.argv) < 4:
        usage()

    ## input parameters
    NPROC = int(sys.argv[1])
    SIM_TYPE = int(sys.argv[2])
    percent = float(sys.argv[3])

    # parameter to update
    if len(sys.argv) == 5:
        param = sys.argv[4]
    else:
        param = "all"

    # checks
    if percent <= 0.0:
        print("Invalid percent {} for step_factor, must be strictly positive".format(percent))
        sys.exit(1)

    model_update(NPROC,SIM_TYPE,percent,param)

    print("")
    print("all done")
    print("")
