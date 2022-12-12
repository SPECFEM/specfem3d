#!/usr/bin/env python
#
# Generates noise distribution and direction binary files
# according to parfile_noise.yaml.
#
# The binary files contain float values for all GLL points
# on the mesh free-surface.
#
# Make sure the free-surface corresponds to the mesh surface
# where you want noise sources to be located.
#
# IMPORTANT: This script needs proc*_free_surface.vtk files in
# DATABASES_MPI. To generate them, set SAVE_MESH_FILES = .true.
# in DATA/Par_file and run xgenerate_databases.
#


import matplotlib.pyplot as plt
import numpy as np
import os
import yaml
from convert_utm2geo import utm_geo


def gaussian_distribution(distribution, mask_noise, xcoor, ycoor):
    if distribution['utm_zone']:
        center_x, center_y = utm_geo(distribution['center_x'],
                                     distribution['center_y'],
                                     distribution['utm_zone'],
                                     2)
    else:
        center_x = distribution['center_x']
        center_y = distribution['center_y']

    dist = np.zeros((xcoor.size))
    a = np.array([center_x, center_y])

    for i in range(0, dist.size):
        b = np.array([xcoor[i], ycoor[i]])
        dist[i] = np.linalg.norm(a-b)

    gaussian = np.exp(- (dist ** 2) / (2 * distribution['sigma_m'] ** 2))
    mask_noise += gaussian * distribution['weight']
    return mask_noise


def uniform_distribution(distribution, mask_noise):
    mask_noise[:] += distribution['weight']
    return mask_noise


def ocean_distribution(distribution, mask_noise, zcoor):
    ocean_idx = np.where(zcoor < 0.0)[0]

    if not ocean_idx.any():
        print('Skipping ocean distribution: There are no gll points on\
               the ocean.')
    else:
        mask_noise[ocean_idx] += distribution['weight']
    return mask_noise


def read_gll_coordinates(ifile):
    X = np.genfromtxt(ifile, skip_header=5)
    xcoor = X[:, 0]
    ycoor = X[:, 1]
    zcoor = X[:, 2]
    return xcoor, ycoor, zcoor


def write_files(proc, mask_noise, normal_x_noise, normal_y_noise,
                normal_z_noise, outdir):
    _write(mask_noise,
           os.path.join(outdir, 'proc{:06}_mask_noise.bin'.format(proc)))
    _write(normal_x_noise,
           os.path.join(outdir, 'proc{:06}_normal_x_noise.bin'.format(proc)))
    _write(normal_y_noise,
           os.path.join(outdir, 'proc{:06}_normal_y_noise.bin'.format(proc)))
    _write(normal_z_noise,
           os.path.join(outdir, 'proc{:06}_normal_z_noise.bin'.format(proc)))
    return


def _write(x, filename):
    n = np.array([4 * len(x)], dtype='int32')
    x = np.array(x, dtype='float32')

    with open(filename, 'wb') as file:
        n.tofile(file)
        x.tofile(file)
        n.tofile(file)


if __name__ == '__main__':
    try:
        with open('parfile_noise.yaml', 'r') as _file:
            par = yaml.load(_file)
    except IOError:
        print('IOError: parfile_noise.yaml not found.')

    ngll_surface = par['NSPEC_SURFACE_EXT_MESH'] * par['NGLLSQUARE']
    proc_ngll_surface = ngll_surface // par['NPROC']

    if par['PLOT_MASK']:
        all_gll_x = np.array([])
        all_gll_y = np.array([])
        all_mask_noise = np.array([])

    for p in range(0, par['NPROC']):
        print('Generating files for process {}'.format(p))

        mask_noise = np.zeros(proc_ngll_surface, dtype='float32')
        normal_x_noise = np.zeros(proc_ngll_surface, dtype='float32')
        normal_y_noise = np.zeros(proc_ngll_surface, dtype='float32')
        normal_z_noise = np.zeros(proc_ngll_surface, dtype='float32')

        # read gll coordinates
        free_surface_file = os.path.join(par['DATABASES_MPI'],
                                         'proc{:06}_free_surface.vtk'.format(p))
        gll_x, gll_y, gll_z = read_gll_coordinates(free_surface_file)

        # set noise direction
        normal_x_noise[:] = par['XDIR']
        normal_y_noise[:] = par['YDIR']
        normal_z_noise[:] = par['ZDIR']

        # set noise distribution
        for d in par['DISTRIBUTIONS']:
            if d['type'] == 'uniform':
                mask_noise = uniform_distribution(d, mask_noise)
            elif d['type'] == 'ocean':
                mask_noise = ocean_distribution(d, mask_noise, gll_z)
            elif d['type'] == 'gaussian':
                mask_noise = gaussian_distribution(d, mask_noise, gll_x, gll_y)
            else:
                print('Undefined noise distribution.')

        if par['WRITE_FILES']:
            write_files(p, mask_noise, normal_x_noise, normal_y_noise,
                        normal_z_noise, par['DATABASES_MPI'])

        if par['PLOT_MASK']:
            all_gll_x = np.append(all_gll_x, gll_x)
            all_gll_y = np.append(all_gll_y, gll_y)
            all_mask_noise = np.append(all_mask_noise, mask_noise)

    if par['PLOT_MASK']:
        fig, ax = plt.subplots()
        im = ax.scatter(all_gll_x, all_gll_y, c=all_mask_noise)
        ax.set_xlabel('x [m]')
        ax.set_ylabel('y [m]')
        ax.set_title('Noise distribution mask')
        cbar = plt.colorbar(im)
        cbar.set_label('Weight')
        plt.show()
        plt.close()
