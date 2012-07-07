/*
 !=====================================================================
 !
 !               S p e c f e m 3 D  V e r s i o n  2 . 1
 !               ---------------------------------------
 !
 !          Main authors: Dimitri Komatitsch and Jeroen Tromp
 !    Princeton University, USA and CNRS / INRIA / University of Pau
 ! (c) Princeton University / California Institute of Technology and CNRS / INRIA / University of Pau
 !                             July 2012
 !
 ! This program is free software; you can redistribute it and/or modify
 ! it under the terms of the GNU General Public License as published by
 ! the Free Software Foundation; either version 2 of the License, or
 ! (at your option) any later version.
 !
 ! This program is distributed in the hope that it will be useful,
 ! but WITHOUT ANY WARRANTY; without even the implied warranty of
 ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 ! GNU General Public License for more details.
 !
 ! You should have received a copy of the GNU General Public License along
 ! with this program; if not, write to the Free Software Foundation, Inc.,
 ! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 !
 !=====================================================================
 */

#ifndef CUDA_HEADER_H
#define CUDA_HEADER_H

typedef float realw;  // type of "working" variables

/* ----------------------------------------------------------------------------------------------- */

// setters for these const arrays (very ugly hack, but will have to do)

// elastic
void setConst_hprime_xx(realw* array,Mesh* mp);
void setConst_hprime_yy(realw* array,Mesh* mp);
void setConst_hprime_zz(realw* array,Mesh* mp);

void setConst_hprimewgll_xx(realw* array,Mesh* mp);
void setConst_hprimewgll_yy(realw* array,Mesh* mp);
void setConst_hprimewgll_zz(realw* array,Mesh* mp);

void setConst_wgllwgll_xy(realw* array,Mesh* mp);
void setConst_wgllwgll_xz(realw* array, Mesh* mp);
void setConst_wgllwgll_yz(realw* array, Mesh* mp);

void setConst_wgll_cube(realw* array, Mesh* mp);

/* ----------------------------------------------------------------------------------------------- */

/* CUDA specific things from specfem3D_kernels.cu */

// older TEXTURE usage. For now just acoustic simulations. See usage
// of USE_TEXTURES_FIELDS elsewhere in code for elastic case
#ifdef USE_TEXTURES

// declaration of textures
texture<realw, 1, cudaReadModeElementType> tex_potential_acoustic;
texture<realw, 1, cudaReadModeElementType> tex_potential_dot_dot_acoustic;


  void bindTexturesPotential(realw* d_potential_acoustic)
  {
    cudaError_t err;

    cudaChannelFormatDesc channelDescFloat = cudaCreateChannelDesc<realw>();

    err = cudaBindTexture(NULL,tex_potential_acoustic, d_potential_acoustic,
                          channelDescFloat, NGLOB*sizeof(realw));
    if (err != cudaSuccess)
    {
      fprintf(stderr, "Error in bindTexturesPotential for potential_acoustic: %s\n", cudaGetErrorString(err));
      exit(1);
    }
  }

  void bindTexturesPotential_dot_dot(realw* d_potential_dot_dot_acoustic)
  {
    cudaError_t err;

    cudaChannelFormatDesc channelDescFloat = cudaCreateChannelDesc<realw>();

    err = cudaBindTexture(NULL,tex_potential_dot_dot_acoustic, d_potential_dot_dot_acoustic,
                          channelDescFloat, NGLOB*sizeof(realw));
    if (err != cudaSuccess)
    {
      fprintf(stderr, "Error in bindTexturesPotential_dot_dot for potential_dot_dot_acoustic: %s\n", cudaGetErrorString(err));
      exit(1);
    }
  }

#endif // USE_TEXTURES


#endif //CUDA_HEADER_H
