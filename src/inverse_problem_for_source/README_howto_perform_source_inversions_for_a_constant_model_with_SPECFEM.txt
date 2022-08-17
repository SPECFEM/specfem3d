
Here is what SPECFEM3D users classically do to relocate earthquake sources (paragraph taken from a
text written by Federica Magnoni and Emanuele Casarotti for Italy, with minor changes by Dimitri Komatitsch):

"As initial source parameters of the chosen events we use Time Domain Moment Tensor solutions
(TDMT, Dreger and Helmberger, 1993; Scognamiglio et al. 2009) that are calculated by inverting full
three-component traces of regional broad-band stations using 1D wave speed models. Before starting
the tomographic inversion, we then exploit the initial 3D model of the structure under study
to recalculate the source parameters of the considered events. The TDMT solutions, based on
1D wave speed models, are considered as initial solutions and we invert for the six moment tensor
components and the event location by using the technique presented in Liu et al. (2004). We use
SPECFEM3D to numerically calculate the Frechet derivatives with respect to the source parameters
based on 3D models and we minimize a waveform misfit function to invert for source solutions. If
the variance reduction of the solution with the 3D model is significant, this is considered as the
new moment tensor solution of the event."

This is performed using one of the two following packages, depending on the language you prefer to use:
CMT3D (written by Qinya Liu et al. in Fortran), or pyCMT3D (same package, but translated to Python by Wenjie Lei).
The two packages are included in this directory.

This source inversion technique comes from the following paper by Qinya Liu et al.:
https://specfem.github.io/komatitsch.free.fr/published_papers/bssa_Qinya_inversion_2004.pdf
and in which the Frechet derivatives are computed using finite-differences.
Here is the full reference of the paper:

@Article{LiPoKoTr04,
  Title                    = {Spectral-element moment tensor inversions for earthquakes in {S}outhern {C}alifornia},
  Author                   = {Qinya Liu and Jascha Polet and Dimitri Komatitsch and Jeroen Tromp},
  Journal                  = {Bulletin of the Seismological Society of America},
  Year                     = {2004},
  Number                   = {5},
  Pages                    = {1748-1761},
  Volume                   = {94},
  Doi                      = {10.1785/012004038}}

