#!/bin/sh
# Author: Hejun Zhu, hejunzhu@princeton.edu
# Princeton University, New Jersey, USA
# Last modified: Fri Mar  4 13:16:22 EST 2011


iter=M00
sigma_h=50
sigma_v=5


for tag in bulk_c_kernel_precond bulk_betav_kernel_precond bulk_betah_kernel_precond eta_kernel_precond
do
  echo $tag

  title="#PBS -N XSMOOTH_$tag"

  sed -e "s/^tag=.*$/tag=$tag/g" \
      -e "s/^iter=.*$/iter=$iter/g" \
      -e "s/^sigma_h=.*$/sigma_h=$sigma_h/g" \
      -e "s/^sigma_v=.*$/sigma_v=$sigma_v/g" \
      -e "s/^#PBS -N.*$/$title/g" \
            XPBS_smooth_kernel.sh > XPBS_smooth_kernel.sh.out

      mv XPBS_smooth_kernel.sh.out XPBS_smooth_kernel.sh
      qsub XPBS_smooth_kernel.sh
      sleep 3
done


