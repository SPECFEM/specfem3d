// updates R_memory in coarse-grained method
__device__  __forceinline__ void compute_element_att_memory_cg(const int tx,const int working_element,
                     const int i,const int j,const int k,const int FULL_ATTENUATION_SOLID,
                     const realw mul,const realw kappal,
                     realw_const_p factor_common,realw_const_p factor_common_kappa,
                     realw_const_p alphaval,realw_const_p betaval,realw_const_p gammaval,
                     realw_p R_trace,realw_p R_xx,realw_p R_yy,
                     realw_p R_xy,realw_p R_xz,realw_p R_yz,
                     const realw epsilondev_trace,const realw epsilondev_xx,const realw epsilondev_yy,
                     const realw epsilondev_xy,const realw epsilondev_xz,const realw epsilondev_yz,
                     const realw epsilondev_trace_loc,
                     const realw epsilondev_xx_loc,const realw epsilondev_yy_loc,const realw epsilondev_xy_loc,
                     const realw epsilondev_xz_loc,const realw epsilondev_yz_loc)
{


  int i_sls=(1-j%2)*(i%2+2*(k%2))+(j%2)*((i+1)%2+2*((k+1)%2));
  // use Runge-Kutta scheme to march in time

  // indices
  int offset = tx + NGLL3*working_element; // (i_sls(i,j,k),i,j,k,ispec)

  realw alphaval_loc = alphaval[i_sls]; // (i_sls)
  realw betaval_loc = betaval[i_sls];
  realw gammaval_loc = gammaval[i_sls];

  if(FULL_ATTENUATION_SOLID) {

    realw factor_loc = kappal * get_global_cr( &factor_common_kappa[offset] ); //kappastore(i,j,k,ispec) * factor_common(i_sls(i,j,k),i,j,k,ispec)
    realw Sn_trace   = factor_loc * epsilondev_trace;
    realw Snp1_trace   = factor_loc * epsilondev_trace_loc;
    realw rl= get_global_cr( &R_trace[offset] );
    R_trace[offset]=alphaval_loc * rl+betaval_loc * Sn_trace + gammaval_loc * Snp1_trace;

  }


  realw factor_loc = mul * get_global_cr( &factor_common[offset] ); //mustore(i,j,k,ispec) * factor_common(i_sls(i,j,k),i,j,k,ispec)
  // term in xx
  realw Sn_xx   = factor_loc * epsilondev_xx;
  realw Snp1_xx   = factor_loc * epsilondev_xx_loc; //(i,j,k)
  realw rxxl = get_global_cr( &R_xx[offset] );
  R_xx[offset]=alphaval_loc * rxxl + betaval_loc * Sn_xx + gammaval_loc * Snp1_xx;

  // term in yy
  realw Sn_yy   = factor_loc * epsilondev_yy;
  realw Snp1_yy   = factor_loc * epsilondev_yy_loc;
  realw ryyl = get_global_cr( &R_yy[offset] );
  R_yy[offset] = alphaval_loc * ryyl + betaval_loc * Sn_yy + gammaval_loc * Snp1_yy;

  // term in zz not computed since zero trace
  // term in xy
  realw Sn_xy   = factor_loc * epsilondev_xy;
  realw Snp1_xy   = factor_loc * epsilondev_xy_loc;
  realw rxyl = get_global_cr( &R_xy[offset] );
  R_xy[offset] = alphaval_loc * rxyl + betaval_loc * Sn_xy + gammaval_loc * Snp1_xy;

  // term in xz
  realw Sn_xz   = factor_loc * epsilondev_xz;
  realw Snp1_xz   = factor_loc * epsilondev_xz_loc;
  realw rxzl = get_global_cr( &R_xz[offset] );
  R_xz[offset] = alphaval_loc * rxzl + betaval_loc * Sn_xz + gammaval_loc * Snp1_xz;

  // term in yz
  realw Sn_yz   = factor_loc * epsilondev_yz;
  realw Snp1_yz   = factor_loc * epsilondev_yz_loc;
  realw ryzl = get_global_cr( &R_yz[offset] );
  R_yz[offset] = alphaval_loc * ryzl + betaval_loc * Sn_yz + gammaval_loc * Snp1_yz;

  return;

}

__device__ __forceinline__ void compute_element_att_stress_cg(const int tx,const int i,const int j,const int k,
                    const int FULL_ATTENUATION_SOLID,
                    const realw* R_trace,const realw* R_xx,const realw* R_yy,
                    const realw* R_xy,const realw* R_xz,const realw* R_yz,
                    realw* sigma_xx,realw* sigma_yy,realw* sigma_zz,
                    realw* sigma_xy,realw* sigma_xz,realw* sigma_yz) {

  realw R_trace_val=0.0f;

  // (Day 1998, Day and Bradley 2001): good
  realw R_xx_val = N_SLS*R_xx[tx];
  realw R_yy_val = N_SLS*R_yy[tx];
  realw R_xy_val = N_SLS*R_xy[tx];
  realw R_xz_val = N_SLS*R_xz[tx];
  realw R_yz_val = N_SLS*R_yz[tx];

  if(FULL_ATTENUATION_SOLID)
    R_trace_val = N_SLS*R_trace[tx];

  /*
  // Kristek and Moczo (2003): not good

  int ip=(i+1)%2;
  int im=(i-1>=0)? ((i-1)%2):ip;
  int jp=(j+1)%2;
  int jm=(j-1>=0)? ((j-1)%2):jp;
  int kp=(k+1)%2;
  int km=(k-1>=0)? ((k-1)%2):kp;

  realw R_xx_val = R_xx[tx]+0.5f*(R_xx[k*NGLL2+j*NGLLX+ip]+
          R_xx[k*NGLL2+j*NGLLX+im]+
          R_xx[k*NGLL2+jp*NGLLX+i]+
          R_xx[k*NGLL2+jm*NGLLX+i]+
          R_xx[kp*NGLL2+j*NGLLX+i]+
          R_xx[km*NGLL2+j*NGLLX+i]);

  realw R_yy_val = R_yy[tx]+0.5f*(R_yy[k*NGLL2+j*NGLLX+ip]+
          R_yy[k*NGLL2+j*NGLLX+im]+
          R_yy[k*NGLL2+jp*NGLLX+i]+
          R_yy[k*NGLL2+jm*NGLLX+i]+
          R_yy[kp*NGLL2+j*NGLLX+i]+
          R_yy[km*NGLL2+j*NGLLX+i]);

  realw R_xy_val = R_xy[tx]+0.5f*(R_xy[k*NGLL2+j*NGLLX+ip]+
          R_xy[k*NGLL2+j*NGLLX+im]+
          R_xy[k*NGLL2+jp*NGLLX+i]+
          R_xy[k*NGLL2+jm*NGLLX+i]+
          R_xy[kp*NGLL2+j*NGLLX+i]+
          R_xy[km*NGLL2+j*NGLLX+i]);

  realw R_xz_val = R_xz[tx]+0.5f*(R_xz[k*NGLL2+j*NGLLX+ip]+
          R_xz[k*NGLL2+j*NGLLX+im]+
          R_xz[k*NGLL2+jp*NGLLX+i]+
          R_xz[k*NGLL2+jm*NGLLX+i]+
          R_xz[kp*NGLL2+j*NGLLX+i]+
          R_xz[km*NGLL2+j*NGLLX+i]);

  realw R_yz_val = R_yz[tx]+0.5f*(R_yz[k*NGLL2+j*NGLLX+ip]+
          R_yz[k*NGLL2+j*NGLLX+im]+
          R_yz[k*NGLL2+jp*NGLLX+i]+
          R_yz[k*NGLL2+jm*NGLLX+i]+
          R_yz[kp*NGLL2+j*NGLLX+i]+
          R_yz[km*NGLL2+j*NGLLX+i]);

  if(FULL_ATTENUATION_SOLID)
    R_trace_val = R_trace[tx]+0.5f*(R_trace[k*NGLL2+j*NGLLX+ip]+
            R_trace[k*NGLL2+j*NGLLX+im]+
            R_trace[k*NGLL2+jp*NGLLX+i]+
            R_trace[k*NGLL2+jm*NGLLX+i]+
            R_trace[kp*NGLL2+j*NGLLX+i]+
            R_trace[km*NGLL2+j*NGLLX+i]);

  */

  *sigma_xx = *sigma_xx - R_xx_val-R_trace_val;
  *sigma_yy = *sigma_yy - R_yy_val-R_trace_val;
  *sigma_zz = *sigma_zz + R_xx_val + R_yy_val-R_trace_val;
  *sigma_xy = *sigma_xy - R_xy_val;
  *sigma_xz = *sigma_xz - R_xz_val;
  *sigma_yz = *sigma_yz - R_yz_val;

  return;
}
