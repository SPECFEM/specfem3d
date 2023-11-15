
c nb de pas de temps
      integer ncycl
      real dtinc
      parameter(ncycl=7500)
      parameter(dtinc=6.5e-3)

c stockage sismogrammes
      integer nseis,isamp
      parameter(isamp=10)
      parameter(nseis=ncycl/isamp)

      integer itaff,itfirstaff
      parameter(itaff=200)
      parameter(itfirstaff=15)

      integer iout
      parameter(iout=6)

