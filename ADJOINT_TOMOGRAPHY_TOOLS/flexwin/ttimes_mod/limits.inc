c   The "j" parameters (1 st line) are intended to be user settable:
c        jmod   Maximum number of (rough) model points.
c        jslo   Maximum number of discrete ray parameters.
c        jdep   Maximum number of discrete model slowness samples.
c               Note:  jdep always need to be larger than jslo by
c               approximately 50% to account for multi-valuedness
c               due to high slowness zones.
c        jsrc   Maximum number of discrete model slowness samples above
c               the maximum source depth of interest.
c        jbrh   Maximum number of tau branches (model discontinuities
c               plus one).
c        jlvz   Maximum number of low velocity zones.
c        jseg   Maximum number of different types of travel-times
c               considered.
c        jbrn   Maximum number of different travel-time branches to be
c               searched.
c        jout   Maximum length of all travel-time branches strung
c               together.
c
      parameter(jmod=160,jslo=350,jdep=450,jsrc=150,jbrh=20,jlvz=5)
      parameter(jseg=30,jbrn=100,jout=2250)
c   The parameters actually used are all derivatives of the "j"
c   parameters and cannot be changed by the user.
      parameter(nmd0=jmod,nsl1=jslo+1,ndp1=jdep+1,nsr0=jsrc)
      parameter(nbr1=jbrh+2,nbr2=jbrh*2,ncp0=2*(jbrh+jlvz),nlvz0=jlvz)
