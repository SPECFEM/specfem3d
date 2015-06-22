      integer mfhdr,mnhdr,mihdr,mlhdr,mkhdr
      parameter (mfhdr=70, mnhdr=15, mihdr=20, mlhdr=5, mkhdr=24)
c common variables
      common / hdr / fhdr,nhdr,ihdr,lhdr,khdr
c header variables (I)
      real fhdr(mfhdr)
      integer nhdr(mnhdr),ihdr(mihdr),lhdr(mlhdr)
      character*8 khdr(mkhdr)
