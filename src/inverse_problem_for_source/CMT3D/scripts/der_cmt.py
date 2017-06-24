#!/usr/bin/env python

import sys
from cmt import *

# derivative cmts
dmoment=1e24; ddepth=2; dlocation=2; npar=9
cmt='CMTSOLUTION'; nevent=1
utm=0; utm_x=0.; utm_y=0.
#ext=['Mrr','Mtt','Mpp','Mrt','Mrp','Mtp','dep','lon','lat']

# generate derivative cmts
if gen_cmt_der(cmt,npar=npar,dmoment=dmoment,ddepth=ddepth,dlocation=dlocation,utm=utm,nevent=nevent) != 0:  # ext has been corrected
  sys.exit('Error generating derivative CMTs')

