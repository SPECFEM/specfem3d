#!/bin/csh -f

# DK DK filter and convert Northridge data to velocity in SAC
# DK DK read sac ALPHA (ascii) format as input

foreach file ($*)
        echo $file

        sac2000 << FIN 
qdp off
read alpha ${file}
bandpass butter corners 0.1 1.0 npoles 4 passes 2
rmean
int
write alpha ${file}.veloc
q
FIN

end


