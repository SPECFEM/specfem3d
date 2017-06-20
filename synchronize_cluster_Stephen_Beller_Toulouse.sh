#!/bin/bash

rsync -e ssh -avz --delete-after --exclude-from=excludeSynchronize ../DSM_FOR_SPECFEM3D beller@testo2.get.obs-mip.fr:~/.


