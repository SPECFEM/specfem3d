#!/bin/sh
# Author: Hejun Zhu, hejunzhu@princeton.edu
# Princeton University, New Jersey, USA
# Last modified: Wed Mar 16 14:57:17 EDT 2011


for dir in CMTSOLUTION_*
do 
	echo $dir 
	tar -czvf $dir.tar.gz $dir 
	rm -rf $dir 

done 
