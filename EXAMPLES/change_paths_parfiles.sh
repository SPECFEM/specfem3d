#!/bin/bash

par_files=`find . -iname \*Par_file\* -exec ls -1 {} \; | grep -v Mesh_Par_file | grep -v Par_file_faults`

for f in $par_files; do
  echo $f
  sed -i 's/\(TOMOGRAPHY_PATH\s*=\s*\)\.\(\.\/DATA\)/\1\2/' $f
  sed -i 's/\(LOCAL_PATH\s*=\s*\)\.\(\.\/OUTPUT\)/\1\2/' $f
done

mesh_par_files=`find . -iname \*Mesh_Par_file\* -exec ls -1 {} \;`

for f in $mesh_par_files; do
  sed -i 's/\(LOCAL_PATH\s*=\s*\)\.\(\.\/OUTPUT\)/\1\2/' $f
done
