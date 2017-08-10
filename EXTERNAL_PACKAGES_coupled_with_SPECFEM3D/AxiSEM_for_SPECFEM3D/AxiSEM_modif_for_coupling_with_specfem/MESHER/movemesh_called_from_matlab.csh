#!/bin/csh -f

if ( ${#argv} < 1 || "$1" == "-h" ) then
    echo "USAGE:"
    echo "  Argument 1:  name for the mesh"
    echo ""
    exit
endif

set homepath = `echo $PWD`
set meshpath = "../SOLVER/MESHES/$1"

if ( ${#argv} == 2 ) then
  set xmeshrundir = $2
else
  set xmeshrundir = "."
endif
echo "Moving mesh from $xmeshrundir to $meshpath"

# test if mesher finished without problems:
if (`grep 'DONE WITH MESHER' $xmeshrundir/OUTPUT | wc -l` != '1') then
  echo "ERROR: MESHER did not finish yet or encountered a problem."
  echo "       Check 'OUTPUT' file for more details."
  exit
else
  echo 'MESHER finished smoothly...'
endif

if ( ! -d ../SOLVER/MESHES ) then
  mkdir ../SOLVER/MESHES
endif

if ( ! -d $meshpath ) then
  mkdir $meshpath
else
  echo "ERROR: the mesh folder " $meshpath " already exists!"
  exit
endif

# convert relative to absolute path:
cd $meshpath
set meshpath = `echo $PWD`
cd $homepath

echo "Moving mesh to" $meshpath

mv $xmeshrundir/meshdb.dat* $meshpath
mv $xmeshrundir/mesh_params.h $meshpath
mv $xmeshrundir/OUTPUT $meshpath
#mv Diags $meshpath
cp -p $xmeshrundir/inparam_mesh $meshpath

set bgmodel = `grep "^BACKGROUND_MODEL" inparam_mesh | awk '{print $2}'`
echo $bgmodel
if ( $bgmodel == 'external') then
  set extmodel = `grep "^EXT_MODEL" inparam_mesh | awk '{print $2}'`
  echo "Copying external model file $extmodel"
  if ( -f $extmodel) then
    cp -p $extmodel $meshpath/external_model.bm
  else
    echo "external model file $extmodel does not exist. Did you delete it? Why?"
    exit
  endif
endif

mkdir $meshpath/Code
cp -p *.[fF]90 $meshpath/Code
cp -p Makefile $meshpath/Code
cp -p makemake.pl $meshpath/Code
cp -p inparam_mesh $meshpath/Code
cp -p xmesh $meshpath/Code
cp -p submit.csh $meshpath/Code
cp -p inparam_mesh $meshpath/Code

mv $xmeshrundir/Diags/* $meshpath

echo "Contents in" $meshpath ":"
ls $meshpath
cd $meshpath
echo "Total size: `du -sh` "
echo "DONE."
