#==========================================================================
function compute_direction ()
{
#
# arguments : $1 = iter
#             $2 = 'alpha', 'beta' or 'gamma'
#             $3 = long_liss
#             $4 = long_lissv
#             $5 = sigma
#
cd bin
$HYBRID_BINNARY/xoptim <<EOF
$1
$2
$3 $3 $4
$5
EOF
cd ..
}
#=========================================================================
function gradient_v_2_gradient_poisson ()
{
#
# argument $1 : iter or -1 when linear search is used
#

cd bin
$HYBRID_BINNARY/gradient_v_2_gradient_poisson<<EOF
$1
EOF
cd ..
}
#==========================================================================
function linear_search_wolf ()
{
# using Wolfe rule to find the best step length for model perturbation

# linear search initialisation
declare -i ssiter nssiter
ssiter=1
nssiter=20
tg=0;
td=0;
pas=1.;

if [ "$iter" -eq 0 ]; then # one may use different step at iteration 0
pas=1.
fi

cd bin
echo $tg "  " $td "  " $pas > t_guess.txt
tail -1 ajsutement.txt > ajustement_0.txt
cd ..

while [ "$ssiter" -le "$nssiter" ] ; do

  echo $(date) >> $flog_file
  echo 'linear search ' $ssiter >> $flog_file

  rm bin/derivee_cout_pas.txt

  # on copie la base de donnee sauvee
  cp ./OUTPUT_FILES/DATABASES_MPI_CURRENT/* ./OUTPUT_FILES/DATABASES_MPI/.

  update_model $pas
  forward_simu ss_iter


  # calcul des sources adjointes
  compute_adjoint_source_FWI

  adjoint_simu ss_iter
  project_grad_in_tomo_grid pas
  project_model_in_tomo_grid pas
  gradient_v_2_gradient_poisson -1


  cd bin
  derivative_step ${iter} 'gamma'
  derivative_step ${iter} 'beta'
  tail -1 ajsutement.txt > ajustement_t.txt
  Wolfe_rule ${iter} ${ssiter}


  pas=$(< pas.txt)
  tg=$(< tg.txt)
  td=$(< td.txt)
  echo $tg "  " $td "  " $pas > t_guess.txt
  fin=$(< fin)

  # sortie
  cd ..

  ssiter="$ssiter+1"

  if [ "$fin" = "oui" ]; then
    ssiter="$nssiter+1"
  fi

done
}
#==============================================================
function Wolfe_rule ()
{
$HYBRID_BINNARY/regle_wolfe<<EOF
$1
$2
EOF
}
#=============================================================
function derivative_step ()
{
$HYBRID_BINNARY/xcalcule_derivee_pas<<EOF
$1
$2
EOF
}
#==============================================================
function update_model ()
{
cd bin

# 1. on projette la direction de descente de la grille tomo a la grille SEM
$MPIRUN $OPTION_MPI $HYBRID_BINNARY/xproject_sem 0 $SLICE
#mpirun -np $NPROC ./xproject_sem 0 $SLICE
# on calcule le nouveau modele
$MPIRUN $OPTION_MPI $HYBRID_BINNARY/xmodel_update $1
#mpirun -np  $NPROC ./xmodel_update $1
cd ..
}
#==================================================================================
function compute_adjoint_source_FWI ()
{
cd bin
cp ../DATA/define_adjoint_sources_for_FWI.par define_adjoint_sources.par
isrc=1
while [ "$isrc" -le "$nsrc" ]; do
   cp ../DATA/data_source_${isrc}.txt .
   $HYBRID_BINNARY/xdefineadj_no_aligne
   isrc="$isrc+1"
done
cd ..
}
#==============================================================================
