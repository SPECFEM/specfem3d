function iteration_fwi ()
{
long_liss=$(< longueur_lissage.txt)
long_lissv=$(< longueur_lissage_vertical.txt)
sigma_beta=$(< sigma_beta.txt)
sigma_gamma=$(< sigma_gamma.txt)


echo ' ' >> $flog_file
echo ' iteration ' $iter >> $flog_file
echo ' '  >> $flog_file
echo $(date) >> $flog_file
echo ' ' >> $flog_file
echo 'longeur de lissage ' $long_liss >> $flog_file
echo 'longeur de lissage v' $long_lissv >> $flog_file
echo 'sigma beta ' $sigma_beta >> $flog_file
echo 'sigma gamma ' $sigma_gamma >> $flog_file

echo ' on sauve la base de donnee SPECFEM '  $iter >> $flog_file
save_cp_current_database
echo $(date) >> $flog_file


echo ' calcul direct SPECFEM3D avec sauvegarde du champ final  ' $iter >> $flog_file
forward_simu iter_${iter}
echo $(date) >> $flog_file

echo ' calcul des sources adjointes  ' $iter >> $flog_file
compute_adjoint_source_FWI
echo $(date) >> $flog_file


echo  'calcul du gradient par etat adjoint de SPECFEM3D ' $iter >> $flog_file
adjoint_simu iter_${iter}
echo $(date) >> $flog_file

echo ' projection du gradient sur la grille tomo  ' $iter  >> $flog_file
project_grad_in_tomo_grid ${iter}
echo $(date) >> $flog_file

echo ' projection du modele sur la grille tomo  ' $iter >> $flog_file
project_model_in_tomo_grid ${iter}
echo $(date) >> $flog_file

cp bin/path_file_for_gradient_0_alpha.par bin/path_file_for_gradient.par  # chemin d'acces aux gradients

echo 'on forme le gradient et le modele gamma ' $iter >> $flog_file
gradient_v_2_gradient_poisson ${iter}

echo ' on calcule la direction de descente gamma  ' $iter >> $flog_file
compute_direction ${iter} 'gamma'  ${long_liss} ${long_lissv} ${sigma_gamma}

echo ' on sauve la direction de descente  ' $iter >> $flog_file
cp bin/direction_gamma.bin bin/direction_gamma_${iter}.bin
echo $(date) >> $flog_file

echo ' on calcule la direction de descente beta  ' $iter >> $flog_file
cp bin/path_file_for_gradient_0_beta.par bin/path_file_for_gradient.par # chemin d'acces au gradients

compute_direction ${iter} 'beta'  ${long_liss} ${long_lissv} ${sigma_beta}

echo ' on sauve la direction de descente  '  $iter >> $flog_file
cp bin/direction_beta.bin bin/direction_beta_${iter}.bin
echo $(date) >> $flog_file

echo  'recherche lineaire du pas optimal  ' $iter >> $flog_file
linear_search_wolf
echo $(date) >> $flog_file

rm bin/derivee_cout.txt  # il faur effacer la valeur de q prime


}
