#!/bin/bash

# DKjob.sh  J.Hertzberg 22/07/2005

# Job pour la soumission à répétition de jobs de calcul :
# - teste la présence de job d'un utilisatuer donné
#   La chaine de recherche est à renseigner dans $CHERCHE
# - relance des jobs lorsque l'utilisateur n'a plus de jobs actifs

# Ne pas oublier de lancer le DKjob en nohup pour qu'il continue
# à s'exécuter après la déconnexion :
# nohup DKjob.sh &

CHERCHE="hertz"

while : ;
      do
      if [[ $(qstat|grep -c $CHERCHE) == 0 ]]; then
# Insérer ici la commande ou le script de soumission. Exemple :
#           /home/toto/script.sh;
#           qsub -q lmatlong.q /home/r/lma/hertz/doc/dev/facto.sh;
# Décommenter le break pour arrêter le DKjob après la soumission
#           break;
      fi
# Fréquence de soumission 24 h
      sleep 86400
      done

