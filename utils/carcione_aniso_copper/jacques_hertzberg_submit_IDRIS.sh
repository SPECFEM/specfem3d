#!/bin/bash

# DKjob.sh  J.Hertzberg 22/07/2005

# Job pour la soumission � r�p�tition de jobs de calcul :
# - teste la pr�sence de job d'un utilisatuer donn�
#   La chaine de recherche est � renseigner dans $CHERCHE
# - relance des jobs lorsque l'utilisateur n'a plus de jobs actifs

# Ne pas oublier de lancer le DKjob en nohup pour qu'il continue
# � s'ex�cuter apr�s la d�connexion :
# nohup DKjob.sh &

CHERCHE="hertz"

while : ;
      do
      if [[ $(qstat|grep -c $CHERCHE) == 0 ]]; then
# Ins�rer ici la commande ou le script de soumission. Exemple :
#           /home/toto/script.sh;
#           qsub -q lmatlong.q /home/r/lma/hertz/doc/dev/facto.sh;
# D�commenter le break pour arr�ter le DKjob apr�s la soumission
#           break;
      fi
# Fr�quence de soumission 24 h
      sleep 86400
      done

