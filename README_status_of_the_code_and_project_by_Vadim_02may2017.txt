
Actuellement, la version de l'inversion couplee avec FK fonctionne + parallelisme sur les sources. Il y a un petit exemple dans  DSM_FOR_SPECFEM3D/Examples/small_exemple_FWI_with_fk/fwi_fk.
J'ai pas encore fini l'exemple avec le couplage AxiSEM, mais ca devrait marcher comme FK.

Actuelement, je suis sur l'UNDO_ATTENUATION + GPU, je ferai l'exemple AxiSEM apres.

Pour passer aux applis, il faut qu'on voie comment vous avez organise les donnees, sans doute il y aura des parametres d'entree supplementaires a prendre en compte que j'ai pas prevu. Il faut aussi tester la
convolution des residus avec l'ondelette source, c'est code mais j'ai pas teste.

Vadim.
02 mai 2017

