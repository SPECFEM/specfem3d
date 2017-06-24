
From Vadim, 05 aug 2016:

Voici ou j'en suis sur l'ecriture  code d'inversion :

sous SVN j'ai mis en place l'arborescence du code d'inversion dans le repertoire ./inverse_problem

j'ai defini un sous repertoire par sous taches :

./input_output

./adjoint_source

./inversion_scheme

./regularization

Le programme principal est program_inverse_problem.f90 qui appelle la subroutine principale inverse_problem_main.f90. A la compilation cette subroutine "voit" tous les modules et subroutines du package specfem3d_Git_devel en plus des subroutines specifiques au probleme inverse. Cela permettra d'appeler Specfem depuis ce code, plusieurs fois en mode forward et adjoint et au cours de chaque iteration. Il faut faire attention a bien reinitialiser les variables et tableaux entre chaque appel. Je suis en train d'ecrire cette routine de reinitialisation.

Pour les I/O j'ai suppose que les donnees seront rangees dans un repertoire par source. Pour l'instant j'ai un fichier par composante, il faut qu'on discute d'un format plus adapte. Aussi, dites-moi quelles sont les infos importantes ; je considere une fenetre d'inversion et un poids par trace, il faudra que l'utilisateur pointe la fenetre et donne une qualite sur la trace. Y a t-il autre chose?

Pour ./adjoint_source, il n'y a pas de difficultes, c'est simplement les residus + filtrage.

./inversion_scheme ce sera L-BFGS + Wolfe. Ensuite on pourra tester d'autres regles de recherche lineaire si besoin. Il suffira de rajouter le fichier dans ce repertoire.

Reste ./regularization, j'ai fait des tests pour definir les matrices de Vandermonde qui permettent de calculer les derivees par DF sur un ensemble de point quelconque. J'ai des routines Matlab et langage C qui fonctionnent (langage C car je m'en servais a l'origine pour un autre projet). Il me reste a traduire ca en Fortran. Le point le plus delicat est lorsqu'on se retrouve a l'interface entre deux tranches MPI, dans ce cas je dois communiquer l'ensemble de l'element voisin et non pas juste la face en contact. J'ai pas encore vu si c'est deja prevu dans Specfem.

D'ici 1 ou 2 semaines je pense avoir une premiere version qui tourne. Je construis aussi un exemple synthetique qui me permet de debugger, des que je code un truc je lance l'exemple avec le code compile en mode debug. A terme, soit je mettrai un exemple synthetique sous SVN pour avoir une illustration de l'utilisation du code soit on peut mettre un petit exemple avec des donnees reelles (deja publiees).

-----------------

Answer from Dimitri:

Pour l'exemple je suggere de mettre les deux : un exemple synthetique et un exemple de donnees reelles deja publiees. Plus y'a d'exemples mieux c'est, ca permet de tester les futures modifs du code aussi en verifiant que tous les exemples tournent toujours.

Tu peux meme mettre des exemples non publies sur SVN, il est entierement prive et y'a acces que pour toi, Seb et moi + Clement, personne d'autre, donc aucun souci, il est a nous !

