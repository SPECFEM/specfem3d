################################################################
#                                                              #
#            define meshfem3D inputs files                     #
#                                                              #
#                                                              #
################################################################



#### DOMAIN ######

size_x_in_km=32
size_y_in_km=32
size_z_in_km=15

#### DISCRETIZATION ######

size_element_x_in_m=1000
size_element_y_in_m=1000
size_element_z_in_m=1000


#### BACKGROUND FK MODEL #####

number_of_layers=2
heigth_of_layer_1_km=10

# the last layer is homgeous half-space thus no heigth is given

rho_1=2000
vp_1=5000
vs_1=3500

rho_2=3600
vp_2=7000
vs_2=4500


######## ADDITIONAL CARTESIAN BOX PERTURBATION INSIDE MODEL ##########

rho_pert=3000
vp_pert=6000
vs_pert=4000

xmin_pert_km=14
xmax_pert_km=18
ymin_pert_km=14
ymax_pert_km=18
zmin_pert_km=3
zmax_pert_km=7


############ TODO CODE TO FILL Mesh_Par_file #######
