#!/usr/bin/env python
# coding: utf-8

# In[1]:


# PURPOSE:
#   This script creates movies from volume files using xcombine_vol_data_vtu 
#   command and corresponding ParaView PVD file for animation.
# DEVELOPER:
#   Hom Nath Gharti, Queen's University.
# RUN:
#   Simply press the 'Run' button or press "Shift+Enter". 
#   This can also be run from the terminal using ipython.
#     - Type and enter 'ipython' in the terminal.
#     - Type and enter 'run create_movie.ipynb'.
#   If necessary, this code can be converted to normal Python script:
#     - Type and enter the command: jupyter nbconvert --to python create_movie.ipynb
# HISTORY:
#   May 25, 2022: First created.
#------------------------------------------------------------------------------
# Import modules.
import subprocess

# Define Input Parameters:
#------------------------------------------------------------------------------
# Input directory.
input_dir = 'OUTPUT_FILES/DATABASES_MPI'
# Output directory.
output_dir = 'movie'
# Movie resolution: 0 = LOW, 1 = HIGH.
res = 1
# VTU file header.
file_head = 'div'
# Start time step.
itstart = 200
# End time step (NTSTEP).
itend = 16000
# Time step increment (DT).
itstep = 200
# Number of processors.
nproc = 40
#------------------------------------------------------------------------------

# Processor indices.
proc_start = 0
proc_end = nproc - 1


# In[2]:


# Set PVD file name and Open to write.
pvd_file = output_dir + '/' + file_head + '.pvd'
pvd = open(pvd_file,'w')
# Write PVD preheader.
pvd.write('<?xml version="1.0"?>\n')
pvd.write('<VTKFile type="Collection" version="0.1" byte_order="LittleEndian" compressor="vtkZLibDataCompressor">\n')
pvd.write('  <Collection>\n')


# In[3]:


# Loop through the time steps.
# itend+1 to include last time step.
for it in range (itstart, itend+1, itstep):
    # File name for xcombine_vol_data_vtu command.
    file_name = file_head + '_it{:06d}'.format(it)  
    # Full xcombine_vol_data_vtu command.
    vtu_command = ['./bin/xcombine_vol_data_vtu',str(proc_start),\
                str(proc_end),file_name,input_dir,\
                output_dir,str(res)]    
    # subprocess.run takes the array as an argument. 
    # The first element is the actual command and rest are the arguements to the command.
    subprocess.run(vtu_command)
    vtu_file = file_head + '_it{:06d}.vtu'.format(it)
    pvd.write('    <DataSet timestep="{:d}" group="" part="0" file="'.format(it) + vtu_file + '"/>\n')


# In[7]:


# Write PVD postfooter.
pvd.write('  </Collection>\n')
pvd.write('</VTKFile>\n')
# Close PVD file.
pvd.close()


# In[ ]:




