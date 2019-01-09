#from scipy.io import netcdf
import numpy as np
import sys
import os
import struct
import scipy.ndimage.filters as filters
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from mpl_toolkits.mplot3d import axes3d



#=====================================
#    Set up parameters             ===
#=====================================
# Please set up the mesh parametes in this section

# Work directory
Work_dir             = "/u/moana/user/weng/Weng/Scripts_dev/temp_output/"
# All the Specfem3D results are saved in the directory named as Model_name
Model_name           = "Kumamoto_semisphere_curvedfault_curvedtopo_4_8_HEX8"
# Grid size used to intepolate the fault nodes
grid_size            = 0.2  # in km
# display the fault horizontally or vertically
display_horizontally = True
#True
# if display_horizontally = False, then display vertically. Rotate=0 means projecting along X-Z plane, and Rotate=90 means projecting along Y-Z plane.
Rotate  =  90.0
# The S wave speed, used to normalized the rupture speed.
Vs = 3.33   # in km/s
# Number of countours of rupture front
ContourT = 5
# Number of countours of final slip
ContourS = 3

#=====================================




## Files check
file_dir  = Work_dir + Model_name
if not os.path.isdir(file_dir):
    print "The directory that contains the data doesn't exist..."
    exit()
file_list = os.listdir(file_dir)
snap_list = []
time_list = []
for name in file_list:
    if(name.find('Snapshot')==0):
        snap_list.append(name)
        time_list.append(int(name.split('Snapshot')[1].split('_')[0]))
num_data   = len(time_list)
time_list  = np.asarray(sorted(time_list))
print 
print "The model name is:", Model_name
print "The model has", num_data, "data files."
print "The time list is:", time_list
print 
if(display_horizontally):
    print "Project the fault horizontally."
else:
    print "Project the fault onto a vertical plane (clockwisely rotated ", Rotate, " degree about X-Z plane)."
print

#==================
#   Functions   ===
#==================
def  FSEM3D_snapshot(filename):
    data = type("",(),{})()
    NDAT    = 14
    length  = 4
    binary_file = open(filename,"rb")
    BinRead = []
    for ii in range(NDAT):
        read_buf = binary_file.read(4)    # the number of bytes of int is 4
        number = struct.unpack('1i',read_buf)[0]
        N = number/length
        read_buf = binary_file.read(number)
        read_d   = struct.unpack(str(N)+'f',read_buf)
        read_d   = np.array(read_d)
        BinRead.append(read_d)
        read_buf = binary_file.read(4)    # the number of bytes of int is 4
        number = struct.unpack('1i',read_buf)[0]
    data.X  = BinRead[0]/1e3   # in km
    data.Y  = BinRead[1]/1e3  # in km
    data.Z  = BinRead[2]/1e3  # in km
    data.Dx = BinRead[3]
    data.Dz = BinRead[4]
    data.Vx = BinRead[5]
    data.Vz = BinRead[6]
    data.Tx = BinRead[7]      # in MPa
    data.Ty = BinRead[8]
    data.Tz = BinRead[9]     # in MPa
    data.S  = BinRead[10]
    data.Sg = BinRead[11]     # in MPa
    data.Trup = BinRead[12]
    data.Tpz = BinRead[13]
    return data

# Load all the data
Data = []
for ii in range(num_data):
    name = file_dir+"/Snapshot"+str(time_list[ii])+"_F1.bin"
    Data.append(FSEM3D_snapshot(name))

# Create the grid arrays of nodes
X    = Data[0].X   
Y    = Data[0].Y
Z    = Data[0].Z   
# if it is a vertical fault (along X or Y axis)
if((np.min(X) == np.max(X) or np.min(Y)==np.max(Y)) and display_horizontally):
    print "Warning: it is a vertical fault along X or Y axis, the fault will be displayed vertically!!!"
    print "Please change display_horizontally to be False and setup the strike of fault to be 0 or 90."
    exit()

if(not display_horizontally):
    X_rot = X * np.cos(Rotate/180.0*np.pi)   + Y * np.sin(Rotate//180.0*np.pi)
    Y_rot = X * - np.sin(Rotate//180.0*np.pi) + Y * np.cos(Rotate//180.0*np.pi)
    X     = X_rot
    Y     = Y_rot

[X_lower,X_upper,Y_lower,Y_upper,Z_lower,Z_upper] = [np.min(X),np.max(X),np.min(Y),np.max(Y),np.min(Z),np.max(Z)]
X_dim = int((X_upper-X_lower)/grid_size+1)
Y_dim = int((Y_upper-Y_lower)/grid_size+1)
Z_dim = int((Z_upper-Z_lower)/grid_size+1)

if(display_horizontally):
    X_Y    = np.column_stack((X, Y))
    Relief = Z
    grid_x, grid_y = np.mgrid[X_lower:X_upper:X_dim*1j, Y_lower:Y_upper:Y_dim*1j]
    fault_range =  [X_lower,X_upper,Y_lower,Y_upper]
    fault_dims  = [X_dim,Y_dim]
else:
    X_Y    = np.column_stack((Y, Z))
    Relief = Y
    grid_x, grid_y = np.mgrid[Y_lower:Y_upper:Y_dim*1j, Z_lower:Z_upper:Z_dim*1j]
    fault_range =  [Y_lower,Y_upper,Z_lower,Z_upper]
    fault_dims  = [Y_dim,Z_dim]


# The relief of the axis that is normal to the projected surface
Z_grid         = griddata(X_Y[:,0:2], Relief, (grid_x,grid_y), method='linear')

X_gradient,Y_gradient  =  np.gradient(Z_grid)
Nor_dir = np.zeros((X_gradient.shape[0],X_gradient.shape[1],3))
for i in range(fault_dims[0]):
    for j in range(fault_dims[1]):
         Nor_dir[i,j,:] = np.cross([1,0,X_gradient[i,j]], [0,1,Y_gradient[i,j]])

# Rupture time
Final_rup_time = Data[len(Data)-1].Trup
Init_t0 = griddata(X_Y[:,0:2], Final_rup_time, (grid_x,grid_y), method='linear')

# Rupture speed and direction
vr      = np.zeros((Init_t0.shape))
vr_dir  = np.zeros((Init_t0.shape))
for x in range(fault_dims[0]):
    for y in range(fault_dims[1]):
        # Skip the boundaries and unbroken nodes
        if(x==0 or x==fault_dims[0]-1 or y==0 or y==fault_dims[1]-1 or np.isnan(Init_t0[x,y])):
            vr[x,y] = np.nan
            continue
        # Calculate the gradient of rupture time
        delta_x = (Init_t0[x+1,y+1]-Init_t0[x-1,y+1]+Init_t0[x+1,y-1]-Init_t0[x-1,y-1]) / 4.0 / ((grid_size*X_gradient[x,y])**2 + grid_size**2)**0.5 
        delta_y = (Init_t0[x+1,y+1]-Init_t0[x+1,y-1]+Init_t0[x-1,y+1]-Init_t0[x-1,y-1]) / 4.0 / ((grid_size*Y_gradient[x,y])**2 + grid_size**2)**0.5 
        # Calculate rupture speed and direction
        if (np.abs(delta_x)<1e-6 and np.abs(delta_x)<1e-6):
            vr[x,y] = np.nan
            vr_dir[x,y]  = np.nan
        else:
            vr[x,y] = 1 / (delta_x**2 + delta_y**2)**0.5 / Vs
            vr_dir[x,y]  =  delta_x / (delta_x**2 + delta_y**2)**0.5


# Final slip
Final_slip_str = Data[len(Data)-1].Dx 
Final_slip_dip = Data[len(Data)-1].Dz
Final_slip     = (Final_slip_str**2 + Final_slip_dip**2)**0.5
Slip_grid      = griddata(X_Y[:,0:2], Final_slip, (grid_x,grid_y), method='linear')


# Stress drop
Stress_drop_str = Data[0].Tx/1e6 - Data[len(Data)-1].Tx/1e6 
Stress_drop_dip = Data[0].Ty/1e6 - Data[len(Data)-1].Ty/1e6 
Stress_drop_nor = Data[len(Data)-1].Tz/1e6
str_grid        = griddata(X_Y[:,0:2], Stress_drop_str, (grid_x,grid_y), method='linear')
dip_grid        = griddata(X_Y[:,0:2], Stress_drop_dip, (grid_x,grid_y), method='linear')
nor_grid        = griddata(X_Y[:,0:2], Stress_drop_nor, (grid_x,grid_y), method='linear')


### Plot figure

# Rupture time
plt.subplot2grid((3,2), (0,0), colspan=1, rowspan=1)
im   = plt.imshow(Init_t0.T, extent=fault_range, origin='lower', cmap='hot_r', aspect='auto')
im.axes.get_xaxis().set_visible(False)
cs   = plt.contour(grid_x, grid_y, Init_t0, ContourT, colors='k', linewidths=0.1)
cbar = plt.colorbar(im, extend='both', shrink=0.5, ticks=[round(np.nanmin(Init_t0),1), round((np.nanmin(Init_t0)+np.nanmax(Init_t0))/2,1), round(np.nanmax(Init_t0),1)])
cbar.set_label('Time (s)')
plt.title('Rupture time')

## 3D view for fault
#ax = plt.add_subplot(111, projection='3d')
#
#fig = plt.subplot2grid((3,2), (0,0), colspan=1, rowspan=1, fig=ax)
#X, Y, Z = axes3d.get_test_data(0.05)
## Plot a basic wireframe.
#ax.plot_wireframe(X, Y, Z, rstride=10, cstride=10)
#
# Rupture speed
plt.subplot2grid((3,2), (1,0), colspan=1, rowspan=1)
im   = plt.imshow(vr.T, extent=fault_range, origin='lower', cmap='hot_r', vmin=0, vmax=1.732, aspect='auto')
im.axes.get_xaxis().set_visible(False)
cs   = plt.contour(grid_x, grid_y, Init_t0, ContourT, colors='k', linewidths=0.1)
cbar = plt.colorbar(im, extend='both', shrink=0.5, ticks=[0,0.5,1,1.5])
cbar.set_label('Vr/Vs')
plt.title('Rupture Speed')

# Final slip
plt.subplot2grid((3,2), (1,1), colspan=1, rowspan=1)
im   = plt.imshow(Slip_grid.T, extent=fault_range, origin='lower', cmap='Reds', aspect='auto')
im.axes.get_xaxis().set_visible(False)
im.axes.get_yaxis().set_visible(False)
cs   = plt.contour(grid_x, grid_y, Slip_grid, ContourS, colors='k', linewidths=0.1)
cbar = plt.colorbar(im, extend='both', shrink=0.5, ticks=[round(np.nanmin(Slip_grid),2), round((np.nanmin(Slip_grid)+np.nanmax(Slip_grid))/2,2), round(np.nanmax(Slip_grid),2)])
cbar.set_label('slip (m)')
plt.title('Final Slip')

# Stress drop at strike direction
plt.subplot2grid((3,2), (2,0), colspan=1, rowspan=1)
im   = plt.imshow(str_grid.T, extent=fault_range, origin='lower', cmap='Reds', aspect='auto')
cbar = plt.colorbar(im, extend='both', shrink=0.5, ticks=[round(np.nanmin(str_grid),1), round((np.nanmin(str_grid)+np.nanmax(str_grid))/2,1), round(np.nanmax(str_grid),1)])
cbar.set_label('stress drop (MPa)')
plt.title('Strike Stress Drop')

# Stress drop at dip direction
plt.subplot2grid((3,2), (2,1), colspan=1, rowspan=1)
im   = plt.imshow(dip_grid.T, extent=fault_range, origin='lower', cmap='Reds', aspect='auto')
im.axes.get_yaxis().set_visible(False)
cbar = plt.colorbar(im, extend='both', shrink=0.5, ticks=[round(np.nanmin(dip_grid),1), round((np.nanmin(dip_grid)+np.nanmax(dip_grid))/2,1), round(np.nanmax(dip_grid),1)])
cbar.set_label('stress drop (MPa)')
plt.title('Dip Stress Drop')

plt.savefig("ps/"+Model_name+"-results.pdf",format="pdf")

print "pdf file is saved."
###

### Save data
output= open("data/"+Model_name+"-results.dat","w")
for x in range(fault_dims[0]):
    for y in range(fault_dims[1]):
        output.writelines(str(grid_x[x,y]))
        output.writelines("  ")
        output.writelines(str(grid_y[x,y]))
        output.writelines("  ")
        output.writelines(str(Init_t0[x,y]))
        output.writelines("  ")
        output.writelines(str(vr[x,y]))
        output.writelines("  ")
        output.writelines(str(Slip_grid[x,y]))
        output.writelines("  ")
        output.writelines(str(str_grid[x,y]))
        output.writelines("  ")
        output.writelines(str(dip_grid[x,y]))
        output.writelines("\n")
output.close()

print "dat file is saved."
