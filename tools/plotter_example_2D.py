import matplotlib
import matplotlib.pyplot as plt
import analysator as alr
import numpy as np

vr = alr.vlsvfile.VlsvReader("state00004000.vlsv")
# simulation box dimensions
[xmin,ymin,zmin,xmax,ymax,zmax] = vr.get_spatial_mesh_extent()
[mx,my,mz] = vr.get_spatial_mesh_size() # how many blocks per direction
[sx,sy,sz] = vr.get_spatial_block_size() # how many cells per block per direction
nx = mx*sx # number of cells along x
ny = my*sy # number of cells along y
nz = mz*sz # number of cells along z
dx = (xmax-xmin)/nx # should be dx = dy = dz in rhybrid

Rp = 2439.7e3 # radius of Mercury
colormapName = "viridis"
#colormapName = "inferno"
#colormapName = "seismic"
#colormapName = "plasma"

# simulation time
fileTime = vr.read_parameter("t")
if fileTime is None:
 fileTime = vr.read_parameter("time")

# read and sort cell ids
cellids_sorted = vr.read_variable("CellID").argsort()

# density
n_i = vr.read_variable_info("n_H+sw")
n = n_i.data[cellids_sorted].reshape(nz,ny,nx)/1e6 # cm-3
n = n[0,:,:] # xy
#nmin = min(n.flatten())
#nmax = max(n.flatten())
nmin = 1
nmax = 2e2

# temperature
T_i = vr.read_variable_info("T_H+sw")
T = T_i.data[cellids_sorted].reshape(nz,ny,nx)/1e6 # MK
T = T[0,:,:] # xy
#Tmin = min(T.flatten())
#Tmax = max(T.flatten())
Tmin = 0.01
Tmax = 10

# speed
V_i = vr.read_variable_info("v_H+sw")
V = np.sqrt(V_i.data[:,0]**2 + V_i.data[:,1]**2 + V_i.data[:,2]**2)/1e3 # Vtot in km/s
V = V[cellids_sorted].reshape(nz,ny,nx)
V = V[0,:,:] # xy
#Vmin = min(V.flatten())
#Vmax = max(V.flatten())
Vmin = 0
Vmax = 500

# magnetic field
B_i = vr.read_variable_info("cellB")
B = np.sqrt(B_i.data[:,0]**2 + B_i.data[:,1]**2 + B_i.data[:,2]**2)/1e-9 # Btot in nT
B = B[cellids_sorted].reshape(nz,ny,nx)
B = B[0,:,:] # xy
#Bmin = min(B.flatten())
#Bmax = max(B.flatten())
Bmin = 1
Bmax = 100

# plot and save figure

plt.figure(figsize=(20,5))

plt.subplot(1,4,1)
plt.imshow(n,vmin=nmin,vmax=nmax,extent=[xmin/Rp,xmax/Rp,ymin/Rp,ymax/Rp],aspect="equal",origin="lower",interpolation="nearest")
plt.set_cmap(colormapName)
plt.xlabel("x [Rp]")
plt.ylabel("y [Rp]")
plt.title("n(H+sw) [cm-3]")
plt.colorbar(fraction=0.046,pad=0.04)

plt.subplot(1,4,2)
plt.imshow(T,vmin=Tmin,vmax=Tmax,extent=[xmin/Rp,xmax/Rp,ymin/Rp,ymax/Rp],aspect="equal",origin="lower",interpolation="nearest")
plt.set_cmap(colormapName)
plt.xlabel("x [Rp]")
#plt.ylabel("y [Rp]")
plt.title("T(H+sw) [MK]")
plt.colorbar(fraction=0.046,pad=0.04)

plt.subplot(1,4,3)
plt.imshow(V,vmin=Vmin,vmax=Vmax,extent=[xmin/Rp,xmax/Rp,ymin/Rp,ymax/Rp],aspect="equal",origin="lower",interpolation="nearest")
plt.set_cmap(colormapName)
plt.xlabel("x [Rp]")
#plt.ylabel("y [Rp]")
plt.title("V(H+sw) [km/s]")
plt.colorbar(fraction=0.046,pad=0.04)

plt.subplot(1,4,4)
plt.imshow(B,vmin=Bmin,vmax=Bmax,extent=[xmin/Rp,xmax/Rp,ymin/Rp,ymax/Rp],aspect="equal",origin="lower",interpolation="nearest")
plt.set_cmap(colormapName)
plt.xlabel("x [Rp]")
#plt.ylabel("y [Rp]")
plt.title("B [nT]")
plt.colorbar(fraction=0.046,pad=0.04)

plt.savefig("example_planet_2D.png")
plt.clf()
plt.close()

