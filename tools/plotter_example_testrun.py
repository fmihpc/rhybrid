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

Rp = 3390e3 # radius of Mars

# simulation time
fileTime = vr.read_parameter("t")
if fileTime is None:
 fileTime = vr.read_parameter("time")

# read and sort cell ids
cellids_sorted = vr.read_variable("CellID").argsort()

# read a vector var: magnetic field
#B_i = vr.read_variable_info("cellBAverage")
#Btot = np.sqrt(B_i.data[:,0]**2 + B_i.data[:,1]**2 + B_i.data[:,2]**2)

# read a scalar var: sw proton number density
D_i = vr.read_variable_info("n_H+sw_ave")
D = D_i.data

# function to choose which plane to plot from 3D array (nx,ny,nz)
def chooseNcol(rPlane,rmin,rmax,nr,rstr):
 if (rPlane > rmax) or (rPlane < rmin):
  rPlane = (rmax + rmin)/2.0
  print(HN + "(chooseNcol) WARNING: " + rstr + "Plane coordinate out of domain, setting: " + rstr + "Plane = " + str(rPlane))
 dr = (rmax-rmin)/nr
 Ncol = int(np.floor((rPlane - rmin)/dr))
 # indices start from zero
 if Ncol >= nr:
  Ncol = nr-1
 if Ncol < 0:
  Ncol = 0
 return Ncol,rPlane

# incices of planes to be plotted
YZ_Ncol,xPlane = chooseNcol(-2.0*Rp,xmin,xmax,nx,"x") # x = -2*Rp plane (yz)
XZ_Ncol,yPlane = chooseNcol(0.0,ymin,ymax,ny,"y") # y = 0 plane (xz)
XY_Ncol,zPlane = chooseNcol(0.0,zmin,zmax,nz,"z") # z = 0 plane (xy)

# reshape variable into 3D array
D = D[cellids_sorted].reshape(nz,ny,nx)

# var unit cm-3
D = D/1e6

# var limits
dmin = min(D.flatten())
dmax = max(D.flatten())

# get 2D slices to plot from 3D data cube
meshD_yz = D[:,:,YZ_Ncol]
meshD_xz = D[:,XZ_Ncol,:]
meshD_xy = D[XY_Ncol,:,:]

# plot and save figure

plt.figure(figsize=(19,10))

plt.subplot(1,3,1)
plt.imshow(meshD_xz,vmin=dmin,vmax=dmax,extent=[xmin/Rp,xmax/Rp,zmin/Rp,zmax/Rp],aspect="equal",origin="lower",interpolation="nearest")
plt.xlabel("x [Rp]")
plt.ylabel("z [Rp]")
plt.title("xz (y = " + str(yPlane/Rp) + " Rp) plane")

plt.subplot(1,3,2)
plt.imshow(meshD_xy,vmin=dmin,vmax=dmax,extent=[xmin/Rp,xmax/Rp,ymin/Rp,ymax/Rp],aspect="equal",origin="lower",interpolation="nearest")
plt.xlabel("x [Rp]")
plt.ylabel("y [Rp]")
plt.title("xy (z = " + str(zPlane/Rp) + " Rp) plane")

plt.subplot(1,3,3)
plt.imshow(meshD_yz,vmin=dmin,vmax=dmax,extent=[ymin/Rp,ymax/Rp,zmin/Rp,zmax/Rp],aspect="equal",origin="lower",interpolation="nearest")
plt.xlabel("y [Rp]")
plt.ylabel("z [Rp]")
plt.title("yz (x = " + str(xPlane/Rp) + " Rp) plane")

clb = plt.colorbar()
clb.ax.set_title("n [cm-3]")

plt.savefig("testrun.png")
plt.clf()
plt.close()

