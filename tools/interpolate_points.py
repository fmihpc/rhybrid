#exec(open("./interpolate_points.py").read())
import sys
import os
import socket
import analysator as alr
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm, colors

matplotlib.rcParams.update({'font.size': 14})
matplotlib.rcParams['lines.linewidth'] = 2
#Rp = 2439.7e3 # Mercury
Rp = 3390e3 # Mars
intpolOrder = 1 # interpolation order


fn = os.getenv("HOME") + "/bin/corsair/testrun/state00004000.vlsv"

# read file
vr = alr.vlsvfile.VlsvReader(fn)
# simulation box dimensions
[xmin,ymin,zmin,xmax,ymax,zmax] = vr.get_spatial_mesh_extent()
[mx,my,mz] = vr.get_spatial_mesh_size() # how many blocks per direction
[sx,sy,sz] = vr.get_spatial_block_size() # how many cells per block per direction
nx = mx*sx # number of cells along x
ny = my*sy # number of cells along y
nz = mz*sz # number of cells along z
dx = (xmax-xmin)/nx # should be dx = dy = dz in rhybrid

# circular orbit around the planet at r = 1.5*Rp
NN = 1000
phi = np.linspace(0,2*np.pi,NN)
xp = 1.5*Rp*np.cos(phi)
yp = 1.5*Rp*np.sin(phi)
zp = np.zeros(len(xp))

#xp = np.linspace(xmin+(1*dx),xmax-(1*dx),NN) # outermost blocks are ghost cells
#yp = np.zeros(len(xp))
#zp = np.zeros(len(xp))
# artifial time parameter along trajectory
tp = np.linspace(0,1,NN)
points = (np.array([xp,yp,zp])).transpose()
varlist = ('cellBAverage','n_H+sw_ave')

# interpolate
[crd,cellids,outVars,header] = alr.calculations.vlsv_intpol_points(vr,points,varlist,interpolation_order=intpolOrder)

# extract variables
n = 0
Bx     = outVars[:,n]; n = n + 1;
By     = outVars[:,n]; n = n + 1;
Bz     = outVars[:,n]; n = n + 1;
nHsw   = outVars[:,n]; n = n + 1;
Btot = np.sqrt(Bx**2 + By**2 + Bz**2)

# plot
figNy = 2
figNx = 1
figNsub = 1

plt.figure(figsize=(12,18))

plt.subplot(figNy,figNx,figNsub); figNsub = figNsub + 1;
plt.plot(tp,Bx/1e-9,'-r',label='Bx')
plt.plot(tp,By/1e-9,'-g',label='By')
plt.plot(tp,Bz/1e-9,'-b',label='Bz')
plt.plot(tp,Btot/1e-9,'-k',label='Btot')
plt.ylabel('B [nT]')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.grid()

plt.subplot(figNy,figNx,figNsub); figNsub = figNsub + 1;
plt.plot(tp,nHsw/1e6,'-b',label='H+sw')
plt.ylabel('n [cm-3]')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.grid()

# plot orbit in 3d
ax = plt.figure(figsize=(12,18)).add_subplot(projection='3d')

# planet (r = Rp)
u,v = np.meshgrid(np.linspace(0,np.pi,50),np.linspace(0,2*np.pi,50))
xs = 1 * np.sin(u) * np.cos(v)
ys = 1 * np.sin(u) * np.sin(v)
zs = 1 * np.cos(u)
cs = np.ones(xs.shape)
cs[xs > 0] = 0
norm=colors.Normalize(vmin=0,vmax=1)
ax.plot_surface(xs,ys,zs,cmap=cm.Greys,facecolors=cm.Greys(norm(cs)),shade=True)

# particle trajectory
ax.plot(xp/Rp,yp/Rp,zp/Rp)
ax.axes.set_xlim3d(xmin=xmin/Rp,xmax=xmax/Rp)
ax.axes.set_ylim3d(ymin=ymin/Rp,ymax=ymax/Rp)
ax.axes.set_zlim3d(zmin=zmin/Rp,zmax=zmax/Rp)
ax.set_xlabel('x [Rp]')
ax.set_ylabel('y [Rp]')
ax.set_zlabel('z [Rp]')

ax.set_box_aspect([1,1,1])

plt.show(block=True)


