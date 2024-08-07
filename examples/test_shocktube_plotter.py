import os
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import pytools as pt
import numpy as np
matplotlib.rcParams.update({'font.size': 14})
matplotlib.rcParams['lines.linewidth'] = 2

runFolder = './'
runFiles = []
for f in sorted(os.listdir(runFolder)):
 if (f.startswith("state") and f.endswith(".vlsv")):
  runFiles.append(f)

vr = pt.vlsvfile.VlsvReader(runFolder + runFiles[0])
# check grid variables
# vr.get_all_variables()

# simulation box dimensions
[xmin,ymin,zmin,xmax,ymax,zmax] = vr.get_spatial_mesh_extent()
[mx,my,mz] = vr.get_spatial_mesh_size() # how many blocks per direction
[sx,sy,sz] = vr.get_spatial_block_size() # how many cells per block per direction
nx = mx*sx # number of cells along x
ny = my*sy # number of cells along y
nz = mz*sz # number of cells along z
dx = (xmax-xmin)/nx # should be dx = dy = dz in rhybrid

# points along x-axis
NN = 1000
xp = np.linspace(xmin+(1*dx),xmax-(1*dx),NN)
yp = np.zeros(len(xp))
zp = np.zeros(len(xp))
points=(np.array([xp,yp,zp])).transpose()
varList = ['cellB','n_H+','v_H+']
Ny = 3
Nx = 1

for f in runFiles:
 # read file
 vr = pt.vlsvfile.VlsvReader(runFolder + f)
 # do interpolation
 [crd,cellids,outVars,header]=pt.calculations.vlsv_intpol_points(vr,points,varList)
 #x = crd[:,0] # should be the same as xp
 #y = crd[:,1] # should be the same as yp
 #z = crd[:,2] # should be the same as zp

 # interpolated variables in points
 n = 0
 Bx     = outVars[:,n]; n = n + 1;
 By     = outVars[:,n]; n = n + 1;
 Bz     = outVars[:,n]; n = n + 1;
 nHsw   = outVars[:,n]; n = n + 1;
 VxHsw  = outVars[:,n]; n = n + 1;
 VyHsw  = outVars[:,n]; n = n + 1;
 VzHsw  = outVars[:,n]; n = n + 1;
 # derived
 Btot = np.sqrt(Bx**2 + By**2 + Bz**2)
 VtotHsw = np.sqrt(VxHsw**2 + VyHsw**2 + VzHsw**2)

 # plotters
 Nsub = 1
 plt.figure(figsize=(12,18))

 plt.subplot(Ny,Nx,Nsub); Nsub = Nsub + 1;
 plt.plot(xp,Bx/1e-9,'-r',label='Bx')
 plt.plot(xp,By/1e-9,'-g',label='By')
 plt.plot(xp,Bz/1e-9,'-b',label='Bz')
 plt.plot(xp,Btot/1e-9,'-k',label='Btot')
 plt.ylabel('B [nT]')
 plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
 plt.grid()

 plt.subplot(Ny,Nx,Nsub); Nsub = Nsub + 1;
 plt.plot(xp,nHsw/1e6,'-b',label='H+sw')
 plt.ylabel('n [cm-3]')
 plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
 plt.grid()
 
 plt.subplot(Ny,Nx,Nsub); Nsub = Nsub + 1;
 plt.plot(xp,VxHsw/1e3,'-r',label='Vx')
 plt.plot(xp,VyHsw/1e3,'-g',label='Vy')
 plt.plot(xp,VzHsw/1e3,'-b',label='Vz')
 plt.plot(xp,VtotHsw/1e3,'-k',label='Vtot')
 plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
 plt.ylabel('v(H+sw) [km/s]')
 plt.grid()
 
 plt.xlabel('x [m]')

 #plt.show()
 plt.savefig(f + '.png',bbox_inches='tight')
 plt.close()

