# spatial interpolation along a spacecraft trajectory in arbitrary points
import os
import sys
from pathlib import Path
import argparse
try:
 import analysator as alr
except ModuleNotFoundError as err:
 print('Analysator not found: ' + str(err))
 sys.exit()
try:
 import numpy as np
except ModuleNotFoundError as err:
 print('NumPy not found: ' + str(err))
 sys.exit()
try:
 from matplotlib import cm, colors
 import matplotlib.pyplot as plt
except ModuleNotFoundError as err:
 print('matplotlib not found: ' + str(err))
 sys.exit()

def list_of_strings(arg):
 return arg.split(',')

# parse input arguments
parser = argparse.ArgumentParser('create_interpolation_along_trajectory.py')
parser.add_argument('input_file',help='Input VLSV file',type=Path)
parser.add_argument('var_list',help='List of variables: var1,var2,var3...',type=list_of_strings)
parser.add_argument('Rp',help='Planet radius',type=float)
args = parser.parse_args()
input_file = str(args.input_file) #input_file = './state00004000.vlsv'
var_list = args.var_list # var_list = ['cellB','n_H+sw_ave','v_H+sw_ave']
Rp = args.Rp #Rp = 3390e3
if Path(input_file).is_file() == False:
 print('ERROR: input file does not exist (' + input_file + ')')
 sys.exit()

save_figures = True # save figures as png files

linear_field_interpolation = True # linear or nearest interpolation

# read file
vr = alr.vlsvfile.VlsvReader(input_file)
# check variables are found
for var in var_list:
 if vr.check_variable(var) == False:
  print('ERROR: variable ' + var + ' not found in ' + input_file)
  print('Available variables: ')
  print(str(vr.get_all_variables()))
  sys.exit()
# simulation box dimensions
[xmin,ymin,zmin,xmax,ymax,zmax] = vr.get_spatial_mesh_extent()
[mx,my,mz] = vr.get_spatial_mesh_size() # how many blocks per direction
[sx,sy,sz] = vr.get_spatial_block_size() # how many cells per block per direction
nx = mx*sx # number of cells along x
ny = my*sy # number of cells along y
nz = mz*sz # number of cells along z
dx = (xmax-xmin)/nx # should be dx = dy = dz in rhybrid

# circular orbit around the planet at r = 1.5*Rp on z = 0 plane
Npoints = 100
phi = np.linspace(0,2*np.pi,Npoints)
xp = 1.5*Rp*np.cos(phi)
yp = 1.5*Rp*np.sin(phi)
zp = np.zeros(len(xp))

# artificial time parameter
tp = np.linspace(0,1,Npoints)

# point array for vlsv_intpol_points
points = (np.array([xp,yp,zp])).transpose()

# interpolate
[crd,cellids,vars_intpol,header] = alr.calculations.vlsv_intpol_points(vr,points,var_list,interpolation_order=linear_field_interpolation)

# plotting

plt.rcParams.update({'font.size': 16})
Nx_sub = 1
Ny_sub = len(var_list)
Nsub = 1
Ncols = vars_intpol.shape[1]
Nc = 0
plt.clf()
plt.close()
plt.figure(figsize=(15,10))
for var in var_list:
 plt.subplot(Ny_sub,Nx_sub,Nsub); Nsub = Nsub + 1;
 dim = vr.read_variable_vectorsize(var)
 if dim == 1:
  vp = vars_intpol[:,Nc]; Nc = Nc + 1
  plt.plot(tp,vp)
 elif dim == 3:
  vpx = vars_intpol[:,Nc]; Nc = Nc + 1
  vpy = vars_intpol[:,Nc]; Nc = Nc + 1
  vpz = vars_intpol[:,Nc]; Nc = Nc + 1
  vp_tot = np.sqrt(np.square(vpx) + np.square(vpy) + np.square(vpz))
  plt.plot(tp,vpx,'-r',label='x comp.')
  plt.plot(tp,vpy,'-g',label='y comp.')
  plt.plot(tp,vpz,'-b',label='z comp.')
  plt.plot(tp,vp_tot,'-k',label='magn.')
  plt.plot(tp,-vp_tot,'--k',label='-magn.')
  plt.legend(bbox_to_anchor=(0.9,1),loc=2,borderaxespad=0.)
 else:
  print('ERROR: only 1 and 3 dimensional variable supporte')
  sys.exit()
 plt.ylabel(var)
 plt.grid()
plt.xlabel('t')

if save_figures == True:
 plt.savefig('create_interpolation_along_trajectory_plot.png')
 plt.clf()
 plt.close()

# plot trajectory points in 3d

ax = plt.figure(figsize=(10,10)).add_subplot(projection='3d')

# planet
u,v = np.meshgrid(np.linspace(0,np.pi,50),np.linspace(0,2*np.pi,50))
xs = 1 * np.sin(u) * np.cos(v)
ys = 1 * np.sin(u) * np.sin(v)
zs = 1 * np.cos(u)
cs = np.ones(xs.shape)
cs[xs > 0] = 0
norm=colors.Normalize(vmin=0,vmax=1)
ax.plot_surface(xs,ys,zs,cmap=cm.Greys,facecolors=cm.Greys(norm(cs)),shade=True)

# interpolation trajectory
ax.plot(xp/Rp,yp/Rp,zp/Rp,'.b')
ax.axes.set_xlim3d(xmin=xmin/Rp,xmax=xmax/Rp)
ax.axes.set_ylim3d(ymin=ymin/Rp,ymax=ymax/Rp)
ax.axes.set_zlim3d(zmin=zmin/Rp,zmax=zmax/Rp)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

ax.set_box_aspect([1,1,1])

if save_figures == True:
 plt.savefig('create_interpolation_along_trajectory_3d.png')
 plt.clf()
 plt.close()
else:
 plt.show(block=True)



