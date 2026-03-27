# test particle tracing
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

# parse input arguments
parser = argparse.ArgumentParser('trace_test_particles.py')
parser.add_argument('input_file',help='Input VLSV file',type=Path)
parser.add_argument('var_name_B',help='Variable name for magnetic field',type=str)
parser.add_argument('var_name_Ue',help='Variable name for electron velocity',type=str)
parser.add_argument('var_name_Ep',help='Variable name for electron pressure term of electric field',type=str)
parser.add_argument('Rp',help='Planet radius',type=float)
args = parser.parse_args()
input_file = str(args.input_file) #input_file = './state00004000.vlsv'
var_name_B = args.var_name_B # var_name_Ue = 'cellB'
var_name_Ue = args.var_name_Ue # var_name_Ue = 'cellUe'
var_name_Ep = args.var_name_Ep # var_name_Ue = 'cellEp'
Rp = args.Rp #Rp = 3390e3
if Path(input_file).is_file() == False:
 print('ERROR: input file does not exist (' + input_file + ')')
 sys.exit()

save_figures = True # save figures as png files

Robstacle = Rp + 200e3 # obstacle radius (particle absorbed below this)
m = 1.672621716e-27 # test particle mass
q = 1.60217653e-19 # test particle charge
linear_field_interpolation = True # linear or nearest interpolation
# timestep [s]
dt = 0.5
# max. number of timesteps
max_timesteps = 500
# constant random seed for reproducibility
np.random.seed(71)

vr = alr.vlsvfile.VlsvReader(input_file)
if vr.check_variable(var_name_B) == False:
 print('ERROR: variable name for B (' + var_name_B + ' / ' + input_file + ')')
 exit()
if vr.check_variable(var_name_Ue) == False:
 print('ERROR: variable name for Ue (' + var_name_Ue + ' / ' + input_file + ')')
 exit()
if vr.check_variable(var_name_Ep) == False:
 print('ERROR: variable name for Ep (' + var_name_Ep + ' / ' + input_file + ')')
 exit()
# simulation box dimensions
[xmin,ymin,zmin,xmax,ymax,zmax] = vr.get_spatial_mesh_extent()
[mx,my,mz] = vr.get_spatial_mesh_size() # how many blocks per direction
[sx,sy,sz] = vr.get_spatial_block_size() # how many cells per block per direction
nx = mx*sx # number of cells along x
ny = my*sy # number of cells along y
nz = mz*sz # number of cells along z
dx = (xmax-xmin)/nx # should be dx = dy = dz in rhybrid
box_eps = 0.1*dx # buffer distance for outer boundary conditions [dx]

# squared obstacle radius
Robstacle2 = Robstacle**2

# INITIALISE TEST PARTICLES

Nparticles = 10

# positions: upstream along a line in y at x = const. and z = 0
x0 = np.ones(Nparticles)*xmax-1.5*dx
y0 = np.linspace(ymin+1.5*dx,ymax-1.5*dx,Nparticles)
z0 = np.zeros(Nparticles)
r0 = np.column_stack((x0,y0,z0))

# velocities [m/s]: bulk velocity + random thermal velocity
Usw = (-430e3,0,0)
vth = 30e3
v0 = np.ones(r0.shape)*Usw + vth*np.random.rand(r0.shape[0],r0.shape[1])

# initial time [s]
t0 = np.zeros(Nparticles)

if (r0.shape != v0.shape) | (len(r0) != len(t0)) | (len(v0) != len(t0)):
 print('ERROR: r0 and v0 should be same shape and same length as t0')
 exit()

# particle ids for indexing results lists
pid = np.arange(Nparticles)

# result lists
r = np.split(r0,Nparticles)
v = np.split(v0,Nparticles)
t = np.split(t0,Nparticles)

# particle positions, velocities and time at a timestep
rnow = np.copy(r0)
vnow = np.copy(v0)
tnow = np.copy(t0)

# Boris-Bunemann algorithm: move and accelerate particles one timestep
def bb_step(E,B,r,v,m,q,dt):
 # move particle
 r = r + dt*v
 vold = v
 #accelerate particle
 qmideltT2= 0.5*q*dt/m
 dv = qmideltT2*E
 t = qmideltT2*B
 t2 = np.sum(t*t,axis=1)
 b2 = 2.0/(1.0 + t2)
 s = b2[:,None]*t
 vm = v + dv
 v0 = vm + np.cross(vm,t)
 vp = vm + np.cross(v0,s)
 v = vp + dv
 return (r,v)

# propagation loop
for ii in range(max_timesteps):
 # check outer boundary conditions
 outer_boundary_condition = (
 (rnow[:,0] < (xmin + box_eps)) |
 (rnow[:,0] > (xmax - box_eps)) |
 (rnow[:,1] < (ymin + box_eps)) |
 (rnow[:,1] > (ymax - box_eps)) |
 (rnow[:,2] < (zmin + box_eps)) |
 (rnow[:,2] > (zmax - box_eps))
 )
 # check inner boundary condition
 r2 = np.sum(rnow*rnow,axis=1)
 inner_boundary_condition = (r2 <= Robstacle2)
 # full boundary conditions
 all_boundary_conditions = (outer_boundary_condition | inner_boundary_condition)
 # remove particles that entered boundaries from propagation
 if any(all_boundary_conditions) == True:
  rnow = rnow[np.where(all_boundary_conditions == False)]
  vnow = vnow[np.where(all_boundary_conditions == False)]
  pid = pid[np.where(all_boundary_conditions == False)]
 # break loop if no more particles
 if len(pid) < 1:
  break
 # interpolate Ue/E and B in particle positions
 var_list = [var_name_B,var_name_Ue,var_name_Ep]
 #varList = [var_name_B,var_name_Ue]
 [crds,cellids,fields,headerstr] = alr.calculations.vlsv_intpol_points(vr,rnow,var_list,interpolation_order=linear_field_interpolation)
 B = fields[:,0:3]
 Ue = fields[:,3:6]
 Ep = fields[:,6:9]
 # E = -Ue x B + Ep
 E = -np.cross(Ue,B) + Ep
 # E = -Ue x B
 #E = -np.cross(Ue,B)

 # Boris-Bunemann step
 (rnow,vnow) = bb_step(E,B,rnow,vnow,m,q,dt)
 # increase time
 tnow += dt
 # split new particle positions and velocities into list of arrays
 rnow_list_split = np.split(rnow,len(rnow))
 vnow_list_split = np.split(vnow,len(vnow))
 # append positions, velocities and times for particles that are still propagated
 for jj in range(len(pid)):
  kk = pid[jj]
  r[kk] = np.vstack((r[kk],rnow_list_split[jj]))
  v[kk] = np.vstack((v[kk],vnow_list_split[jj]))
  t[kk] = np.append(t[kk],tnow[jj])

 # TBD: faster vectorized version
 #r_arr = np.stack(r)
 #rnow_list_split_arr = np.stack(rnow_list_split)
 #r_arr = np.concatenate((r_arr,rnow_list_split_arr),axis=1)
 #r = [a for a in r_arr]

# number of timesteps taken (=number of bb_step/vlsv_intpol_points function calls)
dts_taken = ii + 1
# total number of particle propagations (initial state not counted)
bb_steps_total = sum(tp.size for tp in t) - Nparticles
print('')
print('tracing finished')
print('timesteps taken: ' + str(dts_taken))
print('simulation time: ' + str(dts_taken*dt) + ' s')
print('total number of particle propagations/field interpolations: ' + str(bb_steps_total))

# plotting

Nx_sub = 1
Ny_sub = 6

# time series
plt.figure(figsize=(3*Ny_sub/1.5,3*Ny_sub))
for ii in range(len(r)):
 plt.subplot(Ny_sub,Nx_sub,1)
 plt.plot(t[ii],r[ii][:,0]/Rp)
 plt.ylabel('x')
 plt.grid()
 plt.subplot(Ny_sub,Nx_sub,2)
 plt.plot(t[ii],r[ii][:,1]/Rp)
 plt.ylabel('y')
 plt.grid()
 plt.subplot(Ny_sub,Nx_sub,3)
 plt.plot(t[ii],r[ii][:,2]/Rp)
 plt.ylabel('z')
 plt.grid()
 plt.subplot(Ny_sub,Nx_sub,4)
 plt.plot(t[ii],v[ii][:,0]/1e3)
 plt.ylabel('vx')
 plt.grid()
 plt.subplot(Ny_sub,Nx_sub,5)
 plt.plot(t[ii],v[ii][:,1]/1e3)
 plt.ylabel('vy')
 plt.grid()
 plt.subplot(Ny_sub,Nx_sub,6)
 plt.plot(t[ii],v[ii][:,2]/1e3)
 plt.ylabel('vz')
 plt.grid()
 plt.xlabel('t [s]')

plt.subplot(Ny_sub,Nx_sub,1)
plt.title('variables along test particle trajectories')

if save_figures == True:
 plt.savefig('trace_test_particles_time_series.png')
 plt.clf()
 plt.close()

# particle trajectories in 3D

ax = plt.figure(figsize=(10,10)).add_subplot(projection='3d')

# planet
u,v = np.meshgrid(np.linspace(0,np.pi,50),np.linspace(0,2*np.pi,50))
xs = 1 * np.sin(u) * np.cos(v)
ys = 1 * np.sin(u) * np.sin(v)
zs = 1 * np.cos(u)
cs = np.ones(xs.shape)
cs[xs > 0] = 0
norm = colors.Normalize(vmin=0,vmax=1)
ax.plot_surface(xs,ys,zs,cmap=cm.Greys,facecolors=cm.Greys(norm(cs)),shade=True)
ax.set_xlim(xmin/Rp,xmax/Rp)
ax.set_ylim(ymin/Rp,ymax/Rp)
ax.set_zlim(zmin/Rp,zmax/Rp)
ax.set_box_aspect([1,1,1])
#ax.set_aspect('equal')

# trajectories
for ii in range(len(r)):
 ax.plot(r[ii][:,0]/Rp,r[ii][:,1]/Rp,r[ii][:,2]/Rp)

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
ax.set_title('test particle trajectories')

if save_figures == True:
 plt.savefig('trace_test_particles_trajectories.png')
 plt.clf()
 plt.close()
else:
 plt.show(block=True)

