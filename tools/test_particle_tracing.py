# exec(open("test_particle_tracing.py").read())
# propagete test particles in E and B fields from VLSV files
import os
import pytools as pt
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from bbtracer import loadRun, bbstep

plt.rcParams.update({'font.size': 20})

# radius of planet / inner boundary
Rp = 2439.7e3

# test particle mass and charge
m = 1.672621716e-27
q = 1.60217653e-19
# timestep
dt = 0.5
Ndt = 20

#Ue = np.array((-400e3,0,0))
#B = np.array((0,1e-9,0))
#E = -np.cross(Ue,B)

varList = ['cellB','cellUe'] # variables to interpolate
folderRun = os.path.join(os.getenv("HOME"),'bin/corsair/testrun/'); # run folder
runName = 'testrun01' # run name
intpolOrder = 1 # interpolation order

R = loadRun(folderRun,runName,varList)

# initial particle positionsm velocities and time
rpart = np.array([(R['xmax']-2*R['dx'],R['ymax']-5*R['dx'],0)])
vpart = np.array([(0,0,0)])
tpart = 0

if rpart.shape != vpart.shape:
 print("ERROR: rp and vp should be same size")
 exit

# initialize test particle result arrays
r_res = rpart
v_res = vpart
t_res = np.array((tpart))

Nreader = 10
vr = R['sim']['readers'][Nreader]

# propagate test particles
for ii in range(Ndt):
 vrout = pt.calculations.vlsv_intpol_points(vr,rpart,varList)
 B = vrout[2][0][0:3]
 Ue = vrout[2][0][3:6]
 E = -np.cross(Ue,B)
 (rpart,vpart) = bbstep(E,B,rpart,vpart,m,q,dt)
 # boundary conditions
 # attach particle to result arrays
 r_res = np.vstack((r_res,rpart))
 v_res = np.vstack((v_res,vpart))
 # increase time
 tpart += dt
 t_res = np.vstack((t_res,tpart))

# 1D plot
plt.figure()
plt.subplot(6,1,1);
plt.plot(t_res,r_res[:,0]/Rp)
plt.ylabel('x [Rp]')
plt.subplot(6,1,2);
plt.plot(t_res,r_res[:,1]/Rp)
plt.ylabel('y [Rp]')
plt.subplot(6,1,3);
plt.plot(t_res,r_res[:,2]/Rp)
plt.ylabel('z [Rp]')
plt.subplot(6,1,4);
plt.plot(t_res,v_res[:,0]/1e3)
plt.ylabel('vx [km/s]')
plt.subplot(6,1,5);
plt.plot(t_res,v_res[:,1]/1e3)
plt.ylabel('vy [km/s]')
plt.subplot(6,1,6);
plt.plot(t_res,v_res[:,2]/1e3)
plt.ylabel('vz [km/s]')
plt.xlabel('t [s]')

# 3d plot
ax = plt.figure().add_subplot(projection='3d')

# planet
us = np.linspace(0, 2 * np.pi, 100)
vs = np.linspace(0, np.pi, 100)
xs = 1 * np.outer(np.cos(us), np.sin(vs))
ys = 1 * np.outer(np.sin(us), np.sin(vs))
zs = 1 * np.outer(np.ones(np.size(us)), np.cos(vs))
ax.plot_surface(xs, ys, zs)
#ax.set_aspect('equal')
ax.set_box_aspect([1,1,1])

# particle trajectory
ax.plot(r_res[:,0]/Rp,r_res[:,1]/Rp,r_res[:,2]/Rp)
ax.axes.set_xlim3d(xmin=R['xmin']/Rp,xmax=R['xmax']/Rp)
ax.axes.set_ylim3d(ymin=R['ymin']/Rp,ymax=R['ymax']/Rp)
ax.axes.set_zlim3d(zmin=R['zmin']/Rp,zmax=R['zmax']/Rp)
ax.set_xlabel('x [Rp]')
ax.set_ylabel('y [Rp]')
ax.set_zlabel('z [Rp]')


plt.show(block=True)
