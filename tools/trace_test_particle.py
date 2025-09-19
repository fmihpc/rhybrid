# exec(open("test_particle_tracing.py").read());
# propagete test particles in E and B fields from VLSV files
import os
import pytools as pt
import numpy as np
from matplotlib import cm, colors
import matplotlib.pyplot as plt
from bbtracer import loadRun, bbstep

plt.rcParams.update({'font.size': 20})

RMercury = 2439.7e3
mp = 1.672621716e-27
me = 9.1093825746e-31
qe = 1.60217653e-19

# radius of planet / inner boundary
Rp = RMercury
# test particle mass [kg] and charge [C]
m = mp
q = qe
# timestep [s]
dt = 0.5
# max. number of timesteps
Ndt = 200
# buffer distance for out boundary conditions [dx]
boxEps = 0.1

#Ue = np.array((-400e3,0,0))
#B = np.array((0,1e-9,0))
#E = -np.cross(Ue,B)

varList = ['cellB','cellUe'] # variables to interpolate
folderRun = os.path.join(os.getenv("HOME"),'bin/corsair/testrun/'); # run folder
runName = 'testrun01' # run name
intpolOrder = 1 # interpolation order

R = loadRun(folderRun,runName,varList)
boxEps *= R['dx']

# initial particle positions [m], velocities [m/s] and time [s]
rpart = np.array([(2.5*Rp,0,1.5*Rp)])
vpart = np.array([(-400e3,0,0)])
tpart = 0

if rpart.shape != vpart.shape:
 print("ERROR: rp and vp should be same size")
 exit

# initialize test particle result arrays
r_res = rpart
v_res = vpart
t_res = np.array((tpart))

# select 10th found VLSV file/timestep
Nreader = 10
vr = R['sim']['readers'][Nreader]

# propagate test particles
for ii in range(Ndt):
 # outer boundary conditions (stop if particle reaches box limits)
 if ((rpart[0][0] < (R['xmin'] + boxEps)) or
     (rpart[0][0] > (R['xmax'] - boxEps)) or
     (rpart[0][1] < (R['ymin'] + boxEps)) or
     (rpart[0][1] > (R['ymax'] - boxEps)) or
     (rpart[0][2] < (R['zmin'] + boxEps)) or
     (rpart[0][2] > (R['zmax'] - boxEps))):
  break
 # inner boundary condition (stop if particle reachers inner boundary)
 if np.linalg.norm(rpart) < Rp:
  break
 # interpolate Ue/E and B in particle position
 vrout = pt.calculations.vlsv_intpol_points(vr,rpart,varList)
 B = vrout[2][0][0:3]
 Ue = vrout[2][0][3:6]
 # E = -Ue x B
 E = -np.cross(Ue,B)
 # Boris-Bunemann step
 (rpart,vpart) = bbstep(E,B,rpart,vpart,m,q,dt)
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

# planet (r = Rp)
u,v = np.meshgrid(np.linspace(0,np.pi,50),np.linspace(0,2*np.pi,50))
xs = 1 * np.sin(u) * np.cos(v)
ys = 1 * np.sin(u) * np.sin(v)
zs = 1 * np.cos(u)
cs = np.ones(xs.shape)
cs[xs > 0] = 0
norm=colors.Normalize(vmin=0,vmax=1)
ax.plot_surface(xs,ys,zs,cmap=cm.Greys,facecolors=cm.Greys(norm(cs)),shade=True)
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
