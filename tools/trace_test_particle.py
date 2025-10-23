# exec(open("trace_test_particle.py").read());
# propagete test particles in E and B fields from VLSV files or homogeneous fields
import os
import analysator as alr
import numpy as np
from matplotlib import cm, colors
import matplotlib.pyplot as plt
from bbtracer import loadRun, bbStep

plt.rcParams.update({'font.size': 20})
saveFigures = True # save figures as png files
useHomogeneousFiels = False # use static homogeneous field in tracing

RMe = 2439.7e3 # Mercury
RMa = 3390e3 # Mars
mp = 1.672621716e-27 # proton
me = 9.1093825746e-31 # electron
qe = 1.60217653e-19 # unit charge

# radius of planet / inner boundary, length scale in plots
if useHomogeneousFiels == False:
 Rp = RMa; Rp_str = 'Rp';
 Rp2 = Rp**2
else:
 Rp = 1e3; Rp_str = 'km';

# test particle mass [kg] and charge [C]
m = mp
q = qe
# timestep [s]
dt = 0.5
# max. number of timesteps
Ndt = 200
# constant crandom seed for reproducibility
np.random.seed(71)

if useHomogeneousFiels == True:
 # set homogeneous fields
 UeHomog = np.array((-400e3,0,0))
 BHomog = np.array((0,1e-9,0))
else:
 # set simulation run
 varList = ['cellB','cellUe'] # variables to interpolate
 folderRun = os.path.join(os.getenv("HOME"),'bin/corsair/testrun/')
 intpolOrder = 1 # interpolation order
 simRun = loadRun(folderRun,varList)
 # buffer distance for out boundary conditions [dx]
 boxEps = 0.1*simRun['dx']
 # select a VLSV reader/timestep
 Nreader = 10
 vr = simRun['sim']['readers'][Nreader]

# initial particle positions [m]
r0 = np.array([
(2.5*Rp,0,1.5*Rp),
(2.5*Rp,0,1.2*Rp),
(2.5*Rp,0,-1.5*Rp),
(2.5*Rp,0,-1.3*Rp),
])

Nparticles = len(r0)

# initial particle velocities [m/s]
Usw = (-400e3,0,0)
vth = 30e3
v0 = np.ones(r0.shape)*Usw + vth*np.random.rand(r0.shape[0],r0.shape[1])

if useHomogeneousFiels == True:
 B = np.ones(r0.shape)*BHomog
 Ue = np.ones(r0.shape)*UeHomog
 E = -np.cross(Ue,B)

# initial times [s]
t0 = np.zeros(Nparticles)

if (r0.shape != v0.shape) | (len(r0) != len(t0)) | (len(v0) != len(t0)):
 print("ERROR: r0 and v0 should be same shape and same length as t0")
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
# propagate test particles
for iiDt in range(Ndt):
 # check boundary conditions if not tracing in homogeneous fields
 if useHomogeneousFiels == False:
  # check outer boundary conditions
  oBConds = (
  (rnow[:,0] < (simRun['xmin'] + boxEps)) |
  (rnow[:,0] > (simRun['xmax'] - boxEps)) |
  (rnow[:,1] < (simRun['ymin'] + boxEps)) |
  (rnow[:,1] > (simRun['ymax'] - boxEps)) |
  (rnow[:,2] < (simRun['zmin'] + boxEps)) |
  (rnow[:,2] > (simRun['zmax'] - boxEps))
  )
  # check inner boundary condition
  r2 = np.sum(rnow*rnow,axis=1)
  iBConds = (r2 <= Rp2)
  # full boundary conditions
  BConds = (oBConds | iBConds)
  # remove particles that entered boundaries from propagation
  if any(BConds) == True:
   rnow = rnow[np.where(BConds == False)]
   vnow = vnow[np.where(BConds == False)]
   pid = pid[np.where(BConds == False)]
  # break loop if no more particles
  if len(pid) < 1:
   break

  # interpolate Ue/E and B in particle positions
  [crds,cellids,fields,headerstr] = alr.calculations.vlsv_intpol_points(vr,rnow,varList,interpolation_order=intpolOrder)
  B = fields[:,0:3]
  Ue = fields[:,3:6]
  # E = -Ue x B
  E = -np.cross(Ue,B)

 # Boris-Bunemann step
 (rnow,vnow) = bbStep(E,B,rnow,vnow,m,q,dt)
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

 # faster vectorized version
 #r_arr = np.stack(r)
 #rnow_list_split_arr = np.stack(rnow_list_split)
 #r_arr = np.concatenate((r_arr,rnow_list_split_arr),axis=1)
 #r = [a for a in r_arr]



# number of timesteps taken (=number of bbStep/vlsv_intpol_points function calls)
dtTaken = iiDt + 1
# total number of particle propagations (no initial state)
bbStepsTotal = sum(tp.size for tp in t) - Nparticles
print('')
print('tracing finished')
print('timesteps taken: ' + str(dtTaken))
print('simulation time: ' + str(dtTaken*dt) + ' s')
print('total number of particle propagations/field interpolations: ' + str(bbStepsTotal))

# 1D plot
plt.figure()
for ii in range(len(r)):
 plt.subplot(6,1,1);
 plt.plot(t[ii],r[ii][:,0]/Rp)
 plt.ylabel('x [' + Rp_str +  ']')
 plt.subplot(6,1,2);
 plt.plot(t[ii],r[ii][:,1]/Rp)
 plt.ylabel('y [' + Rp_str + ']')
 plt.subplot(6,1,3);
 plt.plot(t[ii],r[ii][:,2]/Rp)
 plt.ylabel('z [' + Rp_str + ']')
 plt.subplot(6,1,4);
 plt.plot(t[ii],v[ii][:,0]/1e3)
 plt.ylabel('vx [km/s]')
 plt.subplot(6,1,5);
 plt.plot(t[ii],v[ii][:,1]/1e3)
 plt.ylabel('vy [km/s]')
 plt.subplot(6,1,6);
 plt.plot(t[ii],v[ii][:,2]/1e3)
 plt.ylabel('vz [km/s]')
 plt.xlabel('t [s]')

if saveFigures == True:
 plt.savefig("testrun_trace_test_particle_vel.png")
 plt.clf()
 plt.close()

# 3D plot of particle trajectories

ax = plt.figure().add_subplot(projection='3d')

# sphere representing planet (r = Rp)
if useHomogeneousFiels == False:
 u,v = np.meshgrid(np.linspace(0,np.pi,50),np.linspace(0,2*np.pi,50))
 xs = 1 * np.sin(u) * np.cos(v)
 ys = 1 * np.sin(u) * np.sin(v)
 zs = 1 * np.cos(u)
 cs = np.ones(xs.shape)
 cs[xs > 0] = 0
 norm = colors.Normalize(vmin=0,vmax=1)
 ax.plot_surface(xs,ys,zs,cmap=cm.Greys,facecolors=cm.Greys(norm(cs)),shade=True)
 ax.set_box_aspect([1,1,1])

# plot particle trajectory
for ii in range(len(r)):
 ax.plot(r[ii][:,0]/Rp,r[ii][:,1]/Rp,r[ii][:,2]/Rp)

ax.set_xlabel('x [' + Rp_str + ']')
ax.set_ylabel('y [' + Rp_str + ']')
ax.set_zlabel('z [' + Rp_str + ']')
if useHomogeneousFiels == False:
 ax.axes.set_xlim3d(xmin=simRun['xmin']/Rp,xmax=simRun['xmax']/Rp)
 ax.axes.set_ylim3d(ymin=simRun['ymin']/Rp,ymax=simRun['ymax']/Rp)
 ax.axes.set_zlim3d(zmin=simRun['zmin']/Rp,zmax=simRun['zmax']/Rp)

if saveFigures == True:
 plt.savefig("testrun_trace_test_particle_traj.png")
 plt.clf()
 plt.close()
else:
 plt.show(block=True)

