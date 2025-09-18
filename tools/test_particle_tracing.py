# propagete test particles in E and B fields from VLSV files
import os
#import pytools as pt
import numpy as np
import bbtracer

# radius of inner boundary
Rp = 2440e4

# test particle mass and charge
m = 1.672621716e-27
q = 1.60217653e-19
# timestep
dt = 0.1
# initial particle position and velocity
rp = np.array((10e6,-1e5,0))
vp = np.array((0,0,0))

varNameB = 'cellB'
varNameUe = 'cellUe'
#varNameE = 'nodeE'
varListSelected = [varNameB,varNameUe] # variables to interpolate
folderRun = os.path.join(os.getenv("HOME"),'bin/corsair/testrun/'); # run folder
runName = 'testrun01' # run name
intpolOrder = 1 # interpolation order

def getRunInformation(folder,runStr,varList):
 files = [f for f in os.listdir(folder) if f.startswith('state') if f.endswith('.vlsv')]
 # check VLSV files are found in the folder
 if len(files) < 1:
  print('ERROR: no VLSV files found')
  return False
 files.sort()
 Ntimesteps = len(files)
 vr = pt.vlsvfile.VlsvReader(os.path.join(folder,files[0]))
 if set(varList).issubset(vr.get_all_variables()) == False:
  print('ERROR: error variable(s) not found in: ' + str(vr.get_all_variables()))
  return False
 [xmin,ymin,zmin,xmax,ymax,zmax] = vr.get_spatial_mesh_extent()
 [mx,my,mz] = vr.get_spatial_mesh_size() # how many blocks per direction
 [sx,sy,sz] = vr.get_spatial_block_size() # how many cells per block per direction
 nx = mx*sx # number of cells along x
 ny = my*sy # number of cells along y
 nz = mz*sz # number of cells along z
 dx = (xmax-xmin)/nx # should be dx = dy = dz in rhybrid
 res = dict()
 res['folder'] = folder
 res['runName'] = runStr
 res['varList'] = varList
 res['files'] = files
 res['xmin'] = xmin
 res['xmax'] = xmax
 res['ymin'] = ymin
 res['ymax'] = ymax
 res['zmin'] = zmin
 res['zmax'] = zmax
 res['dx'] = dx
 return res
 #cid = np.ravel(vr.get_cellid(rp))[0]
 #vrout = pt.calculations.vlsv_intpol_points(vr,rp,varList)

# Boris-Bunemann algorithm
def bbstep(E,B,r,v,m,q,dt):
 # move particle
 r = r + dt*v;
 vold = v
 #accelerate particle
 qmideltT2= 0.5*q*dt/m
 dv = qmideltT2*E
 t = qmideltT2*B
 t2 = np.linalg.norm(t)
 b2 = 2.0/(1.0 + t2)
 s = b2*t
 vm = v + dv
 v0 = vm + np.cross(vm,t)
 vp = vm + np.cross(v0,s)
 v = vp + dv
 #dv = v-vold

#def propagateTestParticles(vr,):
# bbstep(np.array((1e-9,0,0)),np.array((0,0.1,0)),np.array((0,0,0)),np.array((0,0,0)),1,1,1)

#runInfo = getRunInformation(folderRun,runName,varListSelected)
#print(runInfo)

