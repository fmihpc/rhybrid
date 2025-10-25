import os
import analysator as alr
import numpy as np

def loadRun(folder,varList):
 files = [f for f in os.listdir(folder) if f.startswith('state') if f.endswith('.vlsv')]
 # check VLSV files are found in the folder
 if len(files) < 1:
  print('ERROR: no VLSV files found')
  return False
 files.sort()
 Ntimesteps = len(files)
 # place VLSV readers and simulation time in lists
 vr = list()
 tsim = list()
 for f in files:
  vr.append(alr.vlsvfile.VlsvReader(os.path.join(folder,f)))
  tsim.append(vr[-1].read_parameter('t'))
 # read header information from the first VLSV file
 vr0 = vr[0]
 if set(varList).issubset(vr0.get_all_variables()) == False:
  print('ERROR: error variable(s) not found in: ' + str(vr0.get_all_variables()))
  return False
 [xmin,ymin,zmin,xmax,ymax,zmax] = vr0.get_spatial_mesh_extent()
 [mx,my,mz] = vr0.get_spatial_mesh_size() # how many blocks per direction
 [sx,sy,sz] = vr0.get_spatial_block_size() # how many cells per block per direction
 nx = mx*sx # number of cells along x
 ny = my*sy # number of cells along y
 nz = mz*sz # number of cells along z
 dx = (xmax-xmin)/nx # should be dx = dy = dz in rhybrid
 res = dict()
 res['folder'] = folder
 res['varList'] = varList
 res['xmin'] = xmin
 res['xmax'] = xmax
 res['ymin'] = ymin
 res['ymax'] = ymax
 res['zmin'] = zmin
 res['zmax'] = zmax
 res['dx'] = dx
 res['sim'] = dict()
 res['sim']['readers'] = vr
 res['sim']['t'] = np.array(tsim)
 return res

# Boris-Bunemann algorithm: move and accelerate particles one timestep
def bbStep(E,B,r,v,m,q,dt):
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
