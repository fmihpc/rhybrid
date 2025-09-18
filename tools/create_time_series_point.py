# interpolates variable values in a given point in all VLSV files in the folder and creates a time series ASCII file
import os
import pytools as pt
import numpy as np

pointCoordinates = np.array([[10e6,-1e5,0]]) # point where to create time series
varListSelected = ['cellB','n_H+sw_ave','v_H+sw_ave'] # variables to interpolate
folderRun = os.path.join(os.getenv("HOME"),'bin/corsair/testrun/'); # run folder
runName = 'testrun01' # run name
intpolOrder = 1 # interpolation order

# interpolate time series from a single run
def createTimeSeriesRun(folder,runStr,varList,rp):
 # check only one x,y,z point is given
 if len(rp) > 1:
  print('ERROR: only one point per call is supported')
  return False
 files = [f for f in os.listdir(folder) if f.startswith('state') if f.endswith('.vlsv')]
 # check VLSV files are found in the folder
 if len(files) < 1:
  print('ERROR: not files found')
  return False
 files.sort()
 Ntimesteps = len(files)
 vr = pt.vlsvfile.VlsvReader(os.path.join(folder,files[0]))
 #varList = vr.get_variables() # interpolate all variables
 # check variables are found
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
 cid = np.ravel(vr.get_cellid(rp))[0]
 vrout = pt.calculations.vlsv_intpol_points(vr,rp,varList)
 Nvars = (vrout[2].shape)[1]
 # create header string of output file
 headerStr = 'run name: ' + runStr + '\n'
 headerStr += 'folder: ' + folder + '\n'
 headerStr += 'number of files: ' + str(Ntimesteps) + '\n'
 headerStr += 'coordinates of point [m]: x = ' + str(rp[0][0]) + ', y = ' + str(rp[0][1]) + ', z = ' + str(rp[0][2]) + '\n'
 headerStr += 'cellid: ' + str(cid) + '\n'
 headerStr += 'columns: ' + str(Nvars) + '\n'
 headerStr += 't ' + str(vrout[3][13:-1])
 # create output variables and loop through files (timestep)
 varTimeSeries = np.empty((Ntimesteps,Nvars+1))
 varTimeSeries.fill(np.nan)
 for ii in range(Ntimesteps):
  vr = pt.vlsvfile.VlsvReader(os.path.join(folder,files[ii]))
  t = vr.read_parameter("time")
  vrout = pt.calculations.vlsv_intpol_points(vr,rp,varList,interpolation_order=intpolOrder)
  if np.isnan(np.min(vrout[2])) == 1:
   print('ERROR: interpolation did not succeed')
   return False
  varTimeSeries[ii] = np.append(t,vrout[2])
 # save output file
 np.savetxt('time_series_' + runStr + '.dat',varTimeSeries,header=headerStr)
 return True

# create time series file
runOk = createTimeSeriesRun(folderRun,runName,varListSelected,pointCoordinates)
if runOk == False:
 print('ERROR: failed processing ' + runName)

