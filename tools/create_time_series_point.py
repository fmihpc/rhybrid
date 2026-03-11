# interpolates variable values in a given point in all VLSV files in the folder and creates a time series ASCII file
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

def list_of_strings(arg):
 return arg.split(',')

# parse input arguments
parser = argparse.ArgumentParser('create_time_series_point.py')
parser.add_argument('run_path',help='Input run path',type=Path)
parser.add_argument('var_list',help='List of variables: var1,var2,var3...',type=list_of_strings)
parser.add_argument('x',help='x coordinate',type=float)
parser.add_argument('y',help='y coordinate',type=float)
parser.add_argument('z',help='z coordinate',type=float)
args = parser.parse_args()
run_path = str(args.run_path) # ./test_run/
var_list = args.var_list # ['cellB','n_H+sw_ave','v_H+sw_ave']
x = args.x
y = args.y
z = args.z
if Path(run_path).is_dir() == False:
 print('ERROR: input path does not exist (' + run_path + ')')
 sys.exit()

point_coordinates = np.array([[x,y,z]]) # point where to create time series
linear_field_interpolation = True # linear or nearest interpolation

# interpolate time series in a single point
def create_time_series(folder,var_list,point):
 # check only one x,y,z point is given
 if len(point) > 1:
  print('ERROR: only one point per call is supported')
  return False
 files = [f for f in os.listdir(folder) if f.startswith('state') if f.endswith('.vlsv')]
 # check VLSV files are found in the folder
 if len(files) < 1:
  print('ERROR: no VLSV files found (' + folder + ')')
  return False
 files.sort()
 Ntimesteps = len(files)
 vr = alr.vlsvfile.VlsvReader(os.path.join(folder,files[0]))
 # check variables are found
 for var in var_list:
  if vr.check_variable(var) == False:
   print('ERROR: variable ' + var + ' not found in ' + os.path.join(folder,files[0]))
   print('Available variables: ')
   print(str(vr.get_all_variables()))
   return False
 [xmin,ymin,zmin,xmax,ymax,zmax] = vr.get_spatial_mesh_extent()
 [mx,my,mz] = vr.get_spatial_mesh_size() # how many blocks per direction
 [sx,sy,sz] = vr.get_spatial_block_size() # how many cells per block per direction
 nx = mx*sx # number of cells along x
 ny = my*sy # number of cells along y
 nz = mz*sz # number of cells along z
 dx = (xmax-xmin)/nx # should be dx = dy = dz in rhybrid
 box_eps = 0.1*dx
 xp = point[0][0]
 yp = point[0][1]
 zp = point[0][2]
 if(
 (xp < (xmin + box_eps)) |
 (xp > (xmax - box_eps)) |
 (yp < (ymin + box_eps)) |
 (yp > (ymax - box_eps)) |
 (zp < (zmin + box_eps)) |
 (zp > (zmax - box_eps)) ) == True:
  print('ERROR: point (' + str(point) + ') out of domain (' + str(vr.get_spatial_mesh_extent()) + ')')
  return False
 cid = np.ravel(vr.get_cellid(point))[0]
 vrout = alr.calculations.vlsv_intpol_points(vr,point,var_list)
 Nvars = (vrout[2].shape)[1]
 # create header string of output file
 header_string  = 'folder: ' + folder + '\n'
 header_string += 'number of files: ' + str(Ntimesteps) + '\n'
 header_string += 'coordinates of point [m]: x = ' + str(xp) + ', y = ' + str(yp) + ', z = ' + str(zp) + '\n'
 header_string += 'cellid: ' + str(cid) + '\n'
 header_string += 'columns: ' + str(Nvars) + '\n'
 header_string += 't ' + str(vrout[3][13:-1])
 # create output variables and loop through files (timestep)
 var_time_series = np.empty((Ntimesteps,Nvars+1))
 var_time_series.fill(np.nan)
 for ii in range(Ntimesteps):
  vr = alr.vlsvfile.VlsvReader(os.path.join(folder,files[ii]))
  t = vr.read_parameter('time')
  vrout = alr.calculations.vlsv_intpol_points(vr,point,var_list,interpolation_order=linear_field_interpolation)
  #if np.isnan(np.min(vrout[2])) == 1:
  # print('ERROR: interpolation did not succeed, output contains nans')
  # return False
  var_time_series[ii] = np.append(t,vrout[2])
 # save output file
 np.savetxt('time_series.dat',var_time_series,header=header_string)
 return True

# create time series file
success = create_time_series(run_path,var_list,point_coordinates)
if success == False:
 print('ERROR: failed processing: ' + run_path)

