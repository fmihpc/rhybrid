# convert a VLSV file to a VTK XML file
import os
import sys
from pathlib import Path
import argparse
try:
 import analysator as alr
except ModuleNotFoundError as err:
 print("Analysator not found: " + str(err))
 sys.exit()
import analysator as alr
try:
 import numpy as np
except ModuleNotFoundError as err:
 print("NumPy not found: " + str(err))
 sys.exit()
try:
 import pyvista as pv
except ModuleNotFoundError as err:
 print("PyVista not found: " + str(err))
 sys.exit()

# parse input arguments
parser = argparse.ArgumentParser("vlsv_to_vtk.py")
parser.add_argument("input_file",help="Input VLSV file to be convert to VTK",type=Path)
args = parser.parse_args()
input_file = str(args.input_file)
#input_file = os.path.join(os.getenv('HOME'),'bin/corsair/testrun/state00004000.vlsv')
if Path(input_file).is_file() == False:
 print('ERROR: input file does not exist (' + input_file + ')')
 sys.exit()
output_file_prefix = Path(input_file).stem
output_file_mesh = Path(input_file).stem + ".vti"
output_file_particle = Path(input_file).stem + "_particles" + ".vtp"
if Path(output_file_mesh).is_file() == True:
 print('ERROR: output file already exists (' + output_file_mesh + ')')
 sys.exit()

print('converting: ' + input_file)

vr = alr.vlsvfile.VlsvReader(input_file)
# simulation box dimensions
[xmin,ymin,zmin,xmax,ymax,zmax] = vr.get_spatial_mesh_extent()
[mx,my,mz] = vr.get_spatial_mesh_size() # how many blocks per direction
[sx,sy,sz] = vr.get_spatial_block_size() # how many cells per block per direction
nx = mx*sx # number of cells along x
ny = my*sy # number of cells along y
nz = mz*sz # number of cells along z
dx = (xmax-xmin)/nx # should be dx = dy = dz in rhybrid

print('cells: ' + str(nx) + 'x' + str(ny) + 'x' + str(nz))

# simulation time
file_time = -1
if vr.check_parameter('t') == True:
 file_time = vr.read_parameter('t')
elif vr.check_parameter('time') == True:
 file_time = vr.read_parameter('time')
else:
 print('WARNING: time parameter not found, setting as -1')
file_timestep = -1
#if vr.check_parameter('timestep') == True:
# file_timestep = vr.read_parameter('timestep')
#else:
# print('WARNING: timestep parameter not found, setting as -1')

# read and sort cell ids
if vr.check_variable('CellID') == False:
 print('ERROR: CellID variable not found')
 sys.exit()
cellids = vr.read_variable('CellID')
cellids_sorted = cellids.argsort()

# create pyvista ImageData
# reshape 1d VLSV ordering -> 3d (nz,ny,nx)
cellids = cellids[cellids_sorted].reshape(nz,ny,nx)
# reorder axes for ImageData from (nz,ny,nx) -> (nx,ny,nz)
# old axis 0 → z → new axis 2
# old axis 1 → y → new axis 1
# old axis 2 → x → new axis 0
cellids = np.transpose(cellids,(2,1,0))
grid = pv.ImageData(
 dimensions=np.array(cellids.shape) + 1, # points = cells + 1
 spacing=(dx,dx,dx),
 origin=(xmin,ymin,zmin)
)

# add default variables to display
grid.active_scalars_name = None
grid.active_vectors_name = None

# add simulation time in a VisIt and ParaView friendly way
grid.field_data["TIME"] = np.array([file_time], dtype=np.float64)

# loop through all variables in VLSV file
all_vars = vr.get_all_variables()
for var in all_vars:
 # skip cell flags for detectors
 if "detector_flag_" in var:
  continue
 print('reading variable: ' + var)
 dim = vr.read_variable_vectorsize(var) # get variable dimension
 if dim != 1 and dim != 3:
  print('ERROR: variable should be 1d or 3d: ' + var + ': ' + str(dim))
  sys.exit()
 D = vr.read_variable(var) # get VLSV variable
 if dim == 1:
  # prepare a scalar for ImageData
  D = D[cellids_sorted].reshape(nz,ny,nx)
  D = np.transpose(D,(2,1,0))
  D = D.ravel(order="F") # Fortran order
 elif dim == 3:
  # prepare a vector for ImageData
  D = D[cellids_sorted].reshape(nz,ny,nx,3)
  D = np.transpose(D,(2,1,0,3))
  D = D.reshape(-1,3,order="F")
 grid.cell_data[var] = D # add a mesh variable as cell data to ImageData

# save all mesh variables to VTK XML format
print('saving: ' + output_file_mesh)
grid.save(output_file_mesh,binary=True)

# create separate particle files for testing: works in VisIt
# population 1
#n_particles1 = 100
#points1 = np.random.rand(n_particles1,3)*(xmax-xmin)
#particles1 = pv.PolyData(points1)
#particles1.point_data["id"] = np.arange(n_particles1)
#particles1.save(output_file_particle, binary=True)
# population 1
#n_particles2 = 50
#points2 = np.random.rand(n_particles2, 3)*(xmax-xmin)
#particles2 = pv.PolyData(points2)
#particles2.point_data["charge"] = np.random.randn(n_particles2)
#particles2.save("particles_2.vtp", binary=True)

# create header file for mesh and particle files: does not seem to work in VisIt
#import xml.etree.ElementTree as ET
#files = [output_file_mesh, "particles_1.vtp", "particles_2.vtp"]
#names = ["SpatialGrid", "ParticlesPopulation1", "ParticlesPopulation1"]
#vtkfile = ET.Element("VTKFile", type="Collection", version="0.1", byte_order="LittleEndian")
#collection = ET.SubElement(vtkfile, "Collection")
#for f, n in zip(files, names):
# ET.SubElement(collection, "DataSet", file=f, name=n)
# ET.ElementTree(vtkfile).write("all.pvd", encoding="utf-8", xml_declaration=True)

# create pyvista MultiBlock: but does not show particles in VisIt
#mb = pv.MultiBlock()
#mb["grid"] = grid
#mb["particles_type_1"] = particles1
#mb["particles_type_2"] = particles2
#mb["Simulation_Grid"] = grid
#mb["Particles1"] = particles1
#mb["Particles2"] = particles2
#mb.save("simulation_output.vtm") # save all variables to VTK XML format
