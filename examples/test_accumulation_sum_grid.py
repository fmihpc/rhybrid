import os
import pytools as pt

refvalue = 1e27

runFolder = './'
runFiles = []
for f in sorted(os.listdir(runFolder)):
 if (f.startswith("state") and f.endswith(".vlsv")):
  runFiles.append(f)
for f in runFiles:
 # read file
 vr = pt.vlsvfile.VlsvReader(runFolder + f)
 # simulation box dimensions
 [xmin,ymin,zmin,xmax,ymax,zmax] = vr.get_spatial_mesh_extent()
 [mx,my,mz] = vr.get_spatial_mesh_size() # how many blocks per direction
 [sx,sy,sz] = vr.get_spatial_block_size() # how many cells per block per direction
 nx = mx*sx # number of cells along x
 ny = my*sy # number of cells along y
 nz = mz*sz # number of cells along z
 dx = (xmax-xmin)/nx # should be dx = dy = dz in rhybrid
 dV = dx*dx*dx # cell volume
 nH = vr.read_variable_info("n_H+")  # proton density
 # sum number of protons in every grid cell
 sum_nH = sum(nH.data*dV)
 print(f + ': ' + str(sum_nH) + ' : ' + str(sum_nH/refvalue - 1))

