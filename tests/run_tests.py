# This file is part of the RHybrid simulation.
#
# Copyright 2024- Finnish Meteorological Institute
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# run RHybrid test suite
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import pytools as pt
import numpy as np
import os
import sys
import datetime
import subprocess
matplotlib.rcParams.update({"font.size": 14})
matplotlib.rcParams["lines.linewidth"] = 2

# unique template string for folder and config file names
tmpDateTimeStr = "tmp_" + str(datetime.datetime.now().strftime("%Y%m%d%H%M%S%f"))

# template config file used to create test runs
runCfgFileTemplate = "test_empty_1d.cfg"

# used corsair/rhybrid binary and mpirun commands
runBinary = "~/bin/corsair/corsair_rhybrid"
runCommandPrefix = "mpirun -n 3 " + runBinary + " --runconfig="
runCommandSuffix = " >> std_out_err.txt"

# 1D test parameters

# fixed parameters
Ncells = 40
macroparticlesPerCell = 1000
flowThroughs = 3 # how many times the undisturbed flow crosses the domain before simulation ends
NdataSaves = 10.0 # how many data save steps in total
dx = 100e3
Bx = 0.0
By = 0.0
Bz = 0.0

# varied parameters
listUsw = np.linspace(350,650,2)*1e3
#listUsw = np.array((350e3,))
#listUdir = ("-x","+x","-y","+y","-z","+z")
listUdir = ("-x","+y")
#listUdir = ("-x",)
#listNsw = np.linspace(100,10,2)*1e6
listNsw = np.array((1e6,))
#listTsw = np.array((1e5,1e4))
listTsw = np.array((1e5,))

# estimate of dx/dt (max. signal speed) for time step calculation
maxSignalSpeed = 10.0*max(listUsw)

# round to Nsig significant digits
def roundToSigificantDigits(x,Nsig=2):
 return np.round( x,Nsig - int( np.floor(np.log10(np.abs(x))) )-1 )

# calculate time step
dt = roundToSigificantDigits(dx/maxSignalSpeed,2)

# calculate box dimensions
boxLength = Ncells*dx
boxMin = -boxLength/2.0
boxMax = +boxLength/2.0
if (boxLength%2) != 0:
 print("ERROR: box length not divisible by 2" + str(boxLength))
 exit
if (dx%2) != 0:
 print("ERROR: dx not divisible by 2" + str(dx))
 exit

# calculate other parameters
flowThroughTime = boxLength/max(listUsw)
maxTime = flowThroughTime*flowThroughs
maxTimesteps = int(np.ceil(maxTime/dt))
dataSaveInterval = int(np.floor(maxTimesteps/NdataSaves))

# dictionary for substitution of fixed run parameters to their place holders in templace config file
rp_default = dict()
rp_default["TMP_DT"] = dt
rp_default["TMP_MAXIMUM_TIMESTEPS"] = maxTimesteps
rp_default["TMP_DATA_SAVE_INTERVAL"] = dataSaveInterval
rp_default["TMP_IMF_BX"] = Bx
rp_default["TMP_IMF_BY"] = By
rp_default["TMP_IMF_BZ"] = Bz
rp_default["TMP_MACROPARTICLES_PER_CELL"] = macroparticlesPerCell

# combine varied parameters
comb = list()
[comb.append((nsw,Tsw,Usw,Udir)) for nsw in listNsw for Tsw in listTsw for Usw in listUsw for Udir in listUdir]

# create full run parameter dictionaries
runParams = list()
for nsw,Tsw,Usw,Udir in comb:
 rp = dict(rp_default)
 if "-" in Udir:
  Usw = -np.fabs(Usw)
 elif "+" in Udir:
  Usw = +np.fabs(Usw)
 else:
  print("ERROR: Udir " + str(Udir))
  exit
 if "x" in Udir:
  rp["TMP_X_PERIODIC"] = "no"
  rp["TMP_Y_PERIODIC"] = "yes"
  rp["TMP_Z_PERIODIC"] = "yes"
  rp["TMP_X_MIN"] = boxMin
  rp["TMP_X_MAX"] = boxMax
  rp["TMP_Y_MIN"] = rp["TMP_Z_MIN"] = -dx/2.0
  rp["TMP_Y_MAX"] = rp["TMP_Z_MAX"] = +dx/2.0
  rp["TMP_X_SIZE"] = Ncells
  rp["TMP_Y_SIZE"] = rp["TMP_Z_SIZE"] = 1
  rp["TMP_SW_VELOCITY"] = (Usw,0,0)
 elif "y" in Udir:
  rp["TMP_X_PERIODIC"] = "yes"
  rp["TMP_Y_PERIODIC"] = "no"
  rp["TMP_Z_PERIODIC"] = "yes"
  rp["TMP_Y_MIN"] = boxMin
  rp["TMP_Y_MAX"] = boxMax
  rp["TMP_X_MIN"] = rp["TMP_Z_MIN"] = -dx/2.0
  rp["TMP_X_MAX"] = rp["TMP_Z_MAX"] = +dx/2.0
  rp["TMP_Y_SIZE"] = Ncells
  rp["TMP_X_SIZE"] = rp["TMP_Z_SIZE"] = 1
  rp["TMP_SW_VELOCITY"] = (0,Usw,0)
 elif "z" in Udir:
  rp["TMP_X_PERIODIC"] = "yes"
  rp["TMP_Y_PERIODIC"] = "yes"
  rp["TMP_Z_PERIODIC"] = "no"
  rp["TMP_Z_MIN"] = boxMin
  rp["TMP_Z_MAX"] = boxMax
  rp["TMP_X_MIN"] = rp["TMP_Y_MIN"] = -dx/2.0
  rp["TMP_X_MAX"] = rp["TMP_Y_MAX"] = +dx/2.0
  rp["TMP_Z_SIZE"] = Ncells
  rp["TMP_X_SIZE"] = rp["TMP_Y_SIZE"] = 1
  rp["TMP_SW_VELOCITY"] = (0,0,Usw)
 else:
  print("ERROR Udir = " + Udir)
 rp["TMP_SW_DENSITY"] = nsw
 rp["TMP_SW_TEMPERATURE"] = Tsw
 runParams.append(rp)

# plot all scalar and vector variable from a file along primary axis with interpolation
def doPlotsWithInterpolation(fileName):
 intpolOrder = 0 # 0: NGP interpolation, 1: linear interpolation

 # read file
 vr = pt.vlsvfile.VlsvReader(fileName)
 # simulation box dimensions
 [xmin,ymin,zmin,xmax,ymax,zmax] = vr.get_spatial_mesh_extent()
 [mx,my,mz] = vr.get_spatial_mesh_size() # how many blocks per direction
 [sx,sy,sz] = vr.get_spatial_block_size() # how many cells per block per direction
 nx = mx*sx # number of cells along x
 ny = my*sy # number of cells along y
 nz = mz*sz # number of cells along z
 dx = (xmax-xmin)/nx # should be dx = dy = dz in rhybrid

 # find out run dimensionality
 simDim = ""
 if (nx > 1) & (ny > 1) & (nz > 1):
  simDim = "3d"
 elif (nx > 1) & (ny > 1) & (nz == 1):
  simDim = "2d_xy"
 elif (nx > 1) & (ny == 1) & (nz > 1):
  simDim = "2d_xz"
 elif (nx == 1) & (ny > 1) & (nz > 1):
  simDim = "2d_yz"
 elif (nx > 1) & (ny == 1) & (nz == 1):
  simDim = "1d_x"
 elif (nx == 1) & (ny > 1) & (nz == 1):
  simDim = "1d_y"
 elif (nx == 1) & (ny == 1) & (nz > 1):
  simDim = "1d_z"

 # create interpolation coordinates (x,y,z)
 if (simDim == "3d") | (simDim == "2d_xy") | (simDim == "2d_xz") | (simDim == "1d_x"):
  xp = np.linspace(xmin+(2*dx),xmax-(2*dx),nx-4) # along x-axis
  yp = np.ones(len(xp))*(ymax + ymin)/2.0
  zp = np.ones(len(xp))*(zmax + zmin)/2.0
  pp = xp/1e3; pstr = "x [km]"; # plotting parameter and descr. string
 elif (simDim == "2d_yz") | (simDim == "1d_y"):
  yp = np.linspace(ymin+(2*dx),ymax-(2*dx),ny-4) # along y-axis
  xp = np.ones(len(yp))*(xmax + xmin)/2.0
  zp = np.ones(len(yp))*(zmax + zmin)/2.0
  pp = yp/1e3; pstr = "y [km]"; # plotting parameter and descr. string
 elif (simDim == "1d_z"):
  zp = np.linspace(zmin+(2*dx),zmax-(2*dx),nz-4) # along z-axis
  xp = np.ones(len(zp))*(xmax + xmin)/2.0
  yp = np.ones(len(zp))*(ymax + ymin)/2.0
  pp = zp/1e3; pstr = "z [km]"; # plotting parameter and descr. string
 else:
  print("error: " + simDim)

 Npoints = len(pp)

 #tp = np.linspace(0,1,Npoints) # artifial time parameter along trajectory

 # do interpolation
 points = (np.array([xp,yp,zp])).transpose()
 varList = vr.get_all_variables()
 [crd,cellids,outVars,header]=pt.calculations.vlsv_intpol_points(vr,points,varList,"pass",intpolOrder)
 Nvars = len(varList)
 Nscalars = outVars.shape[1]
 #x = crd[:,0] # should be the same as xp
 #y = crd[:,1] # should be the same as yp
 #z = crd[:,2] # should be the same as zp

 # prepare interpolated variables
 d = {}
 Ndim = np.zeros(Nvars,"int") # variable vector dimension (1: scalar, 3: 3d vector)
 ii = 0
 for jj in range(Nvars):
  varName = varList[jj]
  Ndim[jj] = vr.read_variable_vectorsize(varName)
  if Ndim[jj] == 1:
   d[varName] = outVars[:,ii]; ii = ii+1;
  elif Ndim[jj] == 3:
   for kk in range(Ndim[jj]):
    varNameTmp = varName
    if kk == 0:
     varNameTmp = varNameTmp + "_x"
     d[varName + "_tot"] = np.sqrt(np.square(outVars[:,ii]) + np.square(outVars[:,ii+1]) + np.square(outVars[:,ii+2]))  # magnitude
    elif kk == 1:
     varNameTmp = varNameTmp + "_y"
    elif kk == 2:
     varNameTmp = varNameTmp + "_z"
    d[varNameTmp] = outVars[:,ii]; ii = ii+1;
  else:
   print("error"); exit;

 if ii != Nscalars:
  print("error"); exit;

 Ny = Nvars
 Nx = 1
 Nsub = 1

 plt.figure(figsize=(12,18))

 for ii in range(Nvars):
  varName = varList[ii]
  plt.subplot(Ny,Nx,Nsub); Nsub = Nsub + 1;
  if Ndim[ii] == 1:
   plt.plot(pp/1e3,d[varName])
  elif Ndim[ii] == 3:
   plt.plot(pp/1e3,d[varName + "_x"],"-r",label=varName + "_x")
   plt.plot(pp/1e3,d[varName + "_y"],"-g",label=varName + "_y")
   plt.plot(pp/1e3,d[varName + "_z"],"-b",label=varName + "_z")
   plt.plot(pp/1e3,d[varName + "_tot"],"-k",label=varName + "_tot")
   plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
  else:
   print("error")
   exit

  plt.ylabel(varName)
  plt.grid()
  plt.xlabel(pstr);
 plt.savefig("fig_all_vars_intpol.png",bbox_inches="tight")

# plot all scalar and vector variables from a file along primary axis
def doPlots(fileName):
 # read file
 vr = pt.vlsvfile.VlsvReader(fileName)
 [xmin,ymin,zmin,xmax,ymax,zmax] = vr.get_spatial_mesh_extent() # simulation box dimensions
 [mx,my,mz] = vr.get_spatial_mesh_size() # how many blocks per direction
 [sx,sy,sz] = vr.get_spatial_block_size() # how many cells per block per direction
 nx = mx*sx # number of cells along x
 ny = my*sy # number of cells along y
 nz = mz*sz # number of cells along z
 dx = (xmax-xmin)/nx # should be dx = dy = dz in rhybrid
 # find out run dimensionality
 pstr = ""
 if (nx > 1) & (ny == 1) & (nz == 1):
  pstr = "x"
  pp = np.linspace(xmin + 0.5*dx,xmax - 0.5*dx,nx) # plotting parameter along x axis
 elif (nx == 1) & (ny > 1) & (nz == 1):
  pstr = "y"
  pp = np.linspace(ymin + 0.5*dx,ymax - 0.5*dx,ny) # plotting parameter along y axis
 elif (nx == 1) & (ny == 1) & (nz > 1):
  pstr = "z"
  pp = np.linspace(zmin + 0.5*dx,zmax - 0.5*dx,nz) # plotting parameter along z axis
 else:
  print("ERROR: more than 1D")
  return

 # 1D sorting indices
 iiSorted = vr.read_variable("CellID").argsort()
 #cids = vr.read_variable("CellID")[iiSorted]
 #crds = vr.get_cell_coordinates(cids)
 varList = vr.get_all_variables()

 # create figure and subplots
 plt.figure(figsize=(12,18))
 Ny = len(varList)
 Nx = 1
 Nsub = 1
 for varName in varList:
  plt.subplot(Ny,Nx,Nsub); Nsub = Nsub + 1;
  D = vr.read_variable(varName)[iiSorted]
  Ndim = vr.read_variable_vectorsize(varName)
  if Ndim == 1:
   plt.plot(pp/1e3,D)
  elif Ndim == 3:
   Dtot = np.sqrt(np.square(D[:,0]) + np.square(D[:,1]) + np.square(D[:,2]))
   plt.plot(pp/1e3,D[:,0],"-r",label=varName + "_x")
   plt.plot(pp/1e3,D[:,1],"-g",label=varName + "_y")
   plt.plot(pp/1e3,D[:,2],"-b",label=varName + "_z")
   plt.plot(pp/1e3,Dtot,"--k",label=varName + "_tot")
   plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
  else:
   print("ERROR: dimensionality: " + varName)
  plt.ylabel(varName)
  plt.xlabel(pstr + " [km]")
  plt.grid()
 plt.savefig("fig_all_vars.png",bbox_inches="tight")

def checkFlowConditions(fileName):
 # read file
 vr = pt.vlsvfile.VlsvReader(fileName)
 [xmin,ymin,zmin,xmax,ymax,zmax] = vr.get_spatial_mesh_extent() # simulation box dimensions
 [mx,my,mz] = vr.get_spatial_mesh_size() # how many blocks per direction
 [sx,sy,sz] = vr.get_spatial_block_size() # how many cells per block per direction
 nx = mx*sx # number of cells along x
 ny = my*sy # number of cells along y
 nz = mz*sz # number of cells along z
 dx = (xmax-xmin)/nx # should be dx = dy = dz in rhybrid
 Ncells = nx*ny*nz
 Nbc = round(np.floor(Ncells/10)) # boundary cells not to include in mean

 # 1D sorting indices
 iiSorted = vr.read_variable("CellID").argsort()
 # flow variables
 n = vr.read_variable("n_H+")[iiSorted]
 T = vr.read_variable("T_H+")[iiSorted]
 U = vr.read_variable("v_H+")[iiSorted]; Ux = U[:,0]; Uy = U[:,1]; Uz = U[:,2];
 B = vr.read_variable("cellB")[iiSorted]; Bx = B[:,0]; By = B[:,1]; Bz = B[:,2];
 # average flow conditions in the domain (not including cells near the boundaries of the primary axis)
 n_ = np.mean(n[0+Nbc:-1-Nbc])
 T_ = np.mean(T[0+Nbc:-1-Nbc])
 Ux_ = np.mean(Ux[0+Nbc:-1-Nbc])
 Uy_ = np.mean(Uy[0+Nbc:-1-Nbc])
 Uz_ = np.mean(Uz[0+Nbc:-1-Nbc])
 Bx_ = np.mean(Bx[0+Nbc:-1-Nbc])
 By_ = np.mean(By[0+Nbc:-1-Nbc])
 Bz_ = np.mean(Bz[0+Nbc:-1-Nbc])
 #print("n = " + str(n_) + ", T = " + str(T_) + ", Ux = " + str(Ux_) + ", Uy = " + str(Uy_) + ", Uz = " + str(Uz_) + ", Bx = " + str(Bx_) + ", By = " + str(By_) + ", Bz = " + str(Bz_))

# function: do run analysis
def analyzeRun():
 # check if logfile.txt exists
 if os.path.isfile("logfile.txt") == False:
  print("ERROR: no logfile.txt found")
  return False
 # check if run finished successfully
 if subprocess.run(["grep","-q","(MAIN) Exiting simulation after successful run.","logfile.txt"]).returncode == 1:
  print("RUN FAILED")
 sr = subprocess.run(["grep","ERROR","logfile.txt"],capture_output=True)
 # check if run had errors
 if sr.returncode == 0:
  print("RUN HAD ERRORS:")
  print("\t stdout: " + str(sr.stdout))
  print("\t stderr: " + str(sr.stderr))
 del sr
 # check if run had warnings
 sr = subprocess.run(["grep","WARNING","logfile.txt"],capture_output=True)
 if sr.returncode == 0:
  print("RUN HAD WARNINGS:")
  print("\t stdout: " + str(sr.stdout))
  print("\t stderr: " + str(sr.stderr))
 del sr

 # find VLSV files
 runFiles = []
 for f in sorted(os.listdir("./")):
  if (f.startswith("state") and f.endswith(".vlsv")):
   runFiles.append(f)
 if len(runFiles) < 1:
  print("NO VLSV FILES FOUND")
  return False

 # analyze the second last file
 fileName = runFiles[-2]

 #print(fileName)
 #doPlotsWithInterpolation(fileName)
 doPlots(fileName)
 checkFlowConditions(fileName)
 return True

# create run folder and config file from templace
def createAndPerformRun(rp,runNumber):
 runNumberStr = str(runNumber).zfill(4)
 runFolder = tmpDateTimeStr + "_run" + runNumberStr + "/"
 runCfgFile = tmpDateTimeStr + "_test_empty_1d_run" + runNumberStr + ".cfg"
 print("RUN: " + runFolder)
 # create run folder
 if subprocess.run(["mkdir","-p",runFolder]).returncode != 0:
  print("ERROR creating: " + runFolder)
  quit()
 # cd to run folder
 try:
  os.chdir(runFolder)
 except:
  print("ERROR changing to a folder: " + runFolder,sys.exc_info())
  quit()
 # copy template cfg file to folder
 if subprocess.run(["cp","-p","../test_empty_1d.cfg","./" + runCfgFile]).returncode != 0:
  print("ERROR creating file: " + runCfgFile)
  quit()
 # update variable place holders in template cfg file
 for k in rp.keys():
  #print(k + " : " + str(rp[k]))
  if subprocess.run(["sed","-i","s/" + str(k) + "/" + str(rp[k]) + "/g",runCfgFile]).returncode != 0:
   print("ERROR updating parameter: " + str(k) + " in file: " + runCfgFile)
   quit()
 # check all place holders updated in template cfg file
 if subprocess.run(["grep","-q","TMP_",runCfgFile]).returncode != 1:
  print("ERROR TMP_* place holder variable(s) still exists in : " + runCfgFile)
  quit()

 # perfom simulation run
 runCommandFull = runCommandPrefix + runCfgFile + runCommandSuffix
 print(runCommandFull)
 if subprocess.run([runCommandFull],shell=True).returncode != 0:
  print("RUN FAILED: " + runCommandFull)

 # analyze run
 if analyzeRun() == False:
  print("ANALYSIS: FAILED")
 else:
  print("ANALYSIS: OK")

 del runCommandFull

 # return to parent folder
 try:
  os.chdir("../")
 except:
  print("ERROR changing to folder: ../",sys.exc_info())
  quit()
 # check we really are in the main folder
 if os.path.isfile(runCfgFileTemplate) == False:
  print("ERROR changing to main folder")
  quit()

# perform runs
for ii in range(len(runParams)):
 rp = runParams[ii]
 createAndPerformRun(runParams[ii],ii)




