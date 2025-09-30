# See plotter_rhybrid_2d_slice.py for explanation
import sys
import os
import matplotlib
matplotlib.use("Agg") # this is needed for offline saving of graphics
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import pytools as pt
import numpy as np
import scipy as sp
import operator as oper
import socket
import itertools
import subprocess
from multiprocessing import Pool,Value
from matplotlib.patches import Wedge
plt.switch_backend("agg")

useMultiProcessing = 1
printHeaderOnly = 0
mpFileOpenCnt = Value('i',0) # file opening counter for multiprocesses
Nfiles = -1
Nruns = -1
NfileOpens = -1
Nparams = -1
NfigCols = -1
NfigRows = -1
runDim = -1
runDimAxes = ""
HN = "(hostname = " + socket.gethostname() + ") "
matplotlib.rcParams.update({"font.size": 14})
matplotlib.rcParams["lines.linewidth"] = 2
figDpi = 100
matplotlib.rcParams["figure.dpi"] = figDpi
matplotlib.rcParams["savefig.dpi"] = figDpi
figResolutionX = 1920
figResolutionY = 1080
figureSize = (figResolutionX/figDpi,figResolutionY/figDpi)
tickDir = "out"
tickLength = 2
tickWidth = 1
xTickAngle = 45

outputFileNamePrefix = "./png/" # note: if folder name given here, that folder should exists

xPlane = -2.0 # plotted yz plane at constant x = xPlane*Rp (3D runs, ignored for 2D runs)
yPlane = 0.0 # plotted xz plane at constant y = yPlane*Rp (3D run, ignored for 2D runs)
zPlane = 0.0 # plotted xy plane at constant z = zPlane*Rp (3D run, ignored for 2D runs)
axisLimsZoom = -1 # show full simulation domain
#axisLimsZoom = np.multiply([-4.0,4.0,-4.0,4.0,-4.0,4.0],1.0) # zoom into this domain in Rp
showPlanet = (1,1,0) # if to show a planet in xz, xy, yz plane

colormap1 = "inferno" # basic perceptually uniform colormap, https://matplotlib.org/stable/users/explain/colors/colormaps.html
colormap2 = "seismic" # basic blue-white-red difference colormap

# parameters to be plotted (=found in P_settings below and in VLSV file), filled automatically
P = list()

# setting database for each VLSV file variable known to this script, variables not listed here are not plotted
P_settings = list()

# each P_settings element should have a dictionary the following key:value pairs:
# param: variable name in VLSV file [e.g. "cellB"]
# type: parameter type defined in the plotPanel function ["magnitude"/"scalar"/"xcomp"/"ycomp"/"zcomp"]
# str: variable description string to be show in the figure [e.g. "$|B|$ [nT]"]
# log: use logarihmic color scale [0/1]
# lims: variable values at color scale limit [e.g. (10e-9,1000e-9)]
# unit: scaling factor from variable units in VLSV file to the plot [e.g. 1e-9]
# colormap: name of the colormap [e.g. "inferno"]
# filename: prefix of the saved PNG filename [e.g. "B"]
# sigma: smoothing parameter (Standard deviation for Gaussian kernel) for the sp.ndimage.filters.gaussian_filter [e.g. -1 (not used) or 0.9 (some smoothing)]

P_settings.append({"param":"cellB","type":"magnitude","str":"$|B|$ [nT]","log":0,"lims":(0,40e-9),"unit":1e-9,"colormap":colormap1,"filename":"B","sigma":-1})
P_settings.append({"param":"cellB","type":"xcomp","str":"$B_x$ [nT]","log":0,"lims":(-20e-9,20e-9),"unit":1e-9,"colormap":colormap2,"filename":"Bx","sigma":-1})
P_settings.append({"param":"cellB","type":"ycomp","str":"$B_y$ [nT]","log":0,"lims":(-20e-9,20e-9),"unit":1e-9,"colormap":colormap2,"filename":"By","sigma":-1})
P_settings.append({"param":"cellB","type":"zcomp","str":"$B_z$ [nT]","log":0,"lims":(-20e-9,20e-9),"unit":1e-9,"colormap":colormap2,"filename":"Bz","sigma":-1})

P_settings.append({"param":"cellBAverage","type":"magnitude","str":"$|B|$ [nT]","log":0,"lims":(0,40e-9),"unit":1e-9,"colormap":colormap1,"filename":"B_ave","sigma":-1})
P_settings.append({"param":"cellBAverage","type":"xcomp","str":"$B_x$ [nT]","log":0,"lims":(-20e-9,20e-9),"unit":1e-9,"colormap":colormap2,"filename":"Bx_ave","sigma":-1})
P_settings.append({"param":"cellBAverage","type":"ycomp","str":"$B_y$ [nT]","log":0,"lims":(-20e-9,20e-9),"unit":1e-9,"colormap":colormap2,"filename":"By_ave","sigma":-1})
P_settings.append({"param":"cellBAverage","type":"zcomp","str":"$B_z$ [nT]","log":0,"lims":(-20e-9,20e-9),"unit":1e-9,"colormap":colormap2,"filename":"Bz_ave","sigma":-1})

P_settings.append({"param":"n_H+sw_ave","type":"scalar","str":"$n$(H$^+_\mathrm{sw}$) [m$^{-3}$]","log":1,"lims":(1e4,1e9),"unit":1,"colormap":colormap1,"filename":"Hsw_n","sigma":-1})
P_settings.append({"param":"n_O+_ave","type":"scalar","str":"$n$(O$^{+}_\mathrm{planet}$) [m$^{-3}$]","log":1,"lims":(0.04*1e4,0.04*1e9),"unit":1,"colormap":colormap1,"filename":"Oplanet_n","sigma":-1})

P_settings.append({"param":"nodeE","type":"magnitude","str":"$|E|$ [V/m]","log":1,"lims":(0.0001,0.1),"unit":1,"colormap":colormap1,"filename":"E","sigma":-1})
P_settings.append({"param":"nodeE","type":"xcomp","str":"$E_x$ [mV/m]","log":0,"lims":(-0.02,0.02),"unit":1e-3,"colormap":colormap2,"filename":"Ex","sigma":-1})
P_settings.append({"param":"nodeE","type":"ycomp","str":"$E_y$ [mV/m]","log":0,"lims":(-0.02,0.02),"unit":1e-3,"colormap":colormap2,"filename":"Ey","sigma":-1})
P_settings.append({"param":"nodeE","type":"zcomp","str":"$E_z$ [mV/m]","log":0,"lims":(-0.02,0.02),"unit":1e-3,"colormap":colormap2,"filename":"Ez","sigma":-1})

P_settings.append({"param":"cellUe","type":"magnitude","str":"$|U|(e^-)$ [km/s]","log":0,"lims":(0,700e3),"unit":1e3,"colormap":colormap1,"filename":"e_U","sigma":-1})
P_settings.append({"param":"cellUe","type":"xcomp","str":"$U_x(e^-)$ [km/s]","log":0,"lims":(-600e3,600e3),"unit":1e3,"colormap":colormap2,"filename":"e_Ux","sigma":-1})
P_settings.append({"param":"cellUe","type":"ycomp","str":"$U_y(e^-)$ [km/s]","log":0,"lims":(-600e3,600e3),"unit":1e3,"colormap":colormap2,"filename":"e_Uy","sigma":-1})
P_settings.append({"param":"cellUe","type":"zcomp","str":"$U_z(e^-)$ [km/s]","log":0,"lims":(-600e3,600e3),"unit":1e3,"colormap":colormap2,"filename":"e_Uz","sigma":-1})

P_settings.append({"param":"v_H+sw_ave","type":"magnitude","str":"$|U|$(H$^+_\mathrm{sw}$) [km/s]","log":0,"lims":(0,700e3),"unit":1e3,"colormap":colormap1,"filename":"Hsw_U","sigma":-1})
P_settings.append({"param":"v_H+sw_ave","type":"xcomp","str":"$U_x$(H$^+_\mathrm{sw}$) [km/s]","log":0,"lims":(-600e3,600e3),"unit":1e3,"colormap":colormap2,"filename":"Hsw_Ux","sigma":-1})
P_settings.append({"param":"v_H+sw_ave","type":"ycomp","str":"$U_y$(H$^+_\mathrm{sw}$) [km/s]","log":0,"lims":(-600e3,600e3),"unit":1e3,"colormap":colormap2,"filename":"Hsw_Uy","sigma":-1})
P_settings.append({"param":"v_H+sw_ave","type":"zcomp","str":"$U_z$(H$^+_\mathrm{sw}$) [km/s]","log":0,"lims":(-600e3,600e3),"unit":1e3,"colormap":colormap2,"filename":"Hsw_Uz","sigma":-1})

P_settings.append({"param":"v_O+_ave","type":"magnitude","str":"$|U|$(O$^{+}_\mathrm{planet}$) [km/s]","log":0,"lims":(0,700e3),"unit":1e3,"colormap":colormap1,"filename":"Oplanet_U","sigma":-1})
P_settings.append({"param":"v_O+_ave","type":"xcomp","str":"$U_x$(O$^{+}_\mathrm{planet}$) [km/s]","log":0,"lims":(-600e3,600e3),"unit":1e3,"colormap":colormap2,"filename":"Oplanet_Ux","sigma":-1})
P_settings.append({"param":"v_O+_ave","type":"ycomp","str":"$U_y$(O$^{+}_\mathrm{planet}$) [km/s]","log":0,"lims":(-600e3,600e3),"unit":1e3,"colormap":colormap2,"filename":"Oplanet_Uy","sigma":-1})
P_settings.append({"param":"v_O+_ave","type":"zcomp","str":"$U_z$(O$^{+}_\mathrm{planet}$) [km/s]","log":0,"lims":(-600e3,600e3),"unit":1e3,"colormap":colormap2,"filename":"Oplanet_Uz","sigma":-1})

P_settings.append({"param":"T_H+sw","type":"scalar","str":"$T$(H$^+_\mathrm{sw}$) [K]","log":1,"lims":(1e5,1e8),"unit":1,"colormap":colormap1,"filename":"Hsw_T","sigma":-1})

# check and convert command line arguments
Ncores = -100
runFolder = ""
runDescr = ""
Robject = 1
tStartThisProcess = -100
tEndThisProcess = -100
tStartGlobal = -100
tEndGlobal = -100
Nargs = len(sys.argv)
if (Nargs == 9):
 try:
  Ncores = int(sys.argv[1])
  runFolder = str(sys.argv[2])
  runDescr = str(sys.argv[3])
  Robject = float(sys.argv[4])
  tStartThisProcess = int(sys.argv[5])
  tEndThisProcess = int(sys.argv[6])
  tStartGlobal = int(sys.argv[7])
  tEndGlobal = int(sys.argv[8])
 except ValueError:
  print(HN + "ERROR: bad argument:")
  for ii in range(Nargs):
   print(HN + "arg" + str(ii) + " = " + str(sys.argv[ii]))
  quit()
elif (Nargs == 7):
 try:
  Ncores = int(sys.argv[1])
  runFolder = str(sys.argv[2])
  runDescr = str(sys.argv[3])
  Robject = float(sys.argv[4])
  tStartThisProcess = int(sys.argv[5])
  tEndThisProcess = int(sys.argv[6])
  tStartGlobal = tStartThisProcess
  tEndGlobal = tEndThisProcess
 except ValueError:
  print(HN + "ERROR: bad argument:")
  for ii in range(Nargs):
   print(HN + "arg" + str(ii) + " = " + str(sys.argv[ii]))
  quit()
elif (Nargs == 3):
 try:
  Ncores = Ncores = int(sys.argv[1])
  runFolder = str(sys.argv[2])
 except ValueError:
  print(HN + "ERROR: bad argument:")
  for ii in range(Nargs):
   print(HN + "arg" + str(ii) + " = " + str(sys.argv[ii]))
  quit()
 runDescr = "run"
 Robject = 1
 tStartThisProcess = 0
 tEndThisProcess = 1000000
 tStartGlobal = tStartThisProcess
 tEndGlobal = tEndThisProcess
else:
 print(HN + "ERROR: 2, 6 or 8 command line arguments required")
 for ii in range(Nargs):
  print(HN + "arg" + str(ii) + " = " + str(sys.argv[ii]))
 quit()

# add path separator if not already there
runFolder = os.path.join(runFolder,"")

# do check and print information
if os.path.isdir(runFolder) == False:
 print(HN + "ERROR: cannot read folder: " + runFolder)
 quit()
if not any( (fileName.startswith("state") and fileName.endswith(".vlsv")) for fileName in os.listdir(runFolder)):
 print(HN + "ERROR: no state*.vlsv found in folder: " + runFolder)
 quit()
if (Ncores == -1):
 # get only header information using one core
 printHeaderOnly = 1
 Ncores = 1
elif (Ncores < 1) or (Ncores > 100):
 print(HN + "ERROR: negative or otherwise bad number of cores")
 quit()
if (tStartThisProcess < 0) or (tEndThisProcess < 0) or (tStartThisProcess > tEndThisProcess):
 print(HN + "ERROR: negative or otherwise bad start or end times (this process)")
 quit()
if (tStartGlobal < 0) or (tEndGlobal < 0) or (tStartGlobal > tEndGlobal):
 print(HN + "ERROR: negative or otherwise bad start or end times (whole run)")
 quit()

if printHeaderOnly == 0:
 print("running on " + str(Ncores) + " cores")

# 2D multi-run: add manually description
#runDescr = "run001"
#Robject = 1
#runFolder = "../run001/"

# create simRuns list of tuples
simRuns=[(runFolder,runDescr,Robject,"$R_p$")]

# 2D multi-run: add manually further runs
#runFolder2 = "/path/to/run02/"; runDescr2="run02"
#runFolder3 = "/path/to/run03/"; runDescr3="run03"
#runFolder2 = os.getenv("HOME") + "/run02/"; runDescr2="run02"
#runFolder3 = os.getenv("HOME") + "/run03/"; runDescr3="run03"
#simRuns.append((os.path.join(runFolder2,""),runDescr2,Robject,"$R_p$"))
#simRuns.append((os.path.join(runFolder3,""),runDescr3,Robject,"$R_p$"))

# 2D multi-run: add automatically further runs
#listRunFolders = sorted([a for a in os.listdir("../") if os.path.isdir(os.path.abspath("../" + a)) and a.startswith("run0")])
#for ii in range(len(listRunFolders)):
# if ii > 0:
#  simRuns.append((os.path.join("../" + listRunFolders[ii],""),os.path.basename(listRunFolders[ii]),Robject,"$R_p$"))

# number of runs
Nruns = len(simRuns)

# parse parameters from log file
def parseParameter(logFileName,paramGetCmd,paramUnit):
 res = subprocess.run(paramGetCmd,shell=True,capture_output=True,text=True)
 if res.returncode == 1:
  print("ERROR not found: " + paramGetCmd)
  quit()
 value = float(res.stdout.split()[0])*paramUnit
 #print(value)
 return value

# round number to string
def round2str(x,Ndec=10):
 if Ndec == 1:
  return str(int(round(x*Ndec)/Ndec))
 else:
  return str(round(x*Ndec)/Ndec)

# round number to string with unitstr
def round2strWithUnit(x,unitStr):
 if abs(x) > 0:
  return str(round(x*10)/10) + unitStr
 else:
  return "0"

# configure a subplot panel
def configurePanel(P_ii,cMap,ax,showPlanet,showColorbar):
 if showPlanet == 1:
  w1 = Wedge((0,0),1.0,90,270,fc="dimgray")
  w2 = Wedge((0,0),1.0,270,90,fc="w")
  c1 = plt.Circle((0, 0), 1.0, color="k",fill=False,lw=0.5)
  ax.add_artist(w1)
  ax.add_artist(w2)
  ax.add_artist(c1)
 if P_ii["log"] == 1:
  cMap.set_norm(colors.LogNorm(vmin=P_ii["lims"][0]/P_ii["unit"],vmax=P_ii["lims"][1]/P_ii["unit"]))
 plt.gcf().gca().tick_params(which="both",direction=tickDir,length=tickLength,width=tickWidth)

# find the number of column for given coordinate value r between rmin and rmax with nr cells
def chooseNcol(rPlane,rmin,rmax,nr,rstr,Rp,Rp_str):
#chooseNcol(yPlane*Rp,ymin,ymax,ny)
 if (rPlane > rmax) or (rPlane < rmin):
  rPlane = (rmax + rmin)/2.0
  print(HN + "(chooseNcol) WARNING: " + rstr + "Plane coordinate out of domain, setting: " + rstr + "Plane = " + str(rPlane/Rp) + " Rp")
 dr = (rmax-rmin)/nr
 Ncol = int(np.floor((rPlane - rmin)/dr))
 # indices start from zero
 if Ncol >= nr:
  Ncol = nr-1
 if Ncol < 0:
  Ncol = 0
 return Ncol,rPlane/Rp

# 2D multi-run: check if the file exists, if not choose another file with largest time step from the folder
def checkVlsvFileExists(folder__,fileName__):
 if os.path.isfile(folder__ + fileName__) == False:
  tLatest = -1
  fileNameLastTimestep = ""
  for f in sorted(os.listdir(folder__)):
   if (f.startswith("state") and f.endswith(".vlsv")):
    try:
     t = int(f[5:13])
    except ValueError:
     print(HN + "ERROR: could not parse int from file name" + f)
     continue
    if (t >= tStartGlobal) and (t<= tEndGlobal):
     if t > tLatest:
      tLatest = t
      fileNameLastTimestep = f
  if tLatest > -1:
   return (False,True,fileNameLastTimestep)
  else:
   return (False,False,"")
 else:
  return (True,False,"")

# plotting individual subplot in a figure
def plotPanel(fig,axes,ii_run,P_ii,runFolder,vlsvFileName,runStr,Rp,Rp_str,printHeader=0):
 # counter
 global mpFileOpenCnt
 with mpFileOpenCnt.get_lock():
  if (fig != -1) & (printHeader == 0):
   print(HN + "opening: " + runFolder + " | " + str(vlsvFileName) + " | " + P_ii["param"] + " | " + P_ii["type"] + " (" + str(mpFileOpenCnt.value) + "/" + str(NfileOpens) + ")")
  mpFileOpenCnt.value += 1
 # 2D multi-run: check vlsv file exists, if not try to replace with another file
 if 0:
  fileFound,newFileFound,newFileName = checkVlsvFileExists(runFolder,vlsvFileName)
  if fileFound == False:
   if newFileFound == True:
    vlsvFileName = newFileName
   else:
    print(HN + "ERROR: did not find file from " + runFolder)
    return
 # read file
 vr = pt.vlsvfile.VlsvReader(runFolder + vlsvFileName)
 # simulation box dimensions
 [xmin,ymin,zmin,xmax,ymax,zmax] = vr.get_spatial_mesh_extent()
 [mx,my,mz] = vr.get_spatial_mesh_size() # how many blocks per direction
 [sx,sy,sz] = vr.get_spatial_block_size() # how many cells per block per direction
 nx = mx*sx # number of cells along x
 ny = my*sy # number of cells along y
 nz = mz*sz # number of cells along z
 dx = (xmax-xmin)/nx # should be dx = dy = dz in rhybrid
 # full domain
 axisLims = np.divide([xmin,xmax,ymin,ymax,zmin,zmax],Rp)

 # simulation time
 fileTime = vr.read_parameter("t")
 if fileTime is None:
  fileTime = vr.read_parameter("time")

 # plot title string
 titleStr = runStr
 # 2D multi-run: include time string if vlsv file changed from expected one
 #if newFileFound == True:
 # titleStr += " ($t=$" + round2str(fileTime) + " s)"

 # 2D multi-run: check logfile.txt exists in run folder
 if 0:
  logFileName = runFolder + "logfile.txt"
  Ueobs = Bx = By = Bz = Usw = nsw = Tsw = 0.0;
  if os.path.isfile(logFileName) == False:
   print("WARNING " + logFileName + " not found in " + runFolder)
  else:
   # parse obstacle Uex, IMF components and speed, density and temperature of the first (sw?) population from the run log file
   Ueobs = parseParameter(logFileName,["grep \"Ue(r <= R_fieldObstacle) = (\" " + logFileName + " | cut -f 3 -d '=' |cut -f 2 -d '(' |cut -f 1 -d ','"],1e3)
   Bx = parseParameter(logFileName,["grep -A 3 \"(UPSTREAM IMF)\" " + logFileName + " |grep \"Bx  = \" |cut -f 2 -d '=' |cut -f 2 -d ' '"],1e-9)
   By = parseParameter(logFileName,["grep -A 3 \"(UPSTREAM IMF)\" " + logFileName + " |grep \"By  = \" |cut -f 2 -d '=' |cut -f 2 -d ' '"],1e-9)
   Bz = parseParameter(logFileName,["grep -A 3 \"(UPSTREAM IMF)\" " + logFileName + " |grep \"Bz  = \" |cut -f 2 -d '=' |cut -f 2 -d ' '"],1e-9)
   Usw = parseParameter(logFileName,["grep \") speed\" " + logFileName + " |cut -f 2 -d '=' |cut -f 2 -d ' '"],1e3)
   nsw = parseParameter(logFileName,["grep \") density\" " + logFileName + " |cut -f 2 -d '=' |cut -f 2 -d ' '"],1e6)
   Tsw = parseParameter(logFileName,["grep \") temperature\" " + logFileName + " |cut -f 2 -d '=' |cut -f 2 -d ' '"],1)
   Btot = np.sqrt(Bx*Bx + By*By + Bz*Bz)
   phi = 180 - np.rad2deg(np.arctan2(np.sqrt(By*By + Bz*Bz),Bx)) # IMF cone angle
   titleStr += "\n"
   titleStr += "$B = ($" + round2str(Bx/1e-9,1) + "," + round2str(By/1e-9,1) + "," + round2str(Bz/1e-9,1) + ") nT, "
   titleStr += "$|B|$ = " + round2str(Btot/1e-9,1) + " nT, $\phi$ = " + round2str(phi,1) + "$^\circ$" + "\n"
   titleStr += "$n,U,T$ = " + round2str(nsw/1e6,1) + " cm-3, " + round2str(Usw/1e3,1) + " km/s, " + round2str(Tsw/1e3,1) + " kK"

 # find out dimensionality of the run
 runDimAxes = ""
 runDim = -1
 if nx > 1:
  runDimAxes = "x"
 if ny > 1:
  runDimAxes = runDimAxes + "y"
 if nz > 1:
  runDimAxes = runDimAxes + "z"
 runDim = len(runDimAxes)
 if nx < 1 | ny < 1 | nz < 1:
  print(HN + "ERROR: bad run dimensions, all coordinate axes should have at least one cell")

 # column and row indices
 ii_col = -1
 ii_row = -1
 if runDim == 3:
  ii_col = 0
  ii_row = ii_run
 elif runDim == 2:
  ii_col = ii_run
  ii_row = 0
  # 2D multi-run: use several rows
  #ii_col = ii_run%NfigCols
  #ii_row = np.floor(ii_run/NfigCols).astype(int)
 else:
  print(HN + "ERROR: unsupported run dimensionality (runDim = " + str(runDim) + ")")
  quit()

 # zoomed domain
 global axisLimsZoom
 if len(np.atleast_1d(axisLimsZoom)) < 6:
  axisLimsZoom = axisLims

 # try to automatically decide good coordinate tick marks
 Dx = Dy = Dz = -1
 if nx > 1:
  Dx = axisLimsZoom[1] - axisLimsZoom[0]
 if ny > 1:
  Dy = axisLimsZoom[3] - axisLimsZoom[2]
 if nz > 1:
  Dz = axisLimsZoom[5] - axisLimsZoom[4]
 maxCrdRange = max((Dx,Dy,Dz))
 maxCrd = max(abs(axisLimsZoom)) # maximum absolute coordinate value
 # find maximum absolute coordinate value rounded upward to closest power of ten
 if maxCrd > 0:
  maxCrdTen = pow(10,(np.ceil(np.log10(maxCrd))))
 else:
  maxCrdTen = 1
 # try to find optimal tick step
 tickStep = 1
 for ii,jj,kk in itertools.product(range(1,100),(1,2,5),(-1,1)):
  tickStep =  maxCrdTen/(kk*jj*ii)
  NticksMax = maxCrdRange/tickStep
  if(NticksMax >= 5) & (NticksMax <= 15):
   break # finish loop if good tick step found
 if (NticksMax < 5) | (NticksMax > 15):
  print(HN + "WARNING: no good tick step found: NticksMax = " + str(NticksMax) + ", tickStep = " + str(tickStep))
 # tick marks
 crdTicks = np.arange(-maxCrdTen,+maxCrdTen,step=tickStep)

 # read and sort cell ids
 cellids_sorted = vr.read_variable("CellID").argsort()

 # choose which plane to plot from 3D array (nx,ny,nz)
 XZ_Ncol = XY_Ncol = YZ_Ncol = -1
 if runDim == 3:
  global xPlane, yPlane, zPlane
  YZ_Ncol,xPlane = chooseNcol(xPlane*Rp,xmin,xmax,nx,"x",Rp,Rp_str)
  XZ_Ncol,yPlane = chooseNcol(yPlane*Rp,ymin,ymax,ny,"y",Rp,Rp_str)
  XY_Ncol,zPlane = chooseNcol(zPlane*Rp,zmin,zmax,nz,"z",Rp,Rp_str)
 elif runDim == 2:
  XZ_Ncol = XY_Ncol = YZ_Ncol = -2
 else:
  print(HN + "ERROR: unsupported run dimensionality (runDim = " + str(runDim) + ")")
  quit()

 # function can be called with fig = -1 to get header information and the run dimensionality
 # print header only once
 if printHeader == 1:
  print("\nRun and figure configuration (from the first VLSV file):")
  print("Rp = " + str(Rp/1e3) + " km")
  print("run dimensions: " + str(runDim) + "D, " + runDimAxes)
  print("x = " + str(xmin/Rp) + " ... " + str(xmax/Rp) + " Rp")
  print("y = " + str(ymin/Rp) + " ... " + str(ymax/Rp) + " Rp")
  print("z = " + str(zmin/Rp) + " ... " + str(zmax/Rp) + " Rp")
  print("grid (nx,ny,nz) = (" + str(nx) + "," + str(ny) + "," + str(nz) + ")")
  if runDim == 3:
   print("Plotting a 3D run with:")
   print("xPlane = " + str(xPlane) + " Rp, YZ_Ncol = " + str(YZ_Ncol))
   print("yPlane = " + str(yPlane) + " Rp, XZ_Ncol = " + str(XZ_Ncol))
   print("zPlane = " + str(zPlane) + " Rp, XY_Ncol = " + str(XY_Ncol))
  elif runDim == 2:
   print("Plotting a 2D run: " + runDimAxes)
  else:
   print("ERROR: unsupported run dimensionality (runDim = " + str(runDim) + ")")
   quit()
  print("fig. columns x rows: " + str(NfigCols) + " x " + str(NfigRows) + "\n")

 # return the number of run dimensions
 if fig == -1:
  return runDim,runDimAxes

 # read and sort variable
 if(P_ii["type"] == "magnitude"):
  D_i = vr.read_variable_info(P_ii["param"])
  D = np.sqrt(D_i.data[:,0]**2 + D_i.data[:,1]**2 + D_i.data[:,2]**2)
 elif(P_ii["type"] == "scalar"):
  D_i = vr.read_variable_info(P_ii["param"])
  D = D_i.data
 elif(P_ii["type"] == "xcomp"):
  D_i = vr.read_variable_info(P_ii["param"])
  D = D_i.data[:,0]
 elif(P_ii["type"] == "ycomp"):
  D_i = vr.read_variable_info(P_ii["param"])
  D = D_i.data[:,1]
 elif(P_ii["type"] == "zcomp"):
  D_i = vr.read_variable_info(P_ii["param"])
  D = D_i.data[:,2]
 elif(P_ii["type"] == "nvO"): # example of a custom parameter type
  nO_i = vr.read_variable_info("n_O+_ave")
  VO_i = vr.read_variable_info("v_O+_ave")
  VtotO = np.sqrt(VO_i.data[:,0]**2 + VO_i.data[:,1]**2 + VO_i.data[:,2]**2)
  D = nO_i.data * VtotO
 else:
  print(HN + "ERROR: unknown parameter type: " + P_ii["type"])
  quit()

 # reshape variable into 3D array
 D = D[cellids_sorted].reshape(nz,ny,nx)

 #plt.figure(figsize=figureSize,frameon=False)

 # plot three columns per 3D run
 if runDim == 3:
  # get 2D slices to plot
  meshD_yz = D[:,:,YZ_Ncol]
  meshD_xz = D[:,XZ_Ncol,:]
  meshD_xy = D[XY_Ncol,:,:]
  # apply smoothing
  if P_ii["sigma"] > 0:
   meshD_yz = sp.ndimage.filters.gaussian_filter(D[:,:,YZ_Ncol],sigma=P_ii["sigma"],mode="constant")
   meshD_xz = sp.ndimage.filters.gaussian_filter(D[:,XZ_Ncol,:],sigma=P_ii["sigma"],mode="constant")
   meshD_xy = sp.ndimage.filters.gaussian_filter(D[XY_Ncol,:,:],sigma=P_ii["sigma"],mode="constant")
  #plt.subplot(1,3,1)
  a = axes[ii_row][ii_col].imshow(meshD_xz/P_ii["unit"],vmin=P_ii["lims"][0]/P_ii["unit"],vmax=P_ii["lims"][1]/P_ii["unit"],cmap=P_ii["colormap"],extent=[axisLims[0],axisLims[1],axisLims[4],axisLims[5]],aspect="equal",origin="lower",interpolation="nearest")
  if ii_row == 0:
   axes[ii_row][ii_col].title.set_text("3D: $xz$ ($y=$" + round2strWithUnit(yPlane,Rp_str) + ")")
  if ii_row == (NfigRows-1):
   axes[ii_row][ii_col].set_xlabel("$x$ [" + Rp_str + "]")
  #axes[ii_row][ii_col].set_ylabel(titleStr + "\n\n$z$ [" + Rp_str + "]")
  axes[ii_row][ii_col].set_ylabel("$z$ [" + Rp_str + "]")
  axes[ii_row][ii_col].set_xticks(crdTicks)
  axes[ii_row][ii_col].tick_params("x",labelrotation=xTickAngle)
  axes[ii_row][ii_col].set_yticks(crdTicks)
  axes[ii_row][ii_col].axis("scaled")
  axes[ii_row][ii_col].set_xlim(axisLimsZoom[0:2])
  axes[ii_row][ii_col].set_ylim(axisLimsZoom[4:6])
  configurePanel(P_ii,a,axes[ii_row][ii_col],showPlanet[0],0)
  ii_col += 1

  #plt.subplot(1,3,2)
  a = axes[ii_row][ii_col].imshow(meshD_xy/P_ii["unit"],vmin=P_ii["lims"][0]/P_ii["unit"],vmax=P_ii["lims"][1]/P_ii["unit"],cmap=P_ii["colormap"],extent=[axisLims[0],axisLims[1],axisLims[2],axisLims[3]],aspect="equal",origin="lower",interpolation="nearest")
  if ii_row == 0:
   axes[ii_row][ii_col].title.set_text("$t=$" + round2str(fileTime) + " s\n" "3D: $xy$ ($z=$" + round2strWithUnit(zPlane,Rp_str) + ")")
  if ii_row == (NfigRows-1):
   axes[ii_row][ii_col].set_xlabel("$x$ [" + Rp_str + "]")
  axes[ii_row][ii_col].set_ylabel("$y$ [" + Rp_str + "]")
  axes[ii_row][ii_col].set_xticks(crdTicks)
  axes[ii_row][ii_col].tick_params("x",labelrotation=xTickAngle)
  axes[ii_row][ii_col].set_yticks(crdTicks)
  axes[ii_row][ii_col].axis("scaled")
  axes[ii_row][ii_col].set_xlim(axisLimsZoom[0:2])
  axes[ii_row][ii_col].set_ylim(axisLimsZoom[2:4])
  configurePanel(P_ii,a,axes[ii_row][ii_col],showPlanet[1],0)
  ii_col += 1

  #plt.subplot(1,3,3)
  a = axes[ii_row][ii_col].imshow(meshD_yz/P_ii["unit"],vmin=P_ii["lims"][0]/P_ii["unit"],vmax=P_ii["lims"][1]/P_ii["unit"],cmap=P_ii["colormap"],extent=[axisLims[2],axisLims[3],axisLims[4],axisLims[5]],aspect="equal",origin="lower",interpolation="nearest")
  # first row
  if ii_row == 0:
   axes[ii_row][ii_col].title.set_text("3D: $yz$ ($x=$" + round2strWithUnit(xPlane,Rp_str) + ")")
  # the bottom row
  if ii_row == (NfigRows-1):
   axes[ii_row][ii_col].set_xlabel("$y$ [" + Rp_str + "]")
  axes[ii_row][ii_col].set_ylabel("$z$ [" + Rp_str + "]")
  axes[ii_row][ii_col].set_xticks(crdTicks)
  axes[ii_row][ii_col].tick_params("x",labelrotation=xTickAngle)
  axes[ii_row][ii_col].set_yticks(crdTicks)
  axes[ii_row][ii_col].axis("scaled")
  axes[ii_row][ii_col].set_xlim(axisLimsZoom[2:4])
  axes[ii_row][ii_col].set_ylim(axisLimsZoom[4:6])
  configurePanel(P_ii,a,axes[ii_row][ii_col],showPlanet[2],0)
  ii_col += 1

  # the bottom row
  if ii_row == (NfigRows-1):
   fig.tight_layout()
   if NfigRows == 1:
    clb = fig.colorbar(a,ax=axes.flatten(),shrink=0.5)
   else:
    clb = fig.colorbar(a,ax=axes.flatten())
   clb.ax.set_title(P_ii["str"])
   #plt.savefig(outputFileNamePrefix + P_ii["filename"] + "_" + vlsvFileName + ".png",dpi=figDpi,transparent=False) #,bbox_inches="tight")
   #plt.close()
   #plt.clf()

 # plot one column per 2D run
 elif runDim == 2:
  meshD = -1
  if runDimAxes == "yz":
   meshD = D[:,:,0]
   axisExtend = [axisLims[2],axisLims[3],axisLims[4],axisLims[5]]
   axisLimsZoomX = axisLimsZoom[2:4]
   axisLimsZoomY = axisLimsZoom[4:6]
   xlabelStr = "$y$"
   ylabelStr = "$z$"
  elif runDimAxes == "xz":
   meshD = D[:,0,:]
   axisExtend = [axisLims[0],axisLims[1],axisLims[4],axisLims[5]]
   axisLimsZoomX = axisLimsZoom[0:2]
   axisLimsZoomY = axisLimsZoom[4:6]
   xlabelStr = "$x$"
   ylabelStr = "$z$"
  elif runDimAxes == "xy":
   meshD = D[0,:,:]
   axisExtend = [axisLims[0],axisLims[1],axisLims[2],axisLims[3]]
   axisLimsZoomX = axisLimsZoom[0:2]
   axisLimsZoomY = axisLimsZoom[2:4]
   xlabelStr = "$x$"
   ylabelStr = "$y$"
  else:
   print(HN + "ERROR: unknown run dimensionality (runDimAxes = " + runDimAxes + ")")
   quit()
  #plt.subplot(1,3,2)
  a = axes[ii_row][ii_col].imshow(meshD/P_ii["unit"],vmin=P_ii["lims"][0]/P_ii["unit"],vmax=P_ii["lims"][1]/P_ii["unit"],cmap=P_ii["colormap"],extent=axisExtend,aspect="equal",origin="lower",interpolation="nearest")
  if ii_col == 0:
   axes[ii_row][ii_col].set_ylabel(ylabelStr + " [" + Rp_str + "]")
  if ii_col == 0 and ii_row == 0:
   axes[ii_row][ii_col].title.set_text("$t=$" + round2str(fileTime) + " s, " + str(runDim) + "D: " + runDimAxes + "\n" + titleStr)
  else:
   axes[ii_row][ii_col].title.set_text(titleStr)
  # 2D multi-run: use several row
  #if ii_row == (NfigRows-1):
  axes[ii_row][ii_col].set_xlabel(xlabelStr + " [" + Rp_str + "]")
  axes[ii_row][ii_col].set_xticks(crdTicks)
  #else:
  axes[ii_row][ii_col].set_xticks(crdTicks) #,labels=""
  axes[ii_row][ii_col].tick_params("x",labelrotation=xTickAngle)
  # 2D multi-run:
  #if ii_col == 0:
  # axes[ii_row][ii_col].set_yticks(crdTicks)
  #else:
  # axes[ii_row][ii_col].set_yticks(crdTicks,labels="")
  axes[ii_row][ii_col].axis("scaled")
  axes[ii_row][ii_col].set_xlim(axisLimsZoomX)
  axes[ii_row][ii_col].set_ylim(axisLimsZoomY)
  configurePanel(P_ii,a,axes[ii_row][ii_col],showPlanet[1],0)
  # 2D multi-run: use several rows
  #if ii_col == (NfigCols-1) and ii_row == (NfigRows-1):
  # the last column
  if ii_col == (NfigCols-1):
   fig.tight_layout()
   clb = fig.colorbar(a,ax=axes.flatten(),shrink=0.5)
   clb.ax.set_title(P_ii["str"])

# plot a figure
def plotFigure(vlsvFileNameFullPath,printHeader=0):
 for ii in range(len(P)):
  if printHeader == 0:
   print(HN + "parameter: " + P[ii]["param"] + ", " + P[ii]["type"] + " (" + str(ii+1) + "/" + str(Nparams) + ")")
  fig,axes = plt.subplots(nrows=NfigRows,ncols=NfigCols,figsize=figureSize,frameon=True,squeeze=False)
  jj = 0
  for s in simRuns:
   runFolder = s[0]
   runStr = s[1]
   Rp = s[2]
   Rp_str = s[3]
   vlsvFileName = os.path.basename(vlsvFileNameFullPath)
   plotPanel(fig,axes,jj,P[ii],runFolder,vlsvFileName,runStr,Rp,Rp_str,printHeader)
   # if only priting header, return after reading the first parameter of the first run
   if printHeader == 1:
    return
   jj += 1
  plt.savefig(outputFileNamePrefix + P[ii]["filename"] + "_" + vlsvFileName + ".png",dpi=figDpi,transparent=False) #,bbox_inches="tight")
  plt.clf()
  plt.close()

# find run files
if "runFiles" not in locals():
 runFolder1 = simRuns[0][0]
 print(HN + "processing from: folder = " + runFolder1 + ", time step = " + str(tStartThisProcess) + " ... " + str(tEndThisProcess) + " (this process), whole run = " + str(tStartGlobal) + " ... " + str(tEndGlobal))
 runFiles = []
 for f in sorted(os.listdir(runFolder1)):
  if (f.startswith("state") and f.endswith(".vlsv")):
   try:
    t = int(f[5:13])
   except ValueError:
    print(HN + "ERROR: could not parse int from file name" + f)
    continue
   if (t >= tStartThisProcess) and (t <= tEndThisProcess):
    runFiles.append(runFolder1 + f)

# total number of files (save steps) found to be plotted
Nfiles = len(runFiles)
print(HN + "files found: " + str(Nfiles))
if Nfiles <= 0:
 print(HN + "ERROR: no VLSV files found")
 quit()

# find all variable from the first vlsv file
vr = pt.vlsvfile.VlsvReader(runFiles[0])
all_vars = vr.get_all_variables()
N_vector_vars_found = 0
N_scalar_vars_found = 0
for var_ in all_vars:
 if "flag" in var_ or "detector_" in var_ or "CellID" in var_:
  print(HN + "skipping flag, detector and CellID varibles: " + var_)
  continue
 D_i = vr.read_variable_info(var_)
 isVector = False
 # determined if variable is vector
 if len(D_i.data.shape) > 1:
  if D_i.data.shape[1] == 3:
   isVector = True
   N_vector_vars_found = N_vector_vars_found + 1
  else:
   print(HN + "ERROR: vector size not three: " + str(D_i.data.shape[1]) + " (" + var_ + ")")
   exit()
 else:
  N_scalar_vars_found = N_scalar_vars_found + 1
 var_found = False
 for ii in range(len(P_settings)):
  if P_settings[ii]["param"] == var_:
   var_found = True
   break
 if var_found == False:
  print(HN + "skipping variable not found in P_settings: " + var_)
  continue
 # prepare vector variables
 if isVector == True:
  for ii in range(len(P_settings)):
   if P_settings[ii]["param"] == var_ and P_settings[ii]["type"] == "magnitude":
    P.append({"param":var_,"type":"magnitude","str":P_settings[ii]["str"],"log":P_settings[ii]["log"],"lims":P_settings[ii]["lims"],"unit":P_settings[ii]["unit"],"colormap":P_settings[ii]["colormap"],"filename":P_settings[ii]["filename"],"sigma":P_settings[ii]["sigma"]})
    continue
   elif P_settings[ii]["param"] == var_ and P_settings[ii]["type"] == "xcomp":
    P.append({"param":var_,"type":"xcomp","str":P_settings[ii]["str"],"log":P_settings[ii]["log"],"lims":P_settings[ii]["lims"],"unit":P_settings[ii]["unit"],"colormap":P_settings[ii]["colormap"],"filename":P_settings[ii]["filename"],"sigma":P_settings[ii]["sigma"]})
    continue
   elif P_settings[ii]["param"] == var_ and P_settings[ii]["type"] == "ycomp":
    P.append({"param":var_,"type":"ycomp","str":P_settings[ii]["str"],"log":P_settings[ii]["log"],"lims":P_settings[ii]["lims"],"unit":P_settings[ii]["unit"],"colormap":P_settings[ii]["colormap"],"filename":P_settings[ii]["filename"],"sigma":P_settings[ii]["sigma"]})
    continue
   elif P_settings[ii]["param"] == var_ and P_settings[ii]["type"] == "zcomp":
    P.append({"param":var_,"type":"zcomp","str":P_settings[ii]["str"],"log":P_settings[ii]["log"],"lims":P_settings[ii]["lims"],"unit":P_settings[ii]["unit"],"colormap":P_settings[ii]["colormap"],"filename":P_settings[ii]["filename"],"sigma":P_settings[ii]["sigma"]})
    continue
 # prepare scalar variables
 else:
  for ii in range(len(P_settings)):
   if P_settings[ii]["param"] == var_ and P_settings[ii]["type"] == "scalar":
    P.append({"param":var_,"type":"scalar","str":P_settings[ii]["str"],"log":P_settings[ii]["log"],"lims":P_settings[ii]["lims"],"unit":P_settings[ii]["unit"],"colormap":P_settings[ii]["colormap"],"filename":P_settings[ii]["filename"],"sigma":P_settings[ii]["sigma"]})
    break

# total number of parameters to be plotted (scalar + 3*vector)
Nparams = len(P)
NfileOpens = Nruns*Nfiles*Nparams

if printHeaderOnly == 1:
 print("Plotting variables:")
 for Pii_ in P:
  print(Pii_)
 print("Vector and scalar variables found: " + str(len(all_vars)))
 print("Vector variables found: " + str(N_vector_vars_found))
 print("Scalar variables found: " + str(N_scalar_vars_found))
 print("Plotting fields found: " + str(Nparams))
 print("Runs found: " + str(Nruns))
 print("Total file openings expected: " + str(NfileOpens))

# read header information from the first file of the first run
runDim,runDimAxes = plotPanel(-1,-1,-1,-1,simRuns[0][0],os.path.basename(runFiles[0]),"run0",simRuns[0][2],-1)

# number of figure columns and rows
if runDim == 3:
 NfigCols = 3
 NfigRows = Nruns
elif runDim == 2:
 NfigCols = Nruns
 NfigRows = 1
 # 2D multi-run: use several rows
 #NfigCols = 6
 #NfigRows = 3
else:
 print(HN + "(plotFigure) ERROR: unsupported run dimensionality (runDim = " + str(runDim) + ")")
 quit()

# run plotting with multiple processes
if (useMultiProcessing > 0) & (printHeaderOnly == 0):
 if __name__ == "__main__":
  pool = Pool(Ncores)
  args = [(f,) for f in runFiles]
  pool.starmap(plotFigure,args)
elif (printHeaderOnly == 1):
 plotFigure(runFiles[0],1)
else:
 # without multiple processes (for debugging) or when getting only header information
 print(HN + "DEBUG MODE: multiprocessing not used, serial execution")
 for ff in runFiles:
  plotFigure(ff)

