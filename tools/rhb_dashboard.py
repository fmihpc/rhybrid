# This file is part of the RHybrid simulation platform.
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

import sys
import os
import matplotlib
matplotlib.use("Agg") # this is needed for offline saving of graphics
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import scipy as sp
import operator as oper
import socket
import itertools
plt.switch_backend("agg")

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

runFolder = "./"

# check log files found

if not any( (fileName.startswith("pop") and fileName.endswith(".log")) for fileName in os.listdir(runFolder)):
 print("ERROR: no pop*.log found in folder: " + runFolder)
 quit()
if not any( (fileName.startswith("field") and fileName.endswith(".log")) for fileName in os.listdir(runFolder)):
 print("ERROR: no field.log found in folder: " + runFolder)
 quit()

#a = np.loadtxt( open(runFolder + "field.log","rt").readlines()[:-1],skiprows=0,comments="%",dtype=None)
#t = a[:,0]
