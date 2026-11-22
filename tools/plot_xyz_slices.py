#!/usr/bin/env python3
"""Plot xz, xy and yz slices of a single 3D RHybrid simulation run.

The Slurm job selects state files and writes a NUL-separated manifest.
Usage:
    python plotter_rhybrid_3d_slice_single_node.py Ncores runFolder Robject \
        --file-list FILE --output-dir DIR

Ncores = workers (use -1 to print header information only)
Robject = object radius in metres
FILE = NUL-separated list of paths to selected stateNNN.vlsv files

Each selected state file produces one PNG per active parameter
setting, with three slice panels in each PNG.
"""

import argparse
import itertools
import os
import re
from functools import partial
from multiprocessing import Pool
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np
import rhybridpy as rhb
from matplotlib.patches import Wedge
from scipy.ndimage import gaussian_filter

# Plot presentation: unchanged from the original three-panel 3D plots.
matplotlib.rcParams.update({"font.size": 14})
matplotlib.rcParams["lines.linewidth"] = 2

FIGURE_DPI = 100
FIGURE_WIDTH = 1920
FIGURE_HEIGHT = 1080
FIGURE_SIZE = (FIGURE_WIDTH / FIGURE_DPI, FIGURE_HEIGHT / FIGURE_DPI)
matplotlib.rcParams["figure.dpi"] = FIGURE_DPI
matplotlib.rcParams["savefig.dpi"] = FIGURE_DPI

TICK_DIRECTION = "out"
TICK_LENGTH = 2
TICK_WIDTH = 1
X_TICK_ANGLE = 45

# 3D slice positions in object radii, at constant x, y and z, respectively.
X_PLANE = -2.0  # yz plane
Y_PLANE = 0.0   # xz plane
Z_PLANE = 0.0   # xy plane

# -1: full domain; or [xmin, xmax, ymin, ymax, zmin, zmax] in object radii.
AXIS_LIMITS_ZOOM = -1
# AXIS_LIMITS_ZOOM = np.multiply([-4.0, 4.0, -4.0, 4.0, -4.0, 4.0], 1.0)
SHOW_PLANET = (1, 1, 0)  # xz, xy, yz panels

CMAP_SEQUENTIAL = "inferno"  # perceptually uniform colormap
CMAP_DIVERGING = "seismic"   # blue-white-red colormap

# Physical constants (SI units).
PROTON_MASS = 1.6726217160e-27
HELIUM_MASS = 6.6464761621e-27
OXYGEN_MASS = 2.6567625437e-26
OXYGEN_MOLECULE_MASS = 5.3135250873e-26
MU_0 = 1.25663706e-6

# Only variables enabled in PLOT_SETTINGS and present in a VLSV file are plotted.
# Each dictionary specifies: param (VLSV name), type, str (plot label), log,
# lims (color limits), unit (scale factor), colormap, filename, and sigma
# (Gaussian smoothing; -1 disables it).
PLOT_SETTINGS = []

# Instantaneous parameters

PLOT_SETTINGS.append({
    "param": "cellB",
    "type": "magnitude",
    "str": "$|cellB|$ [nT]",
    "log": 1,
    "lims": (1e-9, 1000e-9),
    "unit": 1e-9,
    "colormap": CMAP_SEQUENTIAL,
    "filename": "B",
    "sigma": -1,
})
# PLOT_SETTINGS.append({"param": "cellB", "type": "xcomp", "str": "$cellB_x$ [nT]", "log": 0, "lims": (-100e-9,100e-9), "unit": 1e-9, "colormap": CMAP_DIVERGING, "filename": "Bx", "sigma": -1})
# PLOT_SETTINGS.append({"param": "cellB", "type": "ycomp", "str": "$cellB_y$ [nT]", "log": 0, "lims": (-100e-9,100e-9), "unit": 1e-9, "colormap": CMAP_DIVERGING, "filename": "By", "sigma": -1})
# PLOT_SETTINGS.append({"param": "cellB", "type": "zcomp", "str": "$cellB_z$ [nT]", "log": 0, "lims": (-100e-9,100e-9), "unit": 1e-9, "colormap": CMAP_DIVERGING, "filename": "Bz", "sigma": -1})
# PLOT_SETTINGS.append({"param": "faceB", "type": "magnitude", "str": "$|faceB|$ [nT]", "log": 1, "lims": (1e-9, 1000e-9), "unit": 1e-9, "colormap": CMAP_SEQUENTIAL, "filename": "faceB", "sigma": -1})
# PLOT_SETTINGS.append({"param": "faceB", "type": "xcomp", "str": "$faceB_x$ [nT]", "log": 0, "lims": (-100e-9,100e-9), "unit": 1e-9, "colormap": CMAP_DIVERGING, "filename": "faceBx", "sigma": -1})
# PLOT_SETTINGS.append({"param": "faceB", "type": "ycomp", "str": "$faceB_y$ [nT]", "log": 0, "lims": (-100e-9,100e-9), "unit": 1e-9, "colormap": CMAP_DIVERGING, "filename": "faceBy", "sigma": -1})
# PLOT_SETTINGS.append({"param": "faceB", "type": "zcomp", "str": "$faceB_z$ [nT]", "log": 0, "lims": (-100e-9,100e-9), "unit": 1e-9, "colormap": CMAP_DIVERGING, "filename": "faceBz", "sigma": -1})
PLOT_SETTINGS.append({
    "param": "n_H+sw",
    "type": "scalar",
    "str": "$n$(H$^{+}_\\mathrm{sw}$) [m$^{-3}$]",
    "log": 1,
    "lims": (1e4, 1e9),
    "unit": 1,
    "colormap": CMAP_SEQUENTIAL,
    "filename": "Hsw_n",
    "sigma": -1,
})
# PLOT_SETTINGS.append({"param": "n_He++sw", "type": "scalar", "str": "$n$(He$^{++}_\mathrm{sw}$) [m$^{-3}$]", "log": 1, "lims": (0.04*1e4,0.04*1e9), "unit": 1, "colormap": CMAP_SEQUENTIAL, "filename": "Hesw_n", "sigma": -1})
# PLOT_SETTINGS.append({"param": "n_Na+", "type": "scalar", "str": "$n$(Na$^{+}$) [m$^{-3}$]", "log": 1, "lims": (10,1e10), "unit": 1, "colormap": CMAP_SEQUENTIAL, "filename": "Na_n", "sigma": -1})
# PLOT_SETTINGS.append({"param": "n_O+", "type": "scalar", "str": "$n$(O$^{+}$) [m$^{-3}$]", "log": 1, "lims": (1,1e7), "unit": 1, "colormap": CMAP_SEQUENTIAL, "filename": "O_n", "sigma": -1})
# PLOT_SETTINGS.append({"param": "n_O2+", "type": "scalar", "str": "$n$(O$^{+}_\mathrm{2}$) [m$^{-3}$]", "log": 1, "lims": (1,1e7), "unit": 1, "colormap": CMAP_SEQUENTIAL, "filename": "O2_n", "sigma": -1})
# PLOT_SETTINGS.append({"param": "n_H+planet", "type": "scalar", "str": "$n$(H$^{+}_\mathrm{planet}$) [m$^{-3}$]", "log": 1, "lims": (1,1e7), "unit": 1, "colormap": CMAP_SEQUENTIAL, "filename": "Hplanet_n", "sigma": -1})
# PLOT_SETTINGS.append({"param": "nodeE", "type": "magnitude", "str": "$|E|$ [V/m]", "log": 1, "lims": (0.0001,0.1), "unit": 1, "colormap": CMAP_SEQUENTIAL, "filename": "nodeE", "sigma": -1})
# PLOT_SETTINGS.append({"param": "nodeE", "type": "xcomp", "str": "$E_x$ [mV/m]", "log": 0, "lims": (-0.02,0.02), "unit": 1e-3, "colormap": CMAP_DIVERGING, "filename": "nodeEx", "sigma": -1})
# PLOT_SETTINGS.append({"param": "nodeE", "type": "ycomp", "str": "$E_y$ [mV/m]", "log": 0, "lims": (-0.02,0.02), "unit": 1e-3, "colormap": CMAP_DIVERGING, "filename": "nodeEy", "sigma": -1})
# PLOT_SETTINGS.append({"param": "nodeE", "type": "zcomp", "str": "$E_z$ [mV/m]", "log": 0, "lims": (-0.02,0.02), "unit": 1e-3, "colormap": CMAP_DIVERGING, "filename": "nodeEz", "sigma": -1})
# PLOT_SETTINGS.append({"param": "cellEp", "type": "magnitude", "str": "$|E_p|$ [V/m]", "log": 1, "lims": (1e-6,1e-3), "unit": 1, "colormap": CMAP_SEQUENTIAL, "filename": "Ep", "sigma": -1})
# PLOT_SETTINGS.append({"param": "cellEp", "type": "xcomp", "str": "$E_{p,x}$ [V/m]", "log": 0, "lims": (-2e-04,2e-04), "unit": 1, "colormap": CMAP_DIVERGING, "filename": "Epx", "sigma": -1})
# PLOT_SETTINGS.append({"param": "cellEp", "type": "ycomp", "str": "$E_{p,y}$ [V/m]", "log": 0, "lims": (-2e-04,2e-04), "unit": 1, "colormap": CMAP_DIVERGING, "filename": "Epy", "sigma": -1})
# PLOT_SETTINGS.append({"param": "cellEp", "type": "zcomp", "str": "$E_{p,z}$ [V/m]", "log": 0, "lims": (-2e-04,2e-04), "unit": 1, "colormap": CMAP_DIVERGING, "filename": "Epz", "sigma": -1})
# PLOT_SETTINGS.append({"param": "nodeJ", "type": "magnitude", "str": "$|J|$ [A/m$^2$]", "log": 1, "lims": (1e-9,1e-6), "unit": 1, "colormap": CMAP_SEQUENTIAL, "filename": "nodeJ", "sigma": -1})
# PLOT_SETTINGS.append({"param": "nodeJ", "type": "xcomp", "str": "$J_x$ [A/m$^2$]", "log": 0, "lims": (-3e-7,3e-7), "unit": 1, "colormap": CMAP_DIVERGING, "filename": "nodeJx", "sigma": -1})
# PLOT_SETTINGS.append({"param": "nodeJ", "type": "ycomp", "str": "$J_y$ [A/m$^2$]", "log": 0, "lims": (-3e-7,3e-7), "unit": 1, "colormap": CMAP_DIVERGING, "filename": "nodeJy", "sigma": -1})
# PLOT_SETTINGS.append({"param": "nodeJ", "type": "zcomp", "str": "$J_z$ [A/m$^2$]", "log": 0, "lims": (-3e-7,3e-7), "unit": 1, "colormap": CMAP_DIVERGING, "filename": "nodeJz", "sigma": -1})
# PLOT_SETTINGS.append({"param": "cellUe", "type": "magnitude", "str": "$|U|(e^-)$ [km/s]", "log": 0, "lims": (0,700e3), "unit": 1e3, "colormap": CMAP_SEQUENTIAL, "filename": "e_U", "sigma": -1})
# PLOT_SETTINGS.append({"param": "cellUe", "type": "xcomp", "str": "$U_x(e^-)$ [km/s]", "log": 0, "lims": (-600e3,600e3), "unit": 1e3, "colormap": CMAP_DIVERGING, "filename": "e_Ux", "sigma": -1})
# PLOT_SETTINGS.append({"param": "cellUe", "type": "ycomp", "str": "$U_y(e^-)$ [km/s]", "log": 0, "lims": (-600e3,600e3), "unit": 1e3, "colormap": CMAP_DIVERGING, "filename": "e_Uy", "sigma": -1})
# PLOT_SETTINGS.append({"param": "cellUe", "type": "zcomp", "str": "$U_z(e^-)$ [km/s]", "log": 0, "lims": (-600e3,600e3), "unit": 1e3, "colormap": CMAP_DIVERGING, "filename": "e_Uz", "sigma": -1})
# PLOT_SETTINGS.append({"param": "v_H+sw", "type": "magnitude", "str": "$|U|$(H$^{+}_\mathrm{sw}$) [km/s]", "log": 0, "lims": (0,700e3), "unit": 1e3, "colormap": CMAP_SEQUENTIAL, "filename": "Hsw_U", "sigma": -1})
# PLOT_SETTINGS.append({"param": "v_H+sw", "type": "xcomp", "str": "$U_x$(H$^{+}_\mathrm{sw}$) [km/s]", "log": 0, "lims": (-600e3,600e3), "unit": 1e3, "colormap": CMAP_DIVERGING, "filename": "Hsw_Ux", "sigma": -1})
# PLOT_SETTINGS.append({"param": "v_H+sw", "type": "ycomp", "str": "$U_y$(H$^{+}_\mathrm{sw}$) [km/s]", "log": 0, "lims": (-600e3,600e3), "unit": 1e3, "colormap": CMAP_DIVERGING, "filename": "Hsw_Uy", "sigma": -1})
# PLOT_SETTINGS.append({"param": "v_H+sw", "type": "zcomp", "str": "$U_z$(H$^{+}_\mathrm{sw}$) [km/s]", "log": 0, "lims": (-600e3,600e3), "unit": 1e3, "colormap": CMAP_DIVERGING, "filename": "Hsw_Uz", "sigma": -1})
# PLOT_SETTINGS.append({"param": "v_He++sw", "type": "magnitude", "str": "$|U|$(He$^{++}_\mathrm{sw}$) [km/s]", "log": 0, "lims": (0,700e3), "unit": 1e3, "colormap": CMAP_SEQUENTIAL, "filename": "Hesw_U", "sigma": -1})
# PLOT_SETTINGS.append({"param": "v_He++sw", "type": "xcomp", "str": "$U_x$(He$^{++}_\mathrm{sw}$) [km/s]", "log": 0, "lims": (-600e3,600e3), "unit": 1e3, "colormap": CMAP_DIVERGING, "filename": "Hesw_Ux", "sigma": -1})
# PLOT_SETTINGS.append({"param": "v_He++sw", "type": "ycomp", "str": "$U_y$(He$^{++}_\mathrm{sw}$) [km/s]", "log": 0, "lims": (-600e3,600e3), "unit": 1e3, "colormap": CMAP_DIVERGING, "filename": "Hesw_Uy", "sigma": -1})
# PLOT_SETTINGS.append({"param": "v_He++sw", "type": "zcomp", "str": "$U_z$(He$^{++}_\mathrm{sw}$) [km/s]", "log": 0, "lims": (-600e3,600e3), "unit": 1e3, "colormap": CMAP_DIVERGING, "filename": "Hesw_Uz", "sigma": -1})
# PLOT_SETTINGS.append({"param": "v_Na+", "type": "magnitude", "str": "$|U|$(Na$^{+}$) [km/s]", "log": 0, "lims": (0,1000e3), "unit": 1e3, "colormap": CMAP_SEQUENTIAL, "filename": "Na_U", "sigma": -1})
# PLOT_SETTINGS.append({"param": "v_Na+", "type": "xcomp", "str": "$U_x$(Na$^{+}$) [km/s]", "log": 0, "lims": (-600e3,600e3), "unit": 1e3, "colormap": CMAP_DIVERGING, "filename": "Na_Ux", "sigma": -1})
# PLOT_SETTINGS.append({"param": "v_Na+", "type": "ycomp", "str": "$U_y$(Na$^{+}$) [km/s]", "log": 0, "lims": (-600e3,600e3), "unit": 1e3, "colormap": CMAP_DIVERGING, "filename": "Na_Uy", "sigma": -1})
# PLOT_SETTINGS.append({"param": "v_Na+", "type": "zcomp", "str": "$U_z$(Na$^{+}$) [km/s]", "log": 0, "lims": (-600e3,600e3), "unit": 1e3, "colormap": CMAP_DIVERGING, "filename": "Na_Uz", "sigma": -1})
# PLOT_SETTINGS.append({"param": "T_H+sw", "type": "scalar", "str": "$T$(H$^{+}_\mathrm{sw}$) [K]", "log": 1, "lims": (1e5,1e8), "unit": 1, "colormap": CMAP_SEQUENTIAL, "filename": "Hsw_T", "sigma": -1})
# PLOT_SETTINGS.append({"param": "T_He++sw", "type": "scalar", "str": "$T$(He$^{++}_\mathrm{sw}$) [K]", "log": 1, "lims": (3.5*1e5,3.5*1e8), "unit": 1, "colormap": CMAP_SEQUENTIAL, "filename": "Hesw_T", "sigma": -1})
# PLOT_SETTINGS.append({"param": "T_Na+", "type": "scalar", "str": "$T$(Na$^{+}$) [K]", "log": 1, "lims": (1e4,1e7), "unit": 1, "colormap": CMAP_SEQUENTIAL, "filename": "Na_T", "sigma": -1})
# PLOT_SETTINGS.append({"param": "T_O+", "type": "scalar", "str": "$T$(O$^{+}$) [K]", "log": 1, "lims": (1e4,1e7), "unit": 1, "colormap": CMAP_SEQUENTIAL, "filename": "O_T", "sigma": -1})
# PLOT_SETTINGS.append({"param": "T_O2+", "type": "scalar", "str": "$T$(O$^{+}_\mathrm{2}$) [K]", "log": 1, "lims": (1e4,1e7), "unit": 1, "colormap": CMAP_SEQUENTIAL, "filename": "O2_T", "sigma": -1})
# PLOT_SETTINGS.append({"param": "T_H+planet", "type": "scalar", "str": "$T$(H$^{+}_\mathrm{planet}$) [K]", "log": 1, "lims": (1e5,10e6), "unit": 1, "colormap": CMAP_SEQUENTIAL, "filename": "Hplanet_T", "sigma": -1})

# Temporally averaged parameters

PLOT_SETTINGS.append({
    "param": "cellBAverage",
    "type": "magnitude",
    "str": "$|B|$ [nT]",
    "log": 1,
    "lims": (1e-9, 1000e-9),
    "unit": 1e-9,
    "colormap": CMAP_SEQUENTIAL,
    "filename": "B_ave_",
    "sigma": -1,
})
# PLOT_SETTINGS.append({"param": "cellBAverage", "type": "xcomp", "str": "$B_x$ [nT]", "log": 0, "lims": (-100e-9,100e-9), "unit": 1e-9, "colormap": CMAP_DIVERGING, "filename": "B_ave_x", "sigma": -1})
# PLOT_SETTINGS.append({"param": "cellBAverage", "type": "ycomp", "str": "$B_y$ [nT]", "log": 0, "lims": (-100e-9,100e-9), "unit": 1e-9, "colormap": CMAP_DIVERGING, "filename": "B_ave_y", "sigma": -1})
# PLOT_SETTINGS.append({"param": "cellBAverage", "type": "zcomp", "str": "$B_z$ [nT]", "log": 0, "lims": (-100e-9,100e-9), "unit": 1e-9, "colormap": CMAP_DIVERGING, "filename": "B_ave_z", "sigma": -1})
PLOT_SETTINGS.append({
    "param": "n_H+sw_ave",
    "type": "scalar",
    "str": "$n$(H$^{+}_\\mathrm{sw,ave}$) [m$^{-3}$]",
    "log": 1,
    "lims": (1e4, 1e9),
    "unit": 1,
    "colormap": CMAP_SEQUENTIAL,
    "filename": "Hsw_ave_n",
    "sigma": -1,
})
# PLOT_SETTINGS.append({"param": "n_He++sw_ave", "type": "scalar", "str": "$n$(He$^{++}_\mathrm{sw,ave}$) [m$^{-3}$]", "log": 1, "lims": (0.04*1e4,0.04*1e9), "unit": 1, "colormap": CMAP_SEQUENTIAL, "filename": "Hesw_ave_n", "sigma": -1})
# PLOT_SETTINGS.append({"param": "n_Na+_ave", "type": "scalar", "str": "$n$(Na$^{+}_\mathrm{ave}$) [m$^{-3}$]", "log": 1, "lims": (10,1e10), "unit": 1, "colormap": CMAP_SEQUENTIAL, "filename": "Na_ave_n", "sigma": -1})
# PLOT_SETTINGS.append({"param": "n_O+_ave", "type": "scalar", "str": "$n$(O$^{+}_\mathrm{sw,ave}$) [m$^{-3}$]", "log": 1, "lims": (1,1e7), "unit": 1, "colormap": CMAP_SEQUENTIAL, "filename": "O_ave_n", "sigma": -1})
# PLOT_SETTINGS.append({"param": "n_O2+_ave", "type": "scalar", "str": "$n$(O$^{+}_\mathrm{2,ave}$) [m$^{-3}$]", "log": 1, "lims": (1,1e7), "unit": 1, "colormap": CMAP_SEQUENTIAL, "filename": "O2_ave_n", "sigma": -1})
# PLOT_SETTINGS.append({"param": "n_H+planet_ave", "type": "scalar", "str": "$n$(H$^{+}_\mathrm{planet,ave}$) [m$^{-3}$]", "log": 1, "lims": (1,1e7), "unit": 1, "colormap": CMAP_SEQUENTIAL, "filename": "Hplanet_ave_n", "sigma": -1})
# PLOT_SETTINGS.append({"param": "v_H+sw_ave", "type": "magnitude", "str": "$|U|$(H$^{+}_\mathrm{sw,ave}$) [km/s]", "log": 0, "lims": (0,700e3), "unit": 1e3, "colormap": CMAP_SEQUENTIAL, "filename": "Hsw_ave_U", "sigma": -1})
# PLOT_SETTINGS.append({"param": "v_H+sw_ave", "type": "xcomp", "str": "$U_x$(H$^{+}_\mathrm{sw,ave}$) [km/s]", "log": 0, "lims": (-600e3,600e3), "unit": 1e3, "colormap": CMAP_DIVERGING, "filename": "Hsw_ave_Ux", "sigma": -1})
# PLOT_SETTINGS.append({"param": "v_H+sw_ave", "type": "ycomp", "str": "$U_y$(H$^{+}_\mathrm{sw,ave}$) [km/s]", "log": 0, "lims": (-600e3,600e3), "unit": 1e3, "colormap": CMAP_DIVERGING, "filename": "Hsw_ave_Uy", "sigma": -1})
# PLOT_SETTINGS.append({"param": "v_H+sw_ave", "type": "zcomp", "str": "$U_z$(H$^{+}_\mathrm{sw,ave}$) [km/s]", "log": 0, "lims": (-600e3,600e3), "unit": 1e3, "colormap": CMAP_DIVERGING, "filename": "Hsw_ave_Uz", "sigma": -1})
# PLOT_SETTINGS.append({"param": "v_He++sw_ave", "type": "magnitude", "str": "$|U|$(He$^{++}_\mathrm{sw,ave}$) [km/s]", "log": 0, "lims": (0,700e3), "unit": 1e3, "colormap": CMAP_SEQUENTIAL, "filename": "Hesw_ave_U", "sigma": -1})
# PLOT_SETTINGS.append({"param": "v_He++sw_ave", "type": "xcomp", "str": "$U_x$(He$^{++}_\mathrm{sw,ave}$) [km/s]", "log": 0, "lims": (-600e3,600e3), "unit": 1e3, "colormap": CMAP_DIVERGING, "filename": "Hesw_ave_Ux", "sigma": -1})
# PLOT_SETTINGS.append({"param": "v_He++sw_ave", "type": "ycomp", "str": "$U_y$(He$^{++}_\mathrm{sw,ave}$) [km/s]", "log": 0, "lims": (-600e3,600e3), "unit": 1e3, "colormap": CMAP_DIVERGING, "filename": "Hesw_ave_Uy", "sigma": -1})
# PLOT_SETTINGS.append({"param": "v_He++sw_ave", "type": "zcomp", "str": "$U_z$(He$^{++}_\mathrm{sw,ave}$) [km/s]", "log": 0, "lims": (-600e3,600e3), "unit": 1e3, "colormap": CMAP_DIVERGING, "filename": "Hesw_ave_Uz", "sigma": -1})
# PLOT_SETTINGS.append({"param": "v_Na+_ave", "type": "magnitude", "str": "$|U|$(Na$^{+}_\mathrm{ave}$) [km/s]", "log": 0, "lims": (0,1000e3), "unit": 1e3, "colormap": CMAP_SEQUENTIAL, "filename": "Na_ave_U", "sigma": -1})
# PLOT_SETTINGS.append({"param": "v_Na+_ave", "type": "xcomp", "str": "$U_x$(Na$^{+}_\mathrm{ave}$) [km/s]", "log": 0, "lims": (-600e3,600e3), "unit": 1e3, "colormap": CMAP_DIVERGING, "filename": "Na_ave_Ux", "sigma": -1})
# PLOT_SETTINGS.append({"param": "v_Na+_ave", "type": "ycomp", "str": "$U_y$(Na$^{+}_\mathrm{ave}$) [km/s]", "log": 0, "lims": (-600e3,600e3), "unit": 1e3, "colormap": CMAP_DIVERGING, "filename": "Na_ave_Uy", "sigma": -1})
# PLOT_SETTINGS.append({"param": "v_Na+_ave", "type": "zcomp", "str": "$U_z$(Na$^{+}_\mathrm{ave}$) [km/s]", "log": 0, "lims": (-600e3,600e3), "unit": 1e3, "colormap": CMAP_DIVERGING, "filename": "Na_ave_Uz", "sigma": -1})


def parse_arguments():
    """Parse CLI arguments and the exact list of files selected by Slurm."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "ncores", metavar="Ncores", type=int,
        help="worker processes (-1: print header only)",
    )
    parser.add_argument(
        "run_folder", metavar="runFolder", type=Path,
        help="directory with selected VLSV files",
    )
    parser.add_argument(
        "object_radius", metavar="Robject", type=float,
        help="object radius in metres",
    )
    parser.add_argument(
        "--file-list", required=True, type=Path,
        help="NUL-separated paths selected by the Slurm script",
    )
    parser.add_argument(
        "--output-dir", type=Path, default=Path("png"),
        help="directory for PNG output",
    )
    args = parser.parse_args()
    args.run_folder = args.run_folder.expanduser().resolve()

    if not args.run_folder.is_dir():
        parser.error(f"cannot read VLSV folder: {args.run_folder}")
    if not np.isfinite(args.object_radius) or args.object_radius <= 0:
        parser.error("Robject must be positive and finite")
    if args.ncores != -1 and not (1 <= args.ncores <= 100):
        parser.error("Ncores must be between 1 and 100, or -1 for header only")

    try:
        raw = args.file_list.read_bytes()
    except OSError as exc:
        parser.error(f"cannot read --file-list: {exc}")
    if not raw or not raw.endswith(b"\0"):
        parser.error("--file-list must contain NUL-terminated paths")

    files = []
    seen = set()
    for raw_path in raw[:-1].split(b"\0"):
        if not raw_path:
            parser.error("empty filename in --file-list")
        path = Path(os.fsdecode(raw_path)).expanduser().resolve()
        if path.parent != args.run_folder:
            parser.error(f"selected VLSV file is outside runFolder: {path}")
        if not re.fullmatch(r"state[0-9]+\.vlsv", path.name):
            parser.error(f"selected file is not stateNNN.vlsv: {path}")
        if not path.is_file():
            parser.error(f"selected VLSV file does not exist: {path}")
        if path in seen:
            parser.error(f"duplicate selected VLSV file: {path}")
        seen.add(path)
        files.append(path)

    if not files:
        parser.error("no VLSV files selected")
    return args, files


def read_metadata(reader, radius):
    """Get 3D mesh dimensions and axes in object radii."""
    xmin, ymin, zmin, xmax, ymax, zmax = reader.get_spatial_mesh_extent()
    mx, my, mz = reader.get_spatial_mesh_size()
    sx, sy, sz = reader.get_spatial_block_size()
    nx, ny, nz = mx * sx, my * sy, mz * sz
    if min(nx, ny, nz) <= 1:
        raise ValueError(f"Only 3D runs are supported: grid = ({nx}, {ny}, {nz})")
    limits = np.array(
        [xmin, xmax, ymin, ymax, zmin, zmax], dtype=float,
    ) / radius
    return (nx, ny, nz), limits, (xmin, ymin, zmin, xmax, ymax, zmax)


def choose_column(plane, rmin, rmax, n, coord, radius):
    """Choose the grid index nearest the requested slice position."""
    # Retain the original selection rule and fallback to the domain midpoint.
    requested = plane * radius
    if requested < rmin or requested > rmax:
        requested = (rmin + rmax) / 2
        print(
            f"WARNING: {coord}Plane outside domain; "
            f"setting to {requested / radius} Rp",
            flush=True,
        )
    index = int(np.floor((requested - rmin) / ((rmax - rmin) / n)))
    return int(np.clip(index, 0, n - 1)), requested / radius


def plane_slices(extent, grid, radius):
    """Return slice indices and actual plane positions (in object radii)."""
    nx, ny, nz = grid
    xmin, ymin, zmin, xmax, ymax, zmax = extent
    yz_index, x_actual = choose_column(X_PLANE, xmin, xmax, nx, "x", radius)
    xz_index, y_actual = choose_column(Y_PLANE, ymin, ymax, ny, "y", radius)
    xy_index, z_actual = choose_column(Z_PLANE, zmin, zmax, nz, "z", radius)
    return (yz_index, xz_index, xy_index), (x_actual, y_actual, z_actual)


def format_number(value, decimals=10):
    """Format values using the original plot-title rounding convention."""
    if decimals == 1:
        return str(int(round(value * decimals) / decimals))
    return str(round(value * decimals) / decimals)


def format_with_unit(value, unit_string):
    """Format a slice position and its unit for a plot title."""
    if abs(value) > 0:
        return str(round(value * 10) / 10) + unit_string
    return "0"


def coordinate_ticks(limits):
    """Use the original automatic tick-step heuristic."""
    ranges = (
        limits[1] - limits[0],
        limits[3] - limits[2],
        limits[5] - limits[4],
    )
    max_range = max(ranges)
    max_coordinate = max(abs(limits))
    max_power = 10 ** np.ceil(np.log10(max_coordinate)) if max_coordinate > 0 else 1
    step = 1
    n_ticks = 0
    for i, j, k in itertools.product(range(1, 100), (1, 2, 5), (-1, 1)):
        step = max_power / (k * j * i)
        n_ticks = max_range / step
        if 5 <= n_ticks <= 15:
            break
    if n_ticks < 5 or n_ticks > 15:
        print(
            f"WARNING: no good tick step: n={n_ticks}, step={step}",
            flush=True,
        )
    return np.arange(-max_power, max_power, step=step)


def detect_settings(reader):
    """Select configured parameters available in the first VLSV file."""
    selected = []
    vector_count = scalar_count = 0
    all_variables = reader.get_all_variables()
    for name in all_variables:
        if "flag" in name or "detector_" in name or "CellID" in name:
            continue
        matching = [setting for setting in PLOT_SETTINGS if setting["param"] == name]
        if not matching:
            print("Skipping variable not in PLOT_SETTINGS:", name)
            continue
        data = reader.read_variable_info(name).data
        is_vector = data.ndim > 1
        if is_vector and data.shape[1] != 3:
            raise ValueError(f"Unsupported vector size for {name}: {data.shape}")
        vector_count += int(is_vector)
        scalar_count += int(not is_vector)
        selected.extend(
            setting for setting in matching
            if setting["type"] in (
                ("magnitude", "xcomp", "ycomp", "zcomp")
                if is_vector else ("scalar",)
            )
        )

    # Add custom derived parameters manually if needed (disabled by default):
    # selected.append({"param": "nvO", "type": "nvO", "str": "$nU$(O$^{+}$) [m$^{-2}$ s$^{-1}$]", "log": 1, "lims": (1e6, 1e11), "unit": 1, "colormap": CMAP_SEQUENTIAL, "filename": "nvO", "sigma": -1})
    # selected.append({"param": "rho_m", "type": "rho_m", "str": "$\\rho_m$ [$\\mu$g/m$^{3}$]", "log": 1, "lims": (1e-21, 1e-18), "unit": 1e-6, "colormap": CMAP_SEQUENTIAL, "filename": "rho_m", "sigma": -1})
    # selected.append({"param": "U_bulk", "type": "U_bulk", "str": "$U_\\mathrm{bulk}$ [km/s]", "log": 0, "lims": (0, 600e3), "unit": 1e3, "colormap": CMAP_SEQUENTIAL, "filename": "U_bulk", "sigma": -1})
    # selected.append({"param": "p_B", "type": "p_B", "str": "$p_B$ [Pa]", "log": 1, "lims": (1e-13, 1e-10), "unit": 1, "colormap": CMAP_SEQUENTIAL, "filename": "p_B", "sigma": -1})
    if not selected:
        raise ValueError("No plottable variables found in first selected VLSV file")
    return selected, len(all_variables), vector_count, scalar_count


def read_plot_data(reader, setting):
    """Compute plot values; scientific calculations are unchanged."""
    kind, name = setting["type"], setting["param"]
    if kind == "magnitude":
        data = reader.read_variable_info(name).data
        return np.sqrt(data[:, 0]**2 + data[:, 1]**2 + data[:, 2]**2)
    if kind == "scalar":
        return reader.read_variable_info(name).data
    if kind in ("xcomp", "ycomp", "zcomp"):
        component = {"xcomp": 0, "ycomp": 1, "zcomp": 2}[kind]
        return reader.read_variable_info(name).data[:, component]
    if kind == "nvO":
        density = reader.read_variable_info("n_O+_ave").data
        velocity = reader.read_variable_info("v_O+_ave").data
        return density * np.sqrt(
            velocity[:, 0]**2 + velocity[:, 1]**2 + velocity[:, 2]**2
        )
    if kind in ("rho_m", "U_bulk"):
        density_names = (
            "n_H+planet_ave", "n_H+sw_ave", "n_He++sw_ave",
            "n_O+_ave", "n_O2+_ave",
        )
        masses = (
            PROTON_MASS, PROTON_MASS, HELIUM_MASS,
            OXYGEN_MASS, OXYGEN_MOLECULE_MASS,
        )
        densities = [reader.read_variable_info(key).data for key in density_names]
        rho = sum(density * mass for density, mass in zip(densities, masses))
        if kind == "rho_m":
            return rho
        velocity_names = (
            "v_H+planet_ave", "v_H+sw_ave", "v_He++sw_ave",
            "v_O+_ave", "v_O2+_ave",
        )
        velocities = [reader.read_variable_info(key).data for key in velocity_names]
        momentum = sum(
            (density * mass)[:, None] * velocity
            for density, mass, velocity in zip(densities, masses, velocities)
        )
        return np.sqrt(np.sum(momentum**2, axis=1)) / rho
    if kind == "p_B":
        field = reader.read_variable_info("cellBAverage").data
        return (
            field[:, 0]**2 + field[:, 1]**2 + field[:, 2]**2
        ) / (2 * MU_0)
    raise ValueError(f"Unknown parameter type: {kind}")


def configure_panel(setting, mappable, axis, draw_planet):
    """Apply planet overlay, color normalization, and tick presentation."""
    if draw_planet:
        axis.add_artist(Wedge((0, 0), 1.0, 90, 270, fc="dimgray"))
        axis.add_artist(Wedge((0, 0), 1.0, 270, 90, fc="w"))
        axis.add_artist(plt.Circle((0, 0), 1.0, color="k", fill=False, lw=0.5))
    if setting["log"]:
        mappable.set_norm(
            colors.LogNorm(
                vmin=setting["lims"][0] / setting["unit"],
                vmax=setting["lims"][1] / setting["unit"],
            )
        )
    axis.tick_params(
        which="both", direction=TICK_DIRECTION,
        length=TICK_LENGTH, width=TICK_WIDTH,
    )


def draw_figure(
    reader, setting, grid, limits, plane_indices, planes,
    tick_marks, time_s, sorted_cell_indices,
):
    """Draw the xz, xy and yz panels for one parameter and timestep."""
    nx, ny, nz = grid
    yz_index, xz_index, xy_index = plane_indices
    x_actual, y_actual, z_actual = planes
    values = read_plot_data(reader, setting)
    if values.size != nx * ny * nz:
        raise ValueError(
            f"{setting['param']}: expected {nx * ny * nz} "
            f"scalar cells, got {values.size}"
        )

    # Reuse the per-file CellID ordering for every plotted variable.
    volume = values[sorted_cell_indices].reshape(nz, ny, nx)
    panels = (
        volume[:, xz_index, :],
        volume[xy_index, :, :],
        volume[:, :, yz_index],
    )
    if setting["sigma"] > 0:
        panels = tuple(
            gaussian_filter(panel, sigma=setting["sigma"], mode="constant")
            for panel in panels
        )

    # Keep the original xz, xy, yz panel order and equal-aspect axes.
    extents = (
        (limits[0], limits[1], limits[4], limits[5]),
        (limits[0], limits[1], limits[2], limits[3]),
        (limits[2], limits[3], limits[4], limits[5]),
    )
    zoom = np.asarray(
        AXIS_LIMITS_ZOOM
        if np.atleast_1d(AXIS_LIMITS_ZOOM).size == 6 else limits
    )
    zooms = (
        (zoom[0:2], zoom[4:6]),
        (zoom[0:2], zoom[2:4]),
        (zoom[2:4], zoom[4:6]),
    )
    titles = (
        "3D: $xz$ ($y=$" + format_with_unit(y_actual, "$R_p$") + ")",
        "$t=$" + format_number(time_s) + " s\n3D: $xy$ ($z=$"
        + format_with_unit(z_actual, "$R_p$") + ")",
        "3D: $yz$ ($x=$" + format_with_unit(x_actual, "$R_p$") + ")",
    )
    labels = (("$x$", "$z$"), ("$x$", "$y$"), ("$y$", "$z$"))
    fig, axes = plt.subplots(1, 3, figsize=FIGURE_SIZE, frameon=True)
    try:
        for index, (axis, panel) in enumerate(zip(axes, panels)):
            mappable = axis.imshow(
                panel / setting["unit"],
                vmin=setting["lims"][0] / setting["unit"],
                vmax=setting["lims"][1] / setting["unit"],
                cmap=setting["colormap"], extent=extents[index],
                aspect="equal", origin="lower", interpolation="nearest",
            )
            axis.set_title(titles[index])
            axis.set_xlabel(labels[index][0] + " [$R_p$]")
            axis.set_ylabel(labels[index][1] + " [$R_p$]")
            axis.set_xticks(tick_marks)
            axis.tick_params("x", labelrotation=X_TICK_ANGLE)
            axis.set_yticks(tick_marks)
            axis.axis("scaled")
            axis.set_xlim(zooms[index][0])
            axis.set_ylim(zooms[index][1])
            configure_panel(setting, mappable, axis, SHOW_PLANET[index])
        fig.tight_layout()
        colorbar = fig.colorbar(mappable, ax=axes, shrink=0.5)
        colorbar.ax.set_title(setting["str"])
        return fig
    except BaseException:
        plt.close(fig)
        raise


def plot_file(path, *, settings, output_dir, radius):
    """Read one VLSV file and save one PNG per enabled parameter."""
    reader = rhb.vlsvfile.VlsvReader(str(path))
    grid, limits, extent = read_metadata(reader, radius)

    # Read and sort CellIDs only once per VLSV file, not once per plot.
    sorted_cell_indices = reader.read_variable("CellID").argsort()
    if sorted_cell_indices.size != np.prod(grid):
        raise ValueError(f"{path.name}: CellID count does not match the grid {grid}")
    plane_indices, planes = plane_slices(extent, grid, radius)
    tick_marks = coordinate_ticks(
        limits if np.atleast_1d(AXIS_LIMITS_ZOOM).size != 6
        else np.asarray(AXIS_LIMITS_ZOOM)
    )
    time_s = reader.read_parameter("t")
    if time_s is None:
        time_s = reader.read_parameter("time")

    for setting in settings:
        print(
            f"Plotting: {path.name} / {setting['param']} / {setting['type']}",
            flush=True,
        )
        fig = draw_figure(
            reader, setting, grid, limits, plane_indices, planes,
            tick_marks, time_s, sorted_cell_indices,
        )
        try:
            filename = f"{setting['filename']}_{path.name}.png"
            fig.savefig(output_dir / filename, dpi=FIGURE_DPI, transparent=False)
        finally:
            plt.close(fig)


def main():
    """Configure one simulation run and plot the explicitly selected files."""
    args, files = parse_arguments()
    print(
        f"Processing {len(files)} selected VLSV files from: {args.run_folder}",
        flush=True,
    )
    first = rhb.vlsvfile.VlsvReader(str(files[0]))
    grid, limits, extent = read_metadata(first, args.object_radius)
    indices, coords = plane_slices(extent, grid, args.object_radius)
    settings, total_vars, vector_vars, scalar_vars = detect_settings(first)

    print("Run and figure configuration (first selected VLSV file):")
    print(
        f"Rp = {args.object_radius / 1e3} km; grid = {grid}; "
        "figure: 3 columns x 1 row"
    )
    print(
        f"x = {limits[0]} ... {limits[1]} Rp; "
        f"y = {limits[2]} ... {limits[3]} Rp; "
        f"z = {limits[4]} ... {limits[5]} Rp"
    )
    print(f"xPlane = {coords[0]} Rp, YZ_Ncol = {indices[0]}")
    print(f"yPlane = {coords[1]} Rp, XZ_Ncol = {indices[1]}")
    print(f"zPlane = {coords[2]} Rp, XY_Ncol = {indices[2]}")
    print(
        f"Found {total_vars} variables; {vector_vars} recognized vectors; "
        f"{scalar_vars} recognized scalars"
    )
    print(
        f"Plotting {len(settings)} parameter(s), "
        f"{len(files) * len(settings)} PNG frames"
    )
    if args.ncores == -1:
        print("Plotting variables:")
        for setting in settings:
            print(setting)
        return

    args.output_dir.mkdir(parents=True, exist_ok=True)
    worker = partial(
        plot_file, settings=settings, output_dir=args.output_dir.resolve(),
        radius=args.object_radius,
    )
    print(f"Running with {args.ncores} Python worker processes", flush=True)
    with Pool(args.ncores) as pool:
        pool.map(worker, files)


if __name__ == "__main__":
    main()
