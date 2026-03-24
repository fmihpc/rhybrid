# create derived variables from a VLSV file and write a new file with them
import os
import sys
from pathlib import Path
import argparse
try:
 import analysator as alr
except ModuleNotFoundError as err:
 print("Analysator not found: " + str(err))
 sys.exit()
try:
 import numpy as np
except ModuleNotFoundError as err:
 print("NumPy not found: " + str(err))
 sys.exit()

# parse input arguments
parser = argparse.ArgumentParser("create_derived_parameters.py")
parser.add_argument("input_file",help="Input VLSV file",type=Path)
parser.add_argument("output_file",help="Output VLSV file",type=Path)
args = parser.parse_args()
input_file = str(args.input_file) #input_file = './state00004000.vlsv'
output_file = str(args.output_file) #output_file = './state00004000_derived_parameters.vlsv'
if Path(input_file).is_file() == False:
 print('ERROR: input file does not exist (' + input_file + ')')
 sys.exit()
if Path(output_file).is_file() == True:
 print('ERROR: output file already exists (' + output_file + ')')
 sys.exit()

# constants
mp = 1.672621716e-27
mHe = 6.6464764e-27
mO = 2.6567625437e-26
mO2 = 5.3135250874e-26
mu0 = 1.25663706e-6

# open and create a VLSV file reader
vr = alr.vlsvfile.VlsvReader(input_file)

# example: list all variables in the file
#print('===== VARIABLES IN ' + fn + '.vlsv')
#print(vr.get_variables())
#print('=====')

# example: read variable info and dimensions (scalar or vector)
#print(vr.read_variable_info('n_H+sw_ave'))
#print(vr.read_variable_vectorsize('n_H+sw_ave'))

# read magnetic field and number densities of all particle populations
B = vr.read_variable('cellBAverage')
Btot = np.sqrt(B[:,0]**2 + B[:,1]**2 + B[:,2]**2)
nHsw = vr.read_variable('n_H+sw_ave')
nHesw = vr.read_variable('n_He++sw_ave')
nO = vr.read_variable('n_O+_ave')
nO2 = vr.read_variable('n_O2+_ave')
nHpla = vr.read_variable('n_H+planet_ave')

# mass density including all populations
rhom = mp*nHsw + mHe*nHesw + mO*nO + mO2*nO2 + mp*nHpla

# Alfven speed
vA = Btot/np.sqrt(mu0*rhom)

# open and create a VLSV file writer
writer = alr.vlsvfile.VlsvWriter(vr,output_file, copy_meshes=['SpatialGrid'])

# write original variables
writer.copy_variables(vr,varlist=['CellID','v_tot','n_tot'])

# write new derived variables
varinfo = alr.calculations.VariableInfo(rhom,name='rhom',units='')
writer.write_variable_info(varinfo,'SpatialGrid',1)
varinfo = alr.calculations.VariableInfo(vA,name='vA',units='')
writer.write_variable_info(varinfo,'SpatialGrid',1)

