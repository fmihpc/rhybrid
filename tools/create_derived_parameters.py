# exec(open("test_particle_tracing.py").read());
# create derived variables from VLSV file and write a new file with them
import analysator as alr
import numpy as np

# constants
mp = 1.672621716e-27
mHe = 6.6464764e-27
mO = 2.6567625437e-26
mO2 = 5.3135250874e-26
mu0 = 1.25663706e-6

vr = alr.vlsvfile.VlsvReader("state00004000.vlsv")

# read magnetic field and number densities of all particle populations
B = vr.read_variable("cellBAverage")
Btot = np.sqrt(B[:,0]**2 + B[:,1]**2 + B[:,2]**2)
nHsw = vr.read_variable("n_H+sw_ave")
nHesw = vr.read_variable("n_He++sw_ave")
nO = vr.read_variable("n_O+_ave")
nO2 = vr.read_variable("n_O2+_ave")
nHpla = vr.read_variable("n_H+planet_ave")

# mass density
rhom = mp*nHsw + mHe*nHesw + mO*nO + mO2*nO2 + mp*nHpla

# Alfven speed
vA = Btot/np.sqrt(mu0*rhom)

# write a new VLSV file with some existing variables
writer = alr.vlsvfile.VlsvWriter(vr,"state00004000_derived_parameters.vlsv", copy_meshes=["SpatialGrid"])
writer.copy_variables(vr,varlist=["CellID","v_tot","n_tot"])

# write new derived variables
varinfo = alr.calculations.VariableInfo(rhom,name="rhom",units="")
writer.write_variable_info(varinfo,"SpatialGrid",1)
varinfo = alr.calculations.VariableInfo(vA,name="vA",units="")
writer.write_variable_info(varinfo,"SpatialGrid",1)



