/** This file is part of the RHybrid simulation.
 *
 *  Copyright 2018- Aalto University
 *  Copyright 2018- Finnish Meteorological Institute
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifdef USE_DETECTORS

#include <cstdlib>
#include <iostream>
#include <cmath>

#include "hybrid.h"
#include "detectors.h"

using namespace std;

bool writeDetectorBulkParams(Simulation& sim,SimulationClasses& simClasses) {
   bool success = true;
   return success;
}

// ascii output writer for gathered particles
bool writeDetectorParticles(Simulation& sim,SimulationClasses& simClasses) {
   // gather all recorded particles on master and write in file
   int N_plesLocal = Hybrid::detParticleOutput.size();
   vector<int> N_plesGlobal(sim.mpiProcesses);
   MPI_Allgather(&N_plesLocal,1,MPI_Type<int>(),&N_plesGlobal[0],1,MPI_Type<int>(),sim.comm);
   int N_plesTotal = accumulate(N_plesGlobal.begin(),N_plesGlobal.end(),0);
   vector<Real> detParticleOutputGlobal(N_plesTotal);
   // create displacement vector for detParticleOutputGlobal
   int displ[sim.mpiProcesses];
   if(sim.mpiRank == sim.MASTER_RANK) {
      int sum = 0.0;
      for (int i = 0; i < sim.mpiProcesses; ++i) {
	 displ[i] = sum;
	 sum += N_plesGlobal[i];
      }
   }
   // gather in master
   MPI_Gatherv(&Hybrid::detParticleOutput[0],N_plesLocal,MPI_Type<Real>(),&detParticleOutputGlobal[0],&N_plesGlobal[0],&displ[0],MPI_Type<Real>(),sim.MASTER_RANK,sim.comm);
   // master writes
   if(sim.mpiRank == sim.MASTER_RANK) {
      if(Hybrid::detFileLineCntParticles <= Hybrid::detMaxRecordedParticles) {
	 if(detParticleOutputGlobal.size() % DETECTOR_PARTICLE_FILE_VARIABLES != 0) {
	    simClasses.logger << "(RHYBRID) DETECTORS: ERROR: error when writing a particle file" << endl << write;
	    return false;
	 }
	 ofstream particleFile;
	 string particleFileName = string("detector_particles_") + int2str(sim.timestep,7) + string(".dat");
	 particleFile.open(particleFileName,ios_base::app);
	 particleFile.precision(6);
	 particleFile << scientific;
	 unsigned long fileLineCnt = 0;
	 for(unsigned int i=0;i<detParticleOutputGlobal.size();i+=DETECTOR_PARTICLE_FILE_VARIABLES) {
	    particleFile
	      << detParticleOutputGlobal[i+0] << " "                             // 01 det: t
	      << static_cast<unsigned int>(detParticleOutputGlobal[i+1]) << " "  // 02 det: popid
	      //<< detParticleOutputGlobal[i+2] << " "                             //  det: weight
	      << static_cast<unsigned int>(detParticleOutputGlobal[i+2]) << " "  // 03 det: block id
	      << detParticleOutputGlobal[i+3] << " "                             // 04 det: vx
	      << detParticleOutputGlobal[i+4] << " "                             // 05 det: vy
	      << detParticleOutputGlobal[i+5] << endl;                           // 06 det: vz
	      /*<< detParticleOutputGlobal[i+7] << " "                             // ini: t
	      << static_cast<unsigned int>(detParticleOutputGlobal[i+8]) << " "  // ini: block id
	      << detParticleOutputGlobal[i+9] << " "                             // ini: x
	      << detParticleOutputGlobal[i+10] << " "                            // ini: y
	      << detParticleOutputGlobal[i+11] << " "                            // ini: z
	      << detParticleOutputGlobal[i+12] << " "                            // ini: vx
	      << detParticleOutputGlobal[i+13] << " "                            // ini: vy
	      << detParticleOutputGlobal[i+14] << endl;                          // ini: vz*/
	    fileLineCnt++;
	 }
	 particleFile << flush;
	 particleFile.close();
	 Hybrid::detFileLineCntParticles += fileLineCnt;
	 simClasses.logger
	   << "(RHYBRID) DETECTORS: written " << particleFileName << " with " << fileLineCnt << " particles" << endl
	   << "(RHYBRID) DETECTORS: maximum particle counter: " << static_cast<long long>(Hybrid::detFileLineCntParticles) << "/" << static_cast<long long>(Hybrid::detMaxRecordedParticles) << endl << write;
      }
   }
   MPI_Barrier(sim.comm);
   MPI_Bcast(&Hybrid::detFileLineCntParticles,1,MPI_Type<Real>(),sim.MASTER_RANK,sim.comm);
   // empty particle output list
   Hybrid::detParticleOutput.clear();
   return true;
}

#endif
