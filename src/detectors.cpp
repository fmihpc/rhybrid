/** This file is part of the RHybrid simulation.
 *
 *  Copyright 2018- Aalto University
 *  Copyright 2015- Finnish Meteorological Institute
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
#include <numeric>

#include "hybrid.h"
#include "detectors.h"

#define DETECTOR_PARTICLE_FILE_VARIABLES 6
#define DETECTOR_BULK_PARAMETER_FILE_VARIABLES 18

using namespace std;

// ascii output writer for particle detector
bool writeDetectorParticle(Simulation& sim,SimulationClasses& simClasses) {
   // gather all recorded particles on master and write in file
   int N_linesLocal = Hybrid::detParticleOutput.size();
   vector<int> N_linesGlobal(sim.mpiProcesses);
   MPI_Allgather(&N_linesLocal,1,MPI_Type<int>(),&N_linesGlobal[0],1,MPI_Type<int>(),sim.comm);
   int N_linesTotal = accumulate(N_linesGlobal.begin(),N_linesGlobal.end(),0);
   vector<Real> detParticleOutputGlobal(N_linesTotal);
   // create displacement vector for detParticleOutputGlobal
   int displ[sim.mpiProcesses];
   if(sim.mpiRank == sim.MASTER_RANK) {
      int sum = 0.0;
      for (int i = 0; i < sim.mpiProcesses; ++i) {
	 displ[i] = sum;
	 sum += N_linesGlobal[i];
      }
   }
   // gather in master
   MPI_Gatherv(&Hybrid::detParticleOutput[0],N_linesLocal,MPI_Type<Real>(),&detParticleOutputGlobal[0],&N_linesGlobal[0],&displ[0],MPI_Type<Real>(),sim.MASTER_RANK,sim.comm);
   // master writes
   if(sim.mpiRank == sim.MASTER_RANK) {
      if(Hybrid::detParticleFileLineCnt <= Hybrid::N_detParticleMaxFileLines) {
	 if(detParticleOutputGlobal.size() % DETECTOR_PARTICLE_FILE_VARIABLES != 0) {
	    simClasses.logger << "(RHYBRID) DETECTORS: ERROR: error when writing a particle detector file" << endl << write;
	    return false;
	 }
	 ofstream particleFile;
	 string particleFileName = string("det_ple_") + int2str(sim.timestep,7) + string(".dat");
	 particleFile.open(particleFileName,ios_base::app);
	 particleFile.precision(6);
	 particleFile << scientific;
	 unsigned long fileLineCnt = 0;
	 particleFile << "% t popid cellid vx vy vz" << endl;
	 for(unsigned int i=0;i<detParticleOutputGlobal.size();i+=DETECTOR_PARTICLE_FILE_VARIABLES) {
	    particleFile.precision(6);
	    particleFile
	      << detParticleOutputGlobal[i+0] << " "                             // 01 t
	      << static_cast<unsigned int>(detParticleOutputGlobal[i+1]) << " "  // 02 popid
	      << static_cast<unsigned int>(detParticleOutputGlobal[i+2]) << " "; // 03 blockid
	    particleFile.precision(4);
	    particleFile
	      << detParticleOutputGlobal[i+3] << " "                             // 04 vx
	      << detParticleOutputGlobal[i+4] << " "                             // 05 vy
	      << detParticleOutputGlobal[i+5] << endl;                           // 06 vz
	      /*<< detParticleOutputGlobal[i+2] << " "                           // weight
	      << detParticleOutputGlobal[i+7] << " "                             // ini: t
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
	 Hybrid::detParticleFileLineCnt += fileLineCnt;
	 simClasses.logger
	   << "(RHYBRID) DETECTORS: written " << fileLineCnt << " particles in " << particleFileName << endl
	   << "(RHYBRID) DETECTORS: particle detector max. counter: " << static_cast<long long>(Hybrid::detParticleFileLineCnt) << "/" << static_cast<long long>(Hybrid::N_detParticleMaxFileLines) << endl << write;
      }
   }
   MPI_Barrier(sim.comm);
   MPI_Bcast(&Hybrid::detParticleFileLineCnt,1,MPI_Type<Real>(),sim.MASTER_RANK,sim.comm);
   // empty particle output list
   Hybrid::detParticleOutput.clear();
   return true;
}

bool recordDetectorBulkParam(Simulation& sim,SimulationClasses& simClasses) {
   bool success = true;
   bool* detBlkFlag = reinterpret_cast<bool*>(simClasses.pargrid.getUserData(Hybrid::dataDetectorBulkParamFlagID));
   Real* cellRhoQi = simClasses.pargrid.getUserDataStatic<Real>(Hybrid::dataCellRhoQiID);
   Real* cellB = simClasses.pargrid.getUserDataStatic<Real>(Hybrid::dataCellBID);
   Real* cellJ = simClasses.pargrid.getUserDataStatic<Real>(Hybrid::dataCellJID);
   Real* cellUe = simClasses.pargrid.getUserDataStatic<Real>(Hybrid::dataCellUeID);
   Real* cellJi = simClasses.pargrid.getUserDataStatic<Real>(Hybrid::dataCellJiID);
   Real* nodeE  = simClasses.pargrid.getUserDataStatic<Real>(Hybrid::dataNodeEID);
   // loop all blocks and cells in this process
   for(pargrid::CellID b=0; b<simClasses.pargrid.getNumberOfLocalCells(); ++b) {
      for(int k=0; k<block::WIDTH_Z; ++k) for(int j=0; j<block::WIDTH_Y; ++j) for(int i=0; i<block::WIDTH_X; ++i) {
	 const pargrid::CellID globalID = simClasses.pargrid.getGlobalIDs()[b]; // Note: we assume here block size = 1 (i.e. blocks == cells)
	 const int n = (b*block::SIZE+block::index(i,j,k));
	 const int n3 = n*3;
	 // record bulk parameters in flagged cells
	 if(detBlkFlag[n] == true && Hybrid::detBulkParamRecording == true) {
	    Hybrid::detBulkParamOutput.push_back( static_cast<Real>(sim.t)    ); // 1. t
	    Hybrid::detBulkParamOutput.push_back( static_cast<Real>(globalID) ); // 2. cellid
	    Hybrid::detBulkParamOutput.push_back( cellRhoQi[n] ); //  3. rhoqi
	    Hybrid::detBulkParamOutput.push_back( cellB[n3+0]  ); //  4. Bx
	    Hybrid::detBulkParamOutput.push_back( cellB[n3+1]  ); //  5. By
	    Hybrid::detBulkParamOutput.push_back( cellB[n3+2]  ); //  6. Bz
	    Hybrid::detBulkParamOutput.push_back( cellJ[n3+0]  ); //  7. Jx
	    Hybrid::detBulkParamOutput.push_back( cellJ[n3+1]  ); //  8. Jy
	    Hybrid::detBulkParamOutput.push_back( cellJ[n3+2]  ); //  9. Jz
	    Hybrid::detBulkParamOutput.push_back( cellUe[n3+0] ); // 10. Uex
	    Hybrid::detBulkParamOutput.push_back( cellUe[n3+1] ); // 11. Uey
	    Hybrid::detBulkParamOutput.push_back( cellUe[n3+2] ); // 12. Uez
	    Hybrid::detBulkParamOutput.push_back( cellJi[n3+0] ); // 13. Jix
	    Hybrid::detBulkParamOutput.push_back( cellJi[n3+1] ); // 14. Jiy
	    Hybrid::detBulkParamOutput.push_back( cellJi[n3+2] ); // 15. Jiz
	    Hybrid::detBulkParamOutput.push_back( nodeE[n3+0]  ); // 16. Ex(node)
	    Hybrid::detBulkParamOutput.push_back( nodeE[n3+1]  ); // 17. Ey(node)
	    Hybrid::detBulkParamOutput.push_back( nodeE[n3+2]  ); // 18. Ez(node)
	 }
      }
   }
   return success;
}

// ascii output writer for bulk parameter detector
bool writeDetectorBulkParam(Simulation& sim,SimulationClasses& simClasses) {
   // gather all recorded bulk values on master and write in file
   int N_linesLocal = Hybrid::detBulkParamOutput.size();
   vector<int> N_linesGlobal(sim.mpiProcesses);
   MPI_Allgather(&N_linesLocal,1,MPI_Type<int>(),&N_linesGlobal[0],1,MPI_Type<int>(),sim.comm);
   int N_linesTotal = accumulate(N_linesGlobal.begin(),N_linesGlobal.end(),0);
   vector<Real> detBulkParamOutputGlobal(N_linesTotal);
   // create displacement vector for detBulkParamOutputGlobal
   int displ[sim.mpiProcesses];
   if(sim.mpiRank == sim.MASTER_RANK) {
      int sum = 0.0;
      for (int i = 0; i < sim.mpiProcesses; ++i) {
	 displ[i] = sum;
	 sum += N_linesGlobal[i];
      }
   }
   // gather in master
   MPI_Gatherv(&Hybrid::detBulkParamOutput[0],N_linesLocal,MPI_Type<Real>(),&detBulkParamOutputGlobal[0],&N_linesGlobal[0],&displ[0],MPI_Type<Real>(),sim.MASTER_RANK,sim.comm);
   // master writes
   if(sim.mpiRank == sim.MASTER_RANK) {
      if(Hybrid::detBulkParamFileLineCnt <= Hybrid::N_detBulkParamMaxFileLines) {
	 if(detBulkParamOutputGlobal.size() % DETECTOR_BULK_PARAMETER_FILE_VARIABLES != 0) {
	    simClasses.logger << "(RHYBRID) DETECTORS: ERROR: error when writing a bulk parameter detector file" << endl << write;
	    return false;
	 }
	 ofstream bulkParamFile;
	 string bulkParamFileName = string("det_blk_") + int2str(sim.timestep,7) + string(".dat");
	 bulkParamFile.open(bulkParamFileName,ios_base::app);
	 bulkParamFile.precision(6);
	 bulkParamFile << scientific;
	 unsigned long fileLineCnt = 0;
	 bulkParamFile << "% t cellid rhoqi Bx By Bz Jx Jy Jz Uex Uey Uez Jix Jiy Jiz Ex(node) Ey(node) Ez(node)" << endl;
	 for(unsigned int i=0;i<detBulkParamOutputGlobal.size();i+=DETECTOR_BULK_PARAMETER_FILE_VARIABLES) {
	    bulkParamFile
	      << detBulkParamOutputGlobal[i+0] << " "
	      << static_cast<unsigned int>(detBulkParamOutputGlobal[i+1]) << " "
	      << detBulkParamOutputGlobal[i+2] << " "
	      << detBulkParamOutputGlobal[i+3] << " "
	      << detBulkParamOutputGlobal[i+4] << " "
	      << detBulkParamOutputGlobal[i+5] << " "
	      << detBulkParamOutputGlobal[i+6] << " "
	      << detBulkParamOutputGlobal[i+7] << " "
	      << detBulkParamOutputGlobal[i+8] << " "
	      << detBulkParamOutputGlobal[i+9] << " "
	      << detBulkParamOutputGlobal[i+10] << " "
	      << detBulkParamOutputGlobal[i+11] << " "
	      << detBulkParamOutputGlobal[i+12] << " "
	      << detBulkParamOutputGlobal[i+13] << " "
	      << detBulkParamOutputGlobal[i+14] << " "
	      << detBulkParamOutputGlobal[i+15] << " "
	      << detBulkParamOutputGlobal[i+16] << " "
	      << detBulkParamOutputGlobal[i+17] << endl;
	    fileLineCnt++;
	 }
	 bulkParamFile << flush;
	 bulkParamFile.close();
	 Hybrid::detBulkParamFileLineCnt += fileLineCnt;
	 simClasses.logger
	   << "(RHYBRID) DETECTORS: written " << fileLineCnt << " entries in " << bulkParamFileName << endl
	   << "(RHYBRID) DETECTORS: bulk parameter detector max. counter: " << static_cast<long long>(Hybrid::detBulkParamFileLineCnt) << "/" << static_cast<long long>(Hybrid::N_detBulkParamMaxFileLines) << endl << write;
      }
   }
   MPI_Barrier(sim.comm);
   MPI_Bcast(&Hybrid::detBulkParamFileLineCnt,1,MPI_Type<Real>(),sim.MASTER_RANK,sim.comm);
   // empty particle output list
   Hybrid::detBulkParamOutput.clear();
   return true;
}

#endif
