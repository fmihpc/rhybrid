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

#define DETECTOR_PARTICLE_FILE_VARIABLES 9
#define DETECTOR_BULK_PARAMETER_FILE_VARIABLES 18

using namespace std;

// ascii output writer for particles from the cell detector and the impact detector
bool writeCellDetectorParticles(Simulation& sim,SimulationClasses& simClasses,const vector<Real>& detParticleData,string fileNamePrefix) {
   // gather all recorded particles on master and write in file
   int N_linesLocal = detParticleData.size();
   vector<int> N_linesGlobal(sim.mpiProcesses);
   MPI_Allgather(&N_linesLocal,1,MPI_INT,&N_linesGlobal[0],1,MPI_INT,sim.comm);
   int N_linesTotal = accumulate(N_linesGlobal.begin(),N_linesGlobal.end(),0);
   vector<Real> detParticleDataGlobal(N_linesTotal);
   // create displacement vector for detParticleDataGlobal
   int displ[sim.mpiProcesses];
   if (sim.mpiRank == sim.MASTER_RANK) {
      int sum = 0;
      for (int i = 0; i < sim.mpiProcesses; ++i) {
	 displ[i] = sum;
	 sum += N_linesGlobal[i];
      }
   }
   // gather in master
   MPI_Gatherv(&detParticleData[0],N_linesLocal,MPI_Type<Real>(),&detParticleDataGlobal[0],&N_linesGlobal[0],&displ[0],MPI_Type<Real>(),sim.MASTER_RANK,sim.comm);
   // master writes
   if (sim.mpiRank == sim.MASTER_RANK) {
      if (Hybrid::detParticleFileLineCnt <= Hybrid::N_detParticleMaxFileLines) {
	 if (detParticleDataGlobal.size() % DETECTOR_PARTICLE_FILE_VARIABLES != 0) {
	    simClasses.logger << "(RHYBRID) DETECTORS: ERROR: error when writing a particle detector file" << endl << write;
	    return false;
	 }
	 if (detParticleDataGlobal.size() < 1) { return true; }
	 ofstream particleFile;
	 string particleFileName = fileNamePrefix + int2str(sim.timestep,7) + string(".dat");
	 particleFile.open(particleFileName,ios_base::app);
	 particleFile.precision(6);
	 particleFile << scientific;
	 unsigned long fileLineCnt = 0;
	 particleFile << "% t popid cellid x y z vx vy vz" << endl;
	 for (unsigned int i=0;i<detParticleDataGlobal.size();i+=DETECTOR_PARTICLE_FILE_VARIABLES) {
	    particleFile.precision(6);
	    particleFile
	      << detParticleDataGlobal[i+0] << " "                             // 01 t
	      << static_cast<unsigned int>(detParticleDataGlobal[i+1]) << " "  // 02 popid
	      << static_cast<unsigned int>(detParticleDataGlobal[i+2]) << " "; // 03 blockid
	    particleFile.precision(4);
	    particleFile
	      << detParticleDataGlobal[i+3] << " "                             // 04 x
	      << detParticleDataGlobal[i+4] << " "                             // 05 y
	      << detParticleDataGlobal[i+5] << " "                             // 06 z
	      << detParticleDataGlobal[i+6] << " "                             // 07 vx
	      << detParticleDataGlobal[i+7] << " "                             // 08 vy
	      << detParticleDataGlobal[i+8] << endl;                           // 09 vz
	      /*<< detParticleDataGlobal[i+2] << " "                           // weight
	      << detParticleDataGlobal[i+7] << " "                             // ini: t
	      << static_cast<unsigned int>(detParticleDataGlobal[i+8]) << " "  // ini: block id
	      << detParticleDataGlobal[i+9] << " "                             // ini: x
	      << detParticleDataGlobal[i+10] << " "                            // ini: y
	      << detParticleDataGlobal[i+11] << " "                            // ini: z
	      << detParticleDataGlobal[i+12] << " "                            // ini: vx
	      << detParticleDataGlobal[i+13] << " "                            // ini: vy
	      << detParticleDataGlobal[i+14] << endl;                          // ini: vz*/
	    fileLineCnt++;
	 }
	 particleFile << flush;
	 particleFile.close();
	 Hybrid::detParticleFileLineCnt += fileLineCnt;
	 simClasses.logger
	   << "(RHYBRID) DETECTORS: written " << fileLineCnt << " particles in " << particleFileName << endl
	   << "(RHYBRID) DETECTORS: particle detector max. counter: " << static_cast<long long>(Hybrid::detParticleFileLineCnt) << "/" << static_cast<long long>(Hybrid::N_detParticleMaxFileLines) << endl;
      }
   }
   MPI_Barrier(sim.comm);
   MPI_Bcast(&Hybrid::detParticleFileLineCnt,1,MPI_Type<Real>(),sim.MASTER_RANK,sim.comm);
   return true;
}

// write particles recorded by all detectors to files
bool writeDetectorParticleOutput(Simulation& sim,SimulationClasses& simClasses) {
   bool success = true;
   if (Hybrid::detCellParticleEnabled == true) {
      // write particles from the cell detector
      success = writeCellDetectorParticles(sim,simClasses,Hybrid::detCellParticleData,"det_cell_ple_");
      // empty particle output list of the cell detector
      Hybrid::detCellParticleData.clear();
   }
   if (Hybrid::detParticleRecordImpacts == true) {
      // write particles from the impact detector
      success = writeCellDetectorParticles(sim,simClasses,Hybrid::detImpactParticleData,"det_impact_ple_");
      // empty particle output list of the impact detector
      Hybrid::detImpactParticleData.clear();
   }
   simClasses.logger << endl << write;
   return success;
}

// store bulk parameter data from cell detectors to be written in files later
bool recordCellDetectorBulkParamData(Simulation& sim,SimulationClasses& simClasses) {
   bool success = true;
   bool* detBlkFlag = reinterpret_cast<bool*>(simClasses.pargrid.getUserData(Hybrid::dataDetectorCellBulkParamFlagID));
   Real* cellRhoQi = simClasses.pargrid.getUserDataStatic<Real>(Hybrid::dataCellRhoQiID);
   Real* cellB = simClasses.pargrid.getUserDataStatic<Real>(Hybrid::dataCellBID);
   Real* cellJ = simClasses.pargrid.getUserDataStatic<Real>(Hybrid::dataCellJID);
   Real* cellUe = simClasses.pargrid.getUserDataStatic<Real>(Hybrid::dataCellUeID);
   Real* cellJi = simClasses.pargrid.getUserDataStatic<Real>(Hybrid::dataCellJiID);
   Real* nodeE  = simClasses.pargrid.getUserDataStatic<Real>(Hybrid::dataNodeEID);
   // loop all blocks and cells in this process
   for (pargrid::CellID b=0; b<simClasses.pargrid.getNumberOfLocalCells(); ++b) {
      for (int k=0; k<block::WIDTH_Z; ++k) for (int j=0; j<block::WIDTH_Y; ++j) for (int i=0; i<block::WIDTH_X; ++i) {
	 const pargrid::CellID globalID = simClasses.pargrid.getGlobalIDs()[b]; // Note: we assume here block size = 1 (i.e. blocks == cells)
	 const int n = (b*block::SIZE+block::index(i,j,k));
	 const int n3 = n*3;
	 // record bulk parameters in flagged cells
	 if (detBlkFlag[n] == true && Hybrid::detBulkParamRecording == true) {
	    Hybrid::detCellBulkParamData.push_back( static_cast<Real>(sim.t)    ); // 1. t
	    Hybrid::detCellBulkParamData.push_back( static_cast<Real>(globalID) ); // 2. cellid
	    Hybrid::detCellBulkParamData.push_back( cellRhoQi[n] ); //  3. rhoqi
	    Hybrid::detCellBulkParamData.push_back( cellB[n3+0]  ); //  4. Bx
	    Hybrid::detCellBulkParamData.push_back( cellB[n3+1]  ); //  5. By
	    Hybrid::detCellBulkParamData.push_back( cellB[n3+2]  ); //  6. Bz
	    Hybrid::detCellBulkParamData.push_back( cellJ[n3+0]  ); //  7. Jx
	    Hybrid::detCellBulkParamData.push_back( cellJ[n3+1]  ); //  8. Jy
	    Hybrid::detCellBulkParamData.push_back( cellJ[n3+2]  ); //  9. Jz
	    Hybrid::detCellBulkParamData.push_back( cellUe[n3+0] ); // 10. Uex
	    Hybrid::detCellBulkParamData.push_back( cellUe[n3+1] ); // 11. Uey
	    Hybrid::detCellBulkParamData.push_back( cellUe[n3+2] ); // 12. Uez
	    Hybrid::detCellBulkParamData.push_back( cellJi[n3+0] ); // 13. Jix
	    Hybrid::detCellBulkParamData.push_back( cellJi[n3+1] ); // 14. Jiy
	    Hybrid::detCellBulkParamData.push_back( cellJi[n3+2] ); // 15. Jiz
	    Hybrid::detCellBulkParamData.push_back( nodeE[n3+0]  ); // 16. Ex(node)
	    Hybrid::detCellBulkParamData.push_back( nodeE[n3+1]  ); // 17. Ey(node)
	    Hybrid::detCellBulkParamData.push_back( nodeE[n3+2]  ); // 18. Ez(node)
	 }
      }
   }
   return success;
}

// ascii output writer for bulk parameter detector
bool writeDetectorBulkParamOutput(Simulation& sim,SimulationClasses& simClasses) {
   // gather all recorded bulk values on master and write in file
   int N_linesLocal = Hybrid::detCellBulkParamData.size();
   vector<int> N_linesGlobal(sim.mpiProcesses);
   MPI_Allgather(&N_linesLocal,1,MPI_INT,&N_linesGlobal[0],1,MPI_INT,sim.comm);
   int N_linesTotal = accumulate(N_linesGlobal.begin(),N_linesGlobal.end(),0);
   vector<Real> detCellBulkParamDataGlobal(N_linesTotal);
   // create displacement vector for detCellBulkParamDataGlobal
   int displ[sim.mpiProcesses];
   if (sim.mpiRank == sim.MASTER_RANK) {
      int sum = 0;
      for (int i = 0; i < sim.mpiProcesses; ++i) {
	 displ[i] = sum;
	 sum += N_linesGlobal[i];
      }
   }
   // gather in master
   MPI_Gatherv(&Hybrid::detCellBulkParamData[0],N_linesLocal,MPI_Type<Real>(),&detCellBulkParamDataGlobal[0],&N_linesGlobal[0],&displ[0],MPI_Type<Real>(),sim.MASTER_RANK,sim.comm);
   // master writes
   if (sim.mpiRank == sim.MASTER_RANK) {
      if (Hybrid::detBulkParamFileLineCnt <= Hybrid::N_detBulkParamMaxFileLines) {
	 if (detCellBulkParamDataGlobal.size() % DETECTOR_BULK_PARAMETER_FILE_VARIABLES != 0) {
	    simClasses.logger << "(RHYBRID) DETECTORS: ERROR: error when writing a bulk parameter detector file" << endl << write;
	    return false;
	 }
	 if (detCellBulkParamDataGlobal.size() < 1) { return true; }
	 ofstream bulkParamFile;
	 string bulkParamFileName = string("det_cell_blk_") + int2str(sim.timestep,7) + string(".dat");
	 bulkParamFile.open(bulkParamFileName,ios_base::app);
	 bulkParamFile.precision(6);
	 bulkParamFile << scientific;
	 unsigned long fileLineCnt = 0;
	 bulkParamFile << "% t cellid rhoqi Bx By Bz Jx Jy Jz Uex Uey Uez Jix Jiy Jiz Ex(node) Ey(node) Ez(node)" << endl;
	 for (unsigned int i=0;i<detCellBulkParamDataGlobal.size();i+=DETECTOR_BULK_PARAMETER_FILE_VARIABLES) {
	    bulkParamFile
	      << detCellBulkParamDataGlobal[i+0] << " "
	      << static_cast<unsigned int>(detCellBulkParamDataGlobal[i+1]) << " "
	      << detCellBulkParamDataGlobal[i+2] << " "
	      << detCellBulkParamDataGlobal[i+3] << " "
	      << detCellBulkParamDataGlobal[i+4] << " "
	      << detCellBulkParamDataGlobal[i+5] << " "
	      << detCellBulkParamDataGlobal[i+6] << " "
	      << detCellBulkParamDataGlobal[i+7] << " "
	      << detCellBulkParamDataGlobal[i+8] << " "
	      << detCellBulkParamDataGlobal[i+9] << " "
	      << detCellBulkParamDataGlobal[i+10] << " "
	      << detCellBulkParamDataGlobal[i+11] << " "
	      << detCellBulkParamDataGlobal[i+12] << " "
	      << detCellBulkParamDataGlobal[i+13] << " "
	      << detCellBulkParamDataGlobal[i+14] << " "
	      << detCellBulkParamDataGlobal[i+15] << " "
	      << detCellBulkParamDataGlobal[i+16] << " "
	      << detCellBulkParamDataGlobal[i+17] << endl;
	    fileLineCnt++;
	 }
	 bulkParamFile << flush;
	 bulkParamFile.close();
	 Hybrid::detBulkParamFileLineCnt += fileLineCnt;
	 simClasses.logger
	   << "(RHYBRID) DETECTORS: written " << fileLineCnt << " entries in " << bulkParamFileName << endl
	   << "(RHYBRID) DETECTORS: bulk parameter detector max. counter: " << static_cast<long long>(Hybrid::detBulkParamFileLineCnt) << "/" << static_cast<long long>(Hybrid::N_detBulkParamMaxFileLines) << endl;
      }
   }
   MPI_Barrier(sim.comm);
   MPI_Bcast(&Hybrid::detBulkParamFileLineCnt,1,MPI_Type<Real>(),sim.MASTER_RANK,sim.comm);
   // empty bulk parameter output list
   Hybrid::detCellBulkParamData.clear();
   simClasses.logger << endl << write;
   return true;
}

#endif
