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

#include <cstdlib>
#include <iostream>
#include <map>
#include <vector>
#include <limits>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <iomanip>

#include "hybrid.h"
#include "operator_userdata.h"
#include "hybrid_propagator.h"
#include "particle_definition.h"
#include "particle_species.h"
#ifdef USE_B_CONSTANT
#include "magnetic_field.h"
#endif

using namespace std;

UserDataOP::UserDataOP(): DataOperator() { 
   profileID = -1;
}

UserDataOP::~UserDataOP() {finalize();}

bool UserDataOP::finalize() {return true;}

std::string UserDataOP::getName() const {return "UserData";}

bool UserDataOP::initialize(ConfigReader& cr,Simulation& sim,SimulationClasses& simClasses) {
   return DataOperator::initialize(cr,sim,simClasses);
}

// log amounts of macroparticles
void logWriteMainMacroparticles(Simulation& sim,SimulationClasses& simClasses,const std::vector<ParticleListBase*>& particleLists) {
    simClasses.logger << "(RHYBRID) Number of macroparticles per population in the simulation:" << endl;
    for(size_t s=0;s<particleLists.size();++s) {
	Real N_macroParticles = 0.0;
	// For now skip particles with invalid data id:
	pargrid::DataID speciesDataID = pargrid::INVALID_DATAID;
	if(particleLists[s]->getParticles(speciesDataID) == true) {
	    pargrid::DataWrapper<Particle<Real> > wrapper = simClasses.pargrid.getUserDataDynamic<Particle<Real> >(speciesDataID);
	    for(pargrid::CellID b=0; b<simClasses.pargrid.getNumberOfLocalCells(); ++b) {
		pargrid::ArraySizetype N_particles = wrapper.size(b);
		N_macroParticles += N_particles;
	    }
	}
	Real N_macroParticlesGlobal = 0.0;
	MPI_Reduce(&N_macroParticles,&N_macroParticlesGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
	const Species* species = reinterpret_cast<const Species*>(particleLists[s]->getSpecies());
	simClasses.logger << "N(" << species->name << ") = " << real2str(N_macroParticlesGlobal,15) << " = " << real2str(N_macroParticlesGlobal/1.0e9,15) << " x 10^9" << endl;
    }
    simClasses.logger << write;
}

bool UserDataOP::writeData(const std::string& spatMeshName,const std::vector<ParticleListBase*>& particleLists) {
   bool success = true;
   if(getInitialized() == false) { return false; }
   profile::start("UserData",profileID);
   // Get the number of local blocks/patches on this process:
   const pargrid::CellID N_blocks = simClasses->pargrid.getNumberOfLocalCells();
   // Number of cells to be written:
   const uint64_t arraySize = N_blocks*block::SIZE;
   // Attribs map for the writer:
   map<string,string> attribs;
   attribs["mesh"] = spatMeshName;
   attribs["type"] = "celldata";
   // write selected ParGrid arrays
   /*for(auto const& p : Hybrid::varReal) {
      if(Hybrid::outputCellParams[p.first] == false) { continue; }
      const uint64_t arraySize = N_blocks*block::SIZE;
      map<string,string> attribs;
      attribs["name"] = p.first;
      attribs["mesh"] = spatMeshName;
      attribs["type"] = "celldata";
      if(simClasses->vlsv.writeArray("VARIABLE",attribs,arraySize,p.second.vectorDim,p.second.ptr) == false) { success = false; }
   }*/
   writeCellDataVariable(spatMeshName,Hybrid::dataFaceBID,                "faceB",               N_blocks,3);
   writeCellDataVariable(spatMeshName,Hybrid::dataFaceJID,                "faceJ",               N_blocks,3);
   writeCellDataVariable(spatMeshName,Hybrid::dataCellRhoQiID,            "cellRhoQi",           N_blocks,1);
#ifdef USE_BACKGROUND_CHARGE_DENSITY
   writeCellDataVariable(spatMeshName,Hybrid::dataCellRhoQiBgID,          "cellRhoQiBg",         N_blocks,1);
#endif
   writeCellDataVariable(spatMeshName,Hybrid::dataCellBID,                "cellB",               N_blocks,3);
   writeCellDataVariable(spatMeshName,Hybrid::dataCellJID,                "cellJ",               N_blocks,3);
   writeCellDataVariable(spatMeshName,Hybrid::dataCellUeID,               "cellUe",              N_blocks,3);
   writeCellDataVariable(spatMeshName,Hybrid::dataCellJiID,               "cellJi",              N_blocks,3);
   if(Hybrid::useElectronPressureElectricField == true) {
      writeCellDataVariable(spatMeshName,Hybrid::dataCellEpID,            "cellEp",              N_blocks,3);
   }
   writeCellDataVariable(spatMeshName,Hybrid::dataNodeRhoQiID,            "nodeRhoQi",           N_blocks,1);
   writeCellDataVariable(spatMeshName,Hybrid::dataNodeEID,                "nodeE",               N_blocks,3);
   writeCellDataVariable(spatMeshName,Hybrid::dataNodeBID,                "nodeB",               N_blocks,3);
   writeCellDataVariable(spatMeshName,Hybrid::dataNodeJID,                "nodeJ",               N_blocks,3);
   writeCellDataVariable(spatMeshName,Hybrid::dataNodeUeID,               "nodeUe",              N_blocks,3);
   writeCellDataVariable(spatMeshName,Hybrid::dataNodeJiID,               "nodeJi",              N_blocks,3);
#ifdef USE_RESISTIVITY
   writeCellDataVariable(spatMeshName,Hybrid::dataNodeEtaID,              "nodeEta",             N_blocks,1);
#endif
#ifdef USE_GRID_CONSTRAINT_COUNTERS
   writeCellDataVariable(spatMeshName,Hybrid::dataGridCounterCellMaxUeID,     "gridCounterCellMaxUe",    N_blocks,1);
   writeCellDataVariable(spatMeshName,Hybrid::dataGridCounterCellMaxViID,     "gridCounterCellMaxVi",    N_blocks,1);
   writeCellDataVariable(spatMeshName,Hybrid::dataGridCounterCellMinRhoQiID,  "gridCounterCellMinRhoQi", N_blocks,1);
   writeCellDataVariable(spatMeshName,Hybrid::dataGridCounterNodeMaxEID,      "gridCounterNodeMaxE",     N_blocks,1);
   writeCellDataVariable(spatMeshName,Hybrid::dataGridCounterNodeMaxVwID,     "gridCounterNodeMaxVw",    N_blocks,1);
#endif
   writeCellDataVariableBool(spatMeshName,Hybrid::dataInnerFlagFieldID,   "innerFlagField",      N_blocks,1);
   writeCellDataVariableBool(spatMeshName,Hybrid::dataInnerFlagNodeID,    "innerFlagNode",       N_blocks,1);
   writeCellDataVariableBool(spatMeshName,Hybrid::dataInnerFlagParticleID,"innerFlagParticle",   N_blocks,1);
   writeCellDataVariableBool(spatMeshName,Hybrid::dataInnerFlagCellEpID,  "innerFlagCellEp",     N_blocks,1);
#ifdef USE_OUTER_BOUNDARY_ZONE
   writeCellDataVariableBool(spatMeshName,Hybrid::dataOuterBoundaryFlagID,"outerBoundaryFlag",   N_blocks,1);
   writeCellDataVariableBool(spatMeshName,Hybrid::dataOuterBoundaryFlagNodeID,"outerBoundaryFlagNode",   N_blocks,1);
#endif
   // write production rates of ionosphere populations
   if(Hybrid::outputCellParams["prod_rate_iono"] == true) {
      Real* const cellIonosphere = reinterpret_cast<Real*>(simClasses->pargrid.getUserData(Hybrid::dataCellIonosphereID));
      for(unsigned int N_ionoPop=0;N_ionoPop<Hybrid::N_ionospherePopulations;++N_ionoPop) {
         vector<Real> iono;
         for(pargrid::CellID b=0; b<simClasses->pargrid.getNumberOfLocalCells(); ++b) {
            for(int k=0; k<block::WIDTH_Z; ++k) for(int j=0; j<block::WIDTH_Y; ++j) for(int i=0; i<block::WIDTH_X; ++i) {
               const int n = (b*block::SIZE+block::index(i,j,k));
               const size_t nIono = n*Hybrid::N_ionospherePopulations + N_ionoPop;
               iono.push_back(cellIonosphere[nIono]);
            }
         }
         attribs["name"] = string("prod_rate_iono") + to_string(N_ionoPop);
         if(simClasses->vlsv.writeArray("VARIABLE",attribs,arraySize,1,&(iono[0])) == false) { success = false; }
      }
   }
   // write production rates of exosphere populations
   if(Hybrid::outputCellParams["prod_rate_exo"] == true) {
      Real* const cellExosphere = reinterpret_cast<Real*>(simClasses->pargrid.getUserData(Hybrid::dataCellExosphereID));
      for(unsigned int N_exoPop=0;N_exoPop<Hybrid::N_exospherePopulations;++N_exoPop) {
         vector<Real> exo;
         for(pargrid::CellID b=0; b<simClasses->pargrid.getNumberOfLocalCells(); ++b) {
            for(int k=0; k<block::WIDTH_Z; ++k) for(int j=0; j<block::WIDTH_Y; ++j) for(int i=0; i<block::WIDTH_X; ++i) {
               const int n = (b*block::SIZE+block::index(i,j,k));
               const size_t nExo = n*Hybrid::N_exospherePopulations + N_exoPop;
               exo.push_back(cellExosphere[nExo]);
            }
         }
         attribs["name"] = string("prod_rate_exo") + to_string(N_exoPop);
         if(simClasses->vlsv.writeArray("VARIABLE",attribs,arraySize,1,&(exo[0])) == false) { success = false; }
      }
   }
#ifdef WRITE_GRID_TEMPORAL_AVERAGES
   Real invAveCnt = 0.0;
   if(Hybrid::gridTemporalAverageCounter > 0) {
      invAveCnt = 1.0/static_cast<Real>(Hybrid::gridTemporalAverageCounter);
   }
   // magnetic field
   if(Hybrid::outputCellParams["cellBAverage"] == true) {
      Real* cellAverageB = reinterpret_cast<Real*>(simClasses->pargrid.getUserData(Hybrid::dataCellAverageBID));
      vector<Real> averageB;
      for(pargrid::CellID b=0; b<simClasses->pargrid.getNumberOfLocalCells(); ++b) {
         for(int k=0; k<block::WIDTH_Z; ++k) for(int j=0; j<block::WIDTH_Y; ++j) for(int i=0; i<block::WIDTH_X; ++i) {
            const int n = (b*block::SIZE+block::index(i,j,k));
            const int n3 = n*3;
            averageB.push_back(cellAverageB[n3+0]*invAveCnt);
            averageB.push_back(cellAverageB[n3+1]*invAveCnt);
            averageB.push_back(cellAverageB[n3+2]*invAveCnt);
            cellAverageB[n3+0] = 0.0;
            cellAverageB[n3+1] = 0.0;
            cellAverageB[n3+2] = 0.0;
         }
      }
      attribs["name"] = "cellBAverage";
      if(simClasses->vlsv.writeArray("VARIABLE",attribs,arraySize,3,&(averageB[0])) == false) { success = false; }
   }
   // particle populations
   if(Hybrid::outputCellParams["n_ave"] == true || Hybrid::outputCellParams["v_ave"] == true || Hybrid::outputCellParams["n_tot_ave"] == true || Hybrid::outputCellParams["v_tot_ave"] == true) {
      vector<Real> averageDensityTot;
      vector<Real> averageDensityTotForVelNorm;
      vector<Real> averageVelocityTot;
      for(pargrid::CellID b=0; b<simClasses->pargrid.getNumberOfLocalCells(); ++b) {
         for(int k=0; k<block::WIDTH_Z; ++k) for(int j=0; j<block::WIDTH_Y; ++j) for(int i=0; i<block::WIDTH_X; ++i) {
            averageDensityTot.push_back(0.0);
            averageDensityTotForVelNorm.push_back(0.0);
            for(int l=0;l<3;l++) { averageVelocityTot.push_back(0.0); }
         }
      }
      for(unsigned int m=0;m<Hybrid::N_outputPopVars;++m) {
         Real* nAveArray  = reinterpret_cast<Real*>(simClasses->pargrid.getUserData(Hybrid::dataCellAverageDensityID[m]));
         Real* vAveArray  = reinterpret_cast<Real*>(simClasses->pargrid.getUserData(Hybrid::dataCellAverageVelocityID[m]));
         // use first particle population in the output particle variable vector to get species
         const Species* species = reinterpret_cast<const Species*>(particleLists[Hybrid::outputPopVarIdVector[m][0]]->getSpecies());
         const Real q = species->q;
         vector<Real> averageDensity;
         vector<Real> averageVelocity;
         for(pargrid::CellID b=0; b<simClasses->pargrid.getNumberOfLocalCells(); ++b) {
            for(int k=0; k<block::WIDTH_Z; ++k) for(int j=0; j<block::WIDTH_Y; ++j) for(int i=0; i<block::WIDTH_X; ++i) {
               const int n = (b*block::SIZE+block::index(i,j,k));
               const int n3 = n*3;
               Real nAve = nAveArray[n]/(q*Hybrid::dV)*invAveCnt;
               averageDensity.push_back(nAve);
               if(nAveArray[n] > 0) {
                  averageVelocity.push_back( vAveArray[n3+0]/nAveArray[n] );
                  averageVelocity.push_back( vAveArray[n3+1]/nAveArray[n] );
                  averageVelocity.push_back( vAveArray[n3+2]/nAveArray[n] );
               }
               else {
                  averageVelocity.push_back(0.0);
                  averageVelocity.push_back(0.0);
                  averageVelocity.push_back(0.0);
               }
               //total average plasma variables
               averageDensityTot[n] += nAve;
               averageDensityTotForVelNorm[n] += nAveArray[n];
               averageVelocityTot[n3+0] += vAveArray[n3+0];
               averageVelocityTot[n3+1] += vAveArray[n3+1];
               averageVelocityTot[n3+2] += vAveArray[n3+2];
               // zero average variables
               nAveArray[n] = 0.0;
               vAveArray[n3+0] = 0.0;
               vAveArray[n3+1] = 0.0;
               vAveArray[n3+2] = 0.0;
            }
         }
         if(Hybrid::outputCellParams["n_ave"] == true) {
            attribs["name"] = string("n_") + Hybrid::outputPopVarStr[m] + "_ave";
            if(simClasses->vlsv.writeArray("VARIABLE",attribs,arraySize,1,&(averageDensity[0])) == false) { success = false; }
         }
         if(Hybrid::outputCellParams["v_ave"] == true) {
            attribs["name"] = string("v_") + Hybrid::outputPopVarStr[m] + "_ave";
            if(simClasses->vlsv.writeArray("VARIABLE",attribs,arraySize,3,&(averageVelocity[0])) == false) { success = false; }
         }
      }
      for(pargrid::CellID b=0; b<simClasses->pargrid.getNumberOfLocalCells(); ++b) {
         for(int k=0; k<block::WIDTH_Z; ++k) for(int j=0; j<block::WIDTH_Y; ++j) for(int i=0; i<block::WIDTH_X; ++i) {
            const int n = (b*block::SIZE+block::index(i,j,k));
            const int n3 = n*3;
            if(averageDensityTotForVelNorm[n] > 0.0) {
               averageVelocityTot[n3+0] /= averageDensityTotForVelNorm[n];
               averageVelocityTot[n3+1] /= averageDensityTotForVelNorm[n];
               averageVelocityTot[n3+2] /= averageDensityTotForVelNorm[n];
            }
         }
      }
      if(Hybrid::outputCellParams["n_tot_ave"] == true) {
         attribs["name"] = string("n_tot_ave");
         if(simClasses->vlsv.writeArray("VARIABLE",attribs,arraySize,1,&(averageDensityTot[0])) == false) { success = false; }
      }
      if(Hybrid::outputCellParams["v_tot_ave"] == true) {
         attribs["name"] = string("v_tot_ave");
         if(simClasses->vlsv.writeArray("VARIABLE",attribs,arraySize,3,&(averageVelocityTot[0])) == false) { success = false; }
      }
   }
   Hybrid::gridTemporalAverageCounter = 0;
#endif

   vector<Real> divB;
   vector<Real> NPles;
#ifdef USE_B_CONSTANT
   vector<Real> B0;
   const Real* crd = getBlockCoordinateArray(*sim,*simClasses);
#endif
   for(pargrid::CellID b=0; b<simClasses->pargrid.getNumberOfLocalCells(); ++b) {
#ifdef USE_B_CONSTANT
      const size_t b3 = 3*b;
#endif
      for(int k=0; k<block::WIDTH_Z; ++k) for(int j=0; j<block::WIDTH_Y; ++j) for(int i=0; i<block::WIDTH_X; ++i) {
	 //const int n = (b*block::SIZE+block::index(i,j,k));
	 divB.push_back(0.0);
	 NPles.push_back(0.0);
#ifdef USE_B_CONSTANT
	 const Real xCellCenter = crd[b3+0] + (i+0.5)*Hybrid::dx;
	 const Real yCellCenter = crd[b3+1] + (j+0.5)*Hybrid::dx;
	 const Real zCellCenter = crd[b3+2] + (k+0.5)*Hybrid::dx;
	 Real B0_temp[3] = {0.0,0.0,0.0};
	 addConstantB(xCellCenter,yCellCenter,zCellCenter,B0_temp);
	 for(int l=0;l<3;l++) { B0.push_back(B0_temp[l]); }
#endif
      }
   }
   if(Hybrid::outputCellParams["cellDivB"] == true) {
      simClasses->pargrid.startNeighbourExchange(pargrid::DEFAULT_STENCIL,Hybrid::dataFaceBID);
      simClasses->pargrid.wait(pargrid::DEFAULT_STENCIL,Hybrid::dataFaceBID);
      Real* const faceB = reinterpret_cast<Real*>(simClasses->pargrid.getUserData(Hybrid::dataFaceBID));
      calcCellDiv(faceB,divB);
      attribs["name"] = "cellDivB";
      if(simClasses->vlsv.writeArray("VARIABLE",attribs,arraySize,1,&(divB[0])) == false) { success = false; }
   }
   if(Hybrid::outputCellParams["cellNPles"] == true) {
      calcCellNPles(NPles,particleLists);
      attribs["name"] = "cellNPles";
      if(simClasses->vlsv.writeArray("VARIABLE",attribs,arraySize,1,&(NPles[0])) == false) { success = false; }
   }
#ifdef USE_B_CONSTANT
   if(Hybrid::outputCellParams["cellB0"] == true) {
      attribs["name"] = "cellB0";
      if(simClasses->vlsv.writeArray("VARIABLE",attribs,arraySize,3,&(B0[0])) == false) { success = false; }
   }
#endif
   // particle bulk parameters (output populations)
   if(Hybrid::outputCellParams["n"] == true || Hybrid::outputCellParams["v"] == true || Hybrid::outputCellParams["T"] == true) {
      for(size_t i=0;i<Hybrid::N_outputPopVars;++i) {
         vector<Real> n,T,U;
         for(pargrid::CellID b=0; b<simClasses->pargrid.getNumberOfLocalCells(); ++b) for(int k=0; k<block::WIDTH_Z; ++k) for(int j=0; j<block::WIDTH_Y; ++j) for(int i=0; i<block::WIDTH_X; ++i) {
            n.push_back(0.0);
            T.push_back(0.0);
            U.push_back(0.0);
            U.push_back(0.0);
            U.push_back(0.0);
         }
         calcCellParticleBulkParameters(n,T,U,particleLists,Hybrid::outputPopVarIdVector[i]);
         const uint64_t arraySize = N_blocks*block::SIZE;
         if(Hybrid::outputCellParams["n"] == true) {
            attribs["name"] = string("n_") + Hybrid::outputPopVarStr[i];
            if(simClasses->vlsv.writeArray("VARIABLE",attribs,arraySize,1,&(n[0])) == false) { success = false; }
         }
         if(Hybrid::outputCellParams["T"] == true) {
            attribs["name"] = string("T_") + Hybrid::outputPopVarStr[i];
            if(simClasses->vlsv.writeArray("VARIABLE",attribs,arraySize,1,&(T[0])) == false) { success = false; }
         }
         if(Hybrid::outputCellParams["v"] == true) {
            attribs["name"] = string("v_") + Hybrid::outputPopVarStr[i];
            if(simClasses->vlsv.writeArray("VARIABLE",attribs,arraySize,3,&(U[0])) == false) { success = false; }
         }
      }
   }
   // particle bulk parameters (total plasma)
   if( (Hybrid::outputPlasmaPopId.size() > 0) &&
       ( (Hybrid::outputCellParams["n_tot"] == true) ||
         (Hybrid::outputCellParams["v_tot"] == true) ||
         (Hybrid::outputCellParams["T_tot"] == true) )) {
      vector<Real> ntot,Ttot,Utot;
      for(pargrid::CellID b=0; b<simClasses->pargrid.getNumberOfLocalCells(); ++b) for(int k=0; k<block::WIDTH_Z; ++k) for(int j=0; j<block::WIDTH_Y; ++j) for(int i=0; i<block::WIDTH_X; ++i) {
         ntot.push_back(0.0);
         Ttot.push_back(0.0);
         Utot.push_back(0.0);
         Utot.push_back(0.0);
         Utot.push_back(0.0);
      }
      calcCellParticleBulkParameters(ntot,Ttot,Utot,particleLists,Hybrid::Hybrid::outputPlasmaPopId);
      const uint64_t arraySize = N_blocks*block::SIZE;
      if(Hybrid::outputCellParams["n_tot"] == true) {
         attribs["name"] = string("n_tot");
         if(simClasses->vlsv.writeArray("VARIABLE",attribs,arraySize,1,&(ntot[0])) == false) { success = false; }
      }
      if(Hybrid::outputCellParams["T_tot"] == true) {
         attribs["name"] = string("T_tot");
         if(simClasses->vlsv.writeArray("VARIABLE",attribs,arraySize,1,&(Ttot[0])) == false) { success = false; }
      }
      if(Hybrid::outputCellParams["v_tot"] == true) {
         attribs["name"] = string("v_tot");
         if(simClasses->vlsv.writeArray("VARIABLE",attribs,arraySize,3,&(Utot[0])) == false) { success = false; }
      }
   }
#ifdef USE_DETECTORS
   // write detector cell masks at the start
   if (sim->timestep <= 0) {
      bool* detPleFlag = reinterpret_cast<bool*>(simClasses->pargrid.getUserData(Hybrid::dataDetectorParticleFlagID));
      attribs["name"] = string("detector_flag_particle");
      if(simClasses->vlsv.writeArray("VARIABLE",attribs,arraySize,1,detPleFlag) == false) { success = false; }
      bool* detBlkFlag = reinterpret_cast<bool*>(simClasses->pargrid.getUserData(Hybrid::dataDetectorBulkParamFlagID));
      attribs["name"] = string("detector_flag_bulk_param");
      if(simClasses->vlsv.writeArray("VARIABLE",attribs,arraySize,1,detBlkFlag) == false) { success = false; }
   }
#endif
   logWriteMainMacroparticles(*sim,*simClasses,particleLists);
   profile::stop();
   return success;
}

bool UserDataOP::writeCellDataVariable(const std::string& spatMeshName,pargrid::DataID& dataVarID,const std::string& dataName,const pargrid::CellID& N_blocks,const uint64_t& vectorDim) {
   if(Hybrid::outputCellParams[dataName] == false) { return true; }
   Real* const d = reinterpret_cast<Real*>(simClasses->pargrid.getUserData(dataVarID));
   if(d != NULL) {
      const uint64_t arraySize = N_blocks*block::SIZE;
      map<string,string> attribs;
      attribs["name"] = dataName;
      attribs["mesh"] = spatMeshName;
      attribs["type"] = "celldata";
      if(simClasses->vlsv.writeArray("VARIABLE",attribs,arraySize,vectorDim,d) == false) { return false; }
   }
   return true;
}

bool UserDataOP::writeCellDataVariableBool(const std::string& spatMeshName,pargrid::DataID& dataVarID,const std::string& dataName,const pargrid::CellID& N_blocks,const uint64_t& vectorDim) {
   if(Hybrid::outputCellParams[dataName] == false) { return true; }
   bool* const d = reinterpret_cast<bool*>(simClasses->pargrid.getUserData(dataVarID));
   if(d != NULL) {
      const uint64_t arraySize = N_blocks*block::SIZE;
      map<string,string> attribs;
      attribs["name"] = dataName;
      attribs["mesh"] = spatMeshName;
      attribs["type"] = "celldata";
      if(simClasses->vlsv.writeArray("VARIABLE",attribs,arraySize,vectorDim,d) == false) { return false; }
   }
   return true;
}

// calculate divergence in cells using face data
void UserDataOP::calcCellDiv(Real* faceData,vector<Real>& cellDiv) {
   for(pargrid::CellID b=0; b<simClasses->pargrid.getNumberOfLocalCells(); ++b) {
      if(simClasses->pargrid.getNeighbourFlags(b) != pargrid::ALL_NEIGHBOURS_EXIST) { continue; }
      const unsigned int size = (block::WIDTH_X+2)*(block::WIDTH_Y+2)*(block::WIDTH_Z+2);
      Real array[size*3];
      fetchData(faceData,array,*simClasses,b,3);
      for(int k=0; k<block::WIDTH_Z; ++k) for(int j=0; j<block::WIDTH_Y; ++j) for(int i=0; i<block::WIDTH_X; ++i) {
	 const int n = (b*block::SIZE+block::index(i,j,k));
	 cellDiv[n] = ((array[(block::arrayIndex(i+1,j+1,k+1))*3+0] - array[(block::arrayIndex(i+0,j+1,k+1))*3+0])
	             + (array[(block::arrayIndex(i+1,j+1,k+1))*3+1] - array[(block::arrayIndex(i+1,j+0,k+1))*3+1])
		     + (array[(block::arrayIndex(i+1,j+1,k+1))*3+2] - array[(block::arrayIndex(i+1,j+1,k+0))*3+2]))/Hybrid::dx;
      }
   }
}

// calculate number of macro particles in cells
void UserDataOP::calcCellNPles(vector<Real>& cellNPles,const vector<ParticleListBase*>& particleLists) {
   for(pargrid::CellID b=0; b<simClasses->pargrid.getNumberOfLocalCells(); ++b) {
      for(size_t s=0;s<particleLists.size();++s) {
	 // For now skip particles with invalid data id:
	 pargrid::DataID speciesDataID = pargrid::INVALID_DATAID;
	 if(particleLists[s]->getParticles(speciesDataID) == false) { continue; }
	 
	 pargrid::DataWrapper<Particle<Real> > wrapper = simClasses->pargrid.getUserDataDynamic<Particle<Real> >(speciesDataID);
	 Particle<Real>** particleList = wrapper.data();
	 Particle<Real>* particles = particleList[b];
	 pargrid::ArraySizetype N_particles = wrapper.size(b);
	 	 
	 for(size_t p=0; p<N_particles; ++p) {
	    const int i = static_cast<int>(floor(particles[p].state[particle::X]/Hybrid::dx));
	    const int j = static_cast<int>(floor(particles[p].state[particle::Y]/Hybrid::dx));
	    const int k = static_cast<int>(floor(particles[p].state[particle::Z]/Hybrid::dx));
	    const int n = (b*block::SIZE+block::index(i,j,k));
	    cellNPles[n]++;
	 }
      }
   }
}

// calculate bulk parameters of particle species s in cells
void UserDataOP::calcCellParticleBulkParameters(vector<Real>& cellDensity,vector<Real>& cellTemperature,vector<Real>& cellVelocity,const vector<ParticleListBase*>& particleLists,vector<unsigned int> S) {
   // density and velocity
   for(pargrid::CellID b=0; b<simClasses->pargrid.getNumberOfLocalCells(); ++b) {
      for(unsigned int i=0;i<S.size();++i) {
         // For now skip particles with invalid data id:
         pargrid::DataID speciesDataID = pargrid::INVALID_DATAID;
         if(particleLists[S[i]]->getParticles(speciesDataID) == false) { continue; }
         pargrid::DataWrapper<Particle<Real> > wrapper = simClasses->pargrid.getUserDataDynamic<Particle<Real> >(speciesDataID);
         Particle<Real>** particleList = wrapper.data();
         Particle<Real>* particles = particleList[b];
         pargrid::ArraySizetype N_particles = wrapper.size(b);
         for(size_t p=0; p<N_particles; ++p) {
            const int i = static_cast<int>(floor(particles[p].state[particle::X]/Hybrid::dx));
            const int j = static_cast<int>(floor(particles[p].state[particle::Y]/Hybrid::dx));
            const int k = static_cast<int>(floor(particles[p].state[particle::Z]/Hybrid::dx));
            const int n = (b*block::SIZE+block::index(i,j,k));
            const int n3 = n*3;
            cellVelocity[n3+0] += particles[p].state[particle::VX]*particles[p].state[particle::WEIGHT];
            cellVelocity[n3+1] += particles[p].state[particle::VY]*particles[p].state[particle::WEIGHT];
            cellVelocity[n3+2] += particles[p].state[particle::VZ]*particles[p].state[particle::WEIGHT];
            cellDensity[n] += particles[p].state[particle::WEIGHT];
         }
      }
   }
   for(size_t i=0; i<cellDensity.size(); ++i) {
      const int i3 = i*3;
      if(cellDensity[i] > 0.0) {
	 cellVelocity[i3+0] /= cellDensity[i];
	 cellVelocity[i3+1] /= cellDensity[i];
	 cellVelocity[i3+2] /= cellDensity[i];
      }
      else {
	 cellVelocity[i3+0] = 0.0;
	 cellVelocity[i3+1] = 0.0;
	 cellVelocity[i3+2] = 0.0;
      }
   }
   // temperature
   for(pargrid::CellID b=0; b<simClasses->pargrid.getNumberOfLocalCells(); ++b) {
      for(unsigned int i=0;i<S.size();++i) {
         // For now skip particles with invalid data id:
         pargrid::DataID speciesDataID = pargrid::INVALID_DATAID;
         if(particleLists[S[i]]->getParticles(speciesDataID) == false) { continue; }
         pargrid::DataWrapper<Particle<Real> > wrapper = simClasses->pargrid.getUserDataDynamic<Particle<Real> >(speciesDataID);
         Particle<Real>** particleList = wrapper.data();
         Particle<Real>* particles = particleList[b];
         pargrid::ArraySizetype N_particles = wrapper.size(b);
         for(size_t p=0; p<N_particles; ++p) {
            const int i = static_cast<int>(floor(particles[p].state[particle::X]/Hybrid::dx));
            const int j = static_cast<int>(floor(particles[p].state[particle::Y]/Hybrid::dx));
            const int k = static_cast<int>(floor(particles[p].state[particle::Z]/Hybrid::dx));
            const int n = (b*block::SIZE+block::index(i,j,k));
            const int n3 = n*3;
            cellTemperature[n] += particles[p].state[particle::WEIGHT]*(sqr(particles[p].state[particle::VX] - cellVelocity[n3+0])
                                                                      + sqr(particles[p].state[particle::VY] - cellVelocity[n3+1])
                                                                      + sqr(particles[p].state[particle::VZ] - cellVelocity[n3+2]));
         }
      }
   }
   const Species* species = reinterpret_cast<const Species*>(particleLists[S[0]]->getSpecies());
   for(size_t i=0; i<cellTemperature.size(); ++i) {
      if(cellDensity[i] > 0.0) {
	 // kinetic temperature
	 // T = 2/3 * 1/kB * sum_i(0.5*w_i*m_i*(v_i - U)^2)/sum_i(w_i)
	 //   = 1/3 * 1/kB * m * sum_i(w_i*(v_i - U)^2)/sum_i(w_i)
	 cellTemperature[i] *= species->m/(3*cellDensity[i]*constants::BOLTZMANN);
      }
      else { cellTemperature[i] = 0.0; }
      cellDensity[i] /= Hybrid::dV;
   }
}

// calculate particle log quantities
void logCalcParticle(Simulation& sim,SimulationClasses& simClasses,vector<LogDataParticle>& logDataParticle,const std::vector<ParticleListBase*>& particleLists,vector<Real>& cellRhoM)
{
   logDataParticle.clear();
   for(size_t s=0;s<particleLists.size();++s) {
      logDataParticle.push_back(LogDataParticle());
      logDataParticle[s].N_macroParticles = 0.0;
      logDataParticle[s].N_realParticles = 0.0;
      logDataParticle[s].sumVx = 0.0;
      logDataParticle[s].sumVy = 0.0;
      logDataParticle[s].sumVz = 0.0;
      logDataParticle[s].sumV = 0.0;
      logDataParticle[s].sumWV2 = 0.0;
      logDataParticle[s].maxVi = 0.0;
   }
   for(pargrid::CellID b=0; b<simClasses.pargrid.getNumberOfLocalCells(); ++b) {
      for(size_t s=0;s<particleLists.size();++s) {
	 pargrid::DataID speciesDataID = pargrid::INVALID_DATAID;
	 if(particleLists[s]->getParticles(speciesDataID) == false) { continue; }
	 pargrid::DataWrapper<Particle<Real> > wrapper = simClasses.pargrid.getUserDataDynamic<Particle<Real> >(speciesDataID);
	 Particle<Real>** particleList = wrapper.data();
	 Particle<Real>* particles = particleList[b];
	 pargrid::ArraySizetype N_particles = wrapper.size(b);
	 logDataParticle[s].N_macroParticles += N_particles;
	 const Species* species = reinterpret_cast<const Species*>(particleLists[s]->getSpecies());
	 for(size_t p=0; p<N_particles; ++p) {
	    logDataParticle[s].N_realParticles += particles[p].state[particle::WEIGHT];
	    logDataParticle[s].sumVx += particles[p].state[particle::WEIGHT]*particles[p].state[particle::VX];
	    logDataParticle[s].sumVy += particles[p].state[particle::WEIGHT]*particles[p].state[particle::VY];
	    logDataParticle[s].sumVz += particles[p].state[particle::WEIGHT]*particles[p].state[particle::VZ];
	    const Real v2 = sqr(particles[p].state[particle::VX]) + sqr(particles[p].state[particle::VY]) + sqr(particles[p].state[particle::VZ]);
	    const Real vtot = sqrt(v2);
	    logDataParticle[s].sumV += particles[p].state[particle::WEIGHT]*vtot;
	    logDataParticle[s].sumWV2 += particles[p].state[particle::WEIGHT]*v2;
	    if(vtot > logDataParticle[s].maxVi) { logDataParticle[s].maxVi = vtot; }
	    // determined alfven speed in each cell
	    const int i = static_cast<int>(floor(particles[p].state[particle::X]/Hybrid::dx));
	    const int j = static_cast<int>(floor(particles[p].state[particle::Y]/Hybrid::dx));
	    const int k = static_cast<int>(floor(particles[p].state[particle::Z]/Hybrid::dx));
	    const int n = (b*block::SIZE+block::index(i,j,k));
	    //const int n3 = n*3;
	    cellRhoM[n] += particles[p].state[particle::WEIGHT]*species->m;
	 }
      }
   }
   for(size_t i=0; i<cellRhoM.size(); ++i) { cellRhoM[i] /= Hybrid::dV; }
   // add background density and implement minimum densities in total mass density
   const Real minRhoMGlobal = Hybrid::minRhoQi*constants::MASS_PROTON/constants::CHARGE_ELEMENTARY; // assume global density minimum is protons
#ifdef USE_BACKGROUND_CHARGE_DENSITY
   Real* cellRhoQiBg = simClasses.pargrid.getUserDataStatic<Real>(Hybrid::dataCellRhoQiBgID);
#endif
#ifdef USE_OUTER_BOUNDARY_ZONE
   bool* outerBoundaryFlag = simClasses.pargrid.getUserDataStatic<bool>(Hybrid::dataOuterBoundaryFlagID);
   const Real minRhoMOuterBoundaryZone = Hybrid::outerBoundaryZone.minRhoQi*constants::MASS_PROTON/constants::CHARGE_ELEMENTARY; // assume outer boundary density minimum is protons
#endif
   for(pargrid::CellID b=0;b<simClasses.pargrid.getNumberOfLocalCells();++b) {
      for(int k=0;k<block::WIDTH_Z;++k) for(int j=0;j<block::WIDTH_Y;++j) for(int i=0;i<block::WIDTH_X;++i) {
	 const int n = (b*block::SIZE+block::index(i,j,k));
#ifdef USE_BACKGROUND_CHARGE_DENSITY
	 // add background density in total mass density
	 cellRhoM[n] += cellRhoQiBg[n]*constants::MASS_PROTON/constants::CHARGE_ELEMENTARY; // assume background density is protons
#endif
	 // implement minimum densities in total mass density
#ifdef USE_OUTER_BOUNDARY_ZONE
	 if(outerBoundaryFlag[n] == true) {
	    if(cellRhoM[n] < minRhoMOuterBoundaryZone) { cellRhoM[n] = minRhoMOuterBoundaryZone; }
	 }
	 else {
	    if(cellRhoM[n] < minRhoMGlobal) { cellRhoM[n] = minRhoMGlobal; }
	 }
#else
	 if(cellRhoM[n] < minRhoMGlobal) { cellRhoM[n] = minRhoMGlobal; }
#endif
      }
   }
}

// calculate field log quantities
void logCalcField(Simulation& sim,SimulationClasses& simClasses,LogDataField& logDataField,vector<Real>& cellRhoM)
{
   // synchronize faceB before using it below
   simClasses.pargrid.startNeighbourExchange(pargrid::DEFAULT_STENCIL,Hybrid::dataFaceBID);
   //profile::start("MPI waits",mpiWaitID);
   simClasses.pargrid.wait(pargrid::DEFAULT_STENCIL,Hybrid::dataFaceBID);
   //profile::stop();
   logDataField.N_cells = 0.0;
   // face magnetic field
   logDataField.sumBx = 0.0;
   logDataField.sumBy = 0.0;
   logDataField.sumBz = 0.0;
   logDataField.sumB = 0.0;
   logDataField.maxB = 0.0;
   logDataField.sumDivB = 0.0;
   logDataField.maxDivB = 0.0;
   logDataField.maxDivBPerB = 0.0;
   logDataField.sumB2 = 0.0;
   // cell ion current density
   logDataField.sumCellJix = 0.0;
   logDataField.sumCellJiy = 0.0;
   logDataField.sumCellJiz = 0.0;
   logDataField.sumCellJi = 0.0;
   logDataField.maxCellJi = 0.0;
   // cell electron pressure electric field
   logDataField.sumCellEpx = 0.0;
   logDataField.sumCellEpy = 0.0;
   logDataField.sumCellEpz = 0.0;
   logDataField.sumCellEp = 0.0;
   logDataField.maxCellEp = 0.0;
   logDataField.sumCellEp2 = 0.0;
   // electric field
   logDataField.sumNodeEx = 0.0;
   logDataField.sumNodeEy = 0.0;
   logDataField.sumNodeEz = 0.0;
   logDataField.sumNodeE = 0.0;
   logDataField.maxNodeE = 0.0;
   logDataField.sumNodeE2 = 0.0;
   logDataField.sumCellE2 = 0.0;
   // spatial, temporal and velocity scales
   logDataField.minInerLengthElectron = numeric_limits<Real>::max();
   logDataField.maxInerLengthElectron = numeric_limits<Real>::min();
   logDataField.minInerLengthProton = numeric_limits<Real>::max();
   logDataField.maxInerLengthProton = numeric_limits<Real>::min();
   logDataField.minTLarmor = numeric_limits<Real>::max();
   logDataField.maxVAlfven = 0.0;
   logDataField.maxUe = 0.0;
   Real* faceB = simClasses.pargrid.getUserDataStatic<Real>(Hybrid::dataFaceBID);
   Real* cellRhoQi = simClasses.pargrid.getUserDataStatic<Real>(Hybrid::dataCellRhoQiID);
   Real* cellJi = simClasses.pargrid.getUserDataStatic<Real>(Hybrid::dataCellJiID);
   Real* cellUe = simClasses.pargrid.getUserDataStatic<Real>(Hybrid::dataCellUeID);
   Real* cellEp = simClasses.pargrid.getUserDataStatic<Real>(Hybrid::dataCellEpID);
   Real* nodeE  = simClasses.pargrid.getUserDataStatic<Real>(Hybrid::dataNodeEID);
   bool* innerFlag = simClasses.pargrid.getUserDataStatic<bool>(Hybrid::dataInnerFlagFieldID);
   for(pargrid::CellID b=0; b<simClasses.pargrid.getNumberOfLocalCells(); ++b) {
      if(simClasses.pargrid.getNeighbourFlags(b) != pargrid::ALL_NEIGHBOURS_EXIST) { continue; }
      const unsigned int size = (block::WIDTH_X+2)*(block::WIDTH_Y+2)*(block::WIDTH_Z+2);
      Real allFaceB[size*3];
      fetchData(faceB,allFaceB,simClasses,b,3);
      for(int k=0; k<block::WIDTH_Z; ++k) for(int j=0; j<block::WIDTH_Y; ++j) for(int i=0; i<block::WIDTH_X; ++i) {
	 const int n = (b*block::SIZE+block::index(i,j,k));
	 const int n3 = n*3;
	 // do not include cells inside the inner boundary
	 if(Hybrid::includeInnerCellsInFieldLog == false) {
	    if(innerFlag[n] == true) { continue; }
	 }
	 // divergence of B in a cell from face magnetic field
	 Real divB = ((allFaceB[(block::arrayIndex(i+1,j+1,k+1))*3+0] - allFaceB[(block::arrayIndex(i+0,j+1,k+1))*3+0])
		    + (allFaceB[(block::arrayIndex(i+1,j+1,k+1))*3+1] - allFaceB[(block::arrayIndex(i+1,j+0,k+1))*3+1])
		    + (allFaceB[(block::arrayIndex(i+1,j+1,k+1))*3+2] - allFaceB[(block::arrayIndex(i+1,j+1,k+0))*3+2]))/Hybrid::dx;
	 divB = fabs(divB);
	 // B in a cell as an average face magnetic field
	 const Real Bx = 0.5*(allFaceB[(block::arrayIndex(i+1,j+1,k+1))*3+0] + allFaceB[(block::arrayIndex(i+0,j+1,k+1))*3+0]);
	 const Real By = 0.5*(allFaceB[(block::arrayIndex(i+1,j+1,k+1))*3+1] + allFaceB[(block::arrayIndex(i+1,j+0,k+1))*3+1]);
	 const Real Bz = 0.5*(allFaceB[(block::arrayIndex(i+1,j+1,k+1))*3+2] + allFaceB[(block::arrayIndex(i+1,j+1,k+0))*3+2]);
	 const Real B2 = sqr(Bx) + sqr(By) + sqr(Bz);
	 const Real Btot = sqrt(B2);
	 Real divBPerB = 0.0;
	 if(Btot > 0.0) { divBPerB = divB/Btot; }
	 logDataField.sumBx += Bx;
	 logDataField.sumBy += By;
	 logDataField.sumBz += Bz;
	 logDataField.sumB += Btot;
	 logDataField.sumDivB += divB;
	 logDataField.sumB2 += B2;
	 if(Btot > logDataField.maxB) { logDataField.maxB = Btot; }
	 if(fabs(divB) > fabs(logDataField.maxDivB)) { logDataField.maxDivB = divB; }
	 if(fabs(divBPerB) > fabs(logDataField.maxDivBPerB)) { logDataField.maxDivBPerB = divBPerB; }

	 // cell ion current density
	 const Real cellJitot = sqrt( sqr(cellJi[n3+0]) + sqr(cellJi[n3+1]) + sqr(cellJi[n3+2]) );
	 logDataField.sumCellJix += cellJi[n3+0];
	 logDataField.sumCellJiy += cellJi[n3+1];
	 logDataField.sumCellJiz += cellJi[n3+2];
	 logDataField.sumCellJi += cellJitot;
	 if(cellJitot > logDataField.maxCellJi) { logDataField.maxCellJi = cellJitot; }

	 // cell electron pressure electric field
	 const Real cellEp2 = sqr(cellEp[n3+0]) + sqr(cellEp[n3+1]) + sqr(cellEp[n3+2]);
	 const Real cellEptot = sqrt(cellEp2);
	 logDataField.sumCellEpx += cellEp[n3+0];
	 logDataField.sumCellEpy += cellEp[n3+1];
	 logDataField.sumCellEpz += cellEp[n3+2];
	 logDataField.sumCellEp += cellEptot;
	 logDataField.sumCellEp2 += cellEp2;
	 if(cellEptot > logDataField.maxCellEp) { logDataField.maxCellEp = cellEptot; }

	 // node electric field -Ue x B + eta*J
	 const Real nodeE2 = sqr(nodeE[n3+0]) + sqr(nodeE[n3+1]) + sqr(nodeE[n3+2]);
	 const Real nodeEtot = sqrt(nodeE2);
	 logDataField.sumNodeEx += nodeE[n3+0];
	 logDataField.sumNodeEy += nodeE[n3+1];
	 logDataField.sumNodeEz += nodeE[n3+2];
	 logDataField.sumNodeE += nodeEtot;
	 logDataField.sumNodeE2 += nodeE2;
	 if(nodeEtot > logDataField.maxNodeE) { logDataField.maxNodeE = nodeEtot; }

	 // cell electric field -Ue x B
	 const Real cellUex = cellUe[n3+0];
	 const Real cellUey = cellUe[n3+1];
	 const Real cellUez = cellUe[n3+2];
	 const Real cellEx = -(cellUey*Bz - cellUez*By);
	 const Real cellEy = -(cellUez*Bx - cellUex*Bz);
	 const Real cellEz = -(cellUex*By - cellUey*Bx);
	 const Real cellE2 = sqr(cellEx) + sqr(cellEy) + sqr(cellEz);
	 logDataField.sumCellE2 += cellE2;

	 // spatial, temporal and velocity scales
	 const Real cellUetot = sqrt(sqr(cellUex) + sqr(cellUey) + sqr(cellUez));
	 if(cellRhoQi[n] > 0.0) {
	    // electron inertial length: sqrt(m_e/(mu0*q_e^2*n_e)) = sqrt(m_e/(mu0*q_e^2*rho_q/q_e)) =  = sqrt(m_e/(mu0*q_e*rho_q))
	    const Real de = sqrt(constants::MASS_ELECTRON/(constants::PERMEABILITY*constants::CHARGE_ELEMENTARY*cellRhoQi[n]));
	    if(de > logDataField.maxInerLengthElectron) { logDataField.maxInerLengthElectron = de; }
	    if(de < logDataField.minInerLengthElectron) { logDataField.minInerLengthElectron = de; }
	 }
	 if(Btot > 0.0) {
	    // Minimum ion (=proton) Larmor period
	    const Real tL = 2.0*M_PI*constants::MASS_PROTON/(constants::CHARGE_ELEMENTARY*Btot);
	    if(tL < logDataField.minTLarmor) { logDataField.minTLarmor = tL; }
	 }
	 if(cellRhoM[n] > 0.0) {
	    // ion inertial length assuming that all particles in total mass density are protons: sqrt(m_p/(mu0*q_e^2*n_i)) = sqrt(m_p/(mu0*q_e^2*(rho_m/m_p))) = sqrt(m_p^2/(mu0*q_e^2*rho_m))
	    const Real di = sqrt(sqr(constants::MASS_PROTON)/(constants::PERMEABILITY*sqr(constants::CHARGE_ELEMENTARY)*cellRhoM[n]));
	    const Real vA = Btot/sqrt(constants::PERMEABILITY*cellRhoM[n]);
	    if(di > logDataField.maxInerLengthProton) { logDataField.maxInerLengthProton = di; }
	    if(di < logDataField.minInerLengthProton) { logDataField.minInerLengthProton = di; }
	    if(vA > logDataField.maxVAlfven) { logDataField.maxVAlfven = vA; }
	 }
	 if(cellUetot > logDataField.maxUe) { logDataField.maxUe = cellUetot; }

	 // update field log cell counter
	 logDataField.N_cells++;
      }
   }
}

// write particle population and field logs
bool logWriteParticleField(Simulation& sim,SimulationClasses& simClasses,const std::vector<ParticleListBase*>& particleLists) {
   bool success = true;
   static int profWriteLogsID = -1;
   //if(getInitialized() == false) { return false; }
   profile::start("logWriteParticleField",profWriteLogsID);
   vector<LogDataParticle> logDataParticle;
   LogDataField logDataField;
   // total mass density for Alfven speed
   vector<Real> cellRhoM;
   for(pargrid::CellID b=0; b<simClasses.pargrid.getNumberOfLocalCells(); ++b) for(int k=0; k<block::WIDTH_Z; ++k) for(int j=0; j<block::WIDTH_Y; ++j) for(int i=0; i<block::WIDTH_X; ++i) {
      cellRhoM.push_back(0.0);
   }
   logCalcParticle(sim,simClasses,logDataParticle,particleLists,cellRhoM);
   logCalcField(sim,simClasses,logDataField,cellRhoM);
   if(sim.mpiRank==sim.MASTER_RANK) {
      for(size_t i=0;i<Hybrid::logParticle.size();++i) {
	 (*Hybrid::logParticle[i]) << sim.t << " ";
      }
      Hybrid::logField << sim.t << " ";
   }
   Real maxViAllPopulations = 0.0;
   const Real Dt = sim.t - Hybrid::logCounterTimeStart;
   // go thru populations and sum particle counters
   for(size_t s=0;s<particleLists.size();++s) {
      Real N_macroParticlesThisProcess = logDataParticle[s].N_macroParticles;
      Real N_macroParticlesGlobal = 0.0;
      Real N_realParticlesThisProcess = logDataParticle[s].N_realParticles;
      Real N_realParticlesGlobal = 0.0;
      Real sumVxThisProcess = logDataParticle[s].sumVx;
      Real sumVxGlobal = 0.0;
      Real sumVyThisProcess = logDataParticle[s].sumVy;
      Real sumVyGlobal = 0.0;
      Real sumVzThisProcess = logDataParticle[s].sumVz;
      Real sumVzGlobal = 0.0;
      Real sumVThisProcess = logDataParticle[s].sumV;
      Real sumVGlobal = 0.0;
      Real sumWV2ThisProcess = logDataParticle[s].sumWV2;
      Real sumWV2Global = 0.0;
      Real maxViThisProcess = logDataParticle[s].maxVi;
      Real maxViGlobal = 0.0;
      Real sumEscapeThisProcess = Hybrid::logCounterParticleEscape[s];
      Real sumEscapeGlobal = 0.0;
      Real sumImpactThisProcess = Hybrid::logCounterParticleImpact[s];
      Real sumImpactGlobal = 0.0;
      Real sumInjectThisProcess = Hybrid::logCounterParticleInject[s];
      Real sumInjectGlobal = 0.0;
      Real sumInjectMacroThisProcess = Hybrid::logCounterParticleInjectMacroparticles[s];
      Real sumInjectMacroGlobal = 0.0;
      Real sumEscapeKineticEnergyThisProcess = Hybrid::logCounterParticleEscapeKineticEnergy[s];
      Real sumEscapeKineticEnergyGlobal = 0.0;
      Real sumImpactKineticEnergyThisProcess = Hybrid::logCounterParticleImpactKineticEnergy[s];
      Real sumImpactKineticEnergyGlobal = 0.0;
      Real sumInjectKineticEnergyThisProcess = Hybrid::logCounterParticleInjectKineticEnergy[s];
      Real sumInjectKineticEnergyGlobal = 0.0;
      Real sumMaxViThisProcess = Hybrid::logCounterParticleMaxVi[s];
      Real sumMaxViGlobal = 0.0;
      MPI_Reduce(&N_macroParticlesThisProcess,&N_macroParticlesGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
      MPI_Reduce(&N_realParticlesThisProcess,&N_realParticlesGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
      MPI_Reduce(&sumVxThisProcess,&sumVxGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
      MPI_Reduce(&sumVyThisProcess,&sumVyGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
      MPI_Reduce(&sumVzThisProcess,&sumVzGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
      MPI_Reduce(&sumVThisProcess,&sumVGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
      MPI_Reduce(&sumWV2ThisProcess,&sumWV2Global,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
      MPI_Reduce(&maxViThisProcess,&maxViGlobal,1,MPI_Type<Real>(),MPI_MAX,sim.MASTER_RANK,sim.comm);
      MPI_Reduce(&sumEscapeThisProcess,&sumEscapeGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
      MPI_Reduce(&sumImpactThisProcess,&sumImpactGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
      MPI_Reduce(&sumInjectThisProcess,&sumInjectGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
      MPI_Reduce(&sumInjectMacroThisProcess,&sumInjectMacroGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
      MPI_Reduce(&sumEscapeKineticEnergyThisProcess,&sumEscapeKineticEnergyGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
      MPI_Reduce(&sumImpactKineticEnergyThisProcess,&sumImpactKineticEnergyGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
      MPI_Reduce(&sumInjectKineticEnergyThisProcess,&sumInjectKineticEnergyGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
      MPI_Reduce(&sumMaxViThisProcess,&sumMaxViGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
      if(sim.mpiRank==sim.MASTER_RANK) {
	 (*Hybrid::logParticle[s]) << N_realParticlesGlobal << " " << N_macroParticlesGlobal << " ";
	 if(N_realParticlesGlobal > 0.0) {
	    (*Hybrid::logParticle[s])
	      << sumVxGlobal/N_realParticlesGlobal << " "
	      << sumVyGlobal/N_realParticlesGlobal << " "
	      << sumVzGlobal/N_realParticlesGlobal << " "
	      << sumVGlobal/N_realParticlesGlobal  << " ";
	 }
	 else {
	    (*Hybrid::logParticle[s]) << 0.0 << " " << 0.0 << " " << 0.0 << " " << 0.0 << " ";
	 }
	 const Species* species = reinterpret_cast<const Species*>(particleLists[s]->getSpecies());
	 (*Hybrid::logParticle[s]) << 0.5*species->m*sumWV2Global << " " << maxViGlobal << " ";
	 if(maxViAllPopulations < maxViGlobal) { maxViAllPopulations = maxViGlobal; }
	 if(Dt > 0) {
	    (*Hybrid::logParticle[s])
	      << sumEscapeGlobal/Dt << " "
	      << sumImpactGlobal/Dt << " "
	      << sumInjectGlobal/Dt << " "
	      << sumInjectMacroGlobal/Dt*sim.dt << " "
	      << 0.5*species->m*sumEscapeKineticEnergyGlobal/Dt << " "
	      << 0.5*species->m*sumImpactKineticEnergyGlobal/Dt << " "
	      << 0.5*species->m*sumInjectKineticEnergyGlobal/Dt << " "
	      << sumMaxViGlobal/Dt*sim.dt << " ";
	 }
	 else {
	    (*Hybrid::logParticle[s]) << 0.0 << " " << 0.0 << " " << 0.0 << " " << 0.0 << " " << 0.0 << " " << 0.0 << " " << 0.0 << " " << 0.0 << " ";
	 }
      }
   }

   // sum global field counters
   Real sumMaxCellUeThisProcess = Hybrid::logCounterFieldMaxCellUe;
   Real sumMaxCellUeGlobal = 0.0;
   Real sumMaxNodeUeThisProcess = Hybrid::logCounterFieldMaxNodeUe;
   Real sumMaxNodeUeGlobal = 0.0;
   Real sumMaxVwThisProcess = Hybrid::logCounterFieldMaxNodeUe;
   Real sumMaxVwGlobal = 0.0;
   Real sumMaxEThisProcess = Hybrid::logCounterFieldMaxE;
   Real sumMaxEGlobal = 0.0;
   Real sumMinCellRhoQiThisProcess = Hybrid::logCounterFieldMinCellRhoQi;
   Real sumMinCellRhoQiGlobal = 0.0;
   Real sumMinNodeRhoQiThisProcess = Hybrid::logCounterFieldMinNodeRhoQi;
   Real sumMinNodeRhoQiGlobal = 0.0;
   MPI_Reduce(&sumMaxCellUeThisProcess,&sumMaxCellUeGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
   MPI_Reduce(&sumMaxNodeUeThisProcess,&sumMaxNodeUeGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
   MPI_Reduce(&sumMaxVwThisProcess,&sumMaxVwGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
   MPI_Reduce(&sumMaxEThisProcess,&sumMaxEGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
   MPI_Reduce(&sumMinCellRhoQiThisProcess,&sumMinCellRhoQiGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
   MPI_Reduce(&sumMinNodeRhoQiThisProcess,&sumMinNodeRhoQiGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
   
   // field log values from logCalcField
   Real N_cellsThisProcess = logDataField.N_cells;
   Real N_cellsGlobal = 0.0;
   MPI_Reduce(&N_cellsThisProcess,&N_cellsGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
   // face magnetic field
   Real sumBxThisProcess = logDataField.sumBx;
   Real sumBxGlobal = 0.0;
   Real sumByThisProcess = logDataField.sumBy;
   Real sumByGlobal = 0.0;
   Real sumBzThisProcess = logDataField.sumBz;
   Real sumBzGlobal = 0.0;
   Real sumBThisProcess = logDataField.sumB;
   Real sumBGlobal = 0.0;
   Real maxBThisProcess = logDataField.maxB;
   Real maxBGlobal = 0.0;
   Real sumDivBThisProcess = logDataField.sumDivB;
   Real sumDivBGlobal = 0.0;
   Real maxDivBThisProcess = logDataField.maxDivB;
   Real maxDivBGlobal = 0.0;
   Real maxDivPerBThisProcess = logDataField.maxDivBPerB;
   Real maxDivPerBGlobal = 0.0;
   Real sumB2ThisProcess = logDataField.sumB2;
   Real sumB2Global = 0.0;
   MPI_Reduce(&sumBxThisProcess,&sumBxGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
   MPI_Reduce(&sumByThisProcess,&sumByGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
   MPI_Reduce(&sumBzThisProcess,&sumBzGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
   MPI_Reduce(&sumBThisProcess,&sumBGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
   MPI_Reduce(&maxBThisProcess,&maxBGlobal,1,MPI_Type<Real>(),MPI_MAX,sim.MASTER_RANK,sim.comm);
   MPI_Reduce(&sumDivBThisProcess,&sumDivBGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
   MPI_Reduce(&maxDivBThisProcess,&maxDivBGlobal,1,MPI_Type<Real>(),MPI_MAX,sim.MASTER_RANK,sim.comm);
   MPI_Reduce(&maxDivPerBThisProcess,&maxDivPerBGlobal,1,MPI_Type<Real>(),MPI_MAX,sim.MASTER_RANK,sim.comm);
   MPI_Reduce(&sumB2ThisProcess,&sumB2Global,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);

   // cell ion current density
   Real sumCellJixThisProcess = logDataField.sumCellJix;
   Real sumCellJixGlobal = 0.0;
   Real sumCellJiyThisProcess = logDataField.sumCellJiy;
   Real sumCellJiyGlobal = 0.0;
   Real sumCellJizThisProcess = logDataField.sumCellJiz;
   Real sumCellJizGlobal = 0.0;
   Real sumCellJiThisProcess = logDataField.sumCellJi;
   Real sumCellJiGlobal = 0.0;
   Real maxCellJiThisProcess = logDataField.maxCellJi;
   Real maxCellJiGlobal = 0.0;
   MPI_Reduce(&sumCellJixThisProcess,&sumCellJixGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
   MPI_Reduce(&sumCellJiyThisProcess,&sumCellJiyGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
   MPI_Reduce(&sumCellJizThisProcess,&sumCellJizGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
   MPI_Reduce(&sumCellJiThisProcess,&sumCellJiGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
   MPI_Reduce(&maxCellJiThisProcess,&maxCellJiGlobal,1,MPI_Type<Real>(),MPI_MAX,sim.MASTER_RANK,sim.comm);

   // cell electron pressure electric field
   Real sumCellEpxThisProcess = logDataField.sumCellEpx;
   Real sumCellEpxGlobal = 0.0;
   Real sumCellEpyThisProcess = logDataField.sumCellEpy;
   Real sumCellEpyGlobal = 0.0;
   Real sumCellEpzThisProcess = logDataField.sumCellEpz;
   Real sumCellEpzGlobal = 0.0;
   Real sumCellEpThisProcess = logDataField.sumCellEp;
   Real sumCellEpGlobal = 0.0;
   Real maxCellEpThisProcess = logDataField.maxCellEp;
   Real maxCellEpGlobal = 0.0;
   Real sumCellEp2ThisProcess = logDataField.sumCellEp2;
   Real sumCellEp2Global = 0.0;
   MPI_Reduce(&sumCellEpxThisProcess,&sumCellEpxGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
   MPI_Reduce(&sumCellEpyThisProcess,&sumCellEpyGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
   MPI_Reduce(&sumCellEpzThisProcess,&sumCellEpzGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
   MPI_Reduce(&sumCellEpThisProcess,&sumCellEpGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
   MPI_Reduce(&maxCellEpThisProcess,&maxCellEpGlobal,1,MPI_Type<Real>(),MPI_MAX,sim.MASTER_RANK,sim.comm);
   MPI_Reduce(&sumCellEp2ThisProcess,&sumCellEp2Global,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);

   // node electric field
   Real sumNodeExThisProcess = logDataField.sumNodeEx;
   Real sumNodeExGlobal = 0.0;
   Real sumNodeEyThisProcess = logDataField.sumNodeEy;
   Real sumNodeEyGlobal = 0.0;
   Real sumNodeEzThisProcess = logDataField.sumNodeEz;
   Real sumNodeEzGlobal = 0.0;
   Real sumNodeEThisProcess = logDataField.sumNodeE;
   Real sumNodeEGlobal = 0.0;
   Real maxNodeEThisProcess = logDataField.maxNodeE;
   Real maxNodeEGlobal = 0.0;
   Real sumNodeE2ThisProcess = logDataField.sumNodeE2;
   Real sumNodeE2Global = 0.0;
   Real sumCellE2ThisProcess = logDataField.sumCellE2;
   Real sumCellE2Global = 0.0;
   MPI_Reduce(&sumNodeExThisProcess,&sumNodeExGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
   MPI_Reduce(&sumNodeEyThisProcess,&sumNodeEyGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
   MPI_Reduce(&sumNodeEzThisProcess,&sumNodeEzGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
   MPI_Reduce(&sumNodeEThisProcess,&sumNodeEGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
   MPI_Reduce(&maxNodeEThisProcess,&maxNodeEGlobal,1,MPI_Type<Real>(),MPI_MAX,sim.MASTER_RANK,sim.comm);
   MPI_Reduce(&sumNodeE2ThisProcess,&sumNodeE2Global,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
   MPI_Reduce(&sumCellE2ThisProcess,&sumCellE2Global,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);

   // spatial, temporal and velocity scales
   Real minInerLengthElectronThisProcess = logDataField.minInerLengthElectron;
   Real minInerLengthElectronGlobal = 0.0;
   Real maxInerLengthElectronThisProcess = logDataField.maxInerLengthElectron;
   Real maxInerLengthElectronGlobal = 0.0;
   Real minInerLengthProtonThisProcess = logDataField.minInerLengthProton;
   Real minInerLengthProtonGlobal = 0.0;
   Real maxInerLengthProtonThisProcess = logDataField.maxInerLengthProton;
   Real maxInerLengthProtonGlobal = 0.0;
   Real minTLarmorThisProcess = logDataField.minTLarmor;
   Real minTLarmorGlobal = 0.0;
   Real maxVAlfvenThisProcess = logDataField.maxVAlfven;
   Real maxVAlfvenGlobal = 0.0;
   Real maxUeThisProcess = logDataField.maxUe;
   Real maxUeGlobal = 0.0;
   MPI_Reduce(&minInerLengthElectronThisProcess,&minInerLengthElectronGlobal,1,MPI_Type<Real>(),MPI_MIN,sim.MASTER_RANK,sim.comm);
   MPI_Reduce(&maxInerLengthElectronThisProcess,&maxInerLengthElectronGlobal,1,MPI_Type<Real>(),MPI_MAX,sim.MASTER_RANK,sim.comm);
   MPI_Reduce(&minInerLengthProtonThisProcess,&minInerLengthProtonGlobal,1,MPI_Type<Real>(),MPI_MIN,sim.MASTER_RANK,sim.comm);
   MPI_Reduce(&maxInerLengthProtonThisProcess,&maxInerLengthProtonGlobal,1,MPI_Type<Real>(),MPI_MAX,sim.MASTER_RANK,sim.comm);
   MPI_Reduce(&minTLarmorThisProcess,&minTLarmorGlobal,1,MPI_Type<Real>(),MPI_MIN,sim.MASTER_RANK,sim.comm);
   MPI_Reduce(&maxVAlfvenThisProcess,&maxVAlfvenGlobal,1,MPI_Type<Real>(),MPI_MAX,sim.MASTER_RANK,sim.comm);
   MPI_Reduce(&maxUeThisProcess,&maxUeGlobal,1,MPI_Type<Real>(),MPI_MAX,sim.MASTER_RANK,sim.comm);

   // write field log
   if(sim.mpiRank==sim.MASTER_RANK) {
      // face magnetic field
      if(N_cellsGlobal > 0) {
	 Hybrid::logField
	   << sumBxGlobal/N_cellsGlobal << " "
	   << sumByGlobal/N_cellsGlobal << " "
	   << sumBzGlobal/N_cellsGlobal << " "
	   << sumBGlobal/N_cellsGlobal << " ";
      }
      else {
	 Hybrid::logField << 0.0 << " " << 0.0 << " " << 0.0 << " " << 0.0 << " ";
      }
      Hybrid::logField << maxBGlobal << " ";
      if(N_cellsGlobal > 0) { Hybrid::logField << sumDivBGlobal/N_cellsGlobal << " "; }
      else { Hybrid::logField << 0.0 << " "; }
      Hybrid::logField << maxDivBGlobal << " " << Hybrid::dx*maxDivPerBGlobal << " " << sumB2Global*Hybrid::dV/(2.0*constants::PERMEABILITY) << " ";

      // cell ion current density
      if(N_cellsGlobal > 0) {
	 Hybrid::logField
	   << sumCellJixGlobal/N_cellsGlobal << " "
	   << sumCellJiyGlobal/N_cellsGlobal << " "
	   << sumCellJizGlobal/N_cellsGlobal << " "
	   << sumCellJiGlobal/N_cellsGlobal << " ";
      }
      else {
	 Hybrid::logField << 0.0 << " " << 0.0 << " " << 0.0 << " " << 0.0 << " ";
      }
      Hybrid::logField << maxCellJiGlobal << " ";

      // cell electron pressure electric field
      if(N_cellsGlobal > 0) {
	 Hybrid::logField
	   << sumCellEpxGlobal/N_cellsGlobal << " "
	   << sumCellEpyGlobal/N_cellsGlobal << " "
	   << sumCellEpzGlobal/N_cellsGlobal << " "
	   << sumCellEpGlobal/N_cellsGlobal << " ";
      }
      else {
	 Hybrid::logField << 0.0 << " " << 0.0 << " " << 0.0 << " " << 0.0 << " ";
      }
      Hybrid::logField << maxCellEpGlobal << " " << sumCellEp2Global*Hybrid::dV*constants::PERMITTIVITY/(2.0) << " ";

      // node electric field
      if(N_cellsGlobal > 0) {
	 Hybrid::logField
	   << sumNodeExGlobal/N_cellsGlobal << " "
	   << sumNodeEyGlobal/N_cellsGlobal << " "
	   << sumNodeEzGlobal/N_cellsGlobal << " "
	   << sumNodeEGlobal/N_cellsGlobal << " ";
      }
      else {
	 Hybrid::logField << 0.0 << " " << 0.0 << " " << 0.0 << " " << 0.0 << " ";
      }
      Hybrid::logField << " " << maxNodeEGlobal << " " << sumNodeE2Global*Hybrid::dV*constants::PERMITTIVITY/(2.0)  << " " << sumCellE2Global*Hybrid::dV*constants::PERMITTIVITY/(2.0)  << " ";

      // global field counters
      if(Dt > 0) {
	 Hybrid::logField
	   << sumMaxCellUeGlobal/Dt*sim.dt << " "
	   << sumMaxNodeUeGlobal/Dt*sim.dt << " "
	   << sumMaxVwGlobal/Dt*sim.dt << " "
	   << sumMaxEGlobal/Dt*sim.dt << " "
	   << sumMinCellRhoQiGlobal/Dt*sim.dt << " "
	   << sumMinNodeRhoQiGlobal/Dt*sim.dt << " ";
      }
      else {
	 Hybrid::logField << 0.0 << " " << 0.0 << " " << 0.0 << " " << 0.0 << " " << 0.0 << " " << 0.0 << " ";
      }

      // spatial, temporal and velocity scales
      Hybrid::logField
	<< minInerLengthElectronGlobal << " "
	<< maxInerLengthElectronGlobal << " "
	<< minInerLengthProtonGlobal << " "
	<< maxInerLengthProtonGlobal << " "
	<< minTLarmorGlobal << " "
	<< maxVAlfvenGlobal << " "
	<< maxUeGlobal << " ";

      // line endings particle logs
      for(size_t i=0;i<Hybrid::logParticle.size();++i) {
	 (*Hybrid::logParticle[i]) << endl;
      }
      // line ending field log
      Hybrid::logField << endl;

      // warnings
      const Real tL_min_dt = minTLarmorGlobal/sim.dt;
      const Real maxVi_dxdt = maxViAllPopulations/(Hybrid::dx/sim.dt);
      const Real maxVA_dxdt = maxVAlfvenGlobal/(Hybrid::dx/sim.dt);
      const Real maxUe_dxdt = maxUeGlobal/(Hybrid::dx/sim.dt);

      // write if on a save step or save step already happened after previous entry
      if(Hybrid::writeMainLogEntriesAfterSaveStep == true && sim.mpiRank == sim.MASTER_RANK) {
	 simClasses.logger
	   << "(RHYBRID) Max. |Vion|             : Vi_max  = " << maxViAllPopulations/1e3 << " km/s = " << maxVi_dxdt << " dx/dt" << endl
	   << "(RHYBRID) Max. |Ue|               : Ue_max  = " << maxUeGlobal/1e3 << " km/s = " << maxUe_dxdt << " dx/dt" << endl
	   << "(RHYBRID) Max. |VAlfven|          : Va_max  = " << maxVAlfvenGlobal/1e3 << " km/s = " << maxVA_dxdt << " dx/dt" << endl
	   << "(RHYBRID) Max. |faceB|            : B_max   = " << maxBGlobal/1e-9 << " nT" << endl
	   << "(RHYBRID) Min. ion Larmor period  : tL_min  = " << minTLarmorGlobal << " s = " << minTLarmorGlobal/sim.dt << " dt" << endl
	   << "(RHYBRID) Max. |cellJi|           : Ji_max  = " << maxCellJiGlobal << " A/m^2" << endl
	   << "(RHYBRID) Max. |cellEp|           : Ep_max  = " << maxCellEpGlobal/1e-3 << " mV/m" << endl
	   << "(RHYBRID) Max. |nodeE|            : E_max   = " << maxNodeEGlobal/1e-3 << " mV/m" << endl
	   << "(RHYBRID) Min. e- inertial length : de_min  = " << minInerLengthElectronGlobal/1e3 << " km = " << minInerLengthElectronGlobal/Hybrid::dx << " dx" << endl
	   << "(RHYBRID) Max. e- inertial length : de_max  = " << maxInerLengthElectronGlobal/1e3 << " km = " << maxInerLengthElectronGlobal/Hybrid::dx << " dx" << endl
	   << "(RHYBRID) Min. H+ inertial length : di_min  = " << minInerLengthProtonGlobal/1e3 << " km = " << minInerLengthProtonGlobal/Hybrid::dx << " dx" << endl
	   << "(RHYBRID) Max. H+ inertial length : di_max  = " << maxInerLengthProtonGlobal/1e3 << " km = " << maxInerLengthProtonGlobal/Hybrid::dx << " dx" << endl;
	 Hybrid::writeMainLogEntriesAfterSaveStep = false;
      }
      if(sim.mpiRank == sim.MASTER_RANK) {
	 if(tL_min_dt < 10)   { simClasses.logger << "(RHYBRID) WARNING: Minimum Larmor period: tL_min/dt < 10 ("      << tL_min_dt  << "), time step = " << sim.timestep << ", time = " << sim.t << " s" << endl; }
	 if(maxVi_dxdt > 0.9) { simClasses.logger << "(RHYBRID) WARNING: Maximum ion speed: Vi_max/(dx/dt) > 0.9 ("    << maxVi_dxdt << "), time step = " << sim.timestep << ", time = " << sim.t << " s" << endl; }
	 if(maxUe_dxdt > 0.9) { simClasses.logger << "(RHYBRID) WARNING: Maximum ion speed: Ue_max/(dx/dt) > 0.9 ("    << maxUe_dxdt << "), time step = " << sim.timestep << ", time = " << sim.t << " s" << endl; }
	 if(maxVA_dxdt > 0.9) { simClasses.logger << "(RHYBRID) WARNING: Maximum Alfven speed: Va_max/(dx/dt) > 0.9 (" << maxVA_dxdt << "), time step = " << sim.timestep << ", time = " << sim.t << " s" << endl; }
      }
      if(maxBGlobal > Hybrid::terminateLimitMaxB) {
         success = false;
	 if(sim.mpiRank == sim.MASTER_RANK) { simClasses.logger << "(RHYBRID) CONSTRAINT: maximum |B| for run termination reached (maxBGlobal = " << maxBGlobal/1e-9 << " nT), time step = " << sim.timestep << ", time = " << sim.t << " s, exiting." << endl << write; }
      }
   }

   // zero particle population counters
   for(size_t s=0;s<particleLists.size();++s) {
      Hybrid::logCounterParticleEscape[s] = 0.0;
      Hybrid::logCounterParticleImpact[s] = 0.0;
      Hybrid::logCounterParticleInject[s] = 0.0;
      Hybrid::logCounterParticleInjectMacroparticles[s] = 0.0;
      Hybrid::logCounterParticleEscapeKineticEnergy[s] = 0.0;
      Hybrid::logCounterParticleImpactKineticEnergy[s] = 0.0;
      Hybrid::logCounterParticleInjectKineticEnergy[s] = 0.0;
      Hybrid::logCounterParticleMaxVi[s] = 0.0;
   }
   
   // zero field counters
   Hybrid::logCounterFieldMaxCellUe = 0.0;
   Hybrid::logCounterFieldMaxNodeUe = 0.0;
   Hybrid::logCounterFieldMaxVw = 0.0;
   Hybrid::logCounterFieldMaxE = 0.0;
   Hybrid::logCounterFieldMinCellRhoQi = 0.0;
   Hybrid::logCounterFieldMinNodeRhoQi = 0.0;

   // reset counter start time
   Hybrid::logCounterTimeStart = sim.t;
   
   /*
   // silo time series (curves.silo)
   Real cnt=0;
   for(pargrid::CellID b=0; b<simClasses.pargrid.getNumberOfLocalCells(); ++b){ cnt+=1.0; }
   
   map<string,string> attribs;
   attribs["name"] = "blah";
   attribs["xlabel"] = "time";
   attribs["ylabel"] = "y laabeli";
   attribs["xunit"] = "s";
   attribs["yunit"] = "y unitti";   
   
   simClasses.vlsv.writeWithReduction("TIMESERIES",attribs,1,&cnt,MPI_SUM);
   */
   
   profile::stop();
   return success;
}


