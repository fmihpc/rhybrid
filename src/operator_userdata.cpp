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
void logMacroparticles(Simulation& sim,SimulationClasses& simClasses,const std::vector<ParticleListBase*>& particleLists) {
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
   writeCellDataVariable(spatMeshName,Hybrid::dataCellRhoQiBgID,          "cellRhoQiBg",         N_blocks,1);
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
   writeCellDataVariable(spatMeshName,Hybrid::dataCounterCellMaxUeID,     "counterCellMaxUe",    N_blocks,1);
   writeCellDataVariable(spatMeshName,Hybrid::dataCounterCellMaxViID,     "counterCellMaxVi",    N_blocks,1);
   writeCellDataVariable(spatMeshName,Hybrid::dataCounterCellMinRhoQiID,  "counterCellMinRhoQi", N_blocks,1);
#ifdef USE_ECUT
   writeCellDataVariable(spatMeshName,Hybrid::dataCounterNodeEcutID,      "counterNodeEcut",     N_blocks,1);
#endif
#ifdef USE_MAXVW
   writeCellDataVariable(spatMeshName,Hybrid::dataCounterNodeMaxVwID,     "counterNodeMaxVw",    N_blocks,1);
#endif
   writeCellDataVariableBool(spatMeshName,Hybrid::dataInnerFlagFieldID,   "innerFlagField",      N_blocks,1);
   writeCellDataVariableBool(spatMeshName,Hybrid::dataInnerFlagNodeID,    "innerFlagNode",       N_blocks,1);
   writeCellDataVariableBool(spatMeshName,Hybrid::dataInnerFlagParticleID,"innerFlagParticle",   N_blocks,1);
   writeCellDataVariableBool(spatMeshName,Hybrid::dataInnerFlagCellEpID,  "innerFlagCellEp",     N_blocks,1);
   writeCellDataVariableBool(spatMeshName,Hybrid::dataOuterBoundaryFlagID,"outerBoundaryFlag",   N_blocks,1);
#ifdef USE_XMIN_BOUNDARY
   writeCellDataVariableBool(spatMeshName,Hybrid::dataXminFlagID,         "xMinFlag",            N_blocks,1);
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
#ifdef WRITE_POPULATION_AVERAGES
   Real invAveCnt = 0.0;
   if(Hybrid::averageCounter > 0) {
      invAveCnt = 1.0/static_cast<Real>(Hybrid::averageCounter);
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
   Hybrid::averageCounter = 0;
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
   logMacroparticles(*sim,*simClasses,particleLists);
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
void calcParticleLog(Simulation& sim,SimulationClasses& simClasses,vector<ParticleLogData>& plogData,const std::vector<ParticleListBase*>& particleLists)
{
   plogData.clear();
   for(size_t s=0;s<particleLists.size();++s) {
      plogData.push_back(ParticleLogData());
      plogData[s].N_macroParticles = 0.0;
      plogData[s].N_realParticles = 0.0;
      plogData[s].sumVx = 0.0;
      plogData[s].sumVy = 0.0;
      plogData[s].sumVz = 0.0;
      plogData[s].sumV = 0.0;
      plogData[s].sumWV2 = 0.0;
   }
#ifdef USE_XMIN_BOUNDARY
   bool* xMinFlag = simClasses.pargrid.getUserDataStatic<bool>(Hybrid::dataXminFlagID);
#endif
   for(pargrid::CellID b=0; b<simClasses.pargrid.getNumberOfLocalCells(); ++b) {
#ifdef USE_XMIN_BOUNDARY
      // do not include particles if this block has a cell with xMinFlag == true
      bool skipBlock = false;
      for(int k=0; k<block::WIDTH_Z; ++k) for(int j=0; j<block::WIDTH_Y; ++j) for(int i=0; i<block::WIDTH_X; ++i) {
	 const int n = (b*block::SIZE+block::index(i,j,k));
         if(xMinFlag[n] == true) { skipBlock = true; }
      }
      if(skipBlock == true) { continue; }
#endif
      for(size_t s=0;s<particleLists.size();++s) {
	 pargrid::DataID speciesDataID = pargrid::INVALID_DATAID;
	 if(particleLists[s]->getParticles(speciesDataID) == false) { continue; }
	 pargrid::DataWrapper<Particle<Real> > wrapper = simClasses.pargrid.getUserDataDynamic<Particle<Real> >(speciesDataID);
	 Particle<Real>** particleList = wrapper.data();
	 Particle<Real>* particles = particleList[b];
	 pargrid::ArraySizetype N_particles = wrapper.size(b);
	 plogData[s].N_macroParticles += N_particles;
	 for(size_t p=0; p<N_particles; ++p) {
	    plogData[s].N_realParticles += particles[p].state[particle::WEIGHT];
	    plogData[s].sumVx += particles[p].state[particle::WEIGHT]*particles[p].state[particle::VX];
	    plogData[s].sumVy += particles[p].state[particle::WEIGHT]*particles[p].state[particle::VY];
	    plogData[s].sumVz += particles[p].state[particle::WEIGHT]*particles[p].state[particle::VZ];
	    plogData[s].sumV += particles[p].state[particle::WEIGHT]*sqrt(sqr(particles[p].state[particle::VX]) +
									  sqr(particles[p].state[particle::VY]) +
									  sqr(particles[p].state[particle::VZ]));
	    plogData[s].sumWV2 += particles[p].state[particle::WEIGHT]*(sqr(particles[p].state[particle::VX]) +
									sqr(particles[p].state[particle::VY]) +
									sqr(particles[p].state[particle::VZ]));
	 }
      }
   }
}

// calculate field log quantities
void calcFieldLog(Simulation& sim,SimulationClasses& simClasses,FieldLogData& flogData)
{
   // synchronize faceB before using it below
   simClasses.pargrid.startNeighbourExchange(pargrid::DEFAULT_STENCIL,Hybrid::dataFaceBID);
   //profile::start("MPI waits",mpiWaitID);
   simClasses.pargrid.wait(pargrid::DEFAULT_STENCIL,Hybrid::dataFaceBID);
   //profile::stop();
   flogData.N_cells = 0.0;
   flogData.sumBx = 0.0;
   flogData.sumBy = 0.0;
   flogData.sumBz = 0.0;
   flogData.sumB = 0.0;
   flogData.maxB = 0.0;
   flogData.sumDivB = 0.0;
   flogData.maxDivB = 0.0;
   flogData.maxDivBPerB = 0.0;
   flogData.sumB2 = 0.0;
   Real* faceB = simClasses.pargrid.getUserDataStatic<Real>(Hybrid::dataFaceBID);
   bool* innerFlag = simClasses.pargrid.getUserDataStatic<bool>(Hybrid::dataInnerFlagFieldID);
#ifdef USE_XMIN_BOUNDARY
   bool* xMinFlag = simClasses.pargrid.getUserDataStatic<bool>(Hybrid::dataXminFlagID);
#endif
   for(pargrid::CellID b=0; b<simClasses.pargrid.getNumberOfLocalCells(); ++b) {
      if(simClasses.pargrid.getNeighbourFlags(b) != pargrid::ALL_NEIGHBOURS_EXIST) { continue; }
      const unsigned int size = (block::WIDTH_X+2)*(block::WIDTH_Y+2)*(block::WIDTH_Z+2);
      Real array[size*3];
      fetchData(faceB,array,simClasses,b,3);
      for(int k=0; k<block::WIDTH_Z; ++k) for(int j=0; j<block::WIDTH_Y; ++j) for(int i=0; i<block::WIDTH_X; ++i) {
	 const int n = (b*block::SIZE+block::index(i,j,k));
	 // do not include cells inside the inner boundary
	 if(Hybrid::includeInnerCellsInFieldLog == false) {
	    if(innerFlag[n] == true) { continue; }
	 }
#ifdef USE_XMIN_BOUNDARY
         // do not include cells at x < xmin
         if(xMinFlag[n] == true) { continue; }
#endif
	 Real divB = ((array[(block::arrayIndex(i+1,j+1,k+1))*3+0] - array[(block::arrayIndex(i+0,j+1,k+1))*3+0])
		    + (array[(block::arrayIndex(i+1,j+1,k+1))*3+1] - array[(block::arrayIndex(i+1,j+0,k+1))*3+1])
		    + (array[(block::arrayIndex(i+1,j+1,k+1))*3+2] - array[(block::arrayIndex(i+1,j+1,k+0))*3+2]))/Hybrid::dx;
	 divB = fabs(divB);
	 Real Bx = 0.5*(array[(block::arrayIndex(i+1,j+1,k+1))*3+0] + array[(block::arrayIndex(i+0,j+1,k+1))*3+0]);
	 Real By = 0.5*(array[(block::arrayIndex(i+1,j+1,k+1))*3+1] + array[(block::arrayIndex(i+1,j+0,k+1))*3+1]);
	 Real Bz = 0.5*(array[(block::arrayIndex(i+1,j+1,k+1))*3+2] + array[(block::arrayIndex(i+1,j+1,k+0))*3+2]);
	 Real B2 = sqr(Bx) + sqr(By) + sqr(Bz);
	 Real Btot = sqrt(B2);
	 Real divBPerB = 0.0;
	 if(Btot > 0.0) { divBPerB = divB/Btot; }
	 flogData.sumBx += Bx;
	 flogData.sumBy += By;
	 flogData.sumBz += Bz;
	 flogData.sumB += Btot;
	 flogData.sumDivB += divB;
	 flogData.sumB2 += B2;
	 if(Btot > flogData.maxB) { flogData.maxB = Btot; }
	 if(fabs(divB) > fabs(flogData.maxDivB)) { flogData.maxDivB = divB; }
	 if(fabs(divBPerB) > fabs(flogData.maxDivBPerB)) { flogData.maxDivBPerB = divBPerB; }
	 flogData.N_cells++;
      }
   }
}

bool writeLogs(Simulation& sim,SimulationClasses& simClasses,const std::vector<ParticleListBase*>& particleLists) {
   bool success = true;
   static int profWriteLogsID = -1;
   //if(getInitialized() == false) { return false; }
   profile::start("writeLogs",profWriteLogsID);
   vector<ParticleLogData> plogData;
   FieldLogData flogData;
   calcParticleLog(sim,simClasses,plogData,particleLists);
   calcFieldLog(sim,simClasses,flogData);
   if(sim.mpiRank==sim.MASTER_RANK) {
      for(size_t i=0;i<Hybrid::plog.size();++i) {
	 (*Hybrid::plog[i]) << sim.t << " ";
      }
      Hybrid::flog << sim.t << " ";
   }
   // go thru populations
   const Real Dt = sim.t - Hybrid::particleCounterTimeStart;
   for(size_t s=0;s<particleLists.size();++s) {
      Real N_macroParticlesThisProcess = plogData[s].N_macroParticles;
      Real N_macroParticlesGlobal = 0.0;
      Real N_realParticlesThisProcess = plogData[s].N_realParticles;
      Real N_realParticlesGlobal = 0.0;
      Real sumVxThisProcess = plogData[s].sumVx;
      Real sumVxGlobal = 0.0;
      Real sumVyThisProcess = plogData[s].sumVy;
      Real sumVyGlobal = 0.0;
      Real sumVzThisProcess = plogData[s].sumVz;
      Real sumVzGlobal = 0.0;
      Real sumVThisProcess = plogData[s].sumV;
      Real sumVGlobal = 0.0;
      Real sumWV2ThisProcess = plogData[s].sumWV2;
      Real sumWV2Global = 0.0;
      Real sumEscapeThisProcess = Hybrid::particleCounterEscape[s];
      Real sumEscapeGlobal = 0.0;
      Real sumImpactThisProcess = Hybrid::particleCounterImpact[s];
      Real sumImpactGlobal = 0.0;
      Real sumInjectThisProcess = Hybrid::particleCounterInject[s];
      Real sumInjectGlobal = 0.0;
      Real sumInjectMacroThisProcess = Hybrid::particleCounterInjectMacroparticles[s];
      Real sumInjectMacroGlobal = 0.0;
      MPI_Reduce(&N_macroParticlesThisProcess,&N_macroParticlesGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
      MPI_Reduce(&N_realParticlesThisProcess,&N_realParticlesGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
      MPI_Reduce(&sumVxThisProcess,&sumVxGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
      MPI_Reduce(&sumVyThisProcess,&sumVyGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
      MPI_Reduce(&sumVzThisProcess,&sumVzGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
      MPI_Reduce(&sumVThisProcess,&sumVGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
      MPI_Reduce(&sumWV2ThisProcess,&sumWV2Global,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
      MPI_Reduce(&sumEscapeThisProcess,&sumEscapeGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
      MPI_Reduce(&sumImpactThisProcess,&sumImpactGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
      MPI_Reduce(&sumInjectThisProcess,&sumInjectGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
      MPI_Reduce(&sumInjectMacroThisProcess,&sumInjectMacroGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
      if(sim.mpiRank==sim.MASTER_RANK) {
	 (*Hybrid::plog[s]) << N_realParticlesGlobal << " " << N_macroParticlesGlobal << " ";
	 if(N_realParticlesGlobal > 0.0) {
	    (*Hybrid::plog[s])
	      << sumVxGlobal/N_realParticlesGlobal << " "
	      << sumVyGlobal/N_realParticlesGlobal << " "
	      << sumVzGlobal/N_realParticlesGlobal << " "
	      << sumVGlobal/N_realParticlesGlobal  << " ";
	 }
	 else {
	    (*Hybrid::plog[s]) << 0.0 << " " << 0.0 << " " << 0.0 << " " << 0.0 << " ";
	 }
	 const Species* species = reinterpret_cast<const Species*>(particleLists[s]->getSpecies());
	 (*Hybrid::plog[s]) << 0.5*species->m*sumWV2Global << " ";
	 if(Dt > 0) {
	    (*Hybrid::plog[s]) << sumEscapeGlobal/Dt << " " << sumImpactGlobal/Dt << " " << sumInjectGlobal/Dt << " " << sumInjectMacroGlobal/Dt*sim.dt << " ";
	 }
	 else {
	    (*Hybrid::plog[s]) << 0.0 << " " << 0.0 << " " << 0.0 << " " << 0.0 << " ";
	 }
      }
   }
   // field log
   Real N_cellsThisProcess = flogData.N_cells;
   Real N_cellsGlobal = 0.0;
   Real sumBxThisProcess = flogData.sumBx;
   Real sumBxGlobal = 0.0;
   Real sumByThisProcess = flogData.sumBy;
   Real sumByGlobal = 0.0;
   Real sumBzThisProcess = flogData.sumBz;
   Real sumBzGlobal = 0.0;
   Real sumBThisProcess = flogData.sumB;
   Real sumBGlobal = 0.0;
   Real maxBThisProcess = flogData.maxB;
   Real maxBGlobal = 0.0;
   Real sumDivBThisProcess = flogData.sumDivB;
   Real sumDivBGlobal = 0.0;
   Real maxDivBThisProcess = flogData.maxDivB;
   Real maxDivBGlobal = 0.0;
   Real maxDivPerBThisProcess = flogData.maxDivBPerB;
   Real maxDivPerBGlobal = 0.0;
   Real sumB2ThisProcess = flogData.sumB2;
   Real sumB2Global = 0.0;
   MPI_Reduce(&N_cellsThisProcess,&N_cellsGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
   MPI_Reduce(&sumBxThisProcess,&sumBxGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
   MPI_Reduce(&sumByThisProcess,&sumByGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
   MPI_Reduce(&sumBzThisProcess,&sumBzGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
   MPI_Reduce(&sumBThisProcess,&sumBGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
   MPI_Reduce(&maxBThisProcess,&maxBGlobal,1,MPI_Type<Real>(),MPI_MAX,sim.MASTER_RANK,sim.comm);
   MPI_Reduce(&sumDivBThisProcess,&sumDivBGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
   MPI_Reduce(&maxDivBThisProcess,&maxDivBGlobal,1,MPI_Type<Real>(),MPI_MAX,sim.MASTER_RANK,sim.comm);
   MPI_Reduce(&maxDivPerBThisProcess,&maxDivPerBGlobal,1,MPI_Type<Real>(),MPI_MAX,sim.MASTER_RANK,sim.comm);
   MPI_Reduce(&sumB2ThisProcess,&sumB2Global,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);

   if(sim.mpiRank==sim.MASTER_RANK) {
      // field
      if(N_cellsGlobal > 0) {
	 Hybrid::flog
	   << sumBxGlobal/N_cellsGlobal << " "
	   << sumByGlobal/N_cellsGlobal << " "
	   << sumBzGlobal/N_cellsGlobal << " "
	   << sumBGlobal/N_cellsGlobal << " ";
      }
      else {
	 Hybrid::flog << 0.0 << " " << 0.0 << " " << 0.0 << " " << 0.0 << " ";
      }
      Hybrid::flog << maxBGlobal << " ";
      if(N_cellsGlobal > 0) { Hybrid::flog << sumDivBGlobal/N_cellsGlobal << " "; }
      else { Hybrid::flog << 0.0 << " "; }
      Hybrid::flog << maxDivBGlobal << " " << Hybrid::dx*maxDivPerBGlobal << " " << sumB2Global*Hybrid::dV/(2.0*constants::PERMEABILITY) << " ";
      // endlines
      for(size_t i=0;i<Hybrid::plog.size();++i) {
	 (*Hybrid::plog[i]) << endl;
      }
      Hybrid::flog << endl;
      if(maxBGlobal > Hybrid::terminateLimitMaxB) {
         success = false;
         simClasses.logger << "(RHYBRID) CONSTRAINT: maximum |B| for run termination reached (maxBGlobal = " << maxBGlobal/1e-9 << " nT), exiting." << endl;
         simClasses.logger << write;
      }
   }

   // zero particle counters
   Hybrid::particleCounterTimeStart = sim.t;
   for(size_t s=0;s<particleLists.size();++s) {
      Hybrid::particleCounterEscape[s] = 0.0;
      Hybrid::particleCounterImpact[s] = 0.0;
      Hybrid::particleCounterInject[s] = 0.0;
      Hybrid::particleCounterInjectMacroparticles[s] = 0.0;
   }
   
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


