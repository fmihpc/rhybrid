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

#ifndef OPERATOR_USER_DATA_H
#define OPERATOR_USER_DATA_H

#include <dataoperator.h>
#include "particle_species.h"

// struct to store some particle population log vales
struct LogDataParticle {
   Real N_macroParticles,N_realParticles,sumVx,sumVy,sumVz,sumV,sumWV2;
};

// struct to store some field log vales
struct LogDataField {
   Real N_cells;
   Real sumBx,sumBy,sumBz,sumB,maxB,sumDivB,maxDivB,maxDivBPerB,sumB2;
   Real sumCellJix,sumCellJiy,sumCellJiz,sumCellJi,maxCellJi;
   Real sumCellEpx,sumCellEpy,sumCellEpz,sumCellEp,maxCellEp;
   Real sumNodeEx,sumNodeEy,sumNodeEz,sumNodeE,maxNodeE;
};

class UserDataOP: public DataOperator {
 public:
   UserDataOP();
   ~UserDataOP();
   bool finalize();
   std::string getName() const;
   bool initialize(ConfigReader& cr,Simulation& sim,SimulationClasses& simClasses);
   virtual bool writeData(const std::string& spatMeshName,const std::vector<ParticleListBase*>& particles);
 private:
   int profileID;
   bool writeCellDataVariable(const std::string& spatMeshName,pargrid::DataID& dataVarID,const std::string& dataName,const pargrid::CellID& N_blocks,const uint64_t& vectorDim);
   bool writeCellDataVariableBool(const std::string& spatMeshName,pargrid::DataID& dataVarID,const std::string& dataName,const pargrid::CellID& N_blocks,const uint64_t& vectorDim);
   void calcCellDiv(Real* faceData,std::vector<Real>& cellDiv);
   void calcCellNPles(std::vector<Real>& cellNPles,const std::vector<ParticleListBase*>& particleLists);
   void calcCellParticleBulkParameters(std::vector<Real>& cellDensity,std::vector<Real>& cellTemperature,std::vector<Real>& cellVelocity,const std::vector<ParticleListBase*>& particleLists,std::vector<unsigned int> s);
};

void logCalcParticle(Simulation& sim,SimulationClasses& simClasses,std::vector<LogDataParticle>& logDataParticle,const std::vector<ParticleListBase*>& particleLists);
void logCalcField(Simulation& sim,SimulationClasses& simClasses,LogDataField& logDataField);
bool logWriteParticleField(Simulation& sim,SimulationClasses& simClasses,const std::vector<ParticleListBase*>& particles);

#endif
