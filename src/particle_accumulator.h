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

#ifndef PARTICLE_ACCUMULATOR_H
#define PARTICLE_ACCUMULATOR_H

#include <cstdlib>
#include <simulation.h>
#include <simulationclasses.h>
#include <base_class_particle_accumulator.h>

#include "particle_definition.h"
#include "particle_species.h"

class Accumulator: public ParticleAccumulatorBase {
 public:
   Accumulator();
   ~Accumulator();
   
   bool accumulateBoundaryCells(pargrid::DataID particleDataID,const unsigned int* N_particles);
   bool accumulateInnerCells(pargrid::DataID particleDataID,const unsigned int* N_particles);
   bool addConfigFileItems(ConfigReader& cr,const std::string& regionName);
   bool addRemoteUpdates();
   bool clearAccumulationArrays();
   bool finalize();
   bool initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,
		   const std::string& regionName,const ParticleListBase* plist);
   bool sendUpdates();
   bool wait();
   
 private:
   static int N_accumulators;        /**< Total number of allocated Accumulators, used together with variable 
				      * accumulatorCounter to determine when MPI transfers should be started.*/
   int myOrderNumber;                /**< Order number of this Accumulator.*/
   const Species* species;
   
   #if PROFILE_LEVEL > 0
      int arrayClearing;
      int dataCopying;
      int particleAccumulation;
   #endif

#ifdef WRITE_POPULATION_AVERAGES
   void accumulateCell(const Species& species,pargrid::CellID blockID,unsigned int N_particles,
		       const Particle<Real>* particles,Real* cellRhoQi,Real* cellJi,
		       Real* nAve,Real* vAve);
#else
   void accumulateCell(const Species& species,pargrid::CellID blockID,unsigned int N_particles,
		       const Particle<Real>* particles,Real* cellRhoQi,Real* cellJi);
#endif
};

inline ParticleAccumulatorBase* AccumulatorMaker() {return new Accumulator();}

#endif
