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

#include <main.h>
#include <gridbuilder.h>
#include <logically_cartesian_builder.h>
#include <particle_accumulator.h>
#include "particle_boundary_cond_hybrid.h"
#include <particle_injector.h>
#include <particle_propagator_boris_buneman.h>
#include <object_factories.h>

#include <dataoperatorcontainer.h>
#include <operator_cellid.h>
#include <operator_mpirank.h>
#include <operator_load.h>
#include <operator_particle.h>
#include "particle_definition.h"
#include "particle_species.h"
#include "operator_userdata.h"

using namespace std;

typedef Particle<Real> PARTICLE;
typedef Species SPECIES;

bool registerObjectMakers(ObjectFactories& objectFactories) {
   bool success = true;
   if (objectFactories.gridBuilders.registerMaker("LogicallyCartesian",LCCreator) == false) { success = false; }

   if (objectFactories.particleAccumulators.registerMaker("HybridAccumulator",AccumulatorMaker) == false) { success = false; }
   if (objectFactories.particleBoundaryConditions.registerMaker("HybridBoundaryCond",HybridBoundaryCondMaker<SPECIES,PARTICLE>) == false) { success = false; }
   if (objectFactories.particleInjectors.registerMaker("UniformInjector",UniformIonCreator) == false) { success = false; }
   if (objectFactories.particleInjectors.registerMaker("SolarWindInjector",SolarWindIonCreator) == false) { success = false; }
   if (objectFactories.particleInjectors.registerMaker("IonosphereInjector",IonosphereIonCreator) == false) { success = false; }
   if (objectFactories.particleInjectors.registerMaker("ChapmanIonosphereInjector",ChapmanIonosphereIonCreator) == false) { success = false; }
   if (objectFactories.particleInjectors.registerMaker("ExosphereInjector",ExosphereIonCreator) == false) { success = false; }
   if (objectFactories.particlePropagators.registerMaker("BorisBuneman",BBMaker<PARTICLE>) == false) {success = false;}
   
   return success;
}

bool registerDataOperators() {
   bool success = true;
   DataOperatorContainer& doc = corsair::getObjectWrapper().dataOperatorContainer;
   //if(doc.registerOperator(new MPIRank) == false)          { success = false; }
   if(doc.registerOperator(new CellIDOP) == false)         { success = false; }
   //if(doc.registerOperator(new LoadOP) == false)           { success = false; }
   if(doc.registerOperator(new UserDataOP) == false)       { success = false; }
   if(doc.registerOperator(new ParticleOperator) == false) { success = false; }
   return success;
}
