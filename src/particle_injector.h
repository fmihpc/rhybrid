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

#ifndef DEFAULT_INJECTOR_H
#define DEFAULT_INJECTOR_H

#include <cstdlib>
#include <climits>

#include <simulation.h>
#include <simulationclasses.h>
#include <configreader.h>
#include <pargrid_userdata_dynamic.h>
#include <base_class_particle_injector.h>
#include "particle_definition.h"
#include "particle_species.h"

class InjectorUniform: public ParticleInjectorBase {
 public:
   InjectorUniform();
   ~InjectorUniform();
   bool addConfigFileItems(ConfigReader& cr,const std::string& regionName);
   bool finalize();
   bool initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,
		   const std::string& regionName,const ParticleListBase* plist);
   bool inject(pargrid::DataID speciesDataID,unsigned int* N_particles);
 private:
   bool initialized;
   Real N_macroParticlesPerCell;
   const Species* species;
   Real U,vth,n,w;
   bool injectParticles(pargrid::CellID blockID,const Species& species,unsigned int* N_particles,
			pargrid::DataWrapper<Particle<Real> >& wrapper);
};

class InjectorSolarWind: public ParticleInjectorBase {
 public:
   InjectorSolarWind();
   ~InjectorSolarWind();
      
   bool addConfigFileItems(ConfigReader& cr,const std::string& regionName);
   bool finalize();
   bool initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,
		   const std::string& regionName,const ParticleListBase* plist);
   bool inject(pargrid::DataID speciesDataID,unsigned int* N_particles);

 private:
   bool initialized;
   Real N_macroParticlesPerCellPerDt;
   Real N_macroParticlesPerCell;
   const Species* species;
   Real U,vth,n,w;
   bool injectParticles(pargrid::CellID blockID,const Species& species,unsigned int* N_particles,
			pargrid::DataWrapper<Particle<Real> >& wrapper);
};

class InjectorIonosphere: public ParticleInjectorBase {
 public:
   InjectorIonosphere();
   ~InjectorIonosphere();
      
   bool addConfigFileItems(ConfigReader& cr,const std::string& regionName);
   bool finalize();
   bool initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,
		   const std::string& regionName,const ParticleListBase* plist); 
   bool inject(pargrid::DataID speciesDataID,unsigned int* N_particles);
  
 private:
   bool initialized;
   const Species* species;
   unsigned int N_ionoPop;
   Real N_macroParticlesPerCell,N_macroParticlesPerDt,vth,w,R;
   bool injectParticles(pargrid::CellID blockID,const Species& species,unsigned int* N_particles,
			pargrid::DataWrapper<Particle<Real> >& wrapper);
};

class InjectorChapmanIonosphere: public ParticleInjectorBase {
 public:
   InjectorChapmanIonosphere();
   ~InjectorChapmanIonosphere();
      
   bool addConfigFileItems(ConfigReader& cr,const std::string& regionName);
   bool finalize();
   bool initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,
		   const std::string& regionName,const ParticleListBase* plist); 
   bool inject(pargrid::DataID speciesDataID,unsigned int* N_particles);
  
 private:
   bool initialized;
   const Species* species;
   unsigned int N_ionoPop;
   Real N_macroParticlesPerCell,N_macroParticlesPerDt,vth,w,R,T,noonFactor, nightFactor;
   bool injectParticles(pargrid::CellID blockID,const Species& species,unsigned int* N_particles,
			pargrid::DataWrapper<Particle<Real> >& wrapper);
};

class InjectorExosphere: public ParticleInjectorBase {
 public:
   InjectorExosphere();
   ~InjectorExosphere();
      
   bool addConfigFileItems(ConfigReader& cr,const std::string& regionName);
   bool finalize();
   bool initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,
		   const std::string& regionName,const ParticleListBase* plist); 
   bool inject(pargrid::DataID speciesDataID,unsigned int* N_particles);
  
 private:
   bool initialized;
   const Species* species;
   unsigned int N_exoPop;
   std::string neutralProfileName;
   Real N_macroParticlesPerCell,N_macroParticlesPerDt,vth,w,r0,R_exobase,R_shadow;
   std::vector<Real> n0,H0,T0,k0;
   bool injectParticles(pargrid::CellID blockID,const Species& species,unsigned int* N_particles,
			pargrid::DataWrapper<Particle<Real> >& wrapper);
};

inline ParticleInjectorBase* UniformIonCreator() {return new InjectorUniform();}
inline ParticleInjectorBase* SolarWindIonCreator() {return new InjectorSolarWind();}
inline ParticleInjectorBase* IonosphereIonCreator() {return new InjectorIonosphere();}
inline ParticleInjectorBase* ChapmanIonosphereIonCreator() {return new InjectorChapmanIonosphere();}
inline ParticleInjectorBase* ExosphereIonCreator() {return new InjectorExosphere();}


#endif
