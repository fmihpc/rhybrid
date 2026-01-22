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

#ifndef PARTICLE_INJECTOR_H
#define PARTICLE_INJECTOR_H

#include <cstdlib>
#include <climits>

#include <simulation.h>
#include <simulationclasses.h>
#include <configreader.h>
#include <pargrid_userdata_dynamic.h>
#include <base_class_particle_injector.h>
#include "particle_definition.h"
#include "particle_species.h"

struct InjectorParameters {
   std::string name,type;
   int popid;
   Real m,q,w,T,vth,n,U,velocity[3];
};

bool getInjectorParameters(ParticleInjectorBase* injBasePtr,InjectorParameters& p);

class InjectorUniform: public ParticleInjectorBase {
 public:
   InjectorUniform();
   ~InjectorUniform();
   bool addConfigFileItems(ConfigReader& cr,const std::string& regionName);
   bool finalize();
   bool initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,const std::string& regionName,const ParticleListBase* plist);
   bool inject(pargrid::DataID speciesDataID,unsigned int* N_particles);
   void getParams(InjectorParameters& p);
 private:
   bool initialized;
   Real N_macroParticlesPerCell;
   const Species* species;
   Real velocity[3];
   Real U,T,vth,n,w,xmin,xmax,ymin,ymax,zmin,zmax;
   bool injectParticles(pargrid::CellID blockID,const Species& species,unsigned int* N_particles,pargrid::DataWrapper<Particle<Real> >& wrapper);
};

class InjectorAmbient: public ParticleInjectorBase {
 public:
   InjectorAmbient();
   ~InjectorAmbient();
   bool addConfigFileItems(ConfigReader& cr,const std::string& regionName);
   bool finalize();
   bool initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,const std::string& regionName,const ParticleListBase* plist);
   bool inject(pargrid::DataID speciesDataID,unsigned int* N_particles);
   void getParams(InjectorParameters& p);
 private:
   bool initialized;
   Real N_macroParticlesPerCellPerDt;
   Real N_macroParticlesPerCell;
   const Species* species;
   Real T,vth,n,w;
   bool injectParticles(pargrid::CellID blockID,const Species& species,unsigned int* N_particles,pargrid::DataWrapper<Particle<Real> >& wrapper,unsigned int wall);
};

class InjectorSolarWind: public ParticleInjectorBase {
 public:
   InjectorSolarWind();
   ~InjectorSolarWind();
   bool addConfigFileItems(ConfigReader& cr,const std::string& regionName);
   bool finalize();
   bool initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,const std::string& regionName,const ParticleListBase* plist);
   bool inject(pargrid::DataID speciesDataID,unsigned int* N_particles);
   void getParams(InjectorParameters& p);
 private:
   bool initialized;
   Real N_macroParticlesPerCellPerDt;
   Real N_macroParticlesPerCell;
   const Species* species;
   Real velocity[3];
   Real U,T,vth,n,w;
   bool (InjectorSolarWind::*checkIfInjectionCellFuncPtr)(const pargrid::CellID b); // pointer to function that checks if particles should be injected in this cell
   void (InjectorSolarWind::*initParticleCrdVelFuncPtr)(Real blockSize[3],Real& x,Real& y,Real& z,Real& vx,Real& vy,Real& vz); // pointer to function that initialize coordinates and velocity of a new particle
   bool checkIfInjectionCellDefault(const pargrid::CellID b);
   bool checkIfInjectionCellXPos(const pargrid::CellID b);
   bool checkIfInjectionCellXNeg(const pargrid::CellID b);
   bool checkIfInjectionCellYPos(const pargrid::CellID b);
   bool checkIfInjectionCellYNeg(const pargrid::CellID b);
   bool checkIfInjectionCellZPos(const pargrid::CellID b);
   bool checkIfInjectionCellZNeg(const pargrid::CellID b);
   void initParticleCrdVelDefault(Real blockSize[3],Real& x,Real& y,Real& z,Real& vx,Real& vy,Real& vz);
   void initParticleCrdVelXPos(Real blockSize[3],Real& x,Real& y,Real& z,Real& vx,Real& vy,Real& vz);
   void initParticleCrdVelXNeg(Real blockSize[3],Real& x,Real& y,Real& z,Real& vx,Real& vy,Real& vz);
   void initParticleCrdVelYPos(Real blockSize[3],Real& x,Real& y,Real& z,Real& vx,Real& vy,Real& vz);
   void initParticleCrdVelYNeg(Real blockSize[3],Real& x,Real& y,Real& z,Real& vx,Real& vy,Real& vz);
   void initParticleCrdVelZPos(Real blockSize[3],Real& x,Real& y,Real& z,Real& vx,Real& vy,Real& vz);
   void initParticleCrdVelZNeg(Real blockSize[3],Real& x,Real& y,Real& z,Real& vx,Real& vy,Real& vz);
   bool injectParticles(pargrid::CellID blockID,const Species& species,unsigned int* N_particles,pargrid::DataWrapper<Particle<Real> >& wrapper);
};

class InjectorFlow: public ParticleInjectorBase {
 public:
   InjectorFlow();
   ~InjectorFlow();
   bool addConfigFileItems(ConfigReader& cr,const std::string& regionName);
   bool finalize();
   bool initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,const std::string& regionName,const ParticleListBase* plist);
   bool inject(pargrid::DataID speciesDataID,unsigned int* N_particles);
   void getParams(InjectorParameters& p);
 private:
   bool initialized;
   uint32_t bitMaskInflowBoundaries;
   Real N_macroParticlesPerCell;
   const Species* species;
   Real velocity[3];
   Real U,T,vth,n,w;
   bool injectParticles(pargrid::CellID blockID,const Species& species,unsigned int* N_particles,pargrid::DataWrapper<Particle<Real> >& wrapper);
};

class InjectorIonosphere: public ParticleInjectorBase {
 public:
   InjectorIonosphere();
   ~InjectorIonosphere();
   bool addConfigFileItems(ConfigReader& cr,const std::string& regionName);
   bool finalize();
   bool initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,const std::string& regionName,const ParticleListBase* plist);
   bool inject(pargrid::DataID speciesDataID,unsigned int* N_particles);
   void getParams(InjectorParameters& p);
 private:
   bool initialized;
   const Species* species;
   unsigned int N_ionoPop;
   Real N_macroParticlesPerCell,N_macroParticlesPerDt,T,vth,w,R;
   bool injectParticles(pargrid::CellID blockID,const Species& species,unsigned int* N_particles,pargrid::DataWrapper<Particle<Real> >& wrapper);
};

class InjectorChapmanIonosphere: public ParticleInjectorBase {
 public:
   InjectorChapmanIonosphere();
   ~InjectorChapmanIonosphere();
   bool addConfigFileItems(ConfigReader& cr,const std::string& regionName);
   bool finalize();
   bool initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,const std::string& regionName,const ParticleListBase* plist);
   bool inject(pargrid::DataID speciesDataID,unsigned int* N_particles);
   void getParams(InjectorParameters& p);
 private:
   bool initialized;
   const Species* species;
   unsigned int N_ionoPop;
   Real N_macroParticlesPerCell,N_macroParticlesPerDt,vth,w,R,T,noonFactor, nightFactor;
   bool injectParticles(pargrid::CellID blockID,const Species& species,unsigned int* N_particles,pargrid::DataWrapper<Particle<Real> >& wrapper);
};

class InjectorExosphere: public ParticleInjectorBase {
 public:
   InjectorExosphere();
   ~InjectorExosphere();
   bool addConfigFileItems(ConfigReader& cr,const std::string& regionName);
   bool finalize();
   bool initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,const std::string& regionName,const ParticleListBase* plist);
   bool inject(pargrid::DataID speciesDataID,unsigned int* N_particles);
   void getParams(InjectorParameters& p);
 private:
   bool initialized;
   const Species* species;
   unsigned int N_exoPop;
   std::string neutralProfileName;
   Real N_macroParticlesPerCell,N_macroParticlesPerDt,T,vth,w,r0,R_exobase,R_shadow;
   std::vector<Real> n0,H0,T0,k0;
   bool injectParticles(pargrid::CellID blockID,const Species& species,unsigned int* N_particles,pargrid::DataWrapper<Particle<Real> >& wrapper);
};

inline ParticleInjectorBase* UniformIonCreator() {return new InjectorUniform();}
inline ParticleInjectorBase* AmbientIonCreator() {return new InjectorAmbient();}
inline ParticleInjectorBase* SolarWindIonCreator() {return new InjectorSolarWind();}
inline ParticleInjectorBase* FlowIonCreator() {return new InjectorFlow();}
inline ParticleInjectorBase* IonosphereIonCreator() {return new InjectorIonosphere();}
inline ParticleInjectorBase* ChapmanIonosphereIonCreator() {return new InjectorChapmanIonosphere();}
inline ParticleInjectorBase* ExosphereIonCreator() {return new InjectorExosphere();}

#endif
