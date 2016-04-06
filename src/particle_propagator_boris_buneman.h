/** This file is part of the RHybrid simulation.
 *
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

#ifndef PARTICLE_PROPAGATOR_BORIS_BUNEMAN_H
#define PARTICLE_PROPAGATOR_BORIS_BUNEMAN_H

#include <configreader.h>
#include <simulation.h>
#include <simulationclasses.h>
#include <base_class_particle_propagator.h>
#include <linear_algebra.h>

#include "hybrid.h"
#include "hybrid_propagator.h"
#include "particle_species.h"
#ifdef USE_B_CONSTANT
#include "magnetic_field.h"
#endif

template<class PARTICLE>
class BorisBuneman: public ParticlePropagatorBase {
 public:
   BorisBuneman();
   
   bool addConfigFileItems(ConfigReader& cr,const std::string& configName);
   bool finalize();
   bool propagateCell(pargrid::CellID block,pargrid::DataID particleDataID,const double* const coordinates,
		      unsigned int N_particles);
   bool initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,
		   const std::string& configName,const ParticleListBase* plist);
   
 private:
   const Species* species;
   
   void propagate(const Real xBlock,const Real yBlock,const Real zBlock,pargrid::CellID blockID,const Species& species,PARTICLE& particle,pargrid::CellID globalID);
};

template<class PARTICLE>
ParticlePropagatorBase* BBMaker() {return new BorisBuneman<PARTICLE>();}

template<class PARTICLE>
BorisBuneman<PARTICLE>::BorisBuneman() { }

template<class PARTICLE>
bool BorisBuneman<PARTICLE>::addConfigFileItems(ConfigReader& cr,const std::string& configName) {
   return true;
}

template<class PARTICLE>
bool BorisBuneman<PARTICLE>::finalize() {return true;}

template<class PARTICLE>
  void BorisBuneman<PARTICLE>::propagate(const Real xBlock,const Real yBlock,const Real zBlock,pargrid::CellID blockID,const Species& species,PARTICLE& particle,pargrid::CellID globalID) {
   bool accelerate = species.accelerate;
   if(simClasses->pargrid.getNeighbourFlags(blockID) != pargrid::ALL_NEIGHBOURS_EXIST) {
      accelerate = false;
   }
   /*if(accelerate == true) {
      Real dU[3] = { particle.state[particle::VX]-Ue[0], particle.state[particle::VY]-Ue[1], particle.state[particle::VZ]-Ue[2] };
      const Real half_alpha = 0.5*species.q*sim->dt/species.m;
      Real b[3] = {half_alpha*B[0], half_alpha*B[1], half_alpha*B[2]};
      const Real b2 = sqr(b[0]) + sqr(b[1]) + sqr(b[2]);
      const Real beta = 2.0/(1.0 + b2);
      Real dUxb[3],dUxbxb[3];
      cross(dU,b,dUxb);
      cross(dUxb,b,dUxbxb);
      particle.state[particle::VX] += beta*( dUxb[0] + dUxbxb[0] );
      particle.state[particle::VY] += beta*( dUxb[1] + dUxbxb[1] );
      particle.state[particle::VZ] += beta*( dUxb[2] + dUxbxb[2] );
      
      const Real v2 = sqr(particle.state[particle::VX]) + sqr(particle.state[particle::VY]) + sqr(particle.state[particle::VZ]);
      if(v2 > Hybrid::maxVi2) {
	 const Real norm = sqrt(Hybrid::maxVi2/v2);
	 particle.state[particle::VX] *= norm;
	 particle.state[particle::VY] *= norm;
	 particle.state[particle::VZ] *= norm;
	 Real* cellMaxVi = reinterpret_cast<Real*>(simClasses->pargrid.getUserData(Hybrid::dataCellMaxViID));
	 cellMaxVi[blockID]++;
      }
   }*/
   
   // accelerate particle
   if(accelerate == true) {
      Real E[3],B[3],Ue[3];
      Real r[3] = {particle.state[particle::X],particle.state[particle::Y],particle.state[particle::Z]};
      getFields(r,B,Ue,*sim,*simClasses,blockID);
#ifdef USE_B_CONSTANT
      addConstantB(r[0]+xBlock,r[1]+yBlock,r[2]+zBlock,B);
#endif
      crossProduct(B,Ue,E);
      Real tx,ty,tz,sx,sy,sz,dvx,dvy,dvz,vmx,vmy,vmz,v0x,v0y,v0z,vpx,vpy,vpz,qmideltT2,t2,b2;
      qmideltT2= 0.5*species.q*sim->dt/species.m;
      dvx=qmideltT2*E[0];
      dvy=qmideltT2*E[1];
      dvz=qmideltT2*E[2];
      tx=qmideltT2*B[0];
      ty=qmideltT2*B[1];
      tz=qmideltT2*B[2];
      t2=tx*tx+ty*ty+tz*tz;
      b2=2./(1.+t2);
      sx=b2*tx;
      sy=b2*ty;
      sz=b2*tz;
      vmx=particle.state[particle::VX]+dvx;
      vmy=particle.state[particle::VY]+dvy;
      vmz=particle.state[particle::VZ]+dvz;
      v0x=vmx+vmy*tz-vmz*ty;
      v0y=vmy+vmz*tx-vmx*tz;
      v0z=vmz+vmx*ty-vmy*tx;
      vpx=vmx+v0y*sz-v0z*sy;
      vpy=vmy+v0z*sx-v0x*sz;
      vpz=vmz+v0x*sy-v0y*sx;
      particle.state[particle::VX]=vpx+dvx;
      particle.state[particle::VY]=vpy+dvy;
      particle.state[particle::VZ]=vpz+dvz;
      const Real v2 = sqr(particle.state[particle::VX]) + sqr(particle.state[particle::VY]) + sqr(particle.state[particle::VZ]);
      if(v2 > Hybrid::maxVi2) {
	 const Real norm = sqrt(Hybrid::maxVi2/v2);
	 particle.state[particle::VX] *= norm;
	 particle.state[particle::VY] *= norm;
	 particle.state[particle::VZ] *= norm;
	 Real* cellMaxVi = reinterpret_cast<Real*>(simClasses->pargrid.getUserData(Hybrid::dataCellMaxViID));
	 cellMaxVi[blockID]++;
      }
   }

#ifdef ION_SPECTRA_ALONG_ORBIT
   bool* spectraFlag = reinterpret_cast<bool*>(simClasses->pargrid.getUserData(Hybrid::dataSpectraFlagID));
     if(spectraFlag[blockID] == true && Hybrid::recordSpectra == true) {
      //if(particle.state[particle::INI_TIME] >= 0.0) {
         Hybrid::spectraParticleOutput.push_back( static_cast<Real>(sim->t) );                               // 1
         Hybrid::spectraParticleOutput.push_back( static_cast<Real>(species.popid) );                        // 2
         Hybrid::spectraParticleOutput.push_back( particle.state[particle::WEIGHT] );                        // 3
         Hybrid::spectraParticleOutput.push_back( static_cast<Real>(globalID) );                             // 4
         Hybrid::spectraParticleOutput.push_back( particle.state[particle::VX] );                            // 5
         Hybrid::spectraParticleOutput.push_back( particle.state[particle::VY] );                            // 6
         Hybrid::spectraParticleOutput.push_back( particle.state[particle::VZ] );                            // 7
         Hybrid::spectraParticleOutput.push_back( particle.state[particle::INI_TIME] );                      // 8
         Hybrid::spectraParticleOutput.push_back( static_cast<Real>(particle.state[particle::INI_CELLID]) ); // 9
         Hybrid::spectraParticleOutput.push_back( particle.state[particle::INI_X] );                         // 10
         Hybrid::spectraParticleOutput.push_back( particle.state[particle::INI_Y] );                         // 11
         Hybrid::spectraParticleOutput.push_back( particle.state[particle::INI_Z] );                         // 12
         Hybrid::spectraParticleOutput.push_back( particle.state[particle::INI_VX] );                        // 13
         Hybrid::spectraParticleOutput.push_back( particle.state[particle::INI_VY] );                        // 14
         Hybrid::spectraParticleOutput.push_back( particle.state[particle::INI_VZ] );                        // 15 = SPECTRA_FILE_VARIABLES
	 //particle.state[particle::INI_TIME] = -100.0; // only detect each particle once
      //}
   }
#endif
   
   // move particle
   particle.state[particle::X] += sim->dt*particle.state[particle::VX]; 
   particle.state[particle::Y] += sim->dt*particle.state[particle::VY];
   particle.state[particle::Z] += sim->dt*particle.state[particle::VZ];
}

template<class PARTICLE>
bool BorisBuneman<PARTICLE>::propagateCell(pargrid::CellID blockID,pargrid::DataID particleDataID,const double* const coordinates,
					   unsigned int N_particles) {
   pargrid::DataWrapper<PARTICLE> wrapper = simClasses->pargrid.getUserDataDynamic<PARTICLE>(particleDataID);
   const Real* crd = getBlockCoordinateArray(*sim,*simClasses);
   const size_t b3 = 3*blockID;
   const Real xBlock = crd[b3+0];
   const Real yBlock = crd[b3+1];
   const Real zBlock = crd[b3+2];
   const pargrid::CellID globalID = simClasses->pargrid.getGlobalIDs()[blockID];
   for(size_t p=0;p<N_particles;++p) { propagate(xBlock,yBlock,zBlock,blockID,*species,wrapper.data()[blockID][p],globalID); }
   return true;
}

template<class PARTICLE>
bool BorisBuneman<PARTICLE>::initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,
					const std::string& regionName,const ParticleListBase* plist) {
   bool success = true;
   if (ParticlePropagatorBase::initialize(sim,simClasses,cr,regionName,plist) == false) success = false;
   species = reinterpret_cast<const Species*>(plist->getSpecies());
   return success;
}

#endif
