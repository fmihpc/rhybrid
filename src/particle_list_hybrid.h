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

#ifndef PARTICLE_LIST_HYBRID_H
#define PARTICLE_LIST_HYBRID_H

#include <particle_list_skeleton.h>
#include <hybrid.h>

namespace hybsave {
   enum VARIABLES {
      XPOS,
      YPOS,
      ZPOS,
      VX,
      VY,
      VZ,
      WEIGHT,
      POPID,
      BLOCKID,
      SIZE
   };
}

template<class SPECIES,class PARTICLE>
class ParticleListHybrid: public ParticleListSkeleton<SPECIES,PARTICLE> {
 public:
   ParticleListHybrid();
   ~ParticleListHybrid();
   bool writeParticles(const std::string& spatMeshName);
};

template<class SPECIES,class PARTICLE> inline
ParticleListHybrid<SPECIES,PARTICLE>::ParticleListHybrid(): ParticleListSkeleton<SPECIES,PARTICLE>() {
   
}

template<class SPECIES,class PARTICLE> inline
ParticleListHybrid<SPECIES,PARTICLE>::~ParticleListHybrid() {
}

template<class SPECIES,class PARTICLE> inline
bool ParticleListHybrid<SPECIES,PARTICLE>::writeParticles(const std::string& spatMeshName) {
   bool success = true;
   if (this->initialized == false) { return false; }

   // check if particles are to be written
   bool writeFlag = 0;
   this->cr->add("Simulation.save_particles","Write particles or not (bool)",0);
   this->cr->parse();
   this->cr->get("Simulation.save_particles",writeFlag);
   if(writeFlag == false) { return success; }
   
   #if PROFILE_LEVEL > 0
      profile::start(this->speciesName+" writing",this->particleWriteID);
   #endif

   const size_t sum_particles = this->size();
   Real* buffer = new Real[sum_particles*hybsave::SIZE];
   const double* crd = getBlockCoordinateArray(*this->sim,*this->simClasses);
   pargrid::DataWrapper<PARTICLE> wrapper = this->simClasses->pargrid.template getUserDataDynamic<PARTICLE>(this->particleDataID);

   size_t counter = 0;
   PARTICLE** particleLists = wrapper.data();
   for (size_t block=0; block<this->simClasses->pargrid.getNumberOfLocalCells(); ++block) {
      for (unsigned int p=0; p<wrapper.size()[block]; ++p) {
#ifdef USE_INDEX_OPERATOR
	 buffer[counter*3+XPOS] = particleLists[block][p][XPOS] + crd[3*block+XPOS];
	 buffer[counter*3+YPOS] = particleLists[block][p][YPOS] + crd[3*block+YPOS];
	 buffer[counter*3+ZPOS] = particleLists[block][p][ZPOS] + crd[3*block+ZPOS];
#else
	 buffer[counter+hybsave::XPOS] = particleLists[block][p].state[particle::X] + crd[3*block+0];
	 buffer[counter+hybsave::YPOS] = particleLists[block][p].state[particle::Y] + crd[3*block+1];
	 buffer[counter+hybsave::ZPOS] = particleLists[block][p].state[particle::Z] + crd[3*block+2];
	 buffer[counter+hybsave::VX] = particleLists[block][p].state[particle::VX];
	 buffer[counter+hybsave::VY] = particleLists[block][p].state[particle::VY];
	 buffer[counter+hybsave::VZ] = particleLists[block][p].state[particle::VZ];
	 buffer[counter+hybsave::WEIGHT] = particleLists[block][p].state[particle::WEIGHT];
	 buffer[counter+hybsave::POPID] = this->species.popid;
	 buffer[counter+hybsave::BLOCKID] = static_cast<double>(this->simClasses->pargrid.getGlobalIDs()[block]);
#endif
	 counter += hybsave::SIZE;
      }
   }
   
   std::map<std::string,std::string> attribs;
   attribs["name"] = this->speciesName;
   attribs["type"] = vlsv::mesh::STRING_POINT;
   if (this->simClasses->vlsv.writeArray("MESH",attribs,sum_particles,hybsave::SIZE,buffer) == false) {
      this->simClasses->logger << "\t ERROR failed to write particle species!" << std::endl;
      success = false;
   }
   delete [] buffer; buffer = NULL;
   
   #if PROFILE_LEVEL > 0
      profile::stop();
   #endif
   return success;
}

#endif
