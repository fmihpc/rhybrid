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

template<class SPECIES,class PARTICLE>
class ParticleListHybrid: public ParticleListSkeleton<SPECIES,PARTICLE> {
 public:
   ParticleListHybrid();
   ~ParticleListHybrid();
   bool writeParticles(const std::string& spatMeshName);
};

template<class SPECIES,class PARTICLE> inline
ParticleListHybrid<SPECIES,PARTICLE>::ParticleListHybrid(): ParticleListSkeleton<SPECIES,PARTICLE>() { }

template<class SPECIES,class PARTICLE> inline
ParticleListHybrid<SPECIES,PARTICLE>::~ParticleListHybrid() { }

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

   // number of particles on this process
   const size_t Nparticles = this->size();

   // particles
   pargrid::DataWrapper<PARTICLE> wrapper = this->simClasses->pargrid.template getUserDataDynamic<PARTICLE>(this->particleDataID);
   PARTICLE** particleLists = wrapper.data();

   // block coordinates
   const double* crd = getBlockCoordinateArray(*this->sim,*this->simClasses);

   std::string particleMeshName = "ParticlePointMesh_" + this->speciesName;

   // write parameters of this particle population
   /*std::map<std::string,std::string> attribsParams;
   if (this->simClasses->pargrid.getRank() == this->sim->MASTER_RANK) {
      attribsParams["mesh"] = particleMeshName;
      attribsParams["name"] = "q";
      if(this->simClasses->vlsv.writeArray("PARAMETER",attribsParams,1,1,&this->species.q) == false) { success = false; }
   }
   else {
      if (this->simClasses->vlsv.writeArray("PARAMETER",attribsParams,0,0,&this->species.q) == false) { success = false; }
   }*/

   // write particle coordinates of this population as a point mesh
   uint64_t vectorSize = 3;
   Real* bufferMeshCrd = new Real[Nparticles*vectorSize];
   size_t cnt = 0;
   for (pargrid::CellID b=0; b<this->simClasses->pargrid.getNumberOfLocalCells(); ++b) {
      const size_t bv = vectorSize*b;
      // global coordinates of the block
      const Real xBlock = crd[bv+0];
      const Real yBlock = crd[bv+1];
      const Real zBlock = crd[bv+2];
      for (unsigned int p=0; p<wrapper.size()[b]; ++p) {
	 // global particle coordinates: local particle coordinates in the block + global coordinates of the block
	 bufferMeshCrd[cnt+0] = particleLists[b][p].state[particle::X] + xBlock;
	 bufferMeshCrd[cnt+1] = particleLists[b][p].state[particle::Y] + yBlock;
	 bufferMeshCrd[cnt+2] = particleLists[b][p].state[particle::Z] + zBlock;
	 cnt += vectorSize;
      }
   }
   std::map<std::string,std::string> attribs;
   attribs["name"] = particleMeshName;
   attribs["type"] = vlsv::mesh::STRING_POINT;
   if (this->simClasses->vlsv.writeArray("MESH",attribs,Nparticles,vectorSize,bufferMeshCrd) == false) {
      this->simClasses->logger << "\t ERROR failed to write particle species!" << std::endl;
      success = false;
   }
   delete [] bufferMeshCrd; bufferMeshCrd = NULL;

   // write particle properties of this population as variables on the point mesh
   vectorSize = 3;
   Real* bufferVar = new Real[Nparticles*vectorSize];
   cnt = 0;
   for (pargrid::CellID b=0; b<this->simClasses->pargrid.getNumberOfLocalCells(); ++b) {
      for (unsigned int p=0; p<wrapper.size()[b]; ++p) {
	 bufferVar[cnt+0] = particleLists[b][p].state[particle::VX];
	 bufferVar[cnt+1] = particleLists[b][p].state[particle::VY];
	 bufferVar[cnt+2] = particleLists[b][p].state[particle::VZ];
	 /*bufferVar[cnt+3] = particleLists[b][p].state[particle::WEIGHT];
	 bufferVar[cnt+4] = this->species.popid;
	 bufferVar[cnt+5] = static_cast<double>(this->simClasses->pargrid.getGlobalIDs()[b]);*/
	 cnt += vectorSize;
      }
   }
   attribs["name"] = "v_" + this->speciesName;
   attribs["mesh"] = particleMeshName;
   attribs["type"] = "pointdata";
   if (this->simClasses->vlsv.writeArray("VARIABLE",attribs,Nparticles,vectorSize,bufferVar) == false) {
      this->simClasses->logger << "\t ERROR failed to write particle species!" << std::endl;
      success = false;
   }
   delete [] bufferVar; bufferVar = NULL;

   // cell ids of particles
   /*vectorSize = 1;
   std::vector<pargrid::CellID> cellIDs; // Note: we assume here block size = 1 (i.e. blocks == cells)
   //bufferVar = new Real[Nparticles*vectorSize];
   //cnt = 0;
   for (pargrid::CellID b=0; b<this->simClasses->pargrid.getNumberOfLocalCells(); ++b) {
      for (unsigned int p=0; p<wrapper.size()[b]; ++p) {
	 cellIDs.push_back(this->simClasses->pargrid.getGlobalIDs()[b]);
	 //bufferVar[cnt+0] = static_cast<double>(this->simClasses->pargrid.getGlobalIDs()[b]);
	 //cnt += vectorSize;
      }
   }
   attribs["name"] = "CellID_" + this->speciesName;
   attribs["mesh"] = particleMeshName;
   attribs["type"] = "pointdata";
   if (this->simClasses->vlsv.writeArray("VARIABLE",attribs,Nparticles,vectorSize,&(cellIDs[0])) == false) {
      this->simClasses->logger << "\t ERROR failed to write particle species!" << std::endl;
      success = false;
   }
   //delete [] bufferVar; bufferVar = NULL;*/

   #if PROFILE_LEVEL > 0
      profile::stop();
   #endif
   return success;
}

#endif
