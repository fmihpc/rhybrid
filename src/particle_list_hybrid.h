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
   if(Hybrid::saveParticles == false) { return success; }

   #if PROFILE_LEVEL > 0
      profile::start(this->speciesName+" writing",this->particleWriteID);
   #endif

   // counter of number of particles on this process to be saved with Nstride (total number: this->size())
   size_t Nparticles = 0;

   // particles
   pargrid::DataWrapper<PARTICLE> wrapper = this->simClasses->pargrid.template getUserDataDynamic<PARTICLE>(this->particleDataID);
   PARTICLE** particleLists = wrapper.data();

   // block coordinates
   const double* crd = getBlockCoordinateArray(*this->sim,*this->simClasses);

   // domain dimensions
   const long nx = this->sim->x_blocks*block::WIDTH_X;
   const long ny = this->sim->y_blocks*block::WIDTH_Y;
   // z-dimension not needed

   // name of the point mesh of this particle population
   std::string particleMeshName = "ParticlePointMesh_" + this->speciesName;

   // global cell ids
   const std::vector<pargrid::CellID>& globalIDs = this->simClasses->pargrid.getGlobalIDs();

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
   uint64_t vectorSize = 3; // mesh has three coordinates: x, y, z
   std::vector<Real> bufferMeshCrd; // buffer variable of data to be written
   for (pargrid::CellID b=0; b<this->simClasses->pargrid.getNumberOfLocalCells(); ++b) {
      pargrid::CellID cid = globalIDs[b];
      if (cid % Hybrid::saveParticlesNstride != 0 ||
          cid % (Hybrid::saveParticlesNstride*nx) >= nx ||
          cid % (Hybrid::saveParticlesNstride*nx*ny) >= nx*ny) { continue; }
      const size_t bv = vectorSize*b;
      // global coordinates of the block
      const Real xBlock = crd[bv+0];
      const Real yBlock = crd[bv+1];
      const Real zBlock = crd[bv+2];
      for (unsigned int p=0; p<wrapper.size()[b]; ++p) {
	 // global particle coordinates: local particle coordinates in the block + global coordinates of the block
	 bufferMeshCrd.push_back(particleLists[b][p].state[particle::X] + xBlock);
	 bufferMeshCrd.push_back(particleLists[b][p].state[particle::Y] + yBlock);
	 bufferMeshCrd.push_back(particleLists[b][p].state[particle::Z] + zBlock);
	 Nparticles++;
      }
   }
   std::map<std::string,std::string> attribs;
   attribs["name"] = particleMeshName;
   attribs["type"] = vlsv::mesh::STRING_POINT;
   if (this->simClasses->vlsv.writeArray("MESH",attribs,Nparticles,vectorSize,&(bufferMeshCrd[0])) == false) {
      this->simClasses->logger << "\t ERROR failed to write particle species!" << std::endl;
      success = false;
   }

   // write particle velocities of this population
   vectorSize = 3; // velocity has three components: vx, vy, vz
   std::vector<Real> bufferVar;
   for (pargrid::CellID b=0; b<this->simClasses->pargrid.getNumberOfLocalCells(); ++b) {
      pargrid::CellID cid = globalIDs[b];
      if (cid % Hybrid::saveParticlesNstride != 0 ||
          cid % (Hybrid::saveParticlesNstride*nx) >= nx ||
          cid % (Hybrid::saveParticlesNstride*nx*ny) >= nx*ny) { continue; }
      for (unsigned int p=0; p<wrapper.size()[b]; ++p) {
	 bufferVar.push_back(particleLists[b][p].state[particle::VX]);
	 bufferVar.push_back(particleLists[b][p].state[particle::VY]);
	 bufferVar.push_back(particleLists[b][p].state[particle::VZ]);
      }
   }
   attribs["name"] = "v_" + this->speciesName;
   attribs["mesh"] = particleMeshName;
   attribs["type"] = "pointdata";
   if (this->simClasses->vlsv.writeArray("VARIABLE",attribs,Nparticles,vectorSize,&(bufferVar[0])) == false) {
      this->simClasses->logger << "\t ERROR failed to write particle species!" << std::endl;
      success = false;
   }

   // write ids of spatial cells where each particle belongs to
   vectorSize = 1; // id is a scalar variable of type pargrid::CellID (unsigned int)
   std::vector<pargrid::CellID> cellIDs; // Note: we assume here block size = 1 (i.e. blocks == cells)
   for (pargrid::CellID b=0; b<this->simClasses->pargrid.getNumberOfLocalCells(); ++b) {
      pargrid::CellID cid = globalIDs[b];
      if (cid % Hybrid::saveParticlesNstride != 0 ||
          cid % (Hybrid::saveParticlesNstride*nx) >= nx ||
          cid % (Hybrid::saveParticlesNstride*nx*ny) >= nx*ny) { continue; }
      for (unsigned int p=0; p<wrapper.size()[b]; ++p) {
	 cellIDs.push_back(globalIDs[b]);
      }
   }
   attribs["name"] = "CellID_" + this->speciesName;
   attribs["mesh"] = particleMeshName;
   attribs["type"] = "pointdata";
   if (this->simClasses->vlsv.writeArray("VARIABLE",attribs,Nparticles,vectorSize,&(cellIDs[0])) == false) {
      this->simClasses->logger << "\t ERROR failed to write particle species!" << std::endl;
      success = false;
   }

   #if PROFILE_LEVEL > 0
      profile::stop();
   #endif
   return success;
}

#endif
