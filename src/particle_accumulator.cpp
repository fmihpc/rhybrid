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

#include "hybrid.h"
#include "particle_accumulator.h"

using namespace std;

int Accumulator::N_accumulators = 0;

Accumulator::Accumulator(): ParticleAccumulatorBase() {
   sim = NULL;
   simClasses = NULL;
   // Increase the total number of Accumulators:
   myOrderNumber = N_accumulators;
   ++N_accumulators;
   #if PROFILE_LEVEL > 0
      dataCopying = -1;
      particleAccumulation = -1;
   #endif
}
   
Accumulator::~Accumulator() {
   finalize();
}

#ifdef WRITE_POPULATION_AVERAGES
void Accumulator::accumulateCell(const Species& species,pargrid::CellID blockID,unsigned int N_particles,
				 const Particle<Real>* particles,Real* cellRhoQi,Real* cellJi,
				 Real* nAve,Real* vAve) {
#else
void Accumulator::accumulateCell(const Species& species,pargrid::CellID blockID,unsigned int N_particles,
				 const Particle<Real>* particles,Real* cellRhoQi,Real* cellJi) {
#endif
   const int accBlockSize  = (block::WIDTH_X+2)*(block::WIDTH_Y+2)*(block::WIDTH_Z+2);
   Real acc1[accBlockSize];
   Real acc2[accBlockSize*3];
   for(int i=0;i<accBlockSize;++i)   { acc1[i] = 0.0; }
   for(int i=0;i<accBlockSize*3;++i) { acc2[i] = 0.0; }
   const Real q = species.q;

#ifdef USE_DETECTORS
   /*bool* detPleFlag = reinterpret_cast<bool*>(simClasses->pargrid.getUserData(Hybrid::dataDetectorParticleFlagID));
   pargrid::DataWrapper<Dist> wrapperSpectra = simClasses->pargrid.getUserDataDynamic<Dist>(Hybrid::dataSpectraID);
   Dist* spectra = wrapperSpectra.data()[blockID];*/
#endif

   #if PROFILE_LEVEL > 1
      profile::start("accumulation",particleAccumulation);
   #endif
   for(unsigned int p=0;p<N_particles;++p) {
      const Real x = particles[p].state[particle::X];
      const Real y = particles[p].state[particle::Y];
      const Real z = particles[p].state[particle::Z];
      const Real w = particles[p].state[particle::WEIGHT];
      const Real wq = w*q;
      
      Real v[3];
      v[0] = particles[p].state[particle::VX];
      v[1] = particles[p].state[particle::VY];
      v[2] = particles[p].state[particle::VZ];
      int i = static_cast<int>(floor(x/Hybrid::dx));
      int j = static_cast<int>(floor(y/Hybrid::dx));
      int k = static_cast<int>(floor(z/Hybrid::dx));
      if(x < (i+0.5)*Hybrid::dx) { --i; }
      if(y < (j+0.5)*Hybrid::dx) { --j; }
      if(z < (k+0.5)*Hybrid::dx) { --k; }
      ++i; ++j; ++k;
      const Real x0 = (i-1+0.5)*Hybrid::dx;
      const Real y0 = (j-1+0.5)*Hybrid::dx;
      const Real z0 = (k-1+0.5)*Hybrid::dx;
      const Real w_x = (x-x0)/Hybrid::dx;
      const Real w_y = (y-y0)/Hybrid::dx;
      const Real w_z = (z-z0)/Hybrid::dx;
      
      // CIC weight factors
      const Real w000 = (1-w_x)*(1-w_y)*(1-w_z) * wq;
      const Real w100 =    w_x *(1-w_y)*(1-w_z) * wq;
      const Real w010 = (1-w_x)*   w_y *(1-w_z) * wq;
      const Real w110 =    w_x *   w_y *(1-w_z) * wq;
      const Real w001 = (1-w_x)*(1-w_y)*   w_z  * wq;
      const Real w101 =    w_x *(1-w_y)*   w_z  * wq;
      const Real w011 = (1-w_x)*   w_y *   w_z  * wq;
      const Real w111 =    w_x *  w_y  *   w_z  * wq;
      
      // indices
      const int ind000 = block::arrayIndex(i+0,j+0,k+0);
      const int ind100 = block::arrayIndex(i+1,j+0,k+0);
      const int ind010 = block::arrayIndex(i+0,j+1,k+0);
      const int ind110 = block::arrayIndex(i+1,j+1,k+0);
      const int ind001 = block::arrayIndex(i+0,j+0,k+1);
      const int ind101 = block::arrayIndex(i+1,j+0,k+1);
      const int ind011 = block::arrayIndex(i+0,j+1,k+1);
      const int ind111 = block::arrayIndex(i+1,j+1,k+1);
      
      // rhoq
      acc1[ind000] += w000;
      acc1[ind100] += w100;
      acc1[ind010] += w010;
      acc1[ind110] += w110;
      acc1[ind001] += w001;
      acc1[ind101] += w101;
      acc1[ind011] += w011;
      acc1[ind111] += w111;
      
      // Ji
      for(int l=0;l<3;++l) {
	 acc2[ind000*3+l] += w000*v[l];
	 acc2[ind100*3+l] += w100*v[l];
	 acc2[ind010*3+l] += w010*v[l];
	 acc2[ind110*3+l] += w110*v[l];
	 acc2[ind001*3+l] += w001*v[l];
	 acc2[ind101*3+l] += w101*v[l];
	 acc2[ind011*3+l] += w011*v[l];
	 acc2[ind111*3+l] += w111*v[l];
      }
#ifdef USE_DETECTORS
      /*if(detPleFlag[blockID] == true && Hybrid::detParticleRecording == true) {
         spectra[species.popid].f[0] += w;
      }*/
#endif
   }
   
   #if PROFILE_LEVEL > 1
      profile::stop();
      profile::start("data copying",dataCopying);
   #endif
   
   block::addValues3D(*simClasses,blockID,acc1,cellRhoQi,1);
   block::addValues3D(*simClasses,blockID,acc2,cellJi,3);

#ifdef WRITE_POPULATION_AVERAGES
   if(nAve != NULL) { block::addValues3D(*simClasses,blockID,acc1,nAve,1); }
   if(vAve != NULL) { block::addValues3D(*simClasses,blockID,acc2,vAve,3); }
#endif
   
   #if PROFILE_LEVEL > 1
      profile::stop();
   #endif
}

/* // NGP accumulator
void Accumulator::accumulateCell(const Species& species,pargrid::CellID blockID,unsigned int N_particles,
				 const Particle<Real>* particles,Real* cellRhoQi,Real* cellJi) {
   const int accBlockSize  = (block::WIDTH_X+2)*(block::WIDTH_Y+2)*(block::WIDTH_Z+2);
   Real acc1[accBlockSize];
   Real acc2[accBlockSize*3];
   for(int i=0;i<accBlockSize;++i)   { acc1[i] = 0.0; }
   for(int i=0;i<accBlockSize*3;++i) { acc2[i] = 0.0; }
   const Real q = species.q;
   
   #if PROFILE_LEVEL > 1
      profile::start("accumulation",particleAccumulation);
   #endif
   for(unsigned int p=0;p<N_particles;++p) {
      const Real x = particles[p].state[particle::X];
      const Real y = particles[p].state[particle::Y];
      const Real z = particles[p].state[particle::Z];
      const Real wq = particles[p].state[particle::WEIGHT]*q;
      
      Real v[3];
      v[0] = particles[p].state[particle::VX];
      v[1] = particles[p].state[particle::VY];
      v[2] = particles[p].state[particle::VZ];
      const int i = 1+static_cast<int>(floor(x/Hybrid::dx));
      const int j = 1+static_cast<int>(floor(y/Hybrid::dx));
      const int k = 1+static_cast<int>(floor(z/Hybrid::dx));
      const int ind000 = block::arrayIndex(i+0,j+0,k+0);
      acc1[ind000] += wq;
      for(int l=0;l<3;++l) { acc2[ind000*3+l] += wq*v[l]; }
   }
   
   #if PROFILE_LEVEL > 1
      profile::stop();
      profile::start("data copying",dataCopying);
   #endif
   
   block::addValues3D(*simClasses,blockID,acc1,cellRhoQi,1);
   block::addValues3D(*simClasses,blockID,acc2,cellJi,3);
   
   #if PROFILE_LEVEL > 1
      profile::stop();
   #endif
}*/

bool Accumulator::accumulateBoundaryCells(pargrid::DataID particleDataID,const unsigned int* N_particles) {
   bool success = true;
   if(species->accumulate == false) { return success; }
   pargrid::DataWrapper<Particle<Real> > wrapper = simClasses->pargrid.getUserDataDynamic<Particle<Real> >(particleDataID);
   Real* cellJi    = simClasses->pargrid.getUserDataStatic<Real>(Hybrid::dataCellJiID);
   Real* cellRhoQi = simClasses->pargrid.getUserDataStatic<Real>(Hybrid::dataCellRhoQiID);
#ifdef WRITE_POPULATION_AVERAGES
   Real* nAve = NULL;
   Real* vAve = NULL;
   const int m = Hybrid::outputPopVarId[species->popid-1];
   if(m >=0 ) {
      nAve = simClasses->pargrid.getUserDataStatic<Real>(Hybrid::dataCellAverageDensityID[m]);
      vAve = simClasses->pargrid.getUserDataStatic<Real>(Hybrid::dataCellAverageVelocityID[m]);
   }
#endif
   const vector<pargrid::CellID>& boundaryBlocks = simClasses->pargrid.getBoundaryCells(Hybrid::accumulationStencilID);
   for(pargrid::CellID b=0; b<boundaryBlocks.size(); ++b) {
      const pargrid::CellID block = boundaryBlocks[b];
#ifdef WRITE_POPULATION_AVERAGES
      accumulateCell(*species,block,N_particles[block],wrapper.data()[block],cellRhoQi,cellJi,nAve,vAve);
#else
      accumulateCell(*species,block,N_particles[block],wrapper.data()[block],cellRhoQi,cellJi);
#endif
   }
   return success;
}

bool Accumulator::accumulateInnerCells(pargrid::DataID particleDataID,const unsigned int* N_particles) {
   bool success = true;
   if(species->accumulate == false) { return success; }
   pargrid::DataWrapper<Particle<Real> > wrapper = simClasses->pargrid.getUserDataDynamic<Particle<Real> >(particleDataID);
   Real* cellJi    = simClasses->pargrid.getUserDataStatic<Real>(Hybrid::dataCellJiID);
   Real* cellRhoQi = simClasses->pargrid.getUserDataStatic<Real>(Hybrid::dataCellRhoQiID);
#ifdef WRITE_POPULATION_AVERAGES
   Real* nAve = NULL;
   Real* vAve = NULL;
   const int m = Hybrid::outputPopVarId[species->popid-1];
   if(m >=0 ) {
      nAve = simClasses->pargrid.getUserDataStatic<Real>(Hybrid::dataCellAverageDensityID[m]);
      vAve = simClasses->pargrid.getUserDataStatic<Real>(Hybrid::dataCellAverageVelocityID[m]);
   }
#endif
   const vector<pargrid::CellID>& innerBlocks = simClasses->pargrid.getInnerCells(Hybrid::accumulationStencilID);
   for(pargrid::CellID b=0; b<innerBlocks.size(); ++b) {
      const pargrid::CellID block = innerBlocks[b];
#ifdef WRITE_POPULATION_AVERAGES
      accumulateCell(*species,block,N_particles[block],wrapper.data()[block],cellRhoQi,cellJi,nAve,vAve);
#else
      accumulateCell(*species,block,N_particles[block],wrapper.data()[block],cellRhoQi,cellJi);
#endif
   }
   return success;
}

bool Accumulator::addConfigFileItems(ConfigReader& cr,const std::string& regionName) {
   return true;
}

bool Accumulator::addRemoteUpdates() {
   bool success = true;
   // Only the very last Accumulator is allowed to add updates from remote processes:
   if (myOrderNumber != N_accumulators-1) { return success; }
   Real* cellJi    = simClasses->pargrid.getUserDataStatic<Real>(Hybrid::dataCellJiID);
   Real* cellRhoQi = simClasses->pargrid.getUserDataStatic<Real>(Hybrid::dataCellRhoQiID);
   Real* cellRhoQiBg = simClasses->pargrid.getUserDataStatic<Real>(Hybrid::dataCellRhoQiBgID);
   Real* counterCellMinRhoQi = simClasses->pargrid.getUserDataStatic<Real>(Hybrid::dataCounterCellMinRhoQiID);
   bool* outerBoundaryFlag   = simClasses->pargrid.getUserDataStatic<bool>(Hybrid::dataOuterBoundaryFlagID);
   unsigned int* offsetsCellRhoQi = NULL;
   Real* buffersCellRhoQi = NULL;
   unsigned int* offsetsCellJi = NULL;
   Real* buffersCellJi = NULL;
   if(simClasses->pargrid.getRemoteUpdates<Real>(Hybrid::accumulationStencilID,Hybrid::dataCellRhoQiID,offsetsCellRhoQi,buffersCellRhoQi) == false) {
      simClasses->logger << "ERROR: Failed to get remote updates of cellRhoQi from ParGrid" << endl << write;
      success = false;
   }
   if(simClasses->pargrid.getRemoteUpdates<Real>(Hybrid::accumulationStencilID,Hybrid::dataCellJiID,offsetsCellJi,buffersCellJi) == false) {
      simClasses->logger << "ERROR: Failed to get remote updates of cellJi from ParGrid" << endl << write;
      success = false;
   }
   const vector<pargrid::CellID>& boundaryBlocks = simClasses->pargrid.getBoundaryCells(Hybrid::accumulationStencilID);
   for(pargrid::CellID b=0;b<boundaryBlocks.size();++b) {
      const pargrid::CellID blockLID = boundaryBlocks[b];
      for(unsigned int i=offsetsCellRhoQi[b];i<offsetsCellRhoQi[b+1];++i) for(int j=0;j<block::SIZE;++j) {
	 cellRhoQi[blockLID*block::SIZE+j] += buffersCellRhoQi[i*block::SIZE + j];
      }
      for(unsigned int i=offsetsCellJi[b];i<offsetsCellJi[b+1];++i) for(int j=0;j<block::SIZE;++j) for(int l=0;l<3;++l) {
	cellJi[(blockLID*block::SIZE+j)*3+l] += buffersCellJi[(i*block::SIZE + j)*3+l];
      }
   }
#ifdef WRITE_POPULATION_AVERAGES
   for(unsigned int m = 0;m<Hybrid::N_outputPopVars;++m) {
      Real* nAve = simClasses->pargrid.getUserDataStatic<Real>(Hybrid::dataCellAverageDensityID[m]);
      Real* vAve = simClasses->pargrid.getUserDataStatic<Real>(Hybrid::dataCellAverageVelocityID[m]);
      unsigned int* offsetsCellAverageDensity = NULL;
      unsigned int* offsetsCellAverageVelocity = NULL;
      Real* buffersCellAverageDensity = NULL;
      Real* buffersCellAverageVelocity = NULL;
      if(simClasses->pargrid.getRemoteUpdates<Real>(Hybrid::accumulationStencilID,Hybrid::dataCellAverageDensityID[m],offsetsCellAverageDensity,buffersCellAverageDensity) == false) {
         simClasses->logger << "ERROR: Failed to get remote updates of cellAverageDensity from ParGrid" << endl << write;
         success = false;
      }
      if(simClasses->pargrid.getRemoteUpdates<Real>(Hybrid::accumulationStencilID,Hybrid::dataCellAverageVelocityID[m],offsetsCellAverageVelocity,buffersCellAverageVelocity) == false) {
         simClasses->logger << "ERROR: Failed to get remote updates of cellAverageVelocity from ParGrid" << endl << write;
         success = false;
      }
      for(pargrid::CellID b=0;b<boundaryBlocks.size();++b) {
         const pargrid::CellID blockLID = boundaryBlocks[b];
         for(unsigned int i=offsetsCellAverageDensity[b];i<offsetsCellAverageDensity[b+1];++i) for(int j=0;j<block::SIZE;++j) {
            nAve[blockLID*block::SIZE+j] += buffersCellAverageDensity[i*block::SIZE + j];
         }
         for(unsigned int i=offsetsCellAverageVelocity[b];i<offsetsCellAverageVelocity[b+1];++i) for(int j=0;j<block::SIZE;++j) for(int l=0;l<3;++l) {
            vAve[(blockLID*block::SIZE+j)*3+l] += buffersCellAverageVelocity[(i*block::SIZE + j)*3+l];
         }
      }
   }
#endif   
   // Finalize accumulate:
   for(pargrid::CellID b=0;b<simClasses->pargrid.getNumberOfLocalCells();++b) {
      for(int k=0;k<block::WIDTH_Z;++k) for(int j=0;j<block::WIDTH_Y;++j) for(int i=0;i<block::WIDTH_X;++i) {
	 const int n = (b*block::SIZE+block::index(i,j,k));
	 const int n3 = n*3;
	 cellRhoQi[n] /= Hybrid::dV;
         cellRhoQi[n] += cellRhoQiBg[n];
         if(outerBoundaryFlag[n] == true) {
            if(cellRhoQi[n] < Hybrid::outerBoundaryZone.minRhoQi) {
               cellRhoQi[n] = Hybrid::outerBoundaryZone.minRhoQi;
               counterCellMinRhoQi[n]++;
            }
         }
	 else if(cellRhoQi[n] < Hybrid::minRhoQi) {
	    cellRhoQi[n] = Hybrid::minRhoQi;
	    counterCellMinRhoQi[n]++;
	 }
	 for(int l=0;l<3;++l) { cellJi[n3+l] /= Hybrid::dV; }
      }
   }
   return success;
}

bool Accumulator::clearAccumulationArrays() {
   bool success = true;
   // Only the first Accumulator to reach this point needs to clear arrays:
   if (myOrderNumber == 0) {
      Real* cellJi    = simClasses->pargrid.getUserDataStatic<Real>(Hybrid::dataCellJiID);
      Real* cellRhoQi = simClasses->pargrid.getUserDataStatic<Real>(Hybrid::dataCellRhoQiID);   
      for (pargrid::CellID b=0; b<simClasses->pargrid.getNumberOfAllCells()*block::SIZE; ++b)   { cellRhoQi[b] = 0.0; }
      for (pargrid::CellID b=0; b<simClasses->pargrid.getNumberOfAllCells()*block::SIZE*3; ++b) { cellJi[b] = 0.0; }
   }
#ifdef WRITE_POPULATION_AVERAGES
   // zero buffer cells for average accumulation arrays
   const int m = Hybrid::outputPopVarId[species->popid-1];
   if(m >= 0) {
      Real* nAve = simClasses->pargrid.getUserDataStatic<Real>(Hybrid::dataCellAverageDensityID[m]);
      Real* vAve = simClasses->pargrid.getUserDataStatic<Real>(Hybrid::dataCellAverageVelocityID[m]);
      size_t startBufferCells = simClasses->pargrid.getNumberOfLocalCells()*block::SIZE;
      size_t endBufferCells   = simClasses->pargrid.getNumberOfAllCells()  *block::SIZE;
      for (pargrid::CellID b=startBufferCells;b<endBufferCells;++b) {
         const size_t b3 = b*3;
         nAve[b] = 0.0;
         for(size_t l=0;l<3;++l) { vAve[b3+l] = 0.0; }
      }
   }
#endif
   return success;
}

bool Accumulator::finalize() {
   bool success = true;
   sim = NULL;
   simClasses = NULL;
   return success;
}

bool Accumulator::initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,
			     const std::string& regionName,const ParticleListBase* plist) {
   bool success = true;
   if (ParticleAccumulatorBase::initialize(sim,simClasses,cr,regionName,plist) == false) {
      simClasses.logger << "(RHYBRID ACCUMULATOR) ERROR: ParticleAccumulatorBase initialization failed" << endl << write;
   }
   this->species = reinterpret_cast<const Species*>(plist->getSpecies());
   return success;
}

bool Accumulator::sendUpdates() {
   bool success = true;
   // Only the very last Accumulator sends updates to remote processes:
   if(myOrderNumber != N_accumulators-1) { return success; }
   if(simClasses->pargrid.startNeighbourExchange(Hybrid::accumulationStencilID,Hybrid::dataCellRhoQiID) == false) { success = false; }
   if(simClasses->pargrid.startNeighbourExchange(Hybrid::accumulationStencilID,Hybrid::dataCellJiID) == false) { success = false; }
#ifdef WRITE_POPULATION_AVERAGES
   for(unsigned int m = 0;m<Hybrid::N_outputPopVars;++m) {
      if(simClasses->pargrid.startNeighbourExchange(Hybrid::accumulationStencilID,Hybrid::dataCellAverageDensityID[m]) == false) { success = false; }
      if(simClasses->pargrid.startNeighbourExchange(Hybrid::accumulationStencilID,Hybrid::dataCellAverageVelocityID[m]) == false) { success = false; }
   }
#endif
   return success;
}

bool Accumulator::wait() {
   bool success = true;
   // Only the very last Accumulator sends updates to remote processes:
   if(myOrderNumber != N_accumulators-1) { return success; }
   if(simClasses->pargrid.wait(Hybrid::accumulationStencilID,Hybrid::dataCellRhoQiID) == false) { success = false; }
   if(simClasses->pargrid.wait(Hybrid::accumulationStencilID,Hybrid::dataCellJiID) == false) { success = false; }
#ifdef WRITE_POPULATION_AVERAGES
   for(unsigned int m = 0;m<Hybrid::N_outputPopVars;++m) {
      if(simClasses->pargrid.wait(Hybrid::accumulationStencilID,Hybrid::dataCellAverageDensityID[m]) == false) { success = false; }
      if(simClasses->pargrid.wait(Hybrid::accumulationStencilID,Hybrid::dataCellAverageVelocityID[m]) == false) { success = false; }
   }
#endif
   return success;
}
