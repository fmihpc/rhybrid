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

#ifndef DIAGNOSTICS_H
#define DIAGNOSTICS_H

#include <cstdlib>
#include <iostream>
#include <map>
#include <vector>
#include <limits>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include "hybrid.h"

namespace diagnostics
{

// plasma bulk, wave, etc, collective parameters
struct PlasmaParametersBulk {
   Real B[3] = {0.0,0.0,0.0};
   Real Ubulk[3] = {0.0,0.0,0.0};
   Real Ec[3] = {0.0,0.0,0.0};
   Real vExB[3] = {0.0,0.0,0.0};
   Real Btot=0.0,Ubulktot=0.0,vExBtot=0.0,Ectot=0.0;
   Real ne=0.0,Te=0.0,rhoq=0.0,rhom=0.0;
   Real vA=0.0,vs=0.0,vms=0.0,vpui=0.0,vw=0.0;
   Real MA=0.0,Ms=0.0,Mms=0.0;
};

// plasma singe particle, per species, per population, etc, parameters
struct PlasmaParametersSingleParticle {
   std::vector<std::string> populationName;
   std::vector<Real> periodPlasma;
   std::vector<Real> periodLarmor;
   std::vector<Real> lengthInertial;
   std::vector<Real> radiusLarmorThermal;
   std::vector<Real> radiusLarmorPickUp;
};

// struct to store some particle population log vales
struct LogDataParticle {
   Real N_macroParticles,N_realParticles,sumVx,sumVy,sumVz,sumV,sumWV2,maxVi;
};

// struct to store some field log vales
struct LogDataField {
   Real N_cells;
   Real sumBx,sumBy,sumBz,sumB,maxB,sumDivB,maxDivB,maxDivBPerB,sumB2;
   Real sumCellJix,sumCellJiy,sumCellJiz,sumCellJi,maxCellJi;
   Real sumCellEpx,sumCellEpy,sumCellEpz,sumCellEp,maxCellEp,sumCellEp2;
   Real sumNodeEx,sumNodeEy,sumNodeEz,sumNodeE,maxNodeE,sumNodeE2,sumCellE2;
   Real minInerLengthElectron,maxInerLengthElectron,minInerLengthProton,maxInerLengthProton,minTLarmor,maxVAlfven,maxUe;
};

// determine average plasma parameters (bulk, waves, etc) from all populations of a selected type ("solarwind" or "uniform")
bool calcPlasmaParametersBulk(SimulationClasses& simClasses,std::string injectorType,std::vector<ParticleListBase*>& particleLists,Real Bx,Real By,Real Bz,Real dx,PlasmaParametersBulk& ppB) {
   // magnetic field
   Real B[3] = {Bx,By,Bz};
   Real Btot = normvec(B);
   // plasma bulk properties
     {
	for(size_t s=0;s<particleLists.size();++s) {
	   // read parameters of this particle population from species and injector
	   InjectorParameters ip;
	   if(getInjectorParameters(particleLists[s]->getInjector(),ip) == false) {
	      simClasses.logger << "(RHYBRID) ERROR: failed to get injector and species parameters (" << ip.name << ")" << std::endl;
	      return false;
	   }
	   if(ip.type.compare(injectorType) != 0) { continue; } // include only chosen type populations
	   ppB.rhoq += ip.n*ip.q; // ion charge density
	   ppB.rhom += ip.n*ip.m; // ion mass density
	   for(size_t ii=0;ii<3;++ii) {
	      ppB.Ubulk[ii] += ip.n*ip.m*ip.velocity[ii]; // bulk velocity (mass density weighted average)
	   }
	   ppB.vs += ip.n*ip.T; // sound velocity
	}
	ppB.ne = ppB.rhoq/constants::CHARGE_ELEMENTARY; // electron number density
	if(ppB.rhom != 0) {
	   for(size_t ii=0;ii<3;++ii) { ppB.Ubulk[ii] /= ppB.rhom; }
	   ppB.vA = Btot/( sqrt(constants::PERMEABILITY*ppB.rhom) ); // alfven velocity
	   ppB.vs = sqrt(5.0/3.0 * constants::BOLTZMANN*ppB.vs/ppB.rhom);
	}
	else {
	   for(size_t ii=0;ii<3;++ii) { ppB.Ubulk[ii] = 0.0; }
	   ppB.vA = 0.0;
	   ppB.vs = 0.0;
	}
	ppB.Ubulktot = normvec(ppB.Ubulk);
	ppB.vms = sqrt( sqr(ppB.vA) + sqr(ppB.vs) ); // magnetosonic velocity
	ppB.MA = (ppB.vA != 0) ? ppB.Ubulktot/ppB.vA : 0.0; // Alfven Mach number
	ppB.Ms = (ppB.vs != 0) ? ppB.Ubulktot/ppB.vs : 0.0; // sonic Mach number
	ppB.Mms = (ppB.vms != 0) ? ppB.Ubulktot/ppB.vms : 0.0; // magnetosonic mach number
	cross(B,ppB.Ubulk,ppB.Ec); // convection electric field Ec = -Ubulk x B
	ppB.Ectot = normvec(ppB.Ec);
	// vExB = E x B/B^2
	if(Btot != 0) {
	   cross(ppB.Ec,B,ppB.vExB);
	   for(size_t ii=0;ii<3;++ii) { ppB.vExB[ii] /= sqr(Btot); }
	}
	ppB.vExBtot = normvec(ppB.vExB);
	ppB.vpui = 2*ppB.vExBtot; // fastest pickup ion velocity
	// fastest whistler signal p. 28 Alho (2016)
	if(ppB.ne != 0 && dx != 0) {
	   ppB.vw = 2.0*Btot*M_PI/( constants::PERMEABILITY*ppB.ne*constants::CHARGE_ELEMENTARY*dx );
	}
     }
   return true;
}

// determine plasma parameters (per population, per species, single particle, etc)
bool calcPlasmaParametersSingleParticle(SimulationClasses& simClasses,std::vector<ParticleListBase*>& particleLists,Real Bx,Real By,Real Bz,Real ne,Real vExBtot,Real Te,PlasmaParametersSingleParticle& ppSP) {
   // magnetic field
   Real B[3] = {Bx,By,Bz};
   Real Btot = normvec(B);

   // electrons (not a particle population in the model but for information purposes included)
     {
	ppSP.populationName.push_back("e-");
	const Real m = constants::MASS_ELECTRON;
	const Real q = constants::CHARGE_ELEMENTARY;
	const Real qB = q*Btot;
	const Real wP = sqrt( ne*sqr(q)/( m*constants::PERMITTIVITY  ) ); // angular plasma frequency
	const Real vth = sqrt(constants::BOLTZMANN*Te/m); // thermal speed
	ppSP.periodPlasma.push_back((wP != 0) ? 2.0*M_PI/wP : 0.0); // plasma period
	ppSP.periodLarmor.push_back((qB != 0) ? 2.0*M_PI*m/qB : 0.0); // Larmor period
	ppSP.lengthInertial.push_back((wP != 0) ? constants::SPEED_LIGHT/wP : 0.0); // inertial length
	ppSP.radiusLarmorThermal.push_back((qB != 0) ? m*vth/qB : 0.0); // thermal Larmor radius
	ppSP.radiusLarmorPickUp.push_back((qB != 0) ? m*vExBtot/qB : 0.0); // PU Larmor radius
     }

   // simulation particle populations
   for(size_t s=0;s<particleLists.size();++s) {
      // read parameters of this particle population from species and injector
      InjectorParameters ip;
      if(getInjectorParameters(particleLists[s]->getInjector(),ip) == false) {
	 simClasses.logger << "(RHYBRID) ERROR: failed to get injector and species parameters (" << ip.name << ")" << std::endl;
	 return false;
      }
	{
	   ppSP.populationName.push_back(ip.name);
	   const Real m = ip.m;
	   const Real q = ip.q;
	   const Real qB = q*Btot;
	   const Real wP = sqrt( ne*sqr(q)/( m*constants::PERMITTIVITY  ) ); // angular plasma frequency (using ne as density)
	   const Real vth = ip.vth; // thermal speed
	   ppSP.periodPlasma.push_back((wP != 0) ? 2.0*M_PI/wP : 0.0); // plasma period
	   ppSP.periodLarmor.push_back((qB != 0) ? 2.0*M_PI*m/qB : 0.0); // Larmor period
	   ppSP.lengthInertial.push_back((wP != 0) ? constants::SPEED_LIGHT/wP : 0.0); // inertial length
	   ppSP.radiusLarmorThermal.push_back((qB != 0) ? m*vth/qB : 0.0); // thermal Larmor radius
	   ppSP.radiusLarmorPickUp.push_back((qB != 0) ? m*vExBtot/qB : 0.0); // PU Larmor radius
	}
   }

   // sanity check just in case
   if(ppSP.populationName.size() != ppSP.periodPlasma.size() ||
      ppSP.populationName.size() != ppSP.periodLarmor.size() ||
      ppSP.populationName.size() != ppSP.lengthInertial.size() ||
      ppSP.populationName.size() != ppSP.radiusLarmorThermal.size() ||
      ppSP.populationName.size() != ppSP.radiusLarmorPickUp.size()
     ) {
      simClasses.logger << "(RHYBRID) ERROR: calcPlasmaParametersSingleParticle function failed" << std::endl;
      return false;
   }
   return true;
}

// log amounts of macroparticles
void logWriteMainMacroparticles(Simulation& sim,SimulationClasses& simClasses,const std::vector<ParticleListBase*>& particleLists) {
    simClasses.logger << "(RHYBRID) Number of macroparticles per population (time step = " << sim.timestep << ", time = " << sim.t << "):" << std::endl;
    for(size_t s=0;s<particleLists.size();++s) {
	Real N_macroParticles = 0.0;
	// For now skip particles with invalid data id:
	pargrid::DataID speciesDataID = pargrid::INVALID_DATAID;
	if(particleLists[s]->getParticles(speciesDataID) == true) {
	    pargrid::DataWrapper<Particle<Real> > wrapper = simClasses.pargrid.getUserDataDynamic<Particle<Real> >(speciesDataID);
	    for(pargrid::CellID b=0; b<simClasses.pargrid.getNumberOfLocalCells(); ++b) {
		pargrid::ArraySizetype N_particles = wrapper.size(b);
		N_macroParticles += N_particles;
	    }
	}
	Real N_macroParticlesGlobal = 0.0;
	MPI_Reduce(&N_macroParticles,&N_macroParticlesGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
	const Species* species = reinterpret_cast<const Species*>(particleLists[s]->getSpecies());
	simClasses.logger << "\t N(" << species->name << ") = " << real2str(N_macroParticlesGlobal,15) << " = " << real2str(N_macroParticlesGlobal/1.0e9,15) << " x 10^9" << std::endl;
    }
   simClasses.logger << write;
}

// calculate particle log quantities
void logCalcParticle(Simulation& sim,SimulationClasses& simClasses,std::vector<LogDataParticle>& logDataParticle,const std::vector<ParticleListBase*>& particleLists,std::vector<Real>& cellRhoM)
{
   logDataParticle.clear();
   for(size_t s=0;s<particleLists.size();++s) {
      logDataParticle.push_back(LogDataParticle());
      logDataParticle[s].N_macroParticles = 0.0;
      logDataParticle[s].N_realParticles = 0.0;
      logDataParticle[s].sumVx = 0.0;
      logDataParticle[s].sumVy = 0.0;
      logDataParticle[s].sumVz = 0.0;
      logDataParticle[s].sumV = 0.0;
      logDataParticle[s].sumWV2 = 0.0;
      logDataParticle[s].maxVi = 0.0;
   }
   for(pargrid::CellID b=0; b<simClasses.pargrid.getNumberOfLocalCells(); ++b) {
      for(size_t s=0;s<particleLists.size();++s) {
	 pargrid::DataID speciesDataID = pargrid::INVALID_DATAID;
	 if(particleLists[s]->getParticles(speciesDataID) == false) { continue; }
	 pargrid::DataWrapper<Particle<Real> > wrapper = simClasses.pargrid.getUserDataDynamic<Particle<Real> >(speciesDataID);
	 Particle<Real>** particleList = wrapper.data();
	 Particle<Real>* particles = particleList[b];
	 pargrid::ArraySizetype N_particles = wrapper.size(b);
	 logDataParticle[s].N_macroParticles += N_particles;
	 const Species* species = reinterpret_cast<const Species*>(particleLists[s]->getSpecies());
	 for(size_t p=0; p<N_particles; ++p) {
	    logDataParticle[s].N_realParticles += particles[p].state[particle::WEIGHT];
	    logDataParticle[s].sumVx += particles[p].state[particle::WEIGHT]*particles[p].state[particle::VX];
	    logDataParticle[s].sumVy += particles[p].state[particle::WEIGHT]*particles[p].state[particle::VY];
	    logDataParticle[s].sumVz += particles[p].state[particle::WEIGHT]*particles[p].state[particle::VZ];
	    const Real v2 = sqr(particles[p].state[particle::VX]) + sqr(particles[p].state[particle::VY]) + sqr(particles[p].state[particle::VZ]);
	    const Real vtot = sqrt(v2);
	    logDataParticle[s].sumV += particles[p].state[particle::WEIGHT]*vtot;
	    logDataParticle[s].sumWV2 += particles[p].state[particle::WEIGHT]*v2;
	    if(vtot > logDataParticle[s].maxVi) { logDataParticle[s].maxVi = vtot; }
	    // determined alfven speed in each cell
	    const int i = static_cast<int>(floor(particles[p].state[particle::X]/Hybrid::dx));
	    const int j = static_cast<int>(floor(particles[p].state[particle::Y]/Hybrid::dx));
	    const int k = static_cast<int>(floor(particles[p].state[particle::Z]/Hybrid::dx));
	    const int n = (b*block::SIZE+block::index(i,j,k));
	    //const int n3 = n*3;
	    cellRhoM[n] += particles[p].state[particle::WEIGHT]*species->m;
	 }
      }
   }
   for(size_t i=0; i<cellRhoM.size(); ++i) { cellRhoM[i] /= Hybrid::dV; }
   // add background density and implement minimum densities in total mass density
   const Real minRhoMGlobal = Hybrid::minRhoQi*constants::MASS_PROTON/constants::CHARGE_ELEMENTARY; // assume global density minimum is protons
#ifdef USE_BACKGROUND_CHARGE_DENSITY
   Real* cellRhoQiBg = simClasses.pargrid.getUserDataStatic<Real>(Hybrid::dataCellRhoQiBgID);
#endif
#ifdef USE_OUTER_BOUNDARY_ZONE
   bool* outerBoundaryFlag = simClasses.pargrid.getUserDataStatic<bool>(Hybrid::dataOuterBoundaryFlagID);
   const Real minRhoMOuterBoundaryZone = Hybrid::outerBoundaryZone.minRhoQi*constants::MASS_PROTON/constants::CHARGE_ELEMENTARY; // assume outer boundary density minimum is protons
#endif
   for(pargrid::CellID b=0;b<simClasses.pargrid.getNumberOfLocalCells();++b) {
      for(int k=0;k<block::WIDTH_Z;++k) for(int j=0;j<block::WIDTH_Y;++j) for(int i=0;i<block::WIDTH_X;++i) {
	 const int n = (b*block::SIZE+block::index(i,j,k));
#ifdef USE_BACKGROUND_CHARGE_DENSITY
	 // add background density in total mass density
	 cellRhoM[n] += cellRhoQiBg[n]*constants::MASS_PROTON/constants::CHARGE_ELEMENTARY; // assume background density is protons
#endif
	 // implement minimum densities in total mass density
#ifdef USE_OUTER_BOUNDARY_ZONE
	 if(outerBoundaryFlag[n] == true) {
	    if(cellRhoM[n] < minRhoMOuterBoundaryZone) { cellRhoM[n] = minRhoMOuterBoundaryZone; }
	 }
	 else {
	    if(cellRhoM[n] < minRhoMGlobal) { cellRhoM[n] = minRhoMGlobal; }
	 }
#else
	 if(cellRhoM[n] < minRhoMGlobal) { cellRhoM[n] = minRhoMGlobal; }
#endif
      }
   }
}

// calculate field log quantities
void logCalcField(Simulation& sim,SimulationClasses& simClasses,LogDataField& logDataField,std::vector<Real>& cellRhoM)
{
   // synchronize faceB before using it below
   simClasses.pargrid.startNeighbourExchange(pargrid::DEFAULT_STENCIL,Hybrid::dataFaceBID);
   //profile::start("MPI waits",mpiWaitID);
   simClasses.pargrid.wait(pargrid::DEFAULT_STENCIL,Hybrid::dataFaceBID);
   //profile::stop();
   logDataField.N_cells = 0.0;
   // face magnetic field
   logDataField.sumBx = 0.0;
   logDataField.sumBy = 0.0;
   logDataField.sumBz = 0.0;
   logDataField.sumB = 0.0;
   logDataField.maxB = 0.0;
   logDataField.sumDivB = 0.0;
   logDataField.maxDivB = 0.0;
   logDataField.maxDivBPerB = 0.0;
   logDataField.sumB2 = 0.0;
   // cell ion current density
   logDataField.sumCellJix = 0.0;
   logDataField.sumCellJiy = 0.0;
   logDataField.sumCellJiz = 0.0;
   logDataField.sumCellJi = 0.0;
   logDataField.maxCellJi = 0.0;
   // cell electron pressure electric field
   logDataField.sumCellEpx = 0.0;
   logDataField.sumCellEpy = 0.0;
   logDataField.sumCellEpz = 0.0;
   logDataField.sumCellEp = 0.0;
   logDataField.maxCellEp = 0.0;
   logDataField.sumCellEp2 = 0.0;
   // electric field
   logDataField.sumNodeEx = 0.0;
   logDataField.sumNodeEy = 0.0;
   logDataField.sumNodeEz = 0.0;
   logDataField.sumNodeE = 0.0;
   logDataField.maxNodeE = 0.0;
   logDataField.sumNodeE2 = 0.0;
   logDataField.sumCellE2 = 0.0;
   // spatial, temporal and velocity scales
   logDataField.minInerLengthElectron = std::numeric_limits<Real>::max();
   logDataField.maxInerLengthElectron = std::numeric_limits<Real>::min();
   logDataField.minInerLengthProton = std::numeric_limits<Real>::max();
   logDataField.maxInerLengthProton = std::numeric_limits<Real>::min();
   logDataField.minTLarmor = std::numeric_limits<Real>::max();
   logDataField.maxVAlfven = 0.0;
   logDataField.maxUe = 0.0;
   Real* faceB = simClasses.pargrid.getUserDataStatic<Real>(Hybrid::dataFaceBID);
   Real* cellRhoQi = simClasses.pargrid.getUserDataStatic<Real>(Hybrid::dataCellRhoQiID);
   Real* cellJi = simClasses.pargrid.getUserDataStatic<Real>(Hybrid::dataCellJiID);
   Real* cellUe = simClasses.pargrid.getUserDataStatic<Real>(Hybrid::dataCellUeID);
   Real* cellEp = simClasses.pargrid.getUserDataStatic<Real>(Hybrid::dataCellEpID);
   Real* nodeE  = simClasses.pargrid.getUserDataStatic<Real>(Hybrid::dataNodeEID);
   bool* innerFlag = simClasses.pargrid.getUserDataStatic<bool>(Hybrid::dataInnerFlagFieldID);
   for(pargrid::CellID b=0; b<simClasses.pargrid.getNumberOfLocalCells(); ++b) {
      if(simClasses.pargrid.getNeighbourFlags(b) != pargrid::ALL_NEIGHBOURS_EXIST) { continue; }
      const size_t vectorDim = 3;
      const size_t s = (block::WIDTH_X+2)*(block::WIDTH_Y+2)*(block::WIDTH_Z+2);
      Real aa[s*vectorDim]; // temp array
      fetchData(faceB,aa,simClasses,b,vectorDim);
      for(int k=0; k<block::WIDTH_Z; ++k) for(int j=0; j<block::WIDTH_Y; ++j) for(int i=0; i<block::WIDTH_X; ++i) {
	 const int n = (b*block::SIZE+block::index(i,j,k));
	 const int n3 = n*3;
	 // do not include cells inside the inner boundary
	 if(Hybrid::includeInnerCellsInFieldLog == false) {
	    if(innerFlag[n] == true) { continue; }
	 }
	 // divergence of B in a cell from face magnetic field
	 Real divB = ((aa[(block::arrayIndex(i+1,j+1,k+1))*vectorDim+0] - aa[(block::arrayIndex(i+0,j+1,k+1))*vectorDim+0])
		    + (aa[(block::arrayIndex(i+1,j+1,k+1))*vectorDim+1] - aa[(block::arrayIndex(i+1,j+0,k+1))*vectorDim+1])
		    + (aa[(block::arrayIndex(i+1,j+1,k+1))*vectorDim+2] - aa[(block::arrayIndex(i+1,j+1,k+0))*vectorDim+2]))/Hybrid::dx;
	 divB = fabs(divB);
	 // B in a cell as an average face magnetic field
	 const Real Bx = 0.5*(aa[(block::arrayIndex(i+1,j+1,k+1))*vectorDim+0] + aa[(block::arrayIndex(i+0,j+1,k+1))*vectorDim+0]);
	 const Real By = 0.5*(aa[(block::arrayIndex(i+1,j+1,k+1))*vectorDim+1] + aa[(block::arrayIndex(i+1,j+0,k+1))*vectorDim+1]);
	 const Real Bz = 0.5*(aa[(block::arrayIndex(i+1,j+1,k+1))*vectorDim+2] + aa[(block::arrayIndex(i+1,j+1,k+0))*vectorDim+2]);
	 const Real B2 = sqr(Bx) + sqr(By) + sqr(Bz);
	 const Real Btot = sqrt(B2);
	 Real divBPerB = 0.0;
	 if(Btot > 0.0) { divBPerB = divB/Btot; }
	 logDataField.sumBx += Bx;
	 logDataField.sumBy += By;
	 logDataField.sumBz += Bz;
	 logDataField.sumB += Btot;
	 logDataField.sumDivB += divB;
	 logDataField.sumB2 += B2;
	 if(Btot > logDataField.maxB) { logDataField.maxB = Btot; }
	 if(fabs(divB) > fabs(logDataField.maxDivB)) { logDataField.maxDivB = divB; }
	 if(fabs(divBPerB) > fabs(logDataField.maxDivBPerB)) { logDataField.maxDivBPerB = divBPerB; }

	 // cell ion current density
	 const Real cellJitot = sqrt( sqr(cellJi[n3+0]) + sqr(cellJi[n3+1]) + sqr(cellJi[n3+2]) );
	 logDataField.sumCellJix += cellJi[n3+0];
	 logDataField.sumCellJiy += cellJi[n3+1];
	 logDataField.sumCellJiz += cellJi[n3+2];
	 logDataField.sumCellJi += cellJitot;
	 if(cellJitot > logDataField.maxCellJi) { logDataField.maxCellJi = cellJitot; }

	 // cell electron pressure electric field
	 const Real cellEp2 = sqr(cellEp[n3+0]) + sqr(cellEp[n3+1]) + sqr(cellEp[n3+2]);
	 const Real cellEptot = sqrt(cellEp2);
	 logDataField.sumCellEpx += cellEp[n3+0];
	 logDataField.sumCellEpy += cellEp[n3+1];
	 logDataField.sumCellEpz += cellEp[n3+2];
	 logDataField.sumCellEp += cellEptot;
	 logDataField.sumCellEp2 += cellEp2;
	 if(cellEptot > logDataField.maxCellEp) { logDataField.maxCellEp = cellEptot; }

	 // node electric field -Ue x B + eta*J
	 const Real nodeE2 = sqr(nodeE[n3+0]) + sqr(nodeE[n3+1]) + sqr(nodeE[n3+2]);
	 const Real nodeEtot = sqrt(nodeE2);
	 logDataField.sumNodeEx += nodeE[n3+0];
	 logDataField.sumNodeEy += nodeE[n3+1];
	 logDataField.sumNodeEz += nodeE[n3+2];
	 logDataField.sumNodeE += nodeEtot;
	 logDataField.sumNodeE2 += nodeE2;
	 if(nodeEtot > logDataField.maxNodeE) { logDataField.maxNodeE = nodeEtot; }

	 // cell electric field -Ue x B
	 const Real cellUex = cellUe[n3+0];
	 const Real cellUey = cellUe[n3+1];
	 const Real cellUez = cellUe[n3+2];
	 const Real cellEx = -(cellUey*Bz - cellUez*By);
	 const Real cellEy = -(cellUez*Bx - cellUex*Bz);
	 const Real cellEz = -(cellUex*By - cellUey*Bx);
	 const Real cellE2 = sqr(cellEx) + sqr(cellEy) + sqr(cellEz);
	 logDataField.sumCellE2 += cellE2;

	 // spatial, temporal and velocity scales
	 const Real cellUetot = sqrt(sqr(cellUex) + sqr(cellUey) + sqr(cellUez));
	 if(cellRhoQi[n] > 0.0) {
	    // electron inertial length: sqrt(m_e/(mu0*q_e^2*n_e)) = sqrt(m_e/(mu0*q_e^2*rho_q/q_e)) =  = sqrt(m_e/(mu0*q_e*rho_q))
	    const Real de = sqrt(constants::MASS_ELECTRON/(constants::PERMEABILITY*constants::CHARGE_ELEMENTARY*cellRhoQi[n]));
	    if(de > logDataField.maxInerLengthElectron) { logDataField.maxInerLengthElectron = de; }
	    if(de < logDataField.minInerLengthElectron) { logDataField.minInerLengthElectron = de; }
	 }
	 if(Btot > 0.0) {
	    // Minimum ion (=proton) Larmor period
	    const Real tL = 2.0*M_PI*constants::MASS_PROTON/(constants::CHARGE_ELEMENTARY*Btot);
	    if(tL < logDataField.minTLarmor) { logDataField.minTLarmor = tL; }
	 }
	 if(cellRhoM[n] > 0.0) {
	    // ion inertial length assuming that all particles in total mass density are protons: sqrt(m_p/(mu0*q_e^2*n_i)) = sqrt(m_p/(mu0*q_e^2*(rho_m/m_p))) = sqrt(m_p^2/(mu0*q_e^2*rho_m))
	    const Real di = sqrt(sqr(constants::MASS_PROTON)/(constants::PERMEABILITY*sqr(constants::CHARGE_ELEMENTARY)*cellRhoM[n]));
	    const Real vA = Btot/sqrt(constants::PERMEABILITY*cellRhoM[n]);
	    if(di > logDataField.maxInerLengthProton) { logDataField.maxInerLengthProton = di; }
	    if(di < logDataField.minInerLengthProton) { logDataField.minInerLengthProton = di; }
	    if(vA > logDataField.maxVAlfven) { logDataField.maxVAlfven = vA; }
	 }
	 if(cellUetot > logDataField.maxUe) { logDataField.maxUe = cellUetot; }

	 // update field log cell counter
	 logDataField.N_cells++;
      }
   }
}

// write particle population and field logs
bool logWriteParticleField(Simulation& sim,SimulationClasses& simClasses,const std::vector<ParticleListBase*>& particleLists) {
   bool success = true;
   static int profWriteLogsID = -1;
   //if(getInitialized() == false) { return false; }
   profile::start("logWriteParticleField",profWriteLogsID);
   std::vector<LogDataParticle> logDataParticle;
   LogDataField logDataField;
   // total mass density for Alfven speed
   std::vector<Real> cellRhoM;
   for(pargrid::CellID b=0; b<simClasses.pargrid.getNumberOfLocalCells(); ++b) for(int k=0; k<block::WIDTH_Z; ++k) for(int j=0; j<block::WIDTH_Y; ++j) for(int i=0; i<block::WIDTH_X; ++i) {
      cellRhoM.push_back(0.0);
   }
   logCalcParticle(sim,simClasses,logDataParticle,particleLists,cellRhoM);
   logCalcField(sim,simClasses,logDataField,cellRhoM);
   if(sim.mpiRank==sim.MASTER_RANK) {
      for(size_t i=0;i<Hybrid::logParticle.size();++i) {
	 (*Hybrid::logParticle[i]) << sim.t << " ";
      }
      Hybrid::logField << sim.t << " ";
   }
   Real maxViAllPopulations = 0.0;
   const Real Dt = sim.t - Hybrid::logCounterTimeStart;
   // go thru populations and sum particle counters
   for(size_t s=0;s<particleLists.size();++s) {
      Real N_macroParticlesThisProcess = logDataParticle[s].N_macroParticles;
      Real N_macroParticlesGlobal = 0.0;
      Real N_realParticlesThisProcess = logDataParticle[s].N_realParticles;
      Real N_realParticlesGlobal = 0.0;
      Real sumVxThisProcess = logDataParticle[s].sumVx;
      Real sumVxGlobal = 0.0;
      Real sumVyThisProcess = logDataParticle[s].sumVy;
      Real sumVyGlobal = 0.0;
      Real sumVzThisProcess = logDataParticle[s].sumVz;
      Real sumVzGlobal = 0.0;
      Real sumVThisProcess = logDataParticle[s].sumV;
      Real sumVGlobal = 0.0;
      Real sumWV2ThisProcess = logDataParticle[s].sumWV2;
      Real sumWV2Global = 0.0;
      Real maxViThisProcess = logDataParticle[s].maxVi;
      Real maxViGlobal = 0.0;
      Real sumEscapeThisProcess = Hybrid::logCounterParticleEscape[s];
      Real sumEscapeGlobal = 0.0;
      Real sumImpactThisProcess = Hybrid::logCounterParticleImpact[s];
      Real sumImpactGlobal = 0.0;
      Real sumInjectThisProcess = Hybrid::logCounterParticleInject[s];
      Real sumInjectGlobal = 0.0;
      Real sumInjectMacroThisProcess = Hybrid::logCounterParticleInjectMacroparticles[s];
      Real sumInjectMacroGlobal = 0.0;
      Real sumEscapeKineticEnergyThisProcess = Hybrid::logCounterParticleEscapeKineticEnergy[s];
      Real sumEscapeKineticEnergyGlobal = 0.0;
      Real sumImpactKineticEnergyThisProcess = Hybrid::logCounterParticleImpactKineticEnergy[s];
      Real sumImpactKineticEnergyGlobal = 0.0;
      Real sumInjectKineticEnergyThisProcess = Hybrid::logCounterParticleInjectKineticEnergy[s];
      Real sumInjectKineticEnergyGlobal = 0.0;
      Real sumMaxViThisProcess = Hybrid::logCounterParticleMaxVi[s];
      Real sumMaxViGlobal = 0.0;
      MPI_Reduce(&N_macroParticlesThisProcess,&N_macroParticlesGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
      MPI_Reduce(&N_realParticlesThisProcess,&N_realParticlesGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
      MPI_Reduce(&sumVxThisProcess,&sumVxGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
      MPI_Reduce(&sumVyThisProcess,&sumVyGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
      MPI_Reduce(&sumVzThisProcess,&sumVzGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
      MPI_Reduce(&sumVThisProcess,&sumVGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
      MPI_Reduce(&sumWV2ThisProcess,&sumWV2Global,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
      MPI_Reduce(&maxViThisProcess,&maxViGlobal,1,MPI_Type<Real>(),MPI_MAX,sim.MASTER_RANK,sim.comm);
      MPI_Reduce(&sumEscapeThisProcess,&sumEscapeGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
      MPI_Reduce(&sumImpactThisProcess,&sumImpactGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
      MPI_Reduce(&sumInjectThisProcess,&sumInjectGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
      MPI_Reduce(&sumInjectMacroThisProcess,&sumInjectMacroGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
      MPI_Reduce(&sumEscapeKineticEnergyThisProcess,&sumEscapeKineticEnergyGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
      MPI_Reduce(&sumImpactKineticEnergyThisProcess,&sumImpactKineticEnergyGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
      MPI_Reduce(&sumInjectKineticEnergyThisProcess,&sumInjectKineticEnergyGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
      MPI_Reduce(&sumMaxViThisProcess,&sumMaxViGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
      if(sim.mpiRank==sim.MASTER_RANK) {
	 (*Hybrid::logParticle[s]) << N_realParticlesGlobal << " " << N_macroParticlesGlobal << " ";
	 if(N_realParticlesGlobal > 0.0) {
	    (*Hybrid::logParticle[s])
	      << sumVxGlobal/N_realParticlesGlobal << " "
	      << sumVyGlobal/N_realParticlesGlobal << " "
	      << sumVzGlobal/N_realParticlesGlobal << " "
	      << sumVGlobal/N_realParticlesGlobal  << " ";
	 }
	 else {
	    (*Hybrid::logParticle[s]) << 0.0 << " " << 0.0 << " " << 0.0 << " " << 0.0 << " ";
	 }
	 const Species* species = reinterpret_cast<const Species*>(particleLists[s]->getSpecies());
	 (*Hybrid::logParticle[s]) << 0.5*species->m*sumWV2Global << " " << maxViGlobal << " ";
	 if(maxViAllPopulations < maxViGlobal) { maxViAllPopulations = maxViGlobal; }
	 if(Dt > 0) {
	    (*Hybrid::logParticle[s])
	      << sumEscapeGlobal/Dt << " "
	      << sumImpactGlobal/Dt << " "
	      << sumInjectGlobal/Dt << " "
	      << sumInjectMacroGlobal/Dt*sim.dt << " "
	      << 0.5*species->m*sumEscapeKineticEnergyGlobal/Dt << " "
	      << 0.5*species->m*sumImpactKineticEnergyGlobal/Dt << " "
	      << 0.5*species->m*sumInjectKineticEnergyGlobal/Dt << " "
	      << sumMaxViGlobal/Dt*sim.dt << " ";
	 }
	 else {
	    (*Hybrid::logParticle[s]) << 0.0 << " " << 0.0 << " " << 0.0 << " " << 0.0 << " " << 0.0 << " " << 0.0 << " " << 0.0 << " " << 0.0 << " ";
	 }
      }
   }

   // sum global field counters
   Real sumMaxCellUeThisProcess = Hybrid::logCounterFieldMaxCellUe;
   Real sumMaxCellUeGlobal = 0.0;
   Real sumMaxNodeUeThisProcess = Hybrid::logCounterFieldMaxNodeUe;
   Real sumMaxNodeUeGlobal = 0.0;
   Real sumMaxVwThisProcess = Hybrid::logCounterFieldMaxNodeUe;
   Real sumMaxVwGlobal = 0.0;
   Real sumMaxEThisProcess = Hybrid::logCounterFieldMaxE;
   Real sumMaxEGlobal = 0.0;
   Real sumMinCellRhoQiThisProcess = Hybrid::logCounterFieldMinCellRhoQi;
   Real sumMinCellRhoQiGlobal = 0.0;
   Real sumMinNodeRhoQiThisProcess = Hybrid::logCounterFieldMinNodeRhoQi;
   Real sumMinNodeRhoQiGlobal = 0.0;
   MPI_Reduce(&sumMaxCellUeThisProcess,&sumMaxCellUeGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
   MPI_Reduce(&sumMaxNodeUeThisProcess,&sumMaxNodeUeGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
   MPI_Reduce(&sumMaxVwThisProcess,&sumMaxVwGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
   MPI_Reduce(&sumMaxEThisProcess,&sumMaxEGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
   MPI_Reduce(&sumMinCellRhoQiThisProcess,&sumMinCellRhoQiGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
   MPI_Reduce(&sumMinNodeRhoQiThisProcess,&sumMinNodeRhoQiGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
   
   // field log values from logCalcField
   Real N_cellsThisProcess = logDataField.N_cells;
   Real N_cellsGlobal = 0.0;
   MPI_Reduce(&N_cellsThisProcess,&N_cellsGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
   // face magnetic field
   Real sumBxThisProcess = logDataField.sumBx;
   Real sumBxGlobal = 0.0;
   Real sumByThisProcess = logDataField.sumBy;
   Real sumByGlobal = 0.0;
   Real sumBzThisProcess = logDataField.sumBz;
   Real sumBzGlobal = 0.0;
   Real sumBThisProcess = logDataField.sumB;
   Real sumBGlobal = 0.0;
   Real maxBThisProcess = logDataField.maxB;
   Real maxBGlobal = 0.0;
   Real sumDivBThisProcess = logDataField.sumDivB;
   Real sumDivBGlobal = 0.0;
   Real maxDivBThisProcess = logDataField.maxDivB;
   Real maxDivBGlobal = 0.0;
   Real maxDivPerBThisProcess = logDataField.maxDivBPerB;
   Real maxDivPerBGlobal = 0.0;
   Real sumB2ThisProcess = logDataField.sumB2;
   Real sumB2Global = 0.0;
   MPI_Reduce(&sumBxThisProcess,&sumBxGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
   MPI_Reduce(&sumByThisProcess,&sumByGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
   MPI_Reduce(&sumBzThisProcess,&sumBzGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
   MPI_Reduce(&sumBThisProcess,&sumBGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
   MPI_Reduce(&maxBThisProcess,&maxBGlobal,1,MPI_Type<Real>(),MPI_MAX,sim.MASTER_RANK,sim.comm);
   MPI_Reduce(&sumDivBThisProcess,&sumDivBGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
   MPI_Reduce(&maxDivBThisProcess,&maxDivBGlobal,1,MPI_Type<Real>(),MPI_MAX,sim.MASTER_RANK,sim.comm);
   MPI_Reduce(&maxDivPerBThisProcess,&maxDivPerBGlobal,1,MPI_Type<Real>(),MPI_MAX,sim.MASTER_RANK,sim.comm);
   MPI_Reduce(&sumB2ThisProcess,&sumB2Global,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);

   // cell ion current density
   Real sumCellJixThisProcess = logDataField.sumCellJix;
   Real sumCellJixGlobal = 0.0;
   Real sumCellJiyThisProcess = logDataField.sumCellJiy;
   Real sumCellJiyGlobal = 0.0;
   Real sumCellJizThisProcess = logDataField.sumCellJiz;
   Real sumCellJizGlobal = 0.0;
   Real sumCellJiThisProcess = logDataField.sumCellJi;
   Real sumCellJiGlobal = 0.0;
   Real maxCellJiThisProcess = logDataField.maxCellJi;
   Real maxCellJiGlobal = 0.0;
   MPI_Reduce(&sumCellJixThisProcess,&sumCellJixGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
   MPI_Reduce(&sumCellJiyThisProcess,&sumCellJiyGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
   MPI_Reduce(&sumCellJizThisProcess,&sumCellJizGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
   MPI_Reduce(&sumCellJiThisProcess,&sumCellJiGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
   MPI_Reduce(&maxCellJiThisProcess,&maxCellJiGlobal,1,MPI_Type<Real>(),MPI_MAX,sim.MASTER_RANK,sim.comm);

   // cell electron pressure electric field
   Real sumCellEpxThisProcess = logDataField.sumCellEpx;
   Real sumCellEpxGlobal = 0.0;
   Real sumCellEpyThisProcess = logDataField.sumCellEpy;
   Real sumCellEpyGlobal = 0.0;
   Real sumCellEpzThisProcess = logDataField.sumCellEpz;
   Real sumCellEpzGlobal = 0.0;
   Real sumCellEpThisProcess = logDataField.sumCellEp;
   Real sumCellEpGlobal = 0.0;
   Real maxCellEpThisProcess = logDataField.maxCellEp;
   Real maxCellEpGlobal = 0.0;
   Real sumCellEp2ThisProcess = logDataField.sumCellEp2;
   Real sumCellEp2Global = 0.0;
   MPI_Reduce(&sumCellEpxThisProcess,&sumCellEpxGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
   MPI_Reduce(&sumCellEpyThisProcess,&sumCellEpyGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
   MPI_Reduce(&sumCellEpzThisProcess,&sumCellEpzGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
   MPI_Reduce(&sumCellEpThisProcess,&sumCellEpGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
   MPI_Reduce(&maxCellEpThisProcess,&maxCellEpGlobal,1,MPI_Type<Real>(),MPI_MAX,sim.MASTER_RANK,sim.comm);
   MPI_Reduce(&sumCellEp2ThisProcess,&sumCellEp2Global,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);

   // node electric field
   Real sumNodeExThisProcess = logDataField.sumNodeEx;
   Real sumNodeExGlobal = 0.0;
   Real sumNodeEyThisProcess = logDataField.sumNodeEy;
   Real sumNodeEyGlobal = 0.0;
   Real sumNodeEzThisProcess = logDataField.sumNodeEz;
   Real sumNodeEzGlobal = 0.0;
   Real sumNodeEThisProcess = logDataField.sumNodeE;
   Real sumNodeEGlobal = 0.0;
   Real maxNodeEThisProcess = logDataField.maxNodeE;
   Real maxNodeEGlobal = 0.0;
   Real sumNodeE2ThisProcess = logDataField.sumNodeE2;
   Real sumNodeE2Global = 0.0;
   Real sumCellE2ThisProcess = logDataField.sumCellE2;
   Real sumCellE2Global = 0.0;
   MPI_Reduce(&sumNodeExThisProcess,&sumNodeExGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
   MPI_Reduce(&sumNodeEyThisProcess,&sumNodeEyGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
   MPI_Reduce(&sumNodeEzThisProcess,&sumNodeEzGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
   MPI_Reduce(&sumNodeEThisProcess,&sumNodeEGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
   MPI_Reduce(&maxNodeEThisProcess,&maxNodeEGlobal,1,MPI_Type<Real>(),MPI_MAX,sim.MASTER_RANK,sim.comm);
   MPI_Reduce(&sumNodeE2ThisProcess,&sumNodeE2Global,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
   MPI_Reduce(&sumCellE2ThisProcess,&sumCellE2Global,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);

   // spatial, temporal and velocity scales
   Real minInerLengthElectronThisProcess = logDataField.minInerLengthElectron;
   Real minInerLengthElectronGlobal = 0.0;
   Real maxInerLengthElectronThisProcess = logDataField.maxInerLengthElectron;
   Real maxInerLengthElectronGlobal = 0.0;
   Real minInerLengthProtonThisProcess = logDataField.minInerLengthProton;
   Real minInerLengthProtonGlobal = 0.0;
   Real maxInerLengthProtonThisProcess = logDataField.maxInerLengthProton;
   Real maxInerLengthProtonGlobal = 0.0;
   Real minTLarmorThisProcess = logDataField.minTLarmor;
   Real minTLarmorGlobal = 0.0;
   Real maxVAlfvenThisProcess = logDataField.maxVAlfven;
   Real maxVAlfvenGlobal = 0.0;
   Real maxUeThisProcess = logDataField.maxUe;
   Real maxUeGlobal = 0.0;
   MPI_Reduce(&minInerLengthElectronThisProcess,&minInerLengthElectronGlobal,1,MPI_Type<Real>(),MPI_MIN,sim.MASTER_RANK,sim.comm);
   MPI_Reduce(&maxInerLengthElectronThisProcess,&maxInerLengthElectronGlobal,1,MPI_Type<Real>(),MPI_MAX,sim.MASTER_RANK,sim.comm);
   MPI_Reduce(&minInerLengthProtonThisProcess,&minInerLengthProtonGlobal,1,MPI_Type<Real>(),MPI_MIN,sim.MASTER_RANK,sim.comm);
   MPI_Reduce(&maxInerLengthProtonThisProcess,&maxInerLengthProtonGlobal,1,MPI_Type<Real>(),MPI_MAX,sim.MASTER_RANK,sim.comm);
   MPI_Reduce(&minTLarmorThisProcess,&minTLarmorGlobal,1,MPI_Type<Real>(),MPI_MIN,sim.MASTER_RANK,sim.comm);
   MPI_Reduce(&maxVAlfvenThisProcess,&maxVAlfvenGlobal,1,MPI_Type<Real>(),MPI_MAX,sim.MASTER_RANK,sim.comm);
   MPI_Reduce(&maxUeThisProcess,&maxUeGlobal,1,MPI_Type<Real>(),MPI_MAX,sim.MASTER_RANK,sim.comm);

   // write field log
   if(sim.mpiRank==sim.MASTER_RANK) {
      // face magnetic field
      if(N_cellsGlobal > 0) {
	 Hybrid::logField
	   << sumBxGlobal/N_cellsGlobal << " "
	   << sumByGlobal/N_cellsGlobal << " "
	   << sumBzGlobal/N_cellsGlobal << " "
	   << sumBGlobal/N_cellsGlobal << " ";
      }
      else {
	 Hybrid::logField << 0.0 << " " << 0.0 << " " << 0.0 << " " << 0.0 << " ";
      }
      Hybrid::logField << maxBGlobal << " ";
      if(N_cellsGlobal > 0) { Hybrid::logField << sumDivBGlobal/N_cellsGlobal << " "; }
      else { Hybrid::logField << 0.0 << " "; }
      Hybrid::logField << maxDivBGlobal << " " << Hybrid::dx*maxDivPerBGlobal << " " << sumB2Global*Hybrid::dV/(2.0*constants::PERMEABILITY) << " ";

      // cell ion current density
      if(N_cellsGlobal > 0) {
	 Hybrid::logField
	   << sumCellJixGlobal/N_cellsGlobal << " "
	   << sumCellJiyGlobal/N_cellsGlobal << " "
	   << sumCellJizGlobal/N_cellsGlobal << " "
	   << sumCellJiGlobal/N_cellsGlobal << " ";
      }
      else {
	 Hybrid::logField << 0.0 << " " << 0.0 << " " << 0.0 << " " << 0.0 << " ";
      }
      Hybrid::logField << maxCellJiGlobal << " ";

      // cell electron pressure electric field
      if(N_cellsGlobal > 0) {
	 Hybrid::logField
	   << sumCellEpxGlobal/N_cellsGlobal << " "
	   << sumCellEpyGlobal/N_cellsGlobal << " "
	   << sumCellEpzGlobal/N_cellsGlobal << " "
	   << sumCellEpGlobal/N_cellsGlobal << " ";
      }
      else {
	 Hybrid::logField << 0.0 << " " << 0.0 << " " << 0.0 << " " << 0.0 << " ";
      }
      Hybrid::logField << maxCellEpGlobal << " " << sumCellEp2Global*Hybrid::dV*constants::PERMITTIVITY/(2.0) << " ";

      // node electric field
      if(N_cellsGlobal > 0) {
	 Hybrid::logField
	   << sumNodeExGlobal/N_cellsGlobal << " "
	   << sumNodeEyGlobal/N_cellsGlobal << " "
	   << sumNodeEzGlobal/N_cellsGlobal << " "
	   << sumNodeEGlobal/N_cellsGlobal << " ";
      }
      else {
	 Hybrid::logField << 0.0 << " " << 0.0 << " " << 0.0 << " " << 0.0 << " ";
      }
      Hybrid::logField << " " << maxNodeEGlobal << " " << sumNodeE2Global*Hybrid::dV*constants::PERMITTIVITY/(2.0)  << " " << sumCellE2Global*Hybrid::dV*constants::PERMITTIVITY/(2.0)  << " ";

      // global field counters
      if(Dt > 0) {
	 Hybrid::logField
	   << sumMaxCellUeGlobal/Dt*sim.dt << " "
	   << sumMaxNodeUeGlobal/Dt*sim.dt << " "
	   << sumMaxVwGlobal/Dt*sim.dt << " "
	   << sumMaxEGlobal/Dt*sim.dt << " "
	   << sumMinCellRhoQiGlobal/Dt*sim.dt << " "
	   << sumMinNodeRhoQiGlobal/Dt*sim.dt << " ";
      }
      else {
	 Hybrid::logField << 0.0 << " " << 0.0 << " " << 0.0 << " " << 0.0 << " " << 0.0 << " " << 0.0 << " ";
      }

      // spatial, temporal and velocity scales
      Hybrid::logField
	<< minInerLengthElectronGlobal << " "
	<< maxInerLengthElectronGlobal << " "
	<< minInerLengthProtonGlobal << " "
	<< maxInerLengthProtonGlobal << " "
	<< minTLarmorGlobal << " "
	<< maxVAlfvenGlobal << " "
	<< maxUeGlobal << " ";

      // line endings particle logs
      for(size_t i=0;i<Hybrid::logParticle.size();++i) {
	 (*Hybrid::logParticle[i]) << std::endl;
      }
      // line ending field log
      Hybrid::logField << std::endl;

      // warnings
      const Real tL_min_dt = minTLarmorGlobal/sim.dt;
      const Real maxVi_dxdt = maxViAllPopulations/(Hybrid::dx/sim.dt);
      const Real maxVA_dxdt = maxVAlfvenGlobal/(Hybrid::dx/sim.dt);
      const Real maxUe_dxdt = maxUeGlobal/(Hybrid::dx/sim.dt);

      // write if on a save step or save step already happened after previous entry
      if(Hybrid::writeMainLogEntriesAfterSaveStep == true && sim.mpiRank == sim.MASTER_RANK) {
	 simClasses.logger
	   << "(RHYBRID) Diagnostics (time step = " << sim.timestep << ", time = " << sim.t << "):" << std::endl
	   << "\t Max. |Vion|             : Vi_max  = " << maxViAllPopulations/1e3 << " km/s = " << maxVi_dxdt << " dx/dt" << std::endl
	   << "\t Max. |Ue|               : Ue_max  = " << maxUeGlobal/1e3 << " km/s = " << maxUe_dxdt << " dx/dt" << std::endl
	   << "\t Max. |VAlfven|          : Va_max  = " << maxVAlfvenGlobal/1e3 << " km/s = " << maxVA_dxdt << " dx/dt" << std::endl
	   << "\t Max. |faceB|            : B_max   = " << maxBGlobal/1e-9 << " nT" << std::endl
	   << "\t Min. ion Larmor period  : tL_min  = " << minTLarmorGlobal << " s = " << minTLarmorGlobal/sim.dt << " dt" << std::endl
	   << "\t Max. |cellJi|           : Ji_max  = " << maxCellJiGlobal << " A/m^2" << std::endl
	   << "\t Max. |cellEp|           : Ep_max  = " << maxCellEpGlobal/1e-3 << " mV/m" << std::endl
	   << "\t Max. |nodeE|            : E_max   = " << maxNodeEGlobal/1e-3 << " mV/m" << std::endl
	   << "\t Min. e- inertial length : de_min  = " << minInerLengthElectronGlobal/1e3 << " km = " << minInerLengthElectronGlobal/Hybrid::dx << " dx" << std::endl
	   << "\t Max. e- inertial length : de_max  = " << maxInerLengthElectronGlobal/1e3 << " km = " << maxInerLengthElectronGlobal/Hybrid::dx << " dx" << std::endl
	   << "\t Min. H+ inertial length : di_min  = " << minInerLengthProtonGlobal/1e3 << " km = " << minInerLengthProtonGlobal/Hybrid::dx << " dx" << std::endl
	   << "\t Max. H+ inertial length : di_max  = " << maxInerLengthProtonGlobal/1e3 << " km = " << maxInerLengthProtonGlobal/Hybrid::dx << " dx" << std::endl << write;
	 Hybrid::writeMainLogEntriesAfterSaveStep = false;
      }
      if(sim.mpiRank == sim.MASTER_RANK) {
	 if(tL_min_dt < 10)   { simClasses.logger << "(RHYBRID) WARNING: Minimum Larmor period: tL_min/dt < 10 ("      << tL_min_dt  << ") (time step = " << sim.timestep << ", time = " << sim.t << ")" << std::endl << write; }
	 if(maxVi_dxdt > 0.9) { simClasses.logger << "(RHYBRID) WARNING: Maximum ion speed: Vi_max/(dx/dt) > 0.9 ("    << maxVi_dxdt << ") (time step = " << sim.timestep << ", time = " << sim.t << ")" << std::endl << write; }
	 if(maxUe_dxdt > 0.9) { simClasses.logger << "(RHYBRID) WARNING: Maximum ion speed: Ue_max/(dx/dt) > 0.9 ("    << maxUe_dxdt << ") (time step = " << sim.timestep << ", time = " << sim.t << ")" << std::endl << write; }
	 if(maxVA_dxdt > 0.9) { simClasses.logger << "(RHYBRID) WARNING: Maximum Alfven speed: Va_max/(dx/dt) > 0.9 (" << maxVA_dxdt << ") (time step = " << sim.timestep << ", time = " << sim.t << ")" << std::endl << write; }
      }
      if(maxBGlobal > Hybrid::terminateLimitMaxB) {
         success = false;
	 if(sim.mpiRank == sim.MASTER_RANK) { simClasses.logger << "(RHYBRID) CONSTRAINT: maximum |B| for run termination reached (maxBGlobal = " << maxBGlobal/1e-9 << " nT) (time step = " << sim.timestep << ", time = " << sim.t << "), exiting." << std::endl << write; }
      }
   }

   // zero particle population counters
   for(size_t s=0;s<particleLists.size();++s) {
      Hybrid::logCounterParticleEscape[s] = 0.0;
      Hybrid::logCounterParticleImpact[s] = 0.0;
      Hybrid::logCounterParticleInject[s] = 0.0;
      Hybrid::logCounterParticleInjectMacroparticles[s] = 0.0;
      Hybrid::logCounterParticleEscapeKineticEnergy[s] = 0.0;
      Hybrid::logCounterParticleImpactKineticEnergy[s] = 0.0;
      Hybrid::logCounterParticleInjectKineticEnergy[s] = 0.0;
      Hybrid::logCounterParticleMaxVi[s] = 0.0;
   }

   // zero field counters
   Hybrid::logCounterFieldMaxCellUe = 0.0;
   Hybrid::logCounterFieldMaxNodeUe = 0.0;
   Hybrid::logCounterFieldMaxVw = 0.0;
   Hybrid::logCounterFieldMaxE = 0.0;
   Hybrid::logCounterFieldMinCellRhoQi = 0.0;
   Hybrid::logCounterFieldMinNodeRhoQi = 0.0;

   // reset counter start time
   Hybrid::logCounterTimeStart = sim.t;

   // silo time series (curves.silo)
   /*Real cnt=0;
   for(pargrid::CellID b=0; b<simClasses.pargrid.getNumberOfLocalCells(); ++b){ cnt+=1.0; }
   
   map<std::string,std::string> attribs;
   attribs["name"] = "blah";
   attribs["xlabel"] = "time";
   attribs["ylabel"] = "y laabeli";
   attribs["xunit"] = "s";
   attribs["yunit"] = "y unitti";   
   
   simClasses.vlsv.writeWithReduction("TIMESERIES",attribs,1,&cnt,MPI_SUM);*/

   profile::stop();
   return success;
}

}

#endif
