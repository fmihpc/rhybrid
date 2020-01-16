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

#include <climits>
#include <string>

#include "particle_injector.h"
#include "hybrid.h"
#include "neutral_profiles.h"

using namespace std;

// probabilistic rounding (used to round non-integer N_macroParticlesPerCellPerDt)
int probround(SimulationClasses& simClasses,Real x) {
   if(x <= 0) { return 0; }
   const int f = static_cast<int>(floor(x));
   const Real r = simClasses.random.uniform();
   if(r < x-f) { return f+1; }
   else { return f; }
}

// gaussian randomness
Real gaussrnd(SimulationClasses& simClasses)
{
   static Real saved;
   static bool is_saved = false;
   Real x,y,r2,fac,result;
   if (is_saved) {
      result = saved;
      is_saved = false;
   } else {
      do {
	 x = 2*simClasses.random.uniform() - 1;
	 y = 2*simClasses.random.uniform() - 1;
	 r2 = x*x + y*y;
      } while (r2 >= 1.0);
      // On average, this do loop is executed 4/pi = 1.27324 times
      fac = sqrt(-2.0*log(r2)/r2);
      result = x*fac;
      saved = y*fac;
      is_saved = true;
   }
   return result;
}

// deriv gaussian randomness
Real derivgaussrnd(Real x0,SimulationClasses& simClasses)
{
   Real x,majorant,pdf;
   const Real invsqrt2 = 1.0/sqrt(2.0);
   const Real sqrt_halfpi = sqrt(0.5*M_PI);
   const Real xm = 0.5*(x0 + sqrt(sqr(x0) + 4.0));
   const Real c = 1.0/(exp(-0.5*x0*x0) + x0*sqrt_halfpi*erfc(-x0*invsqrt2));
   const Real d = sqr(x0-xm);
   const Real cxm = c*xm;
restart:
   x = xm + gaussrnd(simClasses);
   if (x < 0) goto restart;
   majorant = cxm*exp(-0.5*(sqr(x-xm) + d));
   pdf = c*x*exp(-0.5*sqr(x-x0));
   if (simClasses.random.uniform()*majorant > pdf) goto restart;
   return x;
}

string radiusToString(Real R) {
   return
     to_string(R/1e3) + " km = " + 
     to_string(R/Hybrid::R_object) + " R_object = " + 
     to_string((R-Hybrid::R_object)/1e3) + " km + R_object";
}

// UNIFORM INJECTOR

InjectorUniform::InjectorUniform(): ParticleInjectorBase() {
   initialized = false;
   U = vth = n = w = 0.0;
   N_macroParticlesPerCell = -1.0;
}

InjectorUniform::~InjectorUniform() {finalize();}

bool InjectorUniform::finalize() {
   if(initialized == false) { return true; }
   initialized = false;
   return true;
}

bool InjectorUniform::inject(pargrid::DataID speciesDataID,unsigned int* N_particles) {
   if(initialized == false) { return initialized; }
   bool success = true;
   if(sim->timestep != 1) { return success; }
   pargrid::DataWrapper<Particle<Real> > wrapper = simClasses->pargrid.getUserDataDynamic<Particle<Real> >(speciesDataID);
   for(pargrid::CellID b=0;b<simClasses->pargrid.getNumberOfLocalCells();++b) {
      if(injectParticles(b,*species,N_particles,wrapper) == false) { success = false; }
   }
   return success;
}

bool InjectorUniform::injectParticles(pargrid::CellID blockID,const Species& species,unsigned int* N_particles,
				       pargrid::DataWrapper<Particle<Real> >& wrapper) {
/*#ifdef USE_DETECTORS
   const Real* crd = getBlockCoordinateArray(*sim,*simClasses);
   const size_t b3 = 3*blockID;
   const Real xBlock = crd[b3+0];
   const Real yBlock = crd[b3+1];
   const Real zBlock = crd[b3+2];
#endif*/
   vector<Real> xinj,yinj,zinj;
   for(int k=0;k<block::WIDTH_Z;++k) for(int j=0;j<block::WIDTH_Y;++j) for(int i=0;i<block::WIDTH_X;++i) {
      const Real xCell = (i+0.5)*Hybrid::dx;
      const Real yCell = (j+0.5)*Hybrid::dx;
      const Real zCell = (k+0.5)*Hybrid::dx;
      const int N_injectCell = probround(*simClasses,N_macroParticlesPerCell);
      if(N_injectCell <= 0) { continue; }
      for(int s = 0;s<N_injectCell;s++) {
	 const Real eps = 1.0e-2;
	 const Real x = xCell + (1.0-eps)*Hybrid::dx*(simClasses->random.uniform()-0.5);
	 const Real y = yCell + (1.0-eps)*Hybrid::dx*(simClasses->random.uniform()-0.5);
	 const Real z = zCell + (1.0-eps)*Hybrid::dx*(simClasses->random.uniform()-0.5);
	 xinj.push_back(x);
	 yinj.push_back(y);
	 zinj.push_back(z);
      }
   }
   // Make room for new particles:
   const int N_inject = static_cast<int>(xinj.size());
   const pargrid::ArraySizetype oldSize = wrapper.size()[blockID];
   N_particles[blockID] += N_inject;
   wrapper.resize(blockID,oldSize+N_inject);
   Particle<Real>* particles = wrapper.data()[blockID];
   size_t s = 0;
   for(size_t p=oldSize; p<oldSize+N_inject; ++p) {
      particles[p].state[particle::X] = xinj[s];
      particles[p].state[particle::Y] = yinj[s];
      particles[p].state[particle::Z] = zinj[s];
      particles[p].state[particle::VX] = -U + vth*gaussrnd(*simClasses);
      particles[p].state[particle::VY] = vth*gaussrnd(*simClasses);
      particles[p].state[particle::VZ] = vth*gaussrnd(*simClasses);
      particles[p].state[particle::WEIGHT] = w;
/*#ifdef USE_DETECTORS
      particles[p].state[particle::INI_CELLID] = simClasses->pargrid.getGlobalIDs()[blockID];
      particles[p].state[particle::INI_X] = xBlock + particles[p].state[particle::X];
      particles[p].state[particle::INI_Y] = yBlock + particles[p].state[particle::Y];
      particles[p].state[particle::INI_Z] = zBlock + particles[p].state[particle::Z];
      particles[p].state[particle::INI_VX] =  particles[p].state[particle::VX];
      particles[p].state[particle::INI_VY] =  particles[p].state[particle::VY];
      particles[p].state[particle::INI_VZ] =  particles[p].state[particle::VZ];
      particles[p].state[particle::INI_TIME] = sim->t;
#endif*/
      // inject counter
      Hybrid::particleCounterInject[species.popid-1] += w;
      Hybrid::particleCounterInjectMacroparticles[species.popid-1] += 1;
      ++s;
   }
   return true;
}

bool InjectorUniform::addConfigFileItems(ConfigReader& cr,const std::string& configRegionName) {
   cr.add(configRegionName+".speed","Bulk speed [m/s] (float).",(Real)0.0);
   cr.add(configRegionName+".density","Number density [m^-3] (float).",(Real)0.0);
   cr.add(configRegionName+".temperature","Temperature [K] (float).",(Real)0.0);
   cr.add(configRegionName+".macroparticles_per_cell","Number of macroparticles per cell [#] (Real).",(Real)-1.0);
   return true;
}

bool InjectorUniform::initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,
				   const std::string& configRegionName,const ParticleListBase* plist) {
   initialized = ParticleInjectorBase::initialize(sim,simClasses,cr,configRegionName,plist);
   this->species = reinterpret_cast<const Species*>(plist->getSpecies());
   Real T=0;
   cr.parse();
   cr.get(configRegionName+".speed",U);
   cr.get(configRegionName+".density",n);
   cr.get(configRegionName+".temperature",T);
   cr.get(configRegionName+".macroparticles_per_cell",N_macroParticlesPerCell);
   if(U <= 0) { U = 0.0; }
   if(N_macroParticlesPerCell > 0 && n > 0 && Hybrid::dx > 0) {
      w = n*Hybrid::dV/N_macroParticlesPerCell;
   }
   else {
      N_macroParticlesPerCell = n = w = 0.0;
   }
   if(T > 0) { vth = sqrt(constants::BOLTZMANN*T/species->m); }
   else { vth = 0.0; }
   simClasses.logger
     << "(" << species->name << ") speed         = " << U/1e3 << " km/s" << endl
     << "(" << species->name << ") density       = " << n/1e6 << " cm^{-3}" << endl
     << "(" << species->name << ") temperature   = " << T << " K = " << T/constants::EV_TO_KELVIN << " eV" << endl
     << "(" << species->name << ") thermal speed = " << vth/1e3 << " km/s" << endl
     << "(" << species->name << ") macroparticles per cell = " << N_macroParticlesPerCell << endl
     << "(" << species->name << ") macroparticle weight    = " << w << endl << write;
   particlePopulation pp;
   pp.w = w;
   pp.name = species->name;
   Hybrid::allPops.push_back(pp);
   return initialized;
}

// SOLAR WIND INJECTOR

InjectorSolarWind::InjectorSolarWind(): ParticleInjectorBase() { 
   initialized = false;
   U = vth = n = w = 0.0;
   N_macroParticlesPerCellPerDt = -1.0;
   N_macroParticlesPerCell = -1.0;
}

InjectorSolarWind::~InjectorSolarWind() {finalize();}

bool InjectorSolarWind::finalize() {
   if(initialized == false) { return true; }
   initialized = false;
   return true;
}

bool InjectorSolarWind::inject(pargrid::DataID speciesDataID,unsigned int* N_particles) {
   if(initialized == false) { return initialized; }
   bool success = true;
   if(sim->timestep <= 0) { return success; }
   
   pargrid::DataWrapper<Particle<Real> > wrapper = simClasses->pargrid.getUserDataDynamic<Particle<Real> >(speciesDataID);
   
   for(pargrid::CellID b=0; b<simClasses->pargrid.getNumberOfLocalCells(); ++b) {
      if( (simClasses->pargrid.getNeighbourFlags(b) & Hybrid::X_POS_EXISTS) == 0 &&
	  (simClasses->pargrid.getNeighbourFlags(b) & Hybrid::Y_POS_EXISTS) != 0 &&
	  (simClasses->pargrid.getNeighbourFlags(b) & Hybrid::Y_NEG_EXISTS) != 0 &&
	  (simClasses->pargrid.getNeighbourFlags(b) & Hybrid::Z_POS_EXISTS) != 0 &&
	  (simClasses->pargrid.getNeighbourFlags(b) & Hybrid::Z_NEG_EXISTS) != 0) {
	 if(injectParticles(b,*species,N_particles,wrapper) == false) { success = false; }
      }
   }
   return success;
}

bool InjectorSolarWind::injectParticles(pargrid::CellID blockID,const Species& species,unsigned int* N_particles,
				       pargrid::DataWrapper<Particle<Real> >& wrapper) {
/*#ifdef USE_DETECTORS
   const Real* crd = getBlockCoordinateArray(*sim,*simClasses);
   const size_t b3 = 3*blockID;
   const Real xBlock = crd[b3+0];
   const Real yBlock = crd[b3+1];
   const Real zBlock = crd[b3+2];
#endif*/
   Real blockSize[3];
   getBlockSize(*simClasses,*sim,blockID,blockSize);
   // probround
   const int N_inject = probround(*simClasses,N_macroParticlesPerCellPerDt);
   if(N_inject <= 0) { return true; }
   // Make room for new particles:
   const pargrid::ArraySizetype oldSize = wrapper.size()[blockID];
   N_particles[blockID] += N_inject;
   wrapper.resize(blockID,oldSize+N_inject);
   Particle<Real>* particles = wrapper.data()[blockID];
   for(size_t p=oldSize; p<oldSize+N_inject; ++p) {
      particles[p].state[particle::X] = 0;
      particles[p].state[particle::Y] = simClasses->random.uniform()*blockSize[1];
      particles[p].state[particle::Z] = simClasses->random.uniform()*blockSize[2];
      particles[p].state[particle::VX] = -vth*derivgaussrnd(U/vth,*simClasses);
      particles[p].state[particle::VY] = vth*gaussrnd(*simClasses);
      particles[p].state[particle::VZ] = vth*gaussrnd(*simClasses);
      particles[p].state[particle::WEIGHT] = w;
/*#ifdef USE_DETECTORS
      particles[p].state[particle::INI_CELLID] = simClasses->pargrid.getGlobalIDs()[blockID];
      particles[p].state[particle::INI_X] = xBlock + particles[p].state[particle::X];
      particles[p].state[particle::INI_Y] = yBlock + particles[p].state[particle::Y];
      particles[p].state[particle::INI_Z] = zBlock + particles[p].state[particle::Z];
      particles[p].state[particle::INI_VX] =  particles[p].state[particle::VX];
      particles[p].state[particle::INI_VY] =  particles[p].state[particle::VY];
      particles[p].state[particle::INI_VZ] =  particles[p].state[particle::VZ];
      particles[p].state[particle::INI_TIME] = sim->t;
#endif*/
      // inject counter
      Hybrid::particleCounterInject[species.popid-1] += w;
      Hybrid::particleCounterInjectMacroparticles[species.popid-1] += 1;
   }
   return true;
}

bool InjectorSolarWind::addConfigFileItems(ConfigReader& cr,const std::string& configRegionName) {
   cr.add(configRegionName+".speed","Bulk speed [m/s] (float).",(Real)0.0);
   cr.add(configRegionName+".density","Number density [m^-3] (float).",(Real)0.0);
   cr.add(configRegionName+".temperature","Temperature [K] (float).",(Real)0.0);
   cr.add(configRegionName+".macroparticles_per_cell","Number of macroparticles per cell in undisturbed solar wind [#] (Real).",(Real)-1.0);
   return true;
}

bool InjectorSolarWind::initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,
				   const std::string& configRegionName,const ParticleListBase* plist) {
   initialized = ParticleInjectorBase::initialize(sim,simClasses,cr,configRegionName,plist);
   this->species = reinterpret_cast<const Species*>(plist->getSpecies());
   
   Real T=0;
   cr.parse();
   cr.get(configRegionName+".speed",U);
   cr.get(configRegionName+".density",n);
   cr.get(configRegionName+".temperature",T);
   cr.get(configRegionName+".macroparticles_per_cell",N_macroParticlesPerCell);
   if(U <= 0) { U = 0.0; }
   if(N_macroParticlesPerCell > 0 && n > 0 && Hybrid::dx > 0) {
      N_macroParticlesPerCellPerDt = N_macroParticlesPerCell*fabs(U)*sim.dt/Hybrid::dx;
      w = n*Hybrid::dV/N_macroParticlesPerCell;
   }
   else {
      N_macroParticlesPerCellPerDt = N_macroParticlesPerCell = n = w = 0.0;
   }
   if(T > 0) { vth = sqrt(constants::BOLTZMANN*T/species->m); }
   else { vth = 0.0; }
   int N_yz_cells = (sim.y_blocks-2)*(sim.z_blocks-2)*block::WIDTH_Y*block::WIDTH_Z; // blocks != 1 do not work in hybrid
   simClasses.logger
     << "(" << species->name << ") speed         = " << U/1e3 << " km/s" << endl
     << "(" << species->name << ") density       = " << n/1e6 << " cm^{-3}" << endl
     << "(" << species->name << ") temperature   = " << T << " K = " << T/constants::EV_TO_KELVIN << " eV" << endl
     << "(" << species->name << ") thermal speed = " << vth/1e3 << " km/s" << endl
     << "(" << species->name << ") macroparticles per cell         = " << N_macroParticlesPerCell << endl
     << "(" << species->name << ") macroparticles per dt           = " << N_macroParticlesPerCellPerDt*N_yz_cells << endl
     << "(" << species->name << ") macroparticles per cell per dt  = " << N_macroParticlesPerCellPerDt << endl
     << "(" << species->name << ") macroparticle weight            = " << w << endl << write;
   static unsigned int swPopCnt = 0;
   swPopCnt++;
   if(swPopCnt == 1) {
      Hybrid::swMacroParticlesCellPerDt = N_macroParticlesPerCellPerDt*N_yz_cells/N_macroParticlesPerCell;
   }
   particlePopulation pp;
   pp.w = w;
   pp.name = species->name;
   Hybrid::allPops.push_back(pp);
   solarWindPopulation swpop;
   swpop.m = species->m;
   swpop.q = species->q;
   swpop.U = U;
   swpop.n = n;
   swpop.vth = vth;
   swpop.T = T;
   swpop.name = species->name;
   Hybrid::swPops.push_back(swpop);
   return initialized;
}
// CHAPMAN IONOSPHERE EMISSION INJECTOR

InjectorChapmanIonosphere::InjectorChapmanIonosphere(): ParticleInjectorBase() { 
   initialized = false;
   N_ionoPop = -1;
   N_macroParticlesPerCell = -1.0;
   N_macroParticlesPerDt = -1.0;
   vth = w = R = T = noonFactor = nightFactor= 0.0;
}

InjectorChapmanIonosphere::~InjectorChapmanIonosphere() {finalize();}

bool InjectorChapmanIonosphere::finalize() {
   if(initialized == false) { return true; }
   initialized = false;
   return true;
}

bool InjectorChapmanIonosphere::inject(pargrid::DataID speciesDataID,unsigned int* N_particles) {
   if(initialized == false) { return initialized; }
   bool success = true;
   if(sim->timestep <= 0) { return success; }
   pargrid::DataWrapper<Particle<Real> > wrapper = simClasses->pargrid.getUserDataDynamic<Particle<Real> >(speciesDataID);
   for(pargrid::CellID b=0;b<simClasses->pargrid.getNumberOfLocalCells();++b) {
      if(injectParticles(b,*species,N_particles,wrapper) == false) { success = false; }
   }
   return success;
}

bool InjectorChapmanIonosphere::injectParticles(pargrid::CellID blockID,const Species& species,unsigned int* N_particles,
				       pargrid::DataWrapper<Particle<Real> >& wrapper) {
   
   const Real* crd = getBlockCoordinateArray(*sim,*simClasses);
   const size_t b3 = 3*blockID;
   const Real xBlock = crd[b3+0];
   const Real yBlock = crd[b3+1];
   const Real zBlock = crd[b3+2];
   Real* cellIonosphere = simClasses->pargrid.getUserDataStatic<Real>(Hybrid::dataCellIonosphereID);
   vector<Real> xinj,yinj,zinj;
   for(int k=0;k<block::WIDTH_Z;++k) for(int j=0;j<block::WIDTH_Y;++j) for(int i=0;i<block::WIDTH_X;++i) {
      const int n = (blockID*block::SIZE+block::index(i,j,k));
      const size_t nIono = n*Hybrid::N_ionospherePopulations + N_ionoPop;
      const int N_injectCell = probround(*simClasses,cellIonosphere[nIono]);
      if(N_injectCell <= 0) { return true; }
      const Real xCell = crd[b3+0] +(i+0.5)*Hybrid::dx;
      const Real yCell = crd[b3+1] +(j+0.5)*Hybrid::dx;
      const Real zCell = crd[b3+2] +(k+0.5)*Hybrid::dx;
       
      for(int s = 0;s<N_injectCell;s++) {
      	const Real eps = 1.0e-2;
      	const Real x = xCell + (1.0-eps)*Hybrid::dx*(simClasses->random.uniform()-0.5);
      	const Real y = yCell + (1.0-eps)*Hybrid::dx*(simClasses->random.uniform()-0.5);
      	const Real z = zCell + (1.0-eps)*Hybrid::dx*(simClasses->random.uniform()-0.5);
         
         Real r = sqrt(sqr(x)+ sqr(y)+sqr(z));
         Real sza = acos(x/r);
         Real a =0;
         if(sza < M_PI/2) {
             a = noonFactor + (nightFactor - noonFactor) * ( 1-cos(sza) );
         } else {
             a = nightFactor;
         }

         Real g = constants::GRAVITY*constants::MASS_MARS/sqr(constants::DIST_MARS_RADIUS);
         Real H = (constants::BOLTZMANN*T)/(species.m*g);
         Real h_prime = (r-R)/H-log(1.0/a);
         Real a0 = a*exp(1-h_prime-exp(-1*h_prime));
         Real a1 = simClasses->random.uniform();

         //simClasses->logger<<"Injecting:"<<N_injectCell<<"," << a0<<", "<<a<<","<<h_prime <<sza<<endl<<write; 
         if (a1<a0){
	    xinj.push_back(x);
	    yinj.push_back(y);
	    zinj.push_back(z);
         }
         else {
             s--;
         }
      }
   }  

   // make room for new particles:
   const int N_inject = static_cast<int>(xinj.size());
   const pargrid::ArraySizetype oldSize = wrapper.size()[blockID];
   N_particles[blockID] += N_inject;
   wrapper.resize(blockID,oldSize+N_inject);
   Particle<Real>* particles = wrapper.data()[blockID];
   size_t s = 0;
   for(size_t p=oldSize; p<oldSize+N_inject; ++p) {
      particles[p].state[particle::X] = xinj[s]-xBlock;
      particles[p].state[particle::Y] = yinj[s]-yBlock;
      particles[p].state[particle::Z] = zinj[s]-zBlock;
      const Real x = xinj[s];
      const Real y = yinj[s];
      const Real z = zinj[s];
      Real vx = vth*gaussrnd(*simClasses);
      Real vy = vth*gaussrnd(*simClasses);
      Real vz = vth*gaussrnd(*simClasses);
      // make sure that velocity is upwards
      if (vx*x + vy*y + vz*z < 0) {
         vx = -vx;
         vy = -vy;
         vz = -vz;
      }
      particles[p].state[particle::VX] = vx;
      particles[p].state[particle::VY] = vy;
      particles[p].state[particle::VZ] = vz;
      particles[p].state[particle::WEIGHT] = w;
#ifdef ION_SPECTRA_ALONG_ORBIT
      particles[p].state[particle::INI_CELLID] = simClasses->pargrid.getGlobalIDs()[blockID];
      particles[p].state[particle::INI_X] = xBlock + particles[p].state[particle::X];
      particles[p].state[particle::INI_Y] = yBlock + particles[p].state[particle::Y];
      particles[p].state[particle::INI_Z] = zBlock + particles[p].state[particle::Z];
      particles[p].state[particle::INI_VX] =  particles[p].state[particle::VX];
      particles[p].state[particle::INI_VY] =  particles[p].state[particle::VY];
      particles[p].state[particle::INI_VZ] =  particles[p].state[particle::VZ];
      particles[p].state[particle::INI_TIME] = sim->t;
#endif
      // inject counter
      Hybrid::particleCounterInject[species.popid-1] += w;
      Hybrid::particleCounterInjectMacroparticles[species.popid-1] += 1;
      ++s;
   }
   return true;
}

bool InjectorChapmanIonosphere::addConfigFileItems(ConfigReader& cr,const std::string& configRegionName) {
   cr.add(configRegionName+".emission_radius","Radius of the spherical emission shell [m] (Real).",(Real)-1.0);
   cr.add(configRegionName+".noon","Noon (dayside) emission factor [-] (Real).",(Real)-1.0);
   cr.add(configRegionName+".night","Night side emission factor [-] (Real).",(Real)-1.0);
   cr.add(configRegionName+".temperature","Temperature [K] (float).",(Real)0.0);
   cr.add(configRegionName+".total_production_rate","Total production rate of physical particles per second [#/s] (Real).",(Real)-1.0);
   cr.add(configRegionName+".macroparticles_per_cell","Number of macroparticles per cell wtr. to the first solar wind population [#] (Real).",(Real)-1.0);
   return true;
}

bool InjectorChapmanIonosphere::initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,
				    const std::string& configRegionName,const ParticleListBase* plist) {
   initialized = ParticleInjectorBase::initialize(sim,simClasses,cr,configRegionName,plist);
   this->species = reinterpret_cast<const Species*>(plist->getSpecies());
   Real totalRate = 0.0;
   cr.parse();
   cr.get(configRegionName+".emission_radius",R);
   cr.get(configRegionName+".noon",noonFactor);
   cr.get(configRegionName+".night",nightFactor);
   cr.get(configRegionName+".temperature",T);
   cr.get(configRegionName+".total_production_rate",totalRate);
   cr.get(configRegionName+".macroparticles_per_cell",N_macroParticlesPerCell);
   if(Hybrid::swMacroParticlesCellPerDt > 0.0) {
      N_macroParticlesPerDt = N_macroParticlesPerCell*Hybrid::swMacroParticlesCellPerDt;
   }
   else {
      simClasses.logger
        << "(" << species->name << ") WARNING: No solar wind macroparticles cell per dt rate found, assuming 100" << endl << write;
      N_macroParticlesPerDt = N_macroParticlesPerCell*100;
   }
   if(N_macroParticlesPerDt > 0.0 && totalRate > 0.0) {
      w = totalRate*sim.dt/N_macroParticlesPerDt;
   }
   else {
      N_macroParticlesPerDt = totalRate = w = 0.0;
   }
   if(T > 0) { vth = sqrt(constants::BOLTZMANN*T/species->m); }
   else { vth = 0.0; }
   simClasses.logger
     << "(" << species->name << ") emission radius = " << radiusToString(R) << endl
     << "(" << species->name << ") noon emission   = " << noonFactor << endl
     << "(" << species->name << ") night emission  = " << nightFactor << endl
     << "(" << species->name << ") temperature     = " << T << " K = " << T/constants::EV_TO_KELVIN << " eV" << endl
     << "(" << species->name << ") thermal speed   = " << vth/1e3 << " km/s" << endl
     << "(" << species->name << ") density         = " << ( totalRate/( 4*M_PI*sqr(R)*sqrt( constants::BOLTZMANN*T/(2*M_PI*species->m) ) ) )/1e6<< " cm^-3" << endl
     << "(" << species->name << ") total ion production rate = " << totalRate << " 1/s" << endl
     << "(" << species->name << ") macroparticles per cell   = " << N_macroParticlesPerCell << endl
     << "(" << species->name << ") macroparticles per dt     = " << N_macroParticlesPerDt << endl
     << "(" << species->name << ") macroparticle weight      = " << w << endl << write;

   static unsigned int ionoPopCnt = 0;
   N_ionoPop = ionoPopCnt;
   ionoPopCnt++;
   
   const Real* crd = getBlockCoordinateArray(sim,simClasses);
   Real N_insideSum = 0.0;
   Real* cellIonosphere = simClasses.pargrid.getUserDataStatic<Real>(Hybrid::dataCellIonosphereID);
   // go thru local blocks
   for(pargrid::CellID b=0;b<simClasses.pargrid.getNumberOfLocalCells();++b) {
      const size_t b3 = 3*b;
      // go thru cells in a block
      for(int k=0;k<block::WIDTH_Z;++k) for(int j=0;j<block::WIDTH_Y;++j) for(int i=0;i<block::WIDTH_X;++i) {
      	 const int n = (b*block::SIZE+block::index(i,j,k));
      	 const size_t nIono = n*Hybrid::N_ionospherePopulations + N_ionoPop;
      	 // physical (global) coordinates
      	 const Real xCell = crd[b3+0] + (i+0.5)*Hybrid::dx;
      	 const Real yCell = crd[b3+1] + (j+0.5)*Hybrid::dx;
      	 const Real zCell = crd[b3+2] + (k+0.5)*Hybrid::dx;
      	 const Real xCellMin = xCell - 0.5*Hybrid::dx;
      	 const Real yCellMin = yCell - 0.5*Hybrid::dx;
      	 const Real zCellMin = zCell - 0.5*Hybrid::dx;
      	 const Real xCellMax = xCell + 0.5*Hybrid::dx;
      	 const Real yCellMax = yCell + 0.5*Hybrid::dx;
      	 const Real zCellMax = zCell + 0.5*Hybrid::dx;
      	 // generate N random points inside each cell
      	 Real N_inside = 0;
      	 for(int s = 0;s<2000;s++) {
      	   const Real eps = 1.0e-2;
      	   Real x = xCell + (1.0-eps)*Hybrid::dx*(simClasses.random.uniform()-0.5);
      	   Real y = yCell + (1.0-eps)*Hybrid::dx*(simClasses.random.uniform()-0.5);
      	   Real z = zCell + (1.0-eps)*Hybrid::dx*(simClasses.random.uniform()-0.5);



            Real a = 0;
            Real r = sqrt(sqr(x)+ sqr(y)+sqr(z));
            Real sza = acos(x/r);
            if(sza < M_PI/2) {
                a = noonFactor + (nightFactor - noonFactor) * ( 1-cos(sza) );
            } else {
                a = nightFactor;
            }

            Real g = constants::GRAVITY*constants::MASS_MARS/sqr(constants::DIST_MARS_RADIUS);
            Real H = (constants::BOLTZMANN*150)/(species->m*g);
            Real h_prime = (r-R)/H-log(1.0/a);
            Real a0 = a*exp(1-h_prime-exp(-1*h_prime));
            Real a1 = simClasses.random.uniform();
            if (a1<a0){
                N_inside++;
            }
            

      	}
         if (N_inside <10){ N_inside = 0;}
         N_insideSum += N_inside;
         cellIonosphere[nIono] = N_inside;     
      }
   }
   Real N_insideGlobal = 0.0;
   MPI_Reduce(&N_insideSum,&N_insideGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
   MPI_Bcast(&N_insideGlobal,1,MPI_Type<Real>(),sim.MASTER_RANK,sim.comm);
   // determine emission rate in a cell
   for(pargrid::CellID b=0;b<simClasses.pargrid.getNumberOfLocalCells();++b) {
      for(int k=0;k<block::WIDTH_Z;++k) for(int j=0;j<block::WIDTH_Y;++j) for(int i=0;i<block::WIDTH_X;++i) {
	 const int n = (b*block::SIZE+block::index(i,j,k));
	 const size_t nIono = n*Hybrid::N_ionospherePopulations + N_ionoPop;
	 if(N_insideGlobal > 0.0 && N_insideSum > 0.0)  {
	    cellIonosphere[nIono] = N_macroParticlesPerDt*cellIonosphere[nIono]/N_insideGlobal;
	 }
	 else {
	    cellIonosphere[nIono] = 0.0;
	 }
      }
   }

   return initialized;
}

// IONOSPHERE EMISSION INJECTOR

InjectorIonosphere::InjectorIonosphere(): ParticleInjectorBase() { 
   initialized = false;
   N_ionoPop = -1;
   N_macroParticlesPerCell = -1.0;
   N_macroParticlesPerDt = -1.0;
   vth = w = R = 0.0;
}

InjectorIonosphere::~InjectorIonosphere() {finalize();}

bool InjectorIonosphere::finalize() {
   if(initialized == false) { return true; }
   initialized = false;
   return true;
}

bool InjectorIonosphere::inject(pargrid::DataID speciesDataID,unsigned int* N_particles) {
   if(initialized == false) { return initialized; }
   bool success = true;
   if(sim->timestep <= 0) { return success; }
   pargrid::DataWrapper<Particle<Real> > wrapper = simClasses->pargrid.getUserDataDynamic<Particle<Real> >(speciesDataID);
   for(pargrid::CellID b=0;b<simClasses->pargrid.getNumberOfLocalCells();++b) {
      if(injectParticles(b,*species,N_particles,wrapper) == false) { success = false; }
   }
   return success;
}

bool InjectorIonosphere::injectParticles(pargrid::CellID blockID,const Species& species,unsigned int* N_particles,
				       pargrid::DataWrapper<Particle<Real> >& wrapper) {
   const Real* crd = getBlockCoordinateArray(*sim,*simClasses);
   const size_t b3 = 3*blockID;
   const Real xBlock = crd[b3+0];
   const Real yBlock = crd[b3+1];
   const Real zBlock = crd[b3+2];
   Real* cellIonosphere = simClasses->pargrid.getUserDataStatic<Real>(Hybrid::dataCellIonosphereID);
   vector<Real> xinj,yinj,zinj;
   for(int k=0;k<block::WIDTH_Z;++k) for(int j=0;j<block::WIDTH_Y;++j) for(int i=0;i<block::WIDTH_X;++i) {
      const int n = (blockID*block::SIZE+block::index(i,j,k));
      const size_t nIono = n*Hybrid::N_ionospherePopulations + N_ionoPop;
      const int N_injectCell = probround(*simClasses,cellIonosphere[nIono]);
      if(N_injectCell <= 0) { return true; }
      const Real xCell = (i+0.5)*Hybrid::dx;
      const Real yCell = (j+0.5)*Hybrid::dx;
      const Real zCell = (k+0.5)*Hybrid::dx;
      for(int s = 0;s<N_injectCell;s++) {
	 const Real eps = 1.0e-2;
	 const Real x = xCell + (1.0-eps)*Hybrid::dx*(simClasses->random.uniform()-0.5);
	 const Real y = yCell + (1.0-eps)*Hybrid::dx*(simClasses->random.uniform()-0.5);
	 const Real z = zCell + (1.0-eps)*Hybrid::dx*(simClasses->random.uniform()-0.5);
	 xinj.push_back(x);
	 yinj.push_back(y);
	 zinj.push_back(z);
      }
   }   
   // make room for new particles:
   const int N_inject = static_cast<int>(xinj.size());
   const pargrid::ArraySizetype oldSize = wrapper.size()[blockID];
   N_particles[blockID] += N_inject;
   wrapper.resize(blockID,oldSize+N_inject);
   Particle<Real>* particles = wrapper.data()[blockID];
   size_t s = 0;
   for(size_t p=oldSize; p<oldSize+N_inject; ++p) {
      particles[p].state[particle::X] = xinj[s];
      particles[p].state[particle::Y] = yinj[s];
      particles[p].state[particle::Z] = zinj[s];
      const Real x = xBlock + xinj[s];
      const Real y = yBlock + yinj[s];
      const Real z = zBlock + zinj[s];
      Real vx = vth*gaussrnd(*simClasses);
      Real vy = vth*gaussrnd(*simClasses);
      Real vz = vth*gaussrnd(*simClasses);
      // make sure that velocity is upwards
      if (vx*x + vy*y + vz*z < 0) {
         vx = -vx;
         vy = -vy;
         vz = -vz;
      }
      particles[p].state[particle::VX] = vx;
      particles[p].state[particle::VY] = vy;
      particles[p].state[particle::VZ] = vz;
      particles[p].state[particle::WEIGHT] = w;
/*#ifdef USE_DETECTORS
      particles[p].state[particle::INI_CELLID] = simClasses->pargrid.getGlobalIDs()[blockID];
      particles[p].state[particle::INI_X] = xBlock + particles[p].state[particle::X];
      particles[p].state[particle::INI_Y] = yBlock + particles[p].state[particle::Y];
      particles[p].state[particle::INI_Z] = zBlock + particles[p].state[particle::Z];
      particles[p].state[particle::INI_VX] =  particles[p].state[particle::VX];
      particles[p].state[particle::INI_VY] =  particles[p].state[particle::VY];
      particles[p].state[particle::INI_VZ] =  particles[p].state[particle::VZ];
      particles[p].state[particle::INI_TIME] = sim->t;
#endif*/
      // inject counter
      Hybrid::particleCounterInject[species.popid-1] += w;
      Hybrid::particleCounterInjectMacroparticles[species.popid-1] += 1;
      ++s;
   }
   return true;
}

bool InjectorIonosphere::addConfigFileItems(ConfigReader& cr,const std::string& configRegionName) {
   cr.add(configRegionName+".profile_name","Ionospheric emission profile name [-] (string)","");
   cr.add(configRegionName+".emission_radius","Radius of the spherical emission shell [m] (Real)",(Real)-1.0);
   cr.add(configRegionName+".noon","Noon (dayside) emission factor [-] (Real)",(Real)-1.0);
   cr.add(configRegionName+".night","Night side emission factor [-] (Real)",(Real)-1.0);
   cr.add(configRegionName+".temperature","Temperature [K] (float)",(Real)0.0);
   cr.add(configRegionName+".total_production_rate","Total production rate of physical particles per second [#/s] (Real)",(Real)-1.0);
   cr.add(configRegionName+".macroparticles_per_cell","Number of macroparticles per cell wtr. to the first solar wind population [#] (Real)",(Real)-1.0);
   return true;
}

bool InjectorIonosphere::initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,
				    const std::string& configRegionName,const ParticleListBase* plist) {
   initialized = ParticleInjectorBase::initialize(sim,simClasses,cr,configRegionName,plist);
   this->species = reinterpret_cast<const Species*>(plist->getSpecies());
   string profileName = "";
   Real noonFactor = -1.0;
   Real nightFactor = -1.0;
   Real T = 0.0;
   Real totalRate = 0.0;
   cr.parse();
   cr.get(configRegionName+".profile_name",profileName);
   cr.get(configRegionName+".emission_radius",R);
   cr.get(configRegionName+".noon",noonFactor);
   cr.get(configRegionName+".night",nightFactor);
   cr.get(configRegionName+".temperature",T);
   cr.get(configRegionName+".total_production_rate",totalRate);
   cr.get(configRegionName+".macroparticles_per_cell",N_macroParticlesPerCell);
   if(noonFactor < 0.0 || nightFactor < 0.0) {
      simClasses.logger
        << "(" << species->name << ") WARNING: negative noon or night factor" << endl << write;
   }
   if(Hybrid::swMacroParticlesCellPerDt > 0.0) {
      N_macroParticlesPerDt = N_macroParticlesPerCell*Hybrid::swMacroParticlesCellPerDt;
   }
   else {
      simClasses.logger
        << "(" << species->name << ") WARNING: No solar wind macroparticles cell per dt rate found, assuming 100" << endl << write;
      N_macroParticlesPerDt = N_macroParticlesPerCell*100;
   }
   if(N_macroParticlesPerDt > 0.0 && totalRate > 0.0) {
      w = totalRate*sim.dt/N_macroParticlesPerDt;
   }
   else {
      N_macroParticlesPerDt = totalRate = w = 0.0;
   }
   if(T > 0) { vth = sqrt(constants::BOLTZMANN*T/species->m); }
   else { vth = 0.0; }
   simClasses.logger
     << "(" << species->name << ") emission profile = " << profileName << endl
     << "(" << species->name << ") emission radius  = " << radiusToString(R) << endl
     << "(" << species->name << ") noon emission    = " << noonFactor << endl
     << "(" << species->name << ") night emission   = " << nightFactor << endl
     << "(" << species->name << ") temperature      = " << T << " K = " << T/constants::EV_TO_KELVIN << " eV" << endl
     << "(" << species->name << ") thermal speed    = " << vth/1e3 << " km/s" << endl
     << "(" << species->name << ") density          = " << ( totalRate/( 4*M_PI*sqr(R)*sqrt( constants::BOLTZMANN*T/(2*M_PI*species->m) ) ) )/1e6<< " cm^-3" << endl
     << "(" << species->name << ") total ion production rate = " << totalRate << " 1/s" << endl
     << "(" << species->name << ") macroparticles per cell   = " << N_macroParticlesPerCell << endl
     << "(" << species->name << ") macroparticles per dt     = " << N_macroParticlesPerDt << endl
     << "(" << species->name << ") macroparticle weight      = " << w << endl << write;

   static unsigned int ionoPopCnt = 0;
   N_ionoPop = ionoPopCnt;
   ionoPopCnt++;
   
   const Real* crd = getBlockCoordinateArray(sim,simClasses);
   Real N_insideSum = 0.0;
   Real* cellIonosphere = simClasses.pargrid.getUserDataStatic<Real>(Hybrid::dataCellIonosphereID);
   // go thru local blocks
   for(pargrid::CellID b=0;b<simClasses.pargrid.getNumberOfLocalCells();++b) {
      const size_t b3 = 3*b;
      // go thru cells in a block
      for(int k=0;k<block::WIDTH_Z;++k) for(int j=0;j<block::WIDTH_Y;++j) for(int i=0;i<block::WIDTH_X;++i) {
	 const int n = (b*block::SIZE+block::index(i,j,k));
	 const size_t nIono = n*Hybrid::N_ionospherePopulations + N_ionoPop;
	 // physical (global) coordinates
	 const Real xCell = crd[b3+0] + (i+0.5)*Hybrid::dx;
	 const Real yCell = crd[b3+1] + (j+0.5)*Hybrid::dx;
	 const Real zCell = crd[b3+2] + (k+0.5)*Hybrid::dx;
	 const Real xCellMin = xCell - 0.5*Hybrid::dx;
	 const Real yCellMin = yCell - 0.5*Hybrid::dx;
	 const Real zCellMin = zCell - 0.5*Hybrid::dx;
	 const Real xCellMax = xCell + 0.5*Hybrid::dx;
	 const Real yCellMax = yCell + 0.5*Hybrid::dx;
	 const Real zCellMax = zCell + 0.5*Hybrid::dx;
	 // generate N random points inside each cell
	 Real N_inside = 0;
	 for(int s = 0;s<1000;s++) {
	    const Real eps = 1.0e-2;
	    Real x = xCell + (1.0-eps)*Hybrid::dx*(simClasses.random.uniform()-0.5);
	    Real y = yCell + (1.0-eps)*Hybrid::dx*(simClasses.random.uniform()-0.5);
	    Real z = zCell + (1.0-eps)*Hybrid::dx*(simClasses.random.uniform()-0.5);
	    // project points on a shell r = R
	    const Real norm = R/sqrt(sqr(x) + sqr(y) + sqr(z));
	    x *= norm;
	    y *= norm;
	    z *= norm;
	    // check how many of the projected points are inside each cell
	    if(x > xCellMin && x < xCellMax && y > yCellMin && y < yCellMax && z > zCellMin && z < zCellMax) {
	       N_inside++;
	    }
	 }
         // dayside: cos(sza) dependency and nightside: constant
         if(profileName.compare("ionoCosSzaDayConstantNight") == 0) {
            const Real rr = sqrt(sqr(xCell) + sqr(yCell) + sqr(zCell));
            const Real sza = acos(xCell/rr);
            Real a = 0.0;
            if(sza < M_PI/2) {
               a = noonFactor + (nightFactor - noonFactor) * ( 1 - cos(sza) );
            } else {
               a = nightFactor;
            }
            N_inside *= a;
         }
         // dayside: constant and nightside: constant
         else if(profileName.compare("ionoConstantDayConstantNight") == 0) {
            const Real rr = sqrt(sqr(xCell) + sqr(yCell) + sqr(zCell));
            const Real sza = acos(xCell/rr);
            Real a = 0.0;
            if(sza < M_PI/2) {
               a = noonFactor;
            } else {
               a = nightFactor;
            }
            N_inside *= a;
         }
         // polar cap source
         else if(profileName.compare("ionoPolarCap") == 0) {
            // angle towards noon
            const Real t = fabs(noonFactor)*M_PI/180.0;
            // half angle of polar cap in degrees
            const Real cap = fabs(nightFactor)*M_PI/180.0;
            Real x = 0.0;
            Real z = 0.0;
            // rotation around y axis
            if(zCell > 0.0) { // north
               x = xCell*cos(-t) + zCell*sin(-t);
               z = xCell*sin(-t) + zCell*cos(-t);
            }
            else { // south
               x = xCell*cos(t) + zCell*sin(t);
               z = xCell*sin(t) + zCell*cos(t);
            }
            const Real rr = sqrt(sqr(x) + sqr(yCell) + sqr(z));
            const Real sza = acos(x/rr);
            const Real colat = fabs(M_PI/2.0 - sza);
            Real a = 0.0;
            if(colat <= cap) {
               a = 1.0;
            } else {
               a = 0.0;
            }
            N_inside *= a;
         }
         N_insideSum += N_inside;
	 cellIonosphere[nIono] = N_inside;
      }
   }
   Real N_insideGlobal = 0.0;
   MPI_Reduce(&N_insideSum,&N_insideGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
   MPI_Bcast(&N_insideGlobal,1,MPI_Type<Real>(),sim.MASTER_RANK,sim.comm);
   // determine emission rate in a cell
   for(pargrid::CellID b=0;b<simClasses.pargrid.getNumberOfLocalCells();++b) {
      for(int k=0;k<block::WIDTH_Z;++k) for(int j=0;j<block::WIDTH_Y;++j) for(int i=0;i<block::WIDTH_X;++i) {
	 const int n = (b*block::SIZE+block::index(i,j,k));
	 const size_t nIono = n*Hybrid::N_ionospherePopulations + N_ionoPop;
	 if(N_insideGlobal > 0.0 && N_insideSum > 0.0)  {
	    cellIonosphere[nIono] = N_macroParticlesPerDt*cellIonosphere[nIono]/N_insideGlobal;
	 }
	 else {
	    cellIonosphere[nIono] = 0.0;
	 }
      }
   }
   particlePopulation pp;
   pp.w = w;
   pp.name = species->name;
   Hybrid::allPops.push_back(pp);
   return initialized;
}

// EXOSPHERE PHOTOION INJECTOR

InjectorExosphere::InjectorExosphere(): ParticleInjectorBase() { 
   initialized = false;
   N_exoPop = -1;
   neutralProfileName = "";
   N_macroParticlesPerCell = -1.0;
   N_macroParticlesPerDt = -1.0;
   vth = w = r0 = R_exobase = R_shadow = 0.0;
   n0.clear();
   H0.clear();
   T0.clear();
   k0.clear();
}

InjectorExosphere::~InjectorExosphere() {finalize();}

bool InjectorExosphere::finalize() {
   if(initialized == false) { return true; }
   initialized = false;
   return true;
}

bool InjectorExosphere::inject(pargrid::DataID speciesDataID,unsigned int* N_particles) {
   if(initialized == false) { return initialized; }
   bool success = true;
   if(sim->timestep <= 0) { return success; }
   pargrid::DataWrapper<Particle<Real> > wrapper = simClasses->pargrid.getUserDataDynamic<Particle<Real> >(speciesDataID);
   for(pargrid::CellID b=0;b<simClasses->pargrid.getNumberOfLocalCells();++b) {
      if(injectParticles(b,*species,N_particles,wrapper) == false) { success = false; }
   }
   return success;
}

bool InjectorExosphere::injectParticles(pargrid::CellID blockID,const Species& species,unsigned int* N_particles,
				       pargrid::DataWrapper<Particle<Real> >& wrapper) {
/*#ifdef USE_DETECTORS
   const Real* crd = getBlockCoordinateArray(*sim,*simClasses);
   const size_t b3 = 3*blockID;
   const Real xBlock = crd[b3+0];
   const Real yBlock = crd[b3+1];
   const Real zBlock = crd[b3+2];
#endif*/
   Real* cellExosphere = simClasses->pargrid.getUserDataStatic<Real>(Hybrid::dataCellExosphereID);
   vector<Real> xinj,yinj,zinj;
   for(int k=0;k<block::WIDTH_Z;++k) for(int j=0;j<block::WIDTH_Y;++j) for(int i=0;i<block::WIDTH_X;++i) {
      const int n = (blockID*block::SIZE+block::index(i,j,k));
      const size_t nExo = n*Hybrid::N_exospherePopulations + N_exoPop;
      const int N_injectCell = probround(*simClasses,cellExosphere[nExo]);
      if(N_injectCell <= 0) { return true; }
      const Real xCell = (i+0.5)*Hybrid::dx;
      const Real yCell = (j+0.5)*Hybrid::dx;
      const Real zCell = (k+0.5)*Hybrid::dx;
      for(int s = 0;s<N_injectCell;s++) {
	 const Real eps = 1.0e-2;
	 const Real x = xCell + (1.0-eps)*Hybrid::dx*(simClasses->random.uniform()-0.5);
	 const Real y = yCell + (1.0-eps)*Hybrid::dx*(simClasses->random.uniform()-0.5);
	 const Real z = zCell + (1.0-eps)*Hybrid::dx*(simClasses->random.uniform()-0.5);
	 xinj.push_back(x);
	 yinj.push_back(y);
	 zinj.push_back(z);
      }
   }
   // Make room for new particles:
   const int N_inject = static_cast<int>(xinj.size());
   const pargrid::ArraySizetype oldSize = wrapper.size()[blockID];
   N_particles[blockID] += N_inject;
   wrapper.resize(blockID,oldSize+N_inject);
   Particle<Real>* particles = wrapper.data()[blockID];
   size_t s = 0;
   for(size_t p=oldSize; p<oldSize+N_inject; ++p) {
      particles[p].state[particle::X] = xinj[s];
      particles[p].state[particle::Y] = yinj[s];
      particles[p].state[particle::Z] = zinj[s];
      particles[p].state[particle::VX] = vth*gaussrnd(*simClasses);
      particles[p].state[particle::VY] = vth*gaussrnd(*simClasses);
      particles[p].state[particle::VZ] = vth*gaussrnd(*simClasses);
      particles[p].state[particle::WEIGHT] = w;
/*#ifdef USE_DETECTORS
      particles[p].state[particle::INI_CELLID] = simClasses->pargrid.getGlobalIDs()[blockID];
      particles[p].state[particle::INI_X] = xBlock + particles[p].state[particle::X];
      particles[p].state[particle::INI_Y] = yBlock + particles[p].state[particle::Y];
      particles[p].state[particle::INI_Z] = zBlock + particles[p].state[particle::Z];
      particles[p].state[particle::INI_VX] =  particles[p].state[particle::VX];
      particles[p].state[particle::INI_VY] =  particles[p].state[particle::VY];
      particles[p].state[particle::INI_VZ] =  particles[p].state[particle::VZ];
      particles[p].state[particle::INI_TIME] = sim->t;
#endif*/
      // inject counter
      Hybrid::particleCounterInject[species.popid-1] += w;
      Hybrid::particleCounterInjectMacroparticles[species.popid-1] += 1;
      ++s;
   }
   return true;
}

bool InjectorExosphere::addConfigFileItems(ConfigReader& cr,const std::string& configRegionName) {
   cr.add(configRegionName+".neutral_profile","Name of the neutral density profile (string).",string(""));
   cr.add(configRegionName+".neutral_profile.r0","Chamberlain: r0 radius [m] (Real).",(Real)-1.0);
   cr.addComposed(configRegionName+".neutral_profile.n0","Neutral profile: Number densities at r = r0 [m^-3] (Real vector).");
   cr.addComposed(configRegionName+".neutral_profile.H0","Neutral profile: Scale heights at r = r0 [m] (Real vector).");
   cr.addComposed(configRegionName+".neutral_profile.T0","Neutral profile: Temperature at r = r0 [K] (Real vector).");
   cr.addComposed(configRegionName+".neutral_profile.k0","Neutral profile: Exponent at r = r0 [-] (Real vector).");
   cr.add(configRegionName+".temperature","Temperature [K] (Real).",(Real)-1.0);
   cr.add(configRegionName+".exobase_radius","Radius of the exobase [m] (Real).",(Real)-1.0);
   cr.add(configRegionName+".shadow_radius","Radius of the shadow [m] (Real).",(Real)-1.0);
   cr.add(configRegionName+".total_production_rate","Total production rate of physical particles per second [#/s] (Real).",(Real)-1.0);
   cr.add(configRegionName+".macroparticles_per_cell","Number of macroparticles per cell wtr. to the first solar wind population [#] (Real).",(Real)-1.0);
   return true;
}

bool InjectorExosphere::initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,
				    const std::string& configRegionName,const ParticleListBase* plist) {
   initialized = ParticleInjectorBase::initialize(sim,simClasses,cr,configRegionName,plist);
   this->species = reinterpret_cast<const Species*>(plist->getSpecies());
   Real T = 0.0;
   Real totalRate = 0.0;
   cr.parse();
   cr.get(configRegionName+".neutral_profile",neutralProfileName);
   cr.get(configRegionName+".neutral_profile.r0",r0);
   cr.get(configRegionName+".neutral_profile.n0",n0);
   cr.get(configRegionName+".neutral_profile.H0",H0);
   cr.get(configRegionName+".neutral_profile.T0",T0);
   cr.get(configRegionName+".neutral_profile.k0",k0);
   cr.get(configRegionName+".temperature",T);
   cr.get(configRegionName+".exobase_radius",R_exobase);
   cr.get(configRegionName+".shadow_radius",R_shadow);
   cr.get(configRegionName+".total_production_rate",totalRate);
   cr.get(configRegionName+".macroparticles_per_cell",N_macroParticlesPerCell);
   if(Hybrid::swMacroParticlesCellPerDt > 0.0) {
      N_macroParticlesPerDt = N_macroParticlesPerCell*Hybrid::swMacroParticlesCellPerDt;
   }
   else {
      simClasses.logger
        << "(" << species->name << ") WARNING: No solar wind macroparticles cell per dt rate found, assuming 100" << endl << write;
      N_macroParticlesPerDt = N_macroParticlesPerCell*100;
   }
   if(N_macroParticlesPerDt > 0.0 && totalRate > 0.0) {
      w = totalRate*sim.dt/N_macroParticlesPerDt;
   }
   else {
      N_macroParticlesPerDt = totalRate = w = 0.0;
   }
   if(T > 0) { vth = sqrt(constants::BOLTZMANN*T/species->m); }
   else { vth = 0.0; }
   static unsigned int exoPopCnt = 0;
   N_exoPop = exoPopCnt;
   exoPopCnt++;
   
   Real* cellExosphere = simClasses.pargrid.getUserDataStatic<Real>(Hybrid::dataCellExosphereID);
   Real sumThisProcess = 0.0;
   for(pargrid::CellID b=0;b<simClasses.pargrid.getNumberOfLocalCells();++b) {  
      const Real* crd = getBlockCoordinateArray(sim,simClasses);
      const size_t b3 = 3*b;
      for(int k=0;k<block::WIDTH_Z;++k) for(int j=0;j<block::WIDTH_Y;++j) for(int i=0;i<block::WIDTH_X;++i) {
	 const int n = (b*block::SIZE+block::index(i,j,k));
	 const size_t nExo = n*Hybrid::N_exospherePopulations + N_exoPop;
	 const Real xCell = crd[b3+0] + (i+0.5)*Hybrid::dx;
	 const Real yCell = crd[b3+1] + (j+0.5)*Hybrid::dx;
	 const Real zCell = crd[b3+2] + (k+0.5)*Hybrid::dx;
	 NeutralProfileArgs a;
	 a.m = species->m;
	 a.r0 = r0;
	 a.n0 = n0;
	 a.H0 = H0;
	 a.T0 = T0;
	 a.k0 = k0;
	 a.R_exobase = R_exobase;
	 a.R_shadow = R_shadow;
	 cellExosphere[nExo] = getNeutralDensity(simClasses,neutralProfileName,xCell,yCell,zCell,a)*Hybrid::dV;
	 if(cellExosphere[nExo] < 0.0) {
	    simClasses.logger << "(" << species->name << ") ERROR: Neutral profile init failed" << endl << write;
	    initialized = false;
	    goto break_for;
	 }
	 sumThisProcess += cellExosphere[nExo];
      }
   }
break_for:
   Real sumGlobal = 0.0;
   MPI_Reduce(&sumThisProcess,&sumGlobal,1,MPI_Type<Real>(),MPI_SUM,sim.MASTER_RANK,sim.comm);
   MPI_Bcast(&sumGlobal,1,MPI_Type<Real>(),sim.MASTER_RANK,sim.comm);
   const Real Ntot = sumGlobal;
   const Real ionizationRate = totalRate/Ntot;
   for(pargrid::CellID b=0;b<simClasses.pargrid.getNumberOfLocalCells();++b) {
      for(int k=0;k<block::WIDTH_Z;++k) for(int j=0;j<block::WIDTH_Y;++j) for(int i=0;i<block::WIDTH_X;++i) {
	 const int n = (b*block::SIZE+block::index(i,j,k));
	 const size_t nExo = n*Hybrid::N_exospherePopulations + N_exoPop;
	 cellExosphere[nExo] = N_macroParticlesPerDt*cellExosphere[nExo]*ionizationRate/totalRate; }
   }
   simClasses.logger
     << "(" << species->name << ") neutral profile     = " << neutralProfileName << endl
     << "(" << species->name << ") neutral profile: r0 = " << radiusToString(r0) << endl
     << "(" << species->name << ") neutral profile: n0 = ";
   for(size_t i=0;i<n0.size();i++) { simClasses.logger << n0[i]/1e6 << " "; }
   simClasses.logger
     << "cm^-3" << endl
     << "(" << species->name << ") neutral profile: H0 = ";
   for(size_t i=0;i<H0.size();i++) { simClasses.logger << H0[i]/1e3 << " "; }
   simClasses.logger
     << "km" << endl
     << "(" << species->name << ") neutral profile: T0 = ";
   for(size_t i=0;i<T0.size();i++) { simClasses.logger << T0[i] << " "; }
   simClasses.logger << "K = ";
   for(size_t i=0;i<T0.size();i++) { simClasses.logger << T0[i]/constants::EV_TO_KELVIN << " "; }   
   simClasses.logger << "eV" << endl
     << "(" << species->name << ") neutral profile: k0 = ";
   for(size_t i=0;i<k0.size();i++) { simClasses.logger << k0[i] << " "; }
   simClasses.logger
     << endl
     << "(" << species->name << ") temperature    = " << T << " K = " << T/constants::EV_TO_KELVIN << " eV" << endl
     << "(" << species->name << ") thermal speed  = " << vth/1e3 << " km/s" << endl
     << "(" << species->name << ") exobase radius = " << radiusToString(R_exobase) << endl
     << "(" << species->name << ") shadow radius  = " << radiusToString(R_shadow)  << endl
     << "(" << species->name << ") total number of neutrals  = " << Ntot << endl
     << "(" << species->name << ") ionization rate           = " << ionizationRate << " 1/s" << endl
     << "(" << species->name << ") total ion production rate = " << totalRate << " 1/s" << endl
     << "(" << species->name << ") macroparticles per cell   = " << N_macroParticlesPerCell << endl
     << "(" << species->name << ") macroparticles per dt     = " << N_macroParticlesPerDt << endl
     << "(" << species->name << ") macroparticle weight      = " << w << endl << write;
   particlePopulation pp;
   pp.w = w;
   pp.name = species->name;
   Hybrid::allPops.push_back(pp);
   return initialized;
}

