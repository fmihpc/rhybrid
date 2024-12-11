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
#include <algorithm>

#include "particle_injector.h"
#include "hybrid.h"
#include "neutral_profiles.h"

using namespace std;

/* \brief Probabilistic Real2int rounding
 *
 * probround(x) (x >= 0) gives either floor(x) or ceil(x), with probability
 * depending on which one is closer. For example, probround(2.3) gives 2
 * with 70% probability and 3 with 30% probability.
 */
int probround(SimulationClasses& simClasses,Real x) {
   if(x <= 0) { return 0; }
   const int f = static_cast<int>(floor(x));
   const Real r = simClasses.random.uniform();
   if(r < x-f) { return f+1; }
   else { return f; }
}

/** Gaussian randomness
 *
 *  Generate a Gaussian deviate with zero mean and unit
 *  standard deviation: f(x) = (1/(2*pi))*exp(-0.5*x^2), x real.
 *  Algorithm: Generate random pairs (x,y) from unit square
 *  -1 <= x <= 1, -1 <= y <= 1 until (x,y) is within
 *  the unit circle. Compute fac = sqrt(-2.0*log(r2)/r2),
 *  where r2 = x^2 + y^2. Then, x*fac and y*fac are two Gaussian
 *  random numbers.
 */
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

/** Deriv Gaussian randomness
 *
 * Return a random number distributed according to
 * f(x) = c*max(0,x)*exp(-0.5*(x-x0)^2) where the normalization constant c
 * is chosen so that the integrate(f(x),x,-inf,inf)=1 (notice that f(x)=0 for x<=0).
 *
 * Method: F(x)=c*xm*exp(-0.5*(x-xm)^2-0.5*(x0-xm)^2), where xm=0.5*(x0+sqrt(x0^2+4)),
 * is a majorant, i.e. F(x) >= f(x) for all x>=0 and x0. The majorant is Gaussian
 * with unit standard deviation and mean equal to xm. (Note that xm is the abscissa
 * of the maximum of f(x), i.e. f'(xm)=0.) Generate random numbers x from
 * the majorant Gaussian and accept it with probability f(x)/F(x).
 * The area under the majorant curve F(x) is close to unity for x0>=0 so that
 * only a few trials are needed. For x0<0 it is asymptotically proportional
 * to (-x0) so that more and more trials are needed. Therefore, avoid calling
 * the function with x0 < -10.
 */
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

// convert string velocity vector cfg variable to Real vector
bool convertConfigFileVariableVelocity(string velStr,vector<Real>& vel) {
   // check velocity format: (Ux,Uy,Uz)
   bool velStrOk = true;
   if(count(velStr.begin(),velStr.end(),'(') != 1 ||
      count(velStr.begin(),velStr.end(),')') != 1 ||
      count(velStr.begin(),velStr.end(),',') != 2 ||
      velStr.find_first_not_of("(),+-0123456789.e ") != string::npos) {
      velStrOk = false;
   }
   // remove non-numeral characters from string
   string velStrEdit(velStr);
   replace(velStrEdit.begin(),velStrEdit.end(),'(',' ');
   replace(velStrEdit.begin(),velStrEdit.end(),')',' ');
   replace(velStrEdit.begin(),velStrEdit.end(),',',' ');
   // remove multiple whitespace: not needed
   //velStrEdit.erase(unique(velStrEdit.begin(),velStrEdit.end(),[](char a,char b) { return isspace(a) && isspace(b); } ),velStrEdit.end() );
   // convert string to Reals
   vel.clear();
   Real vtmp;
   stringstream ss(velStrEdit);
   while (ss >> vtmp) { vel.push_back(vtmp); }
   // check velocity format: (Ux,Uy,Uz)
   if(vel.size() != 3) { velStrOk = false; }
   return velStrOk;
}

// get injector parameters
bool getInjectorParameters(ParticleInjectorBase* injBasePtr,InjectorParameters& p) {
   p.name = "";
   p.type = "";
   p.m = p.q = p.w = p.T = p.vth = p.n = p.U = p.velocity[0] = p.velocity[1] = p.velocity[2] = 0.0;
   string injectorType = injBasePtr->getType();
   p.type = string(injectorType);
   if(injectorType.compare("uniform") == 0) {
      InjectorUniform* injPtr = dynamic_cast<InjectorUniform*>(injBasePtr);
      injPtr->getParams(p);
   }
   else if(injectorType.compare("ambient") == 0) {
      InjectorAmbient* injPtr = dynamic_cast<InjectorAmbient*>(injBasePtr);
      injPtr->getParams(p);
   }
   else if(injectorType.compare("solarwind") == 0) {
      InjectorSolarWind* injPtr = dynamic_cast<InjectorSolarWind*>(injBasePtr);
      injPtr->getParams(p);
   }
   else if(injectorType.compare("ionosphere") == 0) {
      InjectorIonosphere* injPtr = dynamic_cast<InjectorIonosphere*>(injBasePtr);
      injPtr->getParams(p);
   }
   else if(injectorType.compare("ionosphere_chapman") == 0) {
      InjectorChapmanIonosphere* injPtr = dynamic_cast<InjectorChapmanIonosphere*>(injBasePtr);
      injPtr->getParams(p);
   }
   else if(injectorType.compare("exosphere") == 0) {
      InjectorExosphere* injPtr = dynamic_cast<InjectorExosphere*>(injBasePtr);
      injPtr->getParams(p);
   }
   else { return false; }
   return true;
}

// UNIFORM INJECTOR

InjectorUniform::InjectorUniform(): ParticleInjectorBase() {
   initialized = false;
   velocity[0] = velocity[1] = velocity[2] = U = T = vth = n = w = xmin = xmax = 0.0;
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

bool InjectorUniform::injectParticles(pargrid::CellID blockID,const Species& species,unsigned int* N_particles,pargrid::DataWrapper<Particle<Real> >& wrapper) {
/*#ifdef USE_DETECTORS
   const Real* crd = getBlockCoordinateArray(*sim,*simClasses);
   const size_t b3 = 3*blockID;
   const Real xBlock = crd[b3+0];
   const Real yBlock = crd[b3+1];
   const Real zBlock = crd[b3+2];
#endif*/
   const Real* crd = getBlockCoordinateArray(*sim,*simClasses);
   const size_t b3 = 3*blockID;
   const Real xBlock = crd[b3+0];
   vector<Real> xinj,yinj,zinj;
   for(int k=0;k<block::WIDTH_Z;++k) for(int j=0;j<block::WIDTH_Y;++j) for(int i=0;i<block::WIDTH_X;++i) {
      const Real xCell = (i+0.5)*Hybrid::dx;
      const Real yCell = (j+0.5)*Hybrid::dx;
      const Real zCell = (k+0.5)*Hybrid::dx;
      const Real xCellGlobal = xBlock + xCell;
      if(xCellGlobal < xmin || xCellGlobal > xmax) { continue; }
      const int N_injectCell = probround(*simClasses,N_macroParticlesPerCell);
      //const int N_injectCell = probround(*simClasses,N_macroParticlesPerCell*(1.5 + sin( 10.0*xCellGlobal/(sim->x_max - sim->x_min) ))); // RHBTESTS: init uniform population with sine wave in density
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
      particles[p].state[particle::VX] = velocity[0] + vth*gaussrnd(*simClasses); // not exactly correct
      particles[p].state[particle::VY] = velocity[1] + vth*gaussrnd(*simClasses); // should be replaced with
      particles[p].state[particle::VZ] = velocity[2] + vth*gaussrnd(*simClasses); // derivgaussrnd
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
      Hybrid::logCounterParticleInject[species.popid-1] += w;
      Hybrid::logCounterParticleInjectMacroparticles[species.popid-1] += 1;
      Hybrid::logCounterParticleInjectKineticEnergy[species.popid-1] += w*( sqr(particles[p].state[particle::VX]) + sqr(particles[p].state[particle::VY]) + sqr(particles[p].state[particle::VZ]) );
      ++s;
   }
   return true;
}

bool InjectorUniform::addConfigFileItems(ConfigReader& cr,const std::string& configRegionName) {
   cr.add(configRegionName+".velocity","Bulk velocity vector (Ux,Uy,Uz) [m/s] (float,float,float).",string(""));
   cr.add(configRegionName+".density","Number density [m^-3] (float).",(Real)0.0);
   cr.add(configRegionName+".temperature","Temperature [K] (float).",(Real)0.0);
   cr.add(configRegionName+".macroparticles_per_cell","Number of macroparticles per cell [#] (Real).",(Real)-1.0);
   cr.add(configRegionName+".xmin","Minimum x coordinate [m] (float).",(Real)-1.0e22);
   cr.add(configRegionName+".xmax","Maximum x coordinate [m] (float).",(Real)+1.0e22);
   return true;
}

bool InjectorUniform::initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,const std::string& configRegionName,const ParticleListBase* plist) {
   initialized = ParticleInjectorBase::initialize(sim,simClasses,cr,configRegionName,plist);
   this->type = "uniform";
   this->species = reinterpret_cast<const Species*>(plist->getSpecies());
   string velStr = "";
   cr.parse();
   cr.get(configRegionName+".velocity",velStr);
   cr.get(configRegionName+".density",n);
   cr.get(configRegionName+".temperature",T);
   cr.get(configRegionName+".macroparticles_per_cell",N_macroParticlesPerCell);
   cr.get(configRegionName+".xmin",xmin);
   cr.get(configRegionName+".xmax",xmax);

   // bulk velocity vector
     {
	bool velStrOk = true;
	vector<Real> vel;
	velStrOk = convertConfigFileVariableVelocity(velStr,vel);
	// if not correct format
	if(velStrOk == false) {
	   simClasses.logger << "(" << species->name << ") ERROR: bad format of bulk velocity vector (" << velStr << ")" << endl << write;
	   initialized = false;
	}
	this->velocity[0] = vel[0];
	this->velocity[1] = vel[1];
	this->velocity[2] = vel[2];
     }
   // bulk speed
   U = sqrt(sqr(velocity[0]) + sqr(velocity[1]) + sqr(velocity[2]));

   if(N_macroParticlesPerCell > 0 && n > 0 && Hybrid::dx > 0) {
      w = n*Hybrid::dV/N_macroParticlesPerCell;
   }
   else {
      N_macroParticlesPerCell = n = w = 0.0;
   }
   if(T > 0) { vth = sqrt(constants::BOLTZMANN*T/species->m); }
   else { vth = 0.0; }
   if(xmin >= xmax) {
      simClasses.logger << "(" << species->name << ") ERROR: xmin >= xmax" << endl << write;
      initialized = false;
   }
   simClasses.logger
     << "(" << species->name << ") velocity      = (" << velocity[0]/1e3 << "," << velocity[1]/1e3 << "," << velocity[2]/1e3 << ") km/s" << endl
     << "(" << species->name << ") speed         = " << U/1e3 << " km/s" << endl
     << "(" << species->name << ") density       = " << n/1e6 << " cm^{-3}" << endl
     << "(" << species->name << ") temperature   = " << T << " K = " << T/constants::EV_TO_KELVIN << " eV" << endl
     << "(" << species->name << ") thermal speed = " << vth/1e3 << " km/s" << endl
     << "(" << species->name << ") macroparticles per cell = " << N_macroParticlesPerCell << endl
     << "(" << species->name << ") macroparticle weight    = " << w << endl
     << "(" << species->name << ") xmin = " << xmin/1e3 << " km" << endl
     << "(" << species->name << ") xmax = " << xmax/1e3 << " km" << endl;
   particlePopulationInfo popInfo;
   popInfo.w = w;
   popInfo.name = species->name;
   Hybrid::allPopsInfo.push_back(popInfo);
   return initialized;
}

// get injector parameters
void InjectorUniform::getParams(InjectorParameters& p) {
   p.name = species->name;
   p.m = species->m;
   p.q = species->q;
   p.w = w;
   p.T = T;
   p.vth = vth;
   p.n = n;
   p.U = U;
   p.velocity[0] = velocity[0];
   p.velocity[1] = velocity[1];
   p.velocity[2] = velocity[2];
}

// AMBIENT INJECTOR

InjectorAmbient::InjectorAmbient(): ParticleInjectorBase() {
   initialized = false;
   T = vth = w = 0.0;
   N_macroParticlesPerCellPerDt = -1.0;
   N_macroParticlesPerCell = -1.0;
}

InjectorAmbient::~InjectorAmbient() {finalize();}

bool InjectorAmbient::finalize() {
   if(initialized == false) { return true; }
   initialized = false;
   return true;
}

bool InjectorAmbient::inject(pargrid::DataID speciesDataID,unsigned int* N_particles) {
   if(initialized == false) { return initialized; }
   bool success = true;
   if(sim->timestep <= 0) { return success; }
   pargrid::DataWrapper<Particle<Real> > wrapper = simClasses->pargrid.getUserDataDynamic<Particle<Real> >(speciesDataID);
   for(pargrid::CellID b=0; b<simClasses->pargrid.getNumberOfLocalCells(); ++b) {
      // x walls
      if( (simClasses->pargrid.getNeighbourFlags(b) & Hybrid::Y_POS_EXISTS) != 0 &&
	  (simClasses->pargrid.getNeighbourFlags(b) & Hybrid::Y_NEG_EXISTS) != 0 &&
	  (simClasses->pargrid.getNeighbourFlags(b) & Hybrid::Z_POS_EXISTS) != 0 &&
	  (simClasses->pargrid.getNeighbourFlags(b) & Hybrid::Z_NEG_EXISTS) != 0) {
	 // +x ghost cells
	 if ( (simClasses->pargrid.getNeighbourFlags(b) & Hybrid::X_POS_EXISTS) == 0 ) {
	    if(injectParticles(b,*species,N_particles,wrapper,0) == false) { success = false; }
	 }
	 // -x ghost cells
	 else if ( (simClasses->pargrid.getNeighbourFlags(b) & Hybrid::X_NEG_EXISTS) == 0 ) {
	    if(injectParticles(b,*species,N_particles,wrapper,1) == false) { success = false; }
	 }
      }
      // y walls
      else if( (simClasses->pargrid.getNeighbourFlags(b) & Hybrid::X_POS_EXISTS) != 0 &&
	       (simClasses->pargrid.getNeighbourFlags(b) & Hybrid::X_NEG_EXISTS) != 0 &&
	       (simClasses->pargrid.getNeighbourFlags(b) & Hybrid::Z_POS_EXISTS) != 0 &&
	       (simClasses->pargrid.getNeighbourFlags(b) & Hybrid::Z_NEG_EXISTS) != 0) {
	 // +y ghost cells
	 if ( (simClasses->pargrid.getNeighbourFlags(b) & Hybrid::Y_POS_EXISTS) == 0 ) {
	    if(injectParticles(b,*species,N_particles,wrapper,2) == false) { success = false; }
	 }
	 // -y ghost cells
	 else if ( (simClasses->pargrid.getNeighbourFlags(b) & Hybrid::Y_NEG_EXISTS) == 0 ) {
	    if(injectParticles(b,*species,N_particles,wrapper,3) == false) { success = false; }
	 }
      }
      // z walls
      else if( (simClasses->pargrid.getNeighbourFlags(b) & Hybrid::X_POS_EXISTS) != 0 &&
	       (simClasses->pargrid.getNeighbourFlags(b) & Hybrid::X_NEG_EXISTS) != 0 &&
	       (simClasses->pargrid.getNeighbourFlags(b) & Hybrid::Y_POS_EXISTS) != 0 &&
	       (simClasses->pargrid.getNeighbourFlags(b) & Hybrid::Y_NEG_EXISTS) != 0) {
	 // +z ghost cells
	 if ( (simClasses->pargrid.getNeighbourFlags(b) & Hybrid::Z_POS_EXISTS) == 0 ) {
	    if(injectParticles(b,*species,N_particles,wrapper,4) == false) { success = false; }
	 }
	 // -z ghost cells
	 else if ( (simClasses->pargrid.getNeighbourFlags(b) & Hybrid::Z_NEG_EXISTS) == 0 ) {
	    if(injectParticles(b,*species,N_particles,wrapper,5) == false) { success = false; }
	 }
      }
   }
   return success;
}

// wall: +x = 0 ,-x = 1, +y = 2, -y = 3, +z = 4,-z = 5
bool InjectorAmbient::injectParticles(pargrid::CellID blockID,const Species& species,unsigned int* N_particles,pargrid::DataWrapper<Particle<Real> >& wrapper,unsigned int wall) {
   vector<Real> xinj,yinj,zinj;
   for(int k=0;k<block::WIDTH_Z;++k) for(int j=0;j<block::WIDTH_Y;++j) for(int i=0;i<block::WIDTH_X;++i) {
      const Real xCell = (i+0.5)*Hybrid::dx;
      const Real yCell = (j+0.5)*Hybrid::dx;
      const Real zCell = (k+0.5)*Hybrid::dx;
      const int N_injectCell = probround(*simClasses,N_macroParticlesPerCell);
      if(N_injectCell <= 0) { continue; }
      for(int s = 0;s<N_injectCell;s++) {
	 const Real eps = 1.0e-2;
	 Real x = xCell + (1.0-eps)*Hybrid::dx*(simClasses->random.uniform()-0.5);
	 Real y = yCell + (1.0-eps)*Hybrid::dx*(simClasses->random.uniform()-0.5);
	 Real z = zCell + (1.0-eps)*Hybrid::dx*(simClasses->random.uniform()-0.5);
	 // place particles just at the cell boundary
	 if      (wall == 0) { x = xCell - 0.5*Hybrid::dx; }
	 else if (wall == 1) { x = xCell + 0.5*Hybrid::dx; }
	 else if (wall == 2) { y = yCell - 0.5*Hybrid::dx; }
	 else if (wall == 3) { y = yCell + 0.5*Hybrid::dx; }
	 else if (wall == 4) { z = zCell - 0.5*Hybrid::dx; }
	 else if (wall == 5) { z = zCell + 0.5*Hybrid::dx; }
	 else { return false; }
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
      Real vx = vth*gaussrnd(*simClasses);
      Real vy = vth*gaussrnd(*simClasses);
      Real vz = vth*gaussrnd(*simClasses);
      // make sure directional flux is toward the simulation domain
      if      (wall == 0) { if (vx > 0.0) { vx = -vx; } }
      else if (wall == 1) { if (vx < 0.0) { vx = -vx; } }
      else if (wall == 2) { if (vy > 0.0) { vy = -vy; } }
      else if (wall == 3) { if (vy < 0.0) { vy = -vy; } }
      else if (wall == 4) { if (vz > 0.0) { vz = -vz; } }
      else if (wall == 5) { if (vz < 0.0) { vz = -vz; } }
      else { return false; }
      particles[p].state[particle::VX] = vx;
      particles[p].state[particle::VY] = vy;
      particles[p].state[particle::VZ] = vz;
      particles[p].state[particle::WEIGHT] = w;
      // inject counter
      Hybrid::logCounterParticleInject[species.popid-1] += w;
      Hybrid::logCounterParticleInjectMacroparticles[species.popid-1] += 1;
      Hybrid::logCounterParticleInjectKineticEnergy[species.popid-1] += w*( sqr(particles[p].state[particle::VX]) + sqr(particles[p].state[particle::VY]) + sqr(particles[p].state[particle::VZ]) );
      ++s;
   }
   return true;
}

bool InjectorAmbient::addConfigFileItems(ConfigReader& cr,const std::string& configRegionName) {
   cr.add(configRegionName+".density","Number density [m^-3] (float).",(Real)0.0);
   cr.add(configRegionName+".temperature","Temperature [K] (float).",(Real)0.0);
   cr.add(configRegionName+".macroparticles_per_cell","Number of macroparticles per cell [#] (Real).",(Real)-1.0);
   return true;
}

bool InjectorAmbient::initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,const std::string& configRegionName,const ParticleListBase* plist) {
   initialized = ParticleInjectorBase::initialize(sim,simClasses,cr,configRegionName,plist);
   this->type = "ambient";
   this->species = reinterpret_cast<const Species*>(plist->getSpecies());
   N_macroParticlesPerCell = -100.0;
   cr.parse();
   cr.get(configRegionName+".density",n);
   cr.get(configRegionName+".temperature",T);
   cr.get(configRegionName+".macroparticles_per_cell",N_macroParticlesPerCell);
   if(N_macroParticlesPerCell > 0 && n > 0 && Hybrid::dx > 0) {
      // Random particle flux through a cell face: Bittencourt (2004, p. 180) or Jarvinen et al. (2009)
      const Real N_realParticlesPerCellPerDt = sqr(Hybrid::dx)*n*sqrt( constants::BOLTZMANN*T/(2.0*M_PI*species->m) )*sim.dt;
      w = n*Hybrid::dV/N_macroParticlesPerCell;
      N_macroParticlesPerCellPerDt = N_realParticlesPerCellPerDt/w;
   }
   else {
      N_macroParticlesPerCellPerDt = N_macroParticlesPerCell = n = w = 0.0;
   }
   if(T > 0) { vth = sqrt(constants::BOLTZMANN*T/species->m); }
   else { vth = 0.0; }
   simClasses.logger
     << "(" << species->name << ") density       = " << n/1e6 << " cm^{-3}" << endl
     << "(" << species->name << ") temperature   = " << T << " K = " << T/constants::EV_TO_KELVIN << " eV" << endl
     << "(" << species->name << ") thermal speed = " << vth/1e3 << " km/s" << endl
     << "(" << species->name << ") macroparticles per cell         = " << N_macroParticlesPerCell << endl
     << "(" << species->name << ") macroparticles per cell per dt  = " << N_macroParticlesPerCellPerDt << endl
     << "(" << species->name << ") macroparticle weight    = " << w << endl;
   particlePopulationInfo popInfo;
   popInfo.w = w;
   popInfo.name = species->name;
   Hybrid::allPopsInfo.push_back(popInfo);
   return initialized;
}

// get injector parameters
void InjectorAmbient::getParams(InjectorParameters& p) {
   p.name = species->name;
   p.m = species->m;
   p.q = species->q;
   p.w = w;
   p.T = T;
   p.vth = vth;
}

// SOLAR WIND INJECTOR

InjectorSolarWind::InjectorSolarWind(): ParticleInjectorBase() { 
   initialized = false;
   velocity[0] = velocity[1] = velocity[2] = U = T = vth = n = w = 0.0;
   checkIfInjectionCellFuncPtr = &InjectorSolarWind::checkIfInjectionCellDefault;
   initParticleCrdVelFuncPtr = &InjectorSolarWind::initParticleCrdVelDefault;
   N_macroParticlesPerCellPerDt = -1.0;
   N_macroParticlesPerCell = -1.0;
}

InjectorSolarWind::~InjectorSolarWind() {finalize();}

bool InjectorSolarWind::finalize() {
   if(initialized == false) { return true; }
   initialized = false;
   return true;
}

// default function gives an error and exits if called
bool InjectorSolarWind::checkIfInjectionCellDefault(const pargrid::CellID b) {
   simClasses->logger << "(InjectorSolarWind): ERROR checkIfInjectionCellDefault called " << endl << write;
   exit(1);
   return false;
}

// return true if b is a non-edge cell at positive x outer boundary
bool InjectorSolarWind::checkIfInjectionCellXPos(const pargrid::CellID b) {
   return
    (((simClasses->pargrid.getNeighbourFlags(b) & Hybrid::X_POS_EXISTS) == 0) &&
     ((simClasses->pargrid.getNeighbourFlags(b) & Hybrid::Y_POS_EXISTS) != 0) &&
     ((simClasses->pargrid.getNeighbourFlags(b) & Hybrid::Y_NEG_EXISTS) != 0) &&
     ((simClasses->pargrid.getNeighbourFlags(b) & Hybrid::Z_POS_EXISTS) != 0) &&
     ((simClasses->pargrid.getNeighbourFlags(b) & Hybrid::Z_NEG_EXISTS) != 0));
}

// return true if b is a non-edge cell at negative x outer boundary
bool InjectorSolarWind::checkIfInjectionCellXNeg(const pargrid::CellID b) {
   return
    (((simClasses->pargrid.getNeighbourFlags(b) & Hybrid::X_NEG_EXISTS) == 0) &&
     ((simClasses->pargrid.getNeighbourFlags(b) & Hybrid::Y_POS_EXISTS) != 0) &&
     ((simClasses->pargrid.getNeighbourFlags(b) & Hybrid::Y_NEG_EXISTS) != 0) &&
     ((simClasses->pargrid.getNeighbourFlags(b) & Hybrid::Z_POS_EXISTS) != 0) &&
     ((simClasses->pargrid.getNeighbourFlags(b) & Hybrid::Z_NEG_EXISTS) != 0));
}

// return true if b is a non-edge cell at positive y outer boundary
bool InjectorSolarWind::checkIfInjectionCellYPos(const pargrid::CellID b) {
   return
    (((simClasses->pargrid.getNeighbourFlags(b) & Hybrid::Y_POS_EXISTS) == 0) &&
     ((simClasses->pargrid.getNeighbourFlags(b) & Hybrid::X_POS_EXISTS) != 0) &&
     ((simClasses->pargrid.getNeighbourFlags(b) & Hybrid::X_NEG_EXISTS) != 0) &&
     ((simClasses->pargrid.getNeighbourFlags(b) & Hybrid::Z_POS_EXISTS) != 0) &&
     ((simClasses->pargrid.getNeighbourFlags(b) & Hybrid::Z_NEG_EXISTS) != 0));
}

// return true if b is a non-edge cell at negative y outer boundary
bool InjectorSolarWind::checkIfInjectionCellYNeg(const pargrid::CellID b) {
   return
    (((simClasses->pargrid.getNeighbourFlags(b) & Hybrid::Y_NEG_EXISTS) == 0) &&
     ((simClasses->pargrid.getNeighbourFlags(b) & Hybrid::X_POS_EXISTS) != 0) &&
     ((simClasses->pargrid.getNeighbourFlags(b) & Hybrid::X_NEG_EXISTS) != 0) &&
     ((simClasses->pargrid.getNeighbourFlags(b) & Hybrid::Z_POS_EXISTS) != 0) &&
     ((simClasses->pargrid.getNeighbourFlags(b) & Hybrid::Z_NEG_EXISTS) != 0));
}

// return true if b is a non-edge cell at positive z outer boundary
bool InjectorSolarWind::checkIfInjectionCellZPos(const pargrid::CellID b) {
   return
    (((simClasses->pargrid.getNeighbourFlags(b) & Hybrid::Z_POS_EXISTS) == 0) &&
     ((simClasses->pargrid.getNeighbourFlags(b) & Hybrid::X_POS_EXISTS) != 0) &&
     ((simClasses->pargrid.getNeighbourFlags(b) & Hybrid::X_NEG_EXISTS) != 0) &&
     ((simClasses->pargrid.getNeighbourFlags(b) & Hybrid::Y_POS_EXISTS) != 0) &&
     ((simClasses->pargrid.getNeighbourFlags(b) & Hybrid::Y_NEG_EXISTS) != 0));
}

// return true if b is a non-edge cell at negative z outer boundary
bool InjectorSolarWind::checkIfInjectionCellZNeg(const pargrid::CellID b) {
   return
    (((simClasses->pargrid.getNeighbourFlags(b) & Hybrid::Z_NEG_EXISTS) == 0) &&
     ((simClasses->pargrid.getNeighbourFlags(b) & Hybrid::X_POS_EXISTS) != 0) &&
     ((simClasses->pargrid.getNeighbourFlags(b) & Hybrid::X_NEG_EXISTS) != 0) &&
     ((simClasses->pargrid.getNeighbourFlags(b) & Hybrid::Y_POS_EXISTS) != 0) &&
     ((simClasses->pargrid.getNeighbourFlags(b) & Hybrid::Y_NEG_EXISTS) != 0));
}

// default function gives an error and exits if called
void InjectorSolarWind::initParticleCrdVelDefault(Real blockSize[3],Real& x,Real& y,Real& z,Real& vx,Real& vy,Real& vz) {
   simClasses->logger << "(InjectorSolarWind): ERROR initParticleCrdVelDefault called " << endl << write;
   exit(1);
}

// initialize position and velocity of a new solar wind particle at positive x outer boundary
void InjectorSolarWind::initParticleCrdVelXPos(Real blockSize[3],Real& x,Real& y,Real& z,Real& vx,Real& vy,Real& vz) {
   x = 0;
   y = simClasses->random.uniform()*blockSize[1];
   z = simClasses->random.uniform()*blockSize[2];
   vx = -vth*derivgaussrnd(U/vth,*simClasses);
   vy = vth*gaussrnd(*simClasses);
   vz = vth*gaussrnd(*simClasses);
}

// initialize position and velocity of a new solar wind particle at negative x outer boundary
void InjectorSolarWind::initParticleCrdVelXNeg(Real blockSize[3],Real& x,Real& y,Real& z,Real& vx,Real& vy,Real& vz) {
   x = blockSize[0];
   y = simClasses->random.uniform()*blockSize[1];
   z = simClasses->random.uniform()*blockSize[2];
   vx = +vth*derivgaussrnd(U/vth,*simClasses);
   vy = vth*gaussrnd(*simClasses);
   vz = vth*gaussrnd(*simClasses);
}

// initialize position and velocity of a new solar wind particle at positive y outer boundary
void InjectorSolarWind::initParticleCrdVelYPos(Real blockSize[3],Real& x,Real& y,Real& z,Real& vx,Real& vy,Real& vz) {
   x = simClasses->random.uniform()*blockSize[0];
   y = 0;
   z = simClasses->random.uniform()*blockSize[2];
   vx = vth*gaussrnd(*simClasses);
   vy = -vth*derivgaussrnd(U/vth,*simClasses);
   vz = vth*gaussrnd(*simClasses);
}

// initialize position and velocity of a new solar wind particle at negative y outer boundary
void InjectorSolarWind::initParticleCrdVelYNeg(Real blockSize[3],Real& x,Real& y,Real& z,Real& vx,Real& vy,Real& vz) {
   x = simClasses->random.uniform()*blockSize[0];
   y = blockSize[1];
   z = simClasses->random.uniform()*blockSize[2];
   vx = vth*gaussrnd(*simClasses);
   vy = +vth*derivgaussrnd(U/vth,*simClasses);
   vz = vth*gaussrnd(*simClasses);
}

// initialize position and velocity of a new solar wind particle at positive z outer boundary
void InjectorSolarWind::initParticleCrdVelZPos(Real blockSize[3],Real& x,Real& y,Real& z,Real& vx,Real& vy,Real& vz) {
   x = simClasses->random.uniform()*blockSize[0];
   y = simClasses->random.uniform()*blockSize[1];
   z = 0;
   vx = vth*gaussrnd(*simClasses);
   vy = vth*gaussrnd(*simClasses);
   vz = -vth*derivgaussrnd(U/vth,*simClasses);
}

// initialize position and velocity of a new solar wind particle at negative z outer boundary
void InjectorSolarWind::initParticleCrdVelZNeg(Real blockSize[3],Real& x,Real& y,Real& z,Real& vx,Real& vy,Real& vz) {
   x = simClasses->random.uniform()*blockSize[0];
   y = simClasses->random.uniform()*blockSize[1];
   z = blockSize[2];
   vx = vth*gaussrnd(*simClasses);
   vy = vth*gaussrnd(*simClasses);
   vz = +vth*derivgaussrnd(U/vth,*simClasses);
}

bool InjectorSolarWind::inject(pargrid::DataID speciesDataID,unsigned int* N_particles) {
   if(initialized == false) { return initialized; }
   bool success = true;
   if(sim->timestep <= 0) { return success; }
   
   pargrid::DataWrapper<Particle<Real> > wrapper = simClasses->pargrid.getUserDataDynamic<Particle<Real> >(speciesDataID);
   
   for(pargrid::CellID b=0; b<simClasses->pargrid.getNumberOfLocalCells(); ++b) {
      if((this->*checkIfInjectionCellFuncPtr)(b) == true)  {
	 if(injectParticles(b,*species,N_particles,wrapper) == false) { success = false; }
      }
   }
   return success;
}

bool InjectorSolarWind::injectParticles(pargrid::CellID blockID,const Species& species,unsigned int* N_particles,pargrid::DataWrapper<Particle<Real> >& wrapper) {
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
   //const int N_inject = probround(*simClasses,N_macroParticlesPerCellPerDt*(1.5 + sin(20.0*sim->t/(sim->maximumTimesteps*sim->dt)))); // RHBTESTS: inject sine wave in the solar wind density from the front wall
   if(N_inject <= 0) { return true; }
   // Make room for new particles:
   const pargrid::ArraySizetype oldSize = wrapper.size()[blockID];
   N_particles[blockID] += N_inject;
   wrapper.resize(blockID,oldSize+N_inject);
   Particle<Real>* particles = wrapper.data()[blockID];
   Real xInj,yInj,zInj,vxInj,vyInj,vzInj;
   for(size_t p=oldSize; p<oldSize+N_inject; ++p) {
      // set coordinates and velocity components of a new particle
      (this->*initParticleCrdVelFuncPtr)(blockSize,xInj,yInj,zInj,vxInj,vyInj,vzInj);
      particles[p].state[particle::X] = xInj;
      particles[p].state[particle::Y] = yInj;
      particles[p].state[particle::Z] = zInj;
      particles[p].state[particle::VX] = vxInj;
      particles[p].state[particle::VY] = vyInj;
      particles[p].state[particle::VZ] = vzInj;
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
      Hybrid::logCounterParticleInject[species.popid-1] += w;
      Hybrid::logCounterParticleInjectMacroparticles[species.popid-1] += 1;
      Hybrid::logCounterParticleInjectKineticEnergy[species.popid-1] += w*( sqr(particles[p].state[particle::VX]) + sqr(particles[p].state[particle::VY]) + sqr(particles[p].state[particle::VZ]) );
   }
   return true;
}

bool InjectorSolarWind::addConfigFileItems(ConfigReader& cr,const std::string& configRegionName) {
   cr.add(configRegionName+".velocity","Bulk velocity vector (Ux,Uy,Uz) [m/s] (float,float,float).",string(""));
   cr.add(configRegionName+".density","Number density [m^-3] (float).",(Real)0.0);
   cr.add(configRegionName+".temperature","Temperature [K] (float).",(Real)0.0);
   cr.add(configRegionName+".macroparticles_per_cell","Number of macroparticles per cell in undisturbed solar wind [#] (Real).",(Real)-1.0);
   return true;
}

bool InjectorSolarWind::initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,const std::string& configRegionName,const ParticleListBase* plist) {
   initialized = ParticleInjectorBase::initialize(sim,simClasses,cr,configRegionName,plist);
   this->type = "solarwind";
   this->species = reinterpret_cast<const Species*>(plist->getSpecies());
   string velStr = "";
   cr.parse();
   cr.get(configRegionName+".velocity",velStr);
   cr.get(configRegionName+".density",n);
   cr.get(configRegionName+".temperature",T);
   cr.get(configRegionName+".macroparticles_per_cell",N_macroParticlesPerCell);

   // bulk velocity vector
     {
	bool velStrOk = true;
	vector<Real> vel;
	velStrOk = convertConfigFileVariableVelocity(velStr,vel);
	// if not correct format
	if(velStrOk == false) {
	   simClasses.logger << "(" << species->name << ") ERROR: bad format of bulk velocity vector (" << velStr << ")" << endl << write;
	   initialized = false;
	}
	// check that only one component is non-zero
	if(count(vel.begin(),vel.end(),0) != 2) {
	   simClasses.logger << "(" << species->name << ") ERROR: only bulk velocity vectors aligned or anti-aligned with coordinate axes allowed: (U,0,0), (0,U,0) or (0,0,U) where U < 0 or U > 0 (" << velStr << ")" << endl << write;
	   initialized = false;
	}
	this->velocity[0] = vel[0];
	this->velocity[1] = vel[1];
	this->velocity[2] = vel[2];
     }

   // determine the number of non-ghost cells in the plane perpendicular to undisturbed solar wind flow
   int N_perp_cells = 0;
   int N_x_cells = sim.x_blocks;
   int N_y_cells = sim.y_blocks;
   int N_z_cells = sim.z_blocks;
   if(N_x_cells > 2 && sim.x_periodic == false) { N_x_cells -= 2; } // remove ghost cells
   if(N_y_cells > 2 && sim.y_periodic == false) { N_y_cells -= 2; } // remove ghost cells
   if(N_z_cells > 2 && sim.z_periodic == false) { N_z_cells -= 2; } // remove ghost cells
   const int N_yz_cells = N_y_cells*N_z_cells*block::WIDTH_Y*block::WIDTH_Z; // NOTE: blocks != 1 do not work in RHybrid
   const int N_xz_cells = N_x_cells*N_z_cells*block::WIDTH_X*block::WIDTH_Z; // NOTE: blocks != 1 do not work in RHybrid
   const int N_xy_cells = N_x_cells*N_y_cells*block::WIDTH_X*block::WIDTH_Y; // NOTE: blocks != 1 do not work in RHybrid

   // set function pointers to find injection cells and initializing properties of new particles
   if(velocity[0] < 0) { // bulk velocity is (Ux,0,0), where Ux < 0
      // setup injection from xpos wall
      checkIfInjectionCellFuncPtr = &InjectorSolarWind::checkIfInjectionCellXPos;
      initParticleCrdVelFuncPtr = &InjectorSolarWind::initParticleCrdVelXPos;
      N_perp_cells = N_yz_cells;
      U = fabs(velocity[0]);
   }
   else if(velocity[0] > 0) { // bulk velocity is (Ux,0,0), where Ux > 0
      // setup injection from xneg wall
      checkIfInjectionCellFuncPtr = &InjectorSolarWind::checkIfInjectionCellXNeg;
      initParticleCrdVelFuncPtr = &InjectorSolarWind::initParticleCrdVelXNeg;
      N_perp_cells = N_yz_cells;
      U = fabs(velocity[0]);
   }
   else if(velocity[1] < 0) { // bulk velocity is (0,Uy,0), where Uy < 0
      // setup injection from ypos wall
      checkIfInjectionCellFuncPtr = &InjectorSolarWind::checkIfInjectionCellYPos;
      initParticleCrdVelFuncPtr = &InjectorSolarWind::initParticleCrdVelYPos;
      N_perp_cells = N_xz_cells;
      U = fabs(velocity[1]);
   }
   else if(velocity[1] > 0) { // bulk velocity is (0,Uy,0), where Ux > 0
      // setup injection from yneg wall
      checkIfInjectionCellFuncPtr = &InjectorSolarWind::checkIfInjectionCellYNeg;
      initParticleCrdVelFuncPtr = &InjectorSolarWind::initParticleCrdVelYNeg;
      N_perp_cells = N_xz_cells;
      U = fabs(velocity[1]);
   }
   else if(velocity[2] < 0) { // bulk velocity is (0,0,Uz), where Uz < 0
      // setup injection from zpos wall
      checkIfInjectionCellFuncPtr = &InjectorSolarWind::checkIfInjectionCellZPos;
      initParticleCrdVelFuncPtr = &InjectorSolarWind::initParticleCrdVelZPos;
      N_perp_cells = N_xy_cells;
      U = fabs(velocity[2]);
   }
   else if(velocity[2] > 0) { // bulk velocity is (0,0,Uz), where Uz > 0
      // setup injection from zneg wall
      checkIfInjectionCellFuncPtr = &InjectorSolarWind::checkIfInjectionCellZNeg;
      initParticleCrdVelFuncPtr = &InjectorSolarWind::initParticleCrdVelZNeg;
      N_perp_cells = N_xy_cells;
      U = fabs(velocity[2]);
   }
   else {
      simClasses.logger << "(" << species->name << ") ERROR: zero bulk velocity vector" << endl << write;
      initialized = false;
   }
   if(N_macroParticlesPerCell > 0 && n > 0 && Hybrid::dx > 0) {
      N_macroParticlesPerCellPerDt = N_macroParticlesPerCell*U*sim.dt/Hybrid::dx;
      w = n*Hybrid::dV/N_macroParticlesPerCell;
   }
   else {
      N_macroParticlesPerCellPerDt = N_macroParticlesPerCell = n = w = 0.0;
   }
   if(T > 0) { vth = sqrt(constants::BOLTZMANN*T/species->m); }
   else { vth = 0.0; }
   simClasses.logger
     << "(" << species->name << ") velocity      = (" << velocity[0]/1e3 << "," << velocity[1]/1e3 << "," << velocity[2]/1e3 << ") km/s" << endl
     << "(" << species->name << ") speed         = " << U/1e3 << " km/s" << endl
     << "(" << species->name << ") density       = " << n/1e6 << " cm^{-3}" << endl
     << "(" << species->name << ") temperature   = " << T << " K = " << T/constants::EV_TO_KELVIN << " eV" << endl
     << "(" << species->name << ") thermal speed = " << vth/1e3 << " km/s" << endl
     << "(" << species->name << ") macroparticles per cell         = " << N_macroParticlesPerCell << endl
     << "(" << species->name << ") macroparticles per dt           = " << N_macroParticlesPerCellPerDt*N_perp_cells << endl
     << "(" << species->name << ") macroparticles per cell per dt  = " << N_macroParticlesPerCellPerDt << endl
     << "(" << species->name << ") macroparticle weight            = " << w << endl;
   static unsigned int swPopCnt = 0;
   swPopCnt++;
   if(swPopCnt == 1) {
      Hybrid::swMacroParticlesCellPerDt = N_macroParticlesPerCellPerDt*N_perp_cells/N_macroParticlesPerCell;
   }
   particlePopulationInfo popInfo;
   popInfo.w = w;
   popInfo.name = species->name;
   Hybrid::allPopsInfo.push_back(popInfo);
   solarWindPopulationInfo swPopInfo;
   swPopInfo.m = species->m;
   swPopInfo.q = species->q;
   swPopInfo.U = U;
   swPopInfo.n = n;
   swPopInfo.vth = vth;
   swPopInfo.T = T;
   swPopInfo.name = species->name;
   Hybrid::swPopsInfo.push_back(swPopInfo);
   return initialized;
}

// get injector parameters
void InjectorSolarWind::getParams(InjectorParameters& p) {
   p.name = species->name;
   p.m = species->m;
   p.q = species->q;
   p.w = w;
   p.T = T;
   p.vth = vth;
   p.n = n;
   p.U = U;
   p.velocity[0] = velocity[0];
   p.velocity[1] = velocity[1];
   p.velocity[2] = velocity[2];
}

// IONOSPHERE EMISSION INJECTOR

InjectorIonosphere::InjectorIonosphere(): ParticleInjectorBase() { 
   initialized = false;
   N_ionoPop = -1;
   N_macroParticlesPerCell = -1.0;
   N_macroParticlesPerDt = -1.0;
   T = vth = w = R = 0.0;
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

bool InjectorIonosphere::injectParticles(pargrid::CellID blockID,const Species& species,unsigned int* N_particles,pargrid::DataWrapper<Particle<Real> >& wrapper) {
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
      Hybrid::logCounterParticleInject[species.popid-1] += w;
      Hybrid::logCounterParticleInjectMacroparticles[species.popid-1] += 1;
      Hybrid::logCounterParticleInjectKineticEnergy[species.popid-1] += w*( sqr(particles[p].state[particle::VX]) + sqr(particles[p].state[particle::VY]) + sqr(particles[p].state[particle::VZ]) );
      ++s;
   }
   return true;
}

bool InjectorIonosphere::addConfigFileItems(ConfigReader& cr,const std::string& configRegionName) {
   cr.add(configRegionName+".profile_name","Ionospheric emission profile name [-] (string)",string(""));
   cr.add(configRegionName+".emission_radius","Radius of the spherical emission shell [m] (Real)",(Real)-1.0);
   cr.add(configRegionName+".noon","Noon (dayside) emission factor [-] (Real)",(Real)-1.0);
   cr.add(configRegionName+".night","Night side emission factor [-] (Real)",(Real)-1.0);
   cr.add(configRegionName+".temperature","Temperature [K] (float)",(Real)0.0);
   cr.add(configRegionName+".total_production_rate","Total production rate of physical particles per second [#/s] (Real)",(Real)-1.0);
   cr.add(configRegionName+".macroparticles_per_cell","Number of macroparticles per cell wtr. to the first solar wind population [#] (Real)",(Real)-1.0);
   return true;
}

bool InjectorIonosphere::initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,const std::string& configRegionName,const ParticleListBase* plist) {
   initialized = ParticleInjectorBase::initialize(sim,simClasses,cr,configRegionName,plist);
   this->type = "ionosphere";
   this->species = reinterpret_cast<const Species*>(plist->getSpecies());
   string profileName = "";
   Real noonFactor = -1.0;
   Real nightFactor = -1.0;
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
        << "(" << species->name << ") WARNING: negative noon or night factor" << endl;
   }
   if(Hybrid::swMacroParticlesCellPerDt > 0.0) {
      N_macroParticlesPerDt = N_macroParticlesPerCell*Hybrid::swMacroParticlesCellPerDt;
   }
   else {
      simClasses.logger
        << "(" << species->name << ") WARNING: No solar wind macroparticles cell per dt rate found, assuming 100" << endl;
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
     << "(" << species->name << ") macroparticle weight      = " << w << endl;

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
	 // constant everywhere
         else if(profileName.compare("ionoConstant") == 0) {
            N_inside *= 1.0;
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
	 else {
	    simClasses.logger << "(" << species->name << ") ERROR: unknown name of an ionospheric outflow profile (" << profileName << ")" << endl << write;
	    initialized = false;
	    return initialized;
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
   particlePopulationInfo popInfo;
   popInfo.w = w;
   popInfo.name = species->name;
   Hybrid::allPopsInfo.push_back(popInfo);
   return initialized;
}

// get injector parameters
void InjectorIonosphere::getParams(InjectorParameters& p){
   p.name = species->name;
   p.m = species->m;
   p.q = species->q;
   p.w = w;
   p.T = T;
   p.vth = vth;
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

bool InjectorChapmanIonosphere::injectParticles(pargrid::CellID blockID,const Species& species,unsigned int* N_particles,pargrid::DataWrapper<Particle<Real> >& wrapper) {
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
      Hybrid::logCounterParticleInject[species.popid-1] += w;
      Hybrid::logCounterParticleInjectMacroparticles[species.popid-1] += 1;
      Hybrid::logCounterParticleInjectKineticEnergy[species.popid-1] += w*( sqr(particles[p].state[particle::VX]) + sqr(particles[p].state[particle::VY]) + sqr(particles[p].state[particle::VZ]) );
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

bool InjectorChapmanIonosphere::initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,const std::string& configRegionName,const ParticleListBase* plist) {
   initialized = ParticleInjectorBase::initialize(sim,simClasses,cr,configRegionName,plist);
   this->type = "ionosphere_chapman";
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
        << "(" << species->name << ") WARNING: No solar wind macroparticles cell per dt rate found, assuming 100" << endl;
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
     << "(" << species->name << ") macroparticle weight      = " << w << endl;

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

// get injector parameters
void InjectorChapmanIonosphere::getParams(InjectorParameters& p) {
   p.name = species->name;
   p.m = species->m;
   p.q = species->q;
   p.w = w;
   p.T = T;
   p.vth = vth;
}

// EXOSPHERE PHOTOION INJECTOR

InjectorExosphere::InjectorExosphere(): ParticleInjectorBase() {
   initialized = false;
   N_exoPop = -1;
   neutralProfileName = "";
   N_macroParticlesPerCell = -1.0;
   N_macroParticlesPerDt = -1.0;
   T = vth = w = r0 = R_exobase = R_shadow = 0.0;
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

bool InjectorExosphere::injectParticles(pargrid::CellID blockID,const Species& species,unsigned int* N_particles,pargrid::DataWrapper<Particle<Real> >& wrapper) {
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
      Hybrid::logCounterParticleInject[species.popid-1] += w;
      Hybrid::logCounterParticleInjectMacroparticles[species.popid-1] += 1;
      Hybrid::logCounterParticleInjectKineticEnergy[species.popid-1] += w*( sqr(particles[p].state[particle::VX]) + sqr(particles[p].state[particle::VY]) + sqr(particles[p].state[particle::VZ]) );
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

bool InjectorExosphere::initialize(Simulation& sim,SimulationClasses& simClasses,ConfigReader& cr,const std::string& configRegionName,const ParticleListBase* plist) {
   initialized = ParticleInjectorBase::initialize(sim,simClasses,cr,configRegionName,plist);
   this->type = "exosphere";
   this->species = reinterpret_cast<const Species*>(plist->getSpecies());
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
        << "(" << species->name << ") WARNING: No solar wind macroparticles cell per dt rate found, assuming 100" << endl;
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
     << "(" << species->name << ") macroparticle weight      = " << w << endl;
   particlePopulationInfo popInfo;
   popInfo.w = w;
   popInfo.name = species->name;
   Hybrid::allPopsInfo.push_back(popInfo);
   return initialized;
}

// get injector parameters
void InjectorExosphere::getParams(InjectorParameters& p) {
   p.name = species->name;
   p.m = species->m;
   p.q = species->q;
   p.w = w;
   p.T = T;
   p.vth = vth;
}
