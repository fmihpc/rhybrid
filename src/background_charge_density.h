/** This file is part of the RHybrid simulation.
 *
 *  Copyright 2019- Aalto University
 *  Copyright 2019- Finnish Meteorological Institute
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

#ifndef BACKGROUND_CHARGE_DENSITY_H
#define BACKGROUND_CHARGE_DENSITY_H

#ifdef USE_BACKGROUND_CHARGE_DENSITY

#include <string>
#include <vector>
#include <simulation.h>
#include <simulationclasses.h>
#include "hybrid.h"

struct BackgroundChargeDensityArgs{
   Real R,r0,rhoQi0;
};

// spherical background ion charge density: rhoQi0 at r < R and zero at r >= R
Real bgChargeDensitySphere(SimulationClasses& simClasses,Real x,Real y,Real z,Real R,Real rhoQi0) {
   const Real r = sqrt(sqr(x) + sqr(y) + sqr(z));
   if(r < R) { return rhoQi0; }
   else { return 0.0; }
}

// spherical background ion charge density: rhoQi0 at r < R and rhoQi0*exp(-(r-R)/r0) at r >= R
Real bgChargeDensitySphereSmooth(SimulationClasses& simClasses,Real x,Real y,Real z,Real R,Real r0,Real rhoQi0) {
   const Real r = sqrt(sqr(x) + sqr(y) + sqr(z));
   if(r < R) { return rhoQi0; }
   else { return rhoQi0*exp(-(r-R)/r0); }
}

Real getBackgroundChargeDensity(SimulationClasses& simClasses,std::string name,Real x,Real y,Real z,BackgroundChargeDensityArgs a) {
   if(name.compare("bgChargeDensitySphere") == 0) {
      return bgChargeDensitySphere(simClasses,x,y,z,a.R,a.rhoQi0);
   }
   else if(name.compare("bgChargeDensitySphereSmooth") == 0) {
      return bgChargeDensitySphereSmooth(simClasses,x,y,z,a.R,a.r0,a.rhoQi0);
   }
   else if(name.compare("none") == 0) {
      return 0.0;
   }
   else {
      simClasses.logger << "(getBackgroundChargeDensity) ERROR: unknown name of a background charge density profile (" << name << ")" << std::endl << write;
      MPI_Finalize();
   }
   return 0.0;
}

#endif

#endif
