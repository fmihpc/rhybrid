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

#include <string>
#include <vector>
#include <simulation.h>
#include <simulationclasses.h>
#include "hybrid.h"

struct BackgroundChargeDensityArgs{
   Real R,r0,n0;
};

Real smoothObstacle(SimulationClasses& simClasses,Real x,Real y,Real z,Real R,Real r0,Real n0) {
   const Real r = sqrt(sqr(x) + sqr(y) + sqr(z));
   if(r < R) {
      return n0;
   }
   else {
      return n0*exp(-(r-R)/r0);
   }
}

Real getBackgroundChargeDensity(SimulationClasses& simClasses,std::string name,Real x,Real y,Real z,BackgroundChargeDensityArgs a) {
   if(name.compare("smoothObstacle") == 0) {
      return smoothObstacle(simClasses,x,y,z,a.R,a.r0,a.n0);
   }
   else if(name.compare("none") == 0) {
      return 0.0;
   }
   return 0.0;
}

#endif
