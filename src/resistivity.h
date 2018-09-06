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

#ifndef RESISTIVITY_H
#define RESISTIVITY_H

#ifdef USE_RESISTIVITY

#include <string>
#include <vector>
#include <simulation.h>
#include <simulationclasses.h>
#include "hybrid.h"

Real resistivityConstant(Simulation& sim,SimulationClasses& simClasses,const Real x,const Real y,const Real z) {
   return Hybrid::resistivityEta;
}

Real resistivitySuperConductingSphere(Simulation& sim,SimulationClasses& simClasses,const Real x,const Real y,const Real z) {
   const Real r2 = sqr(x) + sqr(y) + sqr(z);
   if(r2 < Hybrid::resistivityR2) { return 0.0; }
   return Hybrid::resistivityEta;
}

inline bool setResistivityProfile(std::string name) {
   Hybrid::resistivityProfilePtr = NULL;
   if(name.compare("resistivityConstant") == 0) {
      Hybrid::resistivityProfilePtr = &resistivityConstant;
   }
   else if(name.compare("resistivitySuperConductingSphere") == 0){
      Hybrid::resistivityProfilePtr = &resistivitySuperConductingSphere;
   }
   else {
      return false;
   }
   return true;
}

Real getResistivity(Simulation& sim,SimulationClasses& simClasses,const Real x,const Real y,const Real z) {
   Real res = Hybrid::resistivityProfilePtr(sim,simClasses,x,y,z);
   const Real bZone = Hybrid::outerBoundaryZone.sizeEta;
   if(Hybrid::outerBoundaryZone.typeEta == 1) { // all walls
      if(x < (Hybrid::box.xmin + bZone) || x > (Hybrid::box.xmax - bZone) ||
         y < (Hybrid::box.ymin + bZone) || y > (Hybrid::box.ymax - bZone) ||
         z < (Hybrid::box.zmin + bZone) || z > (Hybrid::box.zmax - bZone)) {
         res += Hybrid::outerBoundaryZone.eta;
      }
   }
   else if(Hybrid::outerBoundaryZone.typeEta == 2) { // all edges except +x
      if( (x < (Hybrid::box.xmin + bZone)) && (y < (Hybrid::box.ymin + bZone)) ) {
         // (-x,-y) edge
         res += Hybrid::outerBoundaryZone.eta;
      }
      else if( (x < (Hybrid::box.xmin + bZone)) && (y > (Hybrid::box.ymax - bZone)) ) {
         // (-x,+y) edge
         res += Hybrid::outerBoundaryZone.eta;
      }
      else if( (x < (Hybrid::box.xmin + bZone)) && (z < (Hybrid::box.zmin + bZone)) ) {
         // (-x,-z) edge
         res += Hybrid::outerBoundaryZone.eta;
      }
      else if( (x < (Hybrid::box.xmin + bZone)) && (z > (Hybrid::box.zmax - bZone)) ) {
         // (-x,+z) edge
         res += Hybrid::outerBoundaryZone.eta;
      }
      else if( (y < (Hybrid::box.ymin + bZone)) && (z < (Hybrid::box.zmin + bZone)) ) {
         // (-y,-z) edge
         res += Hybrid::outerBoundaryZone.eta;
      }
      else if( (y < (Hybrid::box.ymin + bZone)) && (z > (Hybrid::box.zmax - bZone)) ) {
         // (-y,+z) edge
         res += Hybrid::outerBoundaryZone.eta;
      }
      else if( (y > (Hybrid::box.ymax - bZone)) && (z < (Hybrid::box.zmin + bZone)) ) {
         // (+y,-z) edge
         res += Hybrid::outerBoundaryZone.eta;
      }
      else if( (y > (Hybrid::box.ymax - bZone)) && (z > (Hybrid::box.zmax - bZone)) ) {
         // (+y,+z) edge
         res += Hybrid::outerBoundaryZone.eta;
      }
   }
   else if(Hybrid::outerBoundaryZone.typeEta == 3) { // -x wall and all edges except +x
      if( x < (Hybrid::box.xmin + bZone) ) {
         // -x wall
         res += Hybrid::outerBoundaryZone.eta;
      }
      else if( (x < (Hybrid::box.xmin + bZone)) && (y < (Hybrid::box.ymin + bZone)) ) {
         // (-x,-y) edge
         res += Hybrid::outerBoundaryZone.eta;
      }
      else if( (x < (Hybrid::box.xmin + bZone)) && (y > (Hybrid::box.ymax - bZone)) ) {
         // (-x,+y) edge
         res += Hybrid::outerBoundaryZone.eta;
      }
      else if( (x < (Hybrid::box.xmin + bZone)) && (z < (Hybrid::box.zmin + bZone)) ) {
         // (-x,-z) edge
         res += Hybrid::outerBoundaryZone.eta;
      }
      else if( (x < (Hybrid::box.xmin + bZone)) && (z > (Hybrid::box.zmax - bZone)) ) {
         // (-x,+z) edge
         res += Hybrid::outerBoundaryZone.eta;
      }
      else if( (y < (Hybrid::box.ymin + bZone)) && (z < (Hybrid::box.zmin + bZone)) ) {
         // (-y,-z) edge
         res += Hybrid::outerBoundaryZone.eta;
      }
      else if( (y < (Hybrid::box.ymin + bZone)) && (z > (Hybrid::box.zmax - bZone)) ) {
         // (-y,+z) edge
         res += Hybrid::outerBoundaryZone.eta;
      }
      else if( (y > (Hybrid::box.ymax - bZone)) && (z < (Hybrid::box.zmin + bZone)) ) {
         // (+y,-z) edge
         res += Hybrid::outerBoundaryZone.eta;
      }
      else if( (y > (Hybrid::box.ymax - bZone)) && (z > (Hybrid::box.zmax - bZone)) ) {
         // (+y,+z) edge
         res += Hybrid::outerBoundaryZone.eta;
      }
   }
   return res;
}

#endif

#endif
