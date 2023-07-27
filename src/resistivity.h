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

Real resistivitySphericalShells(Simulation& sim,SimulationClasses& simClasses,const Real x,const Real y,const Real z) {
   const Real r2 = sqr(x) + sqr(y) + sqr(z);
   const size_t Nsize = Hybrid::resistivitySphericalR2.size();
   if(Nsize <= 0) {
      simClasses.logger << "(resistivitySphericalShells) ERROR: no function parameters given" << std::endl << write;
      MPI_Finalize();
      return 0.0;
   }
   // check if the point is inside or at the first radii: return resistivity of the first shell
   if(r2 <= Hybrid::resistivitySphericalR2[0]) { return Hybrid::resistivitySphericalEta[0]; }
   // loop thru from second to last shell
   for(size_t i=1;i<Hybrid::resistivitySphericalR2.size();i++) {
      // check if the point is between i-1 and i radii
      if(r2 > Hybrid::resistivitySphericalR2[i-1] && r2 <= Hybrid::resistivitySphericalR2[i]) {
	 return Hybrid::resistivitySphericalEta[i];
      }
   }
   // check if the point is above or at the last radii: return resistivity of the last shell
   //if(r2 >= Hybrid::resistivitySphericalR2[Nsize-1]) { return Hybrid::resistivitySphericalEta[Nsize-1]; }
   //simClasses.logger << "(resistivitySphericalShells) ERROR: reached function end, which should never happen" << std::endl << write;
   //MPI_Finalize();
   return 0.0;
}

inline bool setResistivityProfile(std::string name,SimulationClasses& simClasses) {
   Hybrid::resistivityProfilePtr = NULL;
   if(name.compare("resistivityConstant") == 0) {
      Hybrid::resistivityProfilePtr = &resistivityConstant;
      if(Hybrid::resistivitySphericalR2.size() > 0) {
	 simClasses.logger << "(setResistivityProfile) WARNING: when resistivityConstant profile is used, parameters of resistivitySphericalShells are ignored" << std::endl << write;
      }
      if(Hybrid::resistivityR2 > 0) {
	 simClasses.logger << "(setResistivityProfile) WARNING: when resistivityConstant profile is used, R parameter is ignored" << std::endl << write;
      }
      if(Hybrid::resistivityEta <= 0) {
	 simClasses.logger << "(setResistivityProfile) WARNING: eta <= 0 in resistivityConstant" << std::endl << write;
      }
   }
   else if(name.compare("resistivitySuperConductingSphere") == 0){
      Hybrid::resistivityProfilePtr = &resistivitySuperConductingSphere;
      if(Hybrid::resistivitySphericalR2.size() > 0) {
	 simClasses.logger << "(setResistivityProfile) WARNING: when resistivitySuperConductingSphere profile is used, parameters of resistivitySphericalShells are ignored" << std::endl << write;
      }
      if(Hybrid::resistivityR2 <= 0) {
	 simClasses.logger << "(setResistivityProfile) WARNING: R <= 0 in resistivitySuperConductingSphere" << std::endl << write;
      }
      if(Hybrid::resistivityEta <= 0) {
	 simClasses.logger << "(setResistivityProfile) WARNING: eta <= 0 in resistivitySuperConductingSphere" << std::endl << write;
      }
   }
   else if(name.compare("resistivitySphericalShells") == 0){
      Hybrid::resistivityProfilePtr = &resistivitySphericalShells;
      if(Hybrid::resistivityR2 > 0 || Hybrid::resistivityEta > 0) {
	 simClasses.logger << "(setResistivityProfile) WARNING: when resistivitySphericalShells profile is used, eta and R parameters are ignored" << std::endl << write;
      }
   }
   else {
      simClasses.logger << "(setResistivityProfile) ERROR: unknown name of a resistivity profile (" << name << ")" << std::endl << write;
      return false;
   }
   return true;
}

Real getResistivity(Simulation& sim,SimulationClasses& simClasses,const Real x,const Real y,const Real z) {
   Real res = Hybrid::resistivityProfilePtr(sim,simClasses,x,y,z);
#ifdef USE_OUTER_BOUNDARY_ZONE
   const Real bZone = Hybrid::outerBoundaryZone.sizeEta;
   if (Hybrid::outerBoundaryZone.typeEta == 0) { // not used
      res += 0.0;
   }
   else if(Hybrid::outerBoundaryZone.typeEta == 1) { // all walls
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
   else if(Hybrid::outerBoundaryZone.typeEta == 4) { // all walls except -x
      if(                                  x > (Hybrid::box.xmax - bZone) ||
         y < (Hybrid::box.ymin + bZone) || y > (Hybrid::box.ymax - bZone) ||
         z < (Hybrid::box.zmin + bZone) || z > (Hybrid::box.zmax - bZone)) {
         res += Hybrid::outerBoundaryZone.eta;
      }
   }
   else {
      simClasses.logger << "(getResistivity) ERROR: unknown type of an outer boundary zone for eta (" << Hybrid::outerBoundaryZone.typeEta << ")" << std::endl << write;
      MPI_Finalize();
   }
#endif
   return res;
}

#endif

#endif
