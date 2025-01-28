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
   for(size_t i=1;i<Nsize;i++) {
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

inline bool setResistivityProfile(SimulationClasses& simClasses,std::string resProfileName,std::string resValueUnit,Real resValue,std::vector<Real> resSphericalValue,Real resistivityGridUnit,Real Ubulk) {
   Hybrid::resistivityProfilePtr = NULL;

   // convert eta from config file to SI units: resistivityConstant and resistivitySuperConductingSphere profiles
   if(resProfileName.compare("resistivityConstant") == 0 || resProfileName.compare("resistivitySuperConductingSphere") == 0) {
      if(resValueUnit.compare("SI") == 0) {
	 // resValue is already in SI units
	 Hybrid::resistivityEta = (resValue);
      }
      else if(resValueUnit.compare("grid") == 0){
	 // resValue is in grid units, and eta_a = (resValue) * mu0*dx^2/dt, where resValue = eta_c = dimensionless constant
	 Hybrid::resistivityEta = (resValue)*resistivityGridUnit;
      }
      else if(resValueUnit.compare("td") == 0) {
	 if(resValue == 0) {
	    simClasses.logger << "(RHYBRID) ERROR: value to define resistivity cannot be zero in units of td_per_dt (" << resValue << ")" << std::endl << write;
	    return false;
	 }
	 // resValue is diffusion time divided by dt, and eta_a = (1/resValue) * mu0*dx^2/dt, where resValue = td/dt = dimensionless constant
	 Hybrid::resistivityEta = (1.0/resValue) * resistivityGridUnit;
      }
      else if(resValueUnit.compare("Rm") == 0) {
	 if(resValue == 0) {
	    simClasses.logger << "(RHYBRID) ERROR: value to define resistivity cannot be zero in units of Rm" << std::endl << write;
	    return false;
	 }
	 // resValue is minimum magnetic Reynolds number, and eta_a = (1/resValue) * mu0*dx*Ubulk, where resValue = Rm = dimensionless constant
	 Hybrid::resistivityEta = (1.0/resValue) * constants::PERMEABILITY*Hybrid::dx*Ubulk;
      }
      else if(resValueUnit.compare("URm") == 0) {
	 // resValue is velocity divided by minimum magnetic Reynolds number, and eta_a = (resValue) * mu0*dx, where resValue = U/Rm = [m/s]
	 Hybrid::resistivityEta = (resValue) * constants::PERMEABILITY*Hybrid::dx;
      }
      else {
	 simClasses.logger << "(RHYBRID) ERROR: unknown unit and quantity to define resistivity (" << resValueUnit << ")" << std::endl << write;
	 return false;
      }
   }
   // convert eta from config file to SI units: resistivitySphericalShells profile
   else if(resProfileName.compare("resistivitySphericalShells") == 0) {
      if(resValueUnit.compare("SI") == 0) {
	 // resValue is already in SI units
	 for(size_t i=0;i<resSphericalValue.size();i++) {
	    Hybrid::resistivitySphericalEta.push_back( (resSphericalValue[i]) );
	 }
      }
      else if(resValueUnit.compare("grid") == 0){
	 // resValue is in grid units, and eta_a = (resValue) * mu0*dx^2/dt, where resValue = eta_c = dimensionless constant
	 for(size_t i=0;i<resSphericalValue.size();i++) {
	    Hybrid::resistivitySphericalEta.push_back( (resSphericalValue[i])*resistivityGridUnit );
	 }
      }
      else if(resValueUnit.compare("td") == 0) {
	 // resValue is diffusion time divided by dt, and eta_a = (1/resValue) * mu0*dx^2/dt, where resValue = td/dt = dimensionless constant
	 for(size_t i=0;i<resSphericalValue.size();i++) {
	    if(resSphericalValue[i] == 0) {
	       simClasses.logger << "(RHYBRID) ERROR: value_spherical to defined resistivity cannot be zero in units of td_per_dt" << std::endl << write;
	       return false;
	    }
	    Hybrid::resistivitySphericalEta.push_back( (1.0/resSphericalValue[i]) * resistivityGridUnit );
	 }
      }
      else if(resValueUnit.compare("Rm") == 0) {
	 // resValue is minimum magnetic Reynolds number, and eta_a = (1/resValue) * mu0*dx*Ubulk, where resValue = Rm = dimensionless constant
	 for(size_t i=0;i<resSphericalValue.size();i++) {
	    if(resSphericalValue[i] == 0) {
	       simClasses.logger << "(RHYBRID) ERROR: value_spherical to defined resistivity cannot be zero in units of Rm" << std::endl << write;
	       return false;
	    }
	    Hybrid::resistivitySphericalEta.push_back( (1.0/resSphericalValue[i]) * constants::PERMEABILITY*Hybrid::dx*Ubulk );
	 }
      }
      else if(resValueUnit.compare("URm") == 0) {
	 // resValue is velocity divided by minimum magnetic Reynolds number, and eta_a = (resValue) * mu0*dx, where resValue = U/Rm = [m/s]
	 for(size_t i=0;i<resSphericalValue.size();i++) {
	    Hybrid::resistivitySphericalEta.push_back( (resSphericalValue[i]) * constants::PERMEABILITY*Hybrid::dx );
	 }
      }
      else {
	 simClasses.logger << "(RHYBRID) ERROR: unknown unit and quantity to define resistivity (" << resValueUnit << ")" << std::endl << write;
	 return false;
      }
   }
   else {
      simClasses.logger << "(setResistivityProfile) ERROR: unknown name of a resistivity profile (" << resProfileName << ")" << std::endl << write;
      return false;
   }

   if(resProfileName.compare("resistivityConstant") == 0) {
      Hybrid::resistivityProfilePtr = &resistivityConstant;
      // check if config file parameters given not used by this profile
      if(Hybrid::resistivitySphericalR2.size() > 0) {
	 simClasses.logger << "(setResistivityProfile) WARNING: when resistivityConstant profile is used, parameters of resistivitySphericalShells (value_spherical, R_spherical) are ignored" << std::endl;
      }
      if(Hybrid::resistivityR2 > 0) {
	 simClasses.logger << "(setResistivityProfile) WARNING: when resistivityConstant profile is used, R parameter is ignored" << std::endl;
      }
      if(Hybrid::resistivityEta <= 0) {
	 simClasses.logger << "(setResistivityProfile) WARNING: eta <= 0 in resistivityConstant" << std::endl;
      }
   }
   else if(resProfileName.compare("resistivitySuperConductingSphere") == 0){
      Hybrid::resistivityProfilePtr = &resistivitySuperConductingSphere;
      // check if config file parameters given not used by this profile
      if(Hybrid::resistivitySphericalR2.size() > 0) {
	 simClasses.logger << "(setResistivityProfile) WARNING: when resistivitySuperConductingSphere profile is used, parameters of resistivitySphericalShells (value_spherical, R_spherical) are ignored" << std::endl;
      }
      if(Hybrid::resistivityR2 <= 0) {
	 simClasses.logger << "(setResistivityProfile) WARNING: R <= 0 in resistivitySuperConductingSphere" << std::endl;
      }
      if(Hybrid::resistivityEta <= 0) {
	 simClasses.logger << "(setResistivityProfile) WARNING: eta <= 0 in resistivitySuperConductingSphere" << std::endl;
      }
      // set squared radius of super conducting sphere
      Hybrid::resistivityR2 = sqr(Hybrid::resistivityR2);
   }
   else if(resProfileName.compare("resistivitySphericalShells") == 0){
      Hybrid::resistivityProfilePtr = &resistivitySphericalShells;
      // check if config file parameters given not used by this profile
      if(Hybrid::resistivityR2 > 0 || Hybrid::resistivityEta > 0) {
	 simClasses.logger << "(setResistivityProfile) WARNING: when resistivitySphericalShells profile is used, value and R parameters are ignored" << std::endl;
      }
      // check and calculate resistivity radii parameters
      if( ( resSphericalValue.size() != Hybrid::resistivitySphericalR2.size() ) || (resSphericalValue.size() < 1) || (Hybrid::resistivitySphericalR2.size() < 1) ) {
	 simClasses.logger << "(RHYBRID) ERROR: parameter arrays of the spherical shell resistivity model (value_spherical, R_spherical) should be the same non-zero size (" << resSphericalValue.size() << ", " << Hybrid::resistivitySphericalR2.size() << ")" << std::endl << write;
	 return false;;
      }
      for(size_t i=0;i<Hybrid::resistivitySphericalR2.size();i++) {
	 if(Hybrid::resistivitySphericalR2[i] < 0) {
	    simClasses.logger << "(RHYBRID) ERROR: resistivity R_spherical < 0 (" << Hybrid::resistivitySphericalR2[i] << ")" << std::endl << write;
	    return false;
	 }
	 Hybrid::resistivitySphericalR2[i] = sqr(Hybrid::resistivitySphericalR2[i]);
	 // check that the radii are given in monotonically growing order
	 if(i > 0) {
	    if(Hybrid::resistivitySphericalR2[i-1] >= Hybrid::resistivitySphericalR2[i]) {
	       simClasses.logger << "(RHYBRID) ERROR: radii parameters of the spherical shell resistivity model should be given in a monotonically growing order" << std::endl << write;
	       return false;
	    }
	 }
      }
   }
   else {
      simClasses.logger << "(setResistivityProfile) ERROR: unknown name of a resistivity profile (" << resProfileName << ")" << std::endl << write;
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
      if(x < (sim.x_min + bZone) || x > (sim.x_max - bZone) ||
         y < (sim.y_min + bZone) || y > (sim.y_max - bZone) ||
         z < (sim.z_min + bZone) || z > (sim.z_max - bZone)) {
         res += Hybrid::outerBoundaryZone.eta;
      }
   }
   else if(Hybrid::outerBoundaryZone.typeEta == 2) { // all edges except +x
      if( (x < (sim.x_min + bZone)) && (y < (sim.y_min + bZone)) ) {
         // (-x,-y) edge
         res += Hybrid::outerBoundaryZone.eta;
      }
      else if( (x < (sim.x_min + bZone)) && (y > (sim.y_max - bZone)) ) {
         // (-x,+y) edge
         res += Hybrid::outerBoundaryZone.eta;
      }
      else if( (x < (sim.x_min + bZone)) && (z < (sim.z_min + bZone)) ) {
         // (-x,-z) edge
         res += Hybrid::outerBoundaryZone.eta;
      }
      else if( (x < (sim.x_min + bZone)) && (z > (sim.z_max - bZone)) ) {
         // (-x,+z) edge
         res += Hybrid::outerBoundaryZone.eta;
      }
      else if( (y < (sim.y_min + bZone)) && (z < (sim.z_min + bZone)) ) {
         // (-y,-z) edge
         res += Hybrid::outerBoundaryZone.eta;
      }
      else if( (y < (sim.y_min + bZone)) && (z > (sim.z_max - bZone)) ) {
         // (-y,+z) edge
         res += Hybrid::outerBoundaryZone.eta;
      }
      else if( (y > (sim.y_max - bZone)) && (z < (sim.z_min + bZone)) ) {
         // (+y,-z) edge
         res += Hybrid::outerBoundaryZone.eta;
      }
      else if( (y > (sim.y_max - bZone)) && (z > (sim.z_max - bZone)) ) {
         // (+y,+z) edge
         res += Hybrid::outerBoundaryZone.eta;
      }
   }
   else if(Hybrid::outerBoundaryZone.typeEta == 3) { // -x wall and all edges except +x
      if( x < (sim.x_min + bZone) ) {
         // -x wall
         res += Hybrid::outerBoundaryZone.eta;
      }
      else if( (x < (sim.x_min + bZone)) && (y < (sim.y_min + bZone)) ) {
         // (-x,-y) edge
         res += Hybrid::outerBoundaryZone.eta;
      }
      else if( (x < (sim.x_min + bZone)) && (y > (sim.y_max - bZone)) ) {
         // (-x,+y) edge
         res += Hybrid::outerBoundaryZone.eta;
      }
      else if( (x < (sim.x_min + bZone)) && (z < (sim.z_min + bZone)) ) {
         // (-x,-z) edge
         res += Hybrid::outerBoundaryZone.eta;
      }
      else if( (x < (sim.x_min + bZone)) && (z > (sim.z_max - bZone)) ) {
         // (-x,+z) edge
         res += Hybrid::outerBoundaryZone.eta;
      }
      else if( (y < (sim.y_min + bZone)) && (z < (sim.z_min + bZone)) ) {
         // (-y,-z) edge
         res += Hybrid::outerBoundaryZone.eta;
      }
      else if( (y < (sim.y_min + bZone)) && (z > (sim.z_max - bZone)) ) {
         // (-y,+z) edge
         res += Hybrid::outerBoundaryZone.eta;
      }
      else if( (y > (sim.y_max - bZone)) && (z < (sim.z_min + bZone)) ) {
         // (+y,-z) edge
         res += Hybrid::outerBoundaryZone.eta;
      }
      else if( (y > (sim.y_max - bZone)) && (z > (sim.z_max - bZone)) ) {
         // (+y,+z) edge
         res += Hybrid::outerBoundaryZone.eta;
      }
   }
   else if(Hybrid::outerBoundaryZone.typeEta == 4) { // all walls except -x
      if(                                  x > (sim.x_max - bZone) ||
         y < (sim.y_min + bZone) || y > (sim.y_max - bZone) ||
         z < (sim.z_min + bZone) || z > (sim.z_max - bZone)) {
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
