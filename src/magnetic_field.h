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

#ifndef MAGNETIC_FIELD_H
#define MAGNETIC_FIELD_H

#include "hybrid.h"

#if defined(USE_B_INITIAL) || defined(USE_B_CONSTANT)

inline void constantBx(Real x,Real y,Real z,Real B[3]) {
   B[0] += Hybrid::IMFBx;
}

inline void laminarFlowAroundSphereBx(Real x,Real y,Real z,Real B[3]) {
   const Real r2 = sqr(x) + sqr(y) + sqr(z);
   if(r2 < Hybrid::laminarR2) { return; }
   const Real r = sqrt(r2);
   const Real r3 = pow (r,3);
   const Real r5 = pow (r,5);
   const Real coeff = -1.5*Hybrid::IMFBx*Hybrid::laminarR3*x/r5;
   B[0] += coeff*x + Hybrid::IMFBx*( 1 + 0.5 * (Hybrid::laminarR3/r3) );
   B[1] += coeff*y;
   B[2] += coeff*z;
}

inline void hemisphericDipoleB(Real x,Real y,Real z,Real B[3]) {
   // translation
   const Real x1 = x - Hybrid::xDip;
   const Real y1 = y - Hybrid::yDip;
   const Real z1 = z - Hybrid::zDip;
   
   //Rotation matrix R_y(theta) (=theta around y):
   // x'      cos theta     0   sin theta   x
   // y'   =      0         1      0        y
   // z'      -sin theta    0   cos theta   z
   //Rotation matrix R_x(phi) (=phi around x):
   // x'          1         0      0       x
   // y'   =      0   cos phi   -sin phi   y
   // z'          0   sin phi    cos phi   z
   
   // (x2,y2,z2) = R_y(theta) (x1,y1,z1)   
   const Real x2 = x1*cos(Hybrid::thetaDip*M_PI/180.0) + z1*sin(Hybrid::thetaDip*M_PI/180.0);
   const Real y2 = y1;
   const Real z2 = -x1*sin(Hybrid::thetaDip*M_PI/180.0) + z1*cos(Hybrid::thetaDip*M_PI/180.0);
   // (x3,y3,z3) = R_x(phi) (x2,y2,z2)
   const Real x3 = x2;
   const Real y3 = y2*cos(Hybrid::phiDip*M_PI/180.0) - z2*sin(Hybrid::phiDip*M_PI/180.0);
   const Real z3 = y2*sin(Hybrid::phiDip*M_PI/180.0) + z2*cos(Hybrid::phiDip*M_PI/180.0);
   const Real rr = sqr(x3) + sqr(y3) + sqr(z3);
   if(rr < Hybrid::dipMinR2) { return; }
   const Real r = sqrt(rr);
   const Real r5 = sqr(rr)*r;
   const Real coeff = Hybrid::dipMomCoeff/r5;
   // B(x3,y3,z3)
   const Real Bx1 = coeff*x3*(Hybrid::coeffDip*z3 + 0.5*Hybrid::coeffQuad*Hybrid::dipSurfR*(5*sqr(z3)/rr - 1));
   const Real By1 = coeff*y3*(Hybrid::coeffDip*z3 + 0.5*Hybrid::coeffQuad*Hybrid::dipSurfR*(5*sqr(z3)/rr - 1));
   const Real Bz1 = coeff*(Hybrid::coeffDip*(sqr(z3) - rr/3.0) + 0.5*Hybrid::coeffQuad*Hybrid::dipSurfR*z*(5*sqr(z3) - 3*rr)/rr);
   // B2 = R_x(-phi) B1
   const Real Bx2 = Bx1;
   const Real By2 =  By1*cos(Hybrid::phiDip*M_PI/180.0) + Bz1*sin(Hybrid::phiDip*M_PI/180.0);
   const Real Bz2 = -By1*sin(Hybrid::phiDip*M_PI/180.0) + Bz1*cos(Hybrid::phiDip*M_PI/180.0);
   // B = R_y(-theta) B2
   B[0] += Bx2*cos(Hybrid::thetaDip*M_PI/180.0) - Bz2*sin(Hybrid::thetaDip*M_PI/180.0);
   B[1] += By2;
   B[2] += Bx2*sin(Hybrid::thetaDip*M_PI/180.0) + Bz2*cos(Hybrid::thetaDip*M_PI/180.0);
}

inline void translateDipoleB(Real x,Real y,Real z,Real B[3]) {
   // translation
   const Real x1 = x - Hybrid::xDip;
   const Real y1 = y - Hybrid::yDip;
   const Real z1 = z - Hybrid::zDip;
   const Real r2 = sqr(x1) + sqr(y1) + sqr(z1);
   if(r2 < Hybrid::dipMinR2) { return; }
   const Real rr = sqrt(r2);
   const Real r5 = sqr(r2)*rr;
   const Real coeff = Hybrid::dipMomCoeff/r5;
   B[0] += coeff*x1*z1;
   B[1] += coeff*y1*z1;
   B[2] += coeff*(sqr(z1) - r2/3.0);
}

inline void generalDipoleB(Real x,Real y,Real z,Real B[3]) {
   // translation
   const Real x1 = x - Hybrid::xDip;
   const Real y1 = y - Hybrid::yDip;
   const Real z1 = z - Hybrid::zDip;
   
   //Rotation matrix R_y(theta) (=theta around y):
   // x'      cos theta     0   sin theta   x
   // y'   =      0         1      0        y
   // z'      -sin theta    0   cos theta   z
   //Rotation matrix R_x(phi) (=phi around x):
   // x'          1         0      0       x
   // y'   =      0   cos phi   -sin phi   y
   // z'          0   sin phi    cos phi   z
   
   // (x2,y2,z2) = R_y(theta) (x1,y1,z1)   
   const Real x2 = x1*cos(Hybrid::thetaDip*M_PI/180.0) + z1*sin(Hybrid::thetaDip*M_PI/180.0);
   const Real y2 = y1;
   const Real z2 = -x1*sin(Hybrid::thetaDip*M_PI/180.0) + z1*cos(Hybrid::thetaDip*M_PI/180.0);
   // (x3,y3,z3) = R_x(phi) (x2,y2,z2)
   const Real x3 = x2;
   const Real y3 = y2*cos(Hybrid::phiDip*M_PI/180.0) - z2*sin(Hybrid::phiDip*M_PI/180.0);
   const Real z3 = y2*sin(Hybrid::phiDip*M_PI/180.0) + z2*cos(Hybrid::phiDip*M_PI/180.0);
   const Real rr = sqr(x3) + sqr(y3) + sqr(z3);
   if(rr < Hybrid::dipMinR2) { return; }
   const Real r = sqrt(rr);
   const Real r5 = sqr(rr)*r;
   const Real coeff = Hybrid::dipMomCoeff/r5;
   // B(x3,y3,z3)
   const Real Bx1 = coeff*x3*z3;
   const Real By1 = coeff*y3*z3;
   const Real Bz1 = coeff*(sqr(z3) - rr/3.0);
   // B2 = R_x(-phi) B1
   const Real Bx2 = Bx1;
   const Real By2 =  By1*cos(Hybrid::phiDip*M_PI/180.0) + Bz1*sin(Hybrid::phiDip*M_PI/180.0);
   const Real Bz2 = -By1*sin(Hybrid::phiDip*M_PI/180.0) + Bz1*cos(Hybrid::phiDip*M_PI/180.0);
   // B = R_y(-theta) B2
   B[0] += Bx2*cos(Hybrid::thetaDip*M_PI/180.0) - Bz2*sin(Hybrid::thetaDip*M_PI/180.0);
   B[1] += By2;
   B[2] += Bx2*sin(Hybrid::thetaDip*M_PI/180.0) + Bz2*cos(Hybrid::thetaDip*M_PI/180.0);
}

inline void lineDipoleB(Real x,Real y,Real z,Real B[3]) {
   const Real x1 = x - Hybrid::xDip;
   const Real y1 = y - Hybrid::yDip;
   const Real z1 = z - Hybrid::zDip;
   const Real r2 = sqr(x1) + sqr(y1) + sqr(z1);
   if(r2 < Hybrid::dipMinR2) { return; }
   const Real D = -Hybrid::coeffDip;
   B[0] += D*2*x*z/(sqr(r2));
   B[1] += 0.0;
   B[2] += D*( sqr(z) - sqr(x) )/( sqr(r2) );
}

inline void translateDipoleBAndLaminarFlowAroundSphereBx(Real x,Real y,Real z,Real B[3]) {
   translateDipoleB(x,y,z,B);
   laminarFlowAroundSphereBx(x,y,z,B);
}

inline void generalDipoleBAndConstantBx(Real x,Real y,Real z,Real B[3]) {
   generalDipoleB(x,y,z,B);
   constantBx(x,y,z,B);
}

inline bool setMagneticFieldProfile(std::string name) {
   Hybrid::magneticFieldProfilePtr = NULL;
   if(name.compare("constantBx") == 0) {
      Hybrid::magneticFieldProfilePtr = &constantBx;
   }
   else if(name.compare("laminarFlowAroundSphereBx") == 0) {
      Hybrid::magneticFieldProfilePtr = &laminarFlowAroundSphereBx;
   }
   else if(name.compare("hemisphericDipoleB") == 0) {
      Hybrid::magneticFieldProfilePtr = &hemisphericDipoleB;
   }
   else if(name.compare("translateDipoleB") == 0) {
      Hybrid::magneticFieldProfilePtr = &translateDipoleB;
   }
   else if(name.compare("generalDipoleB") == 0) {
      Hybrid::magneticFieldProfilePtr = &generalDipoleB;
   }
   else if(name.compare("lineDipoleB") == 0) {
      Hybrid::magneticFieldProfilePtr = &lineDipoleB;
   }
   else if(name.compare("translateDipoleBAndLaminarFlowAroundSphereBx") == 0) {
      Hybrid::magneticFieldProfilePtr = &translateDipoleBAndLaminarFlowAroundSphereBx;
   }
   else if(name.compare("generalDipoleBAndConstantBx") == 0) {
      Hybrid::magneticFieldProfilePtr = &generalDipoleBAndConstantBx;
   }
   else {
      return false;
   }
   return true;
}

#endif

#ifdef USE_B_INITIAL
inline void setInitialB(const Real x,const Real y,const Real z,Real B[3]) {
   B[0] = B[1] = B[2] = 0.0;
   Hybrid::magneticFieldProfilePtr(x,y,z,B);
}
#endif

#ifdef USE_B_CONSTANT
inline void addConstantB(const Real x,const Real y,const Real z,Real B[3]) {
   Hybrid::magneticFieldProfilePtr(x,y,z,B);
}
#endif

#endif
