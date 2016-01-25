/** This file is part of the RHybrid simulation.
 *
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

#ifndef HYBRID_H
#define HYBRID_H

#include <cstdlib>
#include <vector>

#include <simulationclasses.h>

#define sqr(x) ((x)*(x))
#define cube(x) (x)*(x)*(x)
#define vecsqr(a) (sqr(a[0])+sqr(a[1])+sqr(a[2]))
#define normvec(a) (sqrt(vecsqr(a)))

#ifdef ION_SPECTRA_ALONG_ORBIT
#define EBINS 10
struct Dist {
   Real f[EBINS];
 public:
   void reset() {
      for(unsigned int i=0;i<EBINS;i++){
         f[i] = 0.0;
      }
   }
};
#endif

inline void cross(const Real a[3], const Real b[3], Real result[3]) {
   result[0] = a[1]*b[2] - a[2]*b[1];
   result[1] = a[2]*b[0] - a[0]*b[2];
   result[2] = a[0]*b[1] - a[1]*b[0];
}

struct Hybrid {

   // face data
   static pargrid::DataID dataFaceBID;
   static pargrid::DataID dataFaceJID;

   // cell data
   static pargrid::DataID dataCellRhoQiID;
   static pargrid::DataID dataCellBID;
   static pargrid::DataID dataCellJID;
   static pargrid::DataID dataCellUeID;
   static pargrid::DataID dataCellJiID;
   static pargrid::DataID dataCellMaxUeID;
   static pargrid::DataID dataCellMaxViID;
   static pargrid::DataID dataCellMinRhoQiID;
   static pargrid::DataID dataCellIonosphereID;
   static pargrid::DataID dataCellExosphereID;
   
   // node data
   static pargrid::DataID dataNodeRhoQiID;
   static pargrid::DataID dataNodeEID;
   static pargrid::DataID dataNodeBID;
   static pargrid::DataID dataNodeJID;
   static pargrid::DataID dataNodeUeID;
   static pargrid::DataID dataNodeJiID;

   // stencils
   static pargrid::StencilID accumulationStencilID; /**< ParGrid Stencil used to exchange accumulation array(s).*/
   
   // flags
   static pargrid::DataID dataInnerFlagFieldID;
   static pargrid::DataID dataInnerFlagNodeID;
   static pargrid::DataID dataInnerFlagParticleID;
#ifdef USE_XMIN_BOUNDARY
   static pargrid::DataID dataXminFlagID;
#endif
#ifdef ION_SPECTRA_ALONG_ORBIT
   static pargrid::DataID dataSpectraFlagID;
   static pargrid::DataID dataSpectraID;
#endif
   
   // bit masks
   static uint32_t X_POS_EXISTS;
   static uint32_t X_NEG_EXISTS;
   static uint32_t Y_POS_EXISTS;
   static uint32_t Y_NEG_EXISTS;
   static uint32_t Z_POS_EXISTS;
   static uint32_t Z_NEG_EXISTS;

   static int logInterval;
   static Real dx;
   static Real dV;
   static Real R_object;
   static Real R2_fieldObstacle;
   static Real R2_particleObstacle;
#ifdef USE_XMIN_BOUNDARY
   static Real xMinBoundary;
#endif
   static Real M_object;
   static Real maxUe2;
   static Real maxVi2;
   static Real minRhoQi;
   static Real eta;
   static Real swMacroParticlesCellPerDt;
   static int Efilter;
   static Real IMFBx,IMFBy,IMFBz;
#if defined(USE_B_INITIAL) || defined(USE_B_CONSTANT)
   static Real laminarR2,laminarR3,coeffDip,coeffQuad,dipSurfB,dipSurfR,dipMinR2,dipMomCoeff,xDip,yDip,zDip,thetaDip,phiDip;
#endif
   static unsigned int N_populations;
   static unsigned int N_ionospherePopulations;
   static unsigned int N_exospherePopulations;
   static unsigned int N_outputFields;
   static std::vector<std::string> populationNames;
   static std::vector<std::string> outputFieldStr;
   static std::vector<unsigned int> outputFieldId;
   static std::vector< std::vector<unsigned int> > outputFieldIdVector;

   static std::vector<std::ofstream*> plog;
   static std::ofstream flog;
   static std::vector<Real> particleCounterEscape;
   static std::vector<Real> particleCounterImpact;
   static std::vector<Real> particleCounterInject;
   static std::vector<Real> particleCounterInjectMacroparticles;
   static Real particleCounterTimeStart;
   
#ifdef WRITE_POPULATION_AVERAGES
   static pargrid::DataID dataCellAverageBID;
   static std::vector<pargrid::DataID> dataCellAverageDensityID;
   static std::vector<pargrid::DataID> dataCellAverageVelocityID;
   static int averageCounter;
#endif
   
   static std::map<std::string,bool> outParams;
};

#endif
