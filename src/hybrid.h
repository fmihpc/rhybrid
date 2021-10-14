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

#ifndef HYBRID_H
#define HYBRID_H

#include <cstdlib>
#include <vector>
#include <map>
#include <sstream>
#include <simulationclasses.h>
#include <detectors.h>

#define sqr(x) ((x)*(x))
#define cube(x) (x)*(x)*(x)
#define vecsqr(a) (sqr(a[0])+sqr(a[1])+sqr(a[2]))
#define normvec(a) (sqrt(vecsqr(a)))

inline std::string real2str(Real x,unsigned int prec) {
    std::stringstream ss;
    ss.precision(prec);
    ss << x;
    return ss.str();
}

inline std::string int2str(int x,unsigned int N) {
    std::stringstream ss;
    ss.width(N);
    ss.fill('0');
    ss << x;
    return ss.str();
}

template<typename T>
struct HybridVariable {
   std::string name = "";
   std::string type = "";
   std::size_t vectorDim = 0;
   pargrid::DataID id = pargrid::INVALID_DATAID;
   T* ptr = NULL;
};

inline void cross(const Real a[3], const Real b[3], Real result[3]) {
   result[0] = a[1]*b[2] - a[2]*b[1];
   result[1] = a[2]*b[0] - a[0]*b[2];
   result[2] = a[0]*b[1] - a[1]*b[0];
}

struct particlePopulation {
    Real w;
    std::string name;
};

struct solarWindPopulation {
   Real m,q,U,n,vth,T;
   std::string name;
};

struct Box {
   Real xmin=0.0,xmax=0.0,ymin=0.0,ymax=0.0,zmin=0.0,zmax=0.0;
};

struct OuterBoundaryZone {
   int typeEta=0,typeMinRhoQi=0;
   Real sizeEta=0.0,sizeMinRhoQi=0.0,minRhoQi=0.0,eta=0.0;
};

struct Hybrid {
   static std::map< std::string, HybridVariable<Real> > varReal;
   static std::map< std::string, HybridVariable<bool> > varBool;

   // face data
   static pargrid::DataID dataFaceBID;
   static pargrid::DataID dataFaceJID;

   // cell data
   static pargrid::DataID dataCellRhoQiID;
   static pargrid::DataID dataCellRhoQiBgID;
   static pargrid::DataID dataCellBID;
   static pargrid::DataID dataCellJID;
   static pargrid::DataID dataCellUeID;
   static pargrid::DataID dataCellJiID;
   static pargrid::DataID dataCellEpID;
   static pargrid::DataID dataCellIonosphereID;
   static pargrid::DataID dataCellExosphereID;
   
   // node data
   static pargrid::DataID dataNodeRhoQiID;
   static pargrid::DataID dataNodeEID;
   static pargrid::DataID dataNodeBID;
   static pargrid::DataID dataNodeJID;
   static pargrid::DataID dataNodeUeID;
   static pargrid::DataID dataNodeJiID;

#ifdef USE_RESISTIVITY
   static pargrid::DataID dataNodeEtaID;
#endif

   // counters
   static pargrid::DataID dataCounterCellMaxUeID;
   static pargrid::DataID dataCounterCellMaxViID;
   static pargrid::DataID dataCounterCellMinRhoQiID;
#ifdef USE_ECUT
   static pargrid::DataID dataCounterNodeEcutID;
#endif
#ifdef USE_MAXVW
   static pargrid::DataID dataCounterNodeMaxVwID;
#endif

   // stencils
   static pargrid::StencilID accumulationStencilID; /**< ParGrid Stencil used to exchange accumulation array(s).*/
   
   // flags
   static pargrid::DataID dataInnerFlagFieldID;
   static pargrid::DataID dataInnerFlagNodeID;
   static pargrid::DataID dataInnerFlagParticleID;
   static pargrid::DataID dataInnerFlagCellEpID;
   static pargrid::DataID dataOuterBoundaryFlagID;
#ifdef USE_XMIN_BOUNDARY
   static pargrid::DataID dataXminFlagID;
#endif
#ifdef USE_DETECTORS
   // detector: particles
   static pargrid::DataID dataDetectorParticleFlagID;
   static Real detParticleStartTime;
   static Real detParticleEndTime;
   static Real N_detParticleMaxFileLines;
   static Real detParticleWriteInterval;
   static Real detParticleTimestepCnt;
   static Real detParticleFileLineCnt;
   static bool detParticleRecording;
   static std::vector<Real> detParticleOutput;
   // detector: bulk parameters
   static pargrid::DataID dataDetectorBulkParamFlagID;
   static Real detBulkParamStartTime;
   static Real detBulkParamEndTime;
   static Real N_detBulkParamMaxFileLines;
   static Real detBulkParamWriteInterval;
   static Real detBulkParamTimestepCnt;
   static Real detBulkParamFileLineCnt;
   static bool detBulkParamRecording;
   static std::vector<Real> detBulkParamOutput;
#endif
   
   // bit masks
   static uint32_t X_POS_EXISTS;
   static uint32_t X_NEG_EXISTS;
   static uint32_t Y_POS_EXISTS;
   static uint32_t Y_NEG_EXISTS;
   static uint32_t Z_POS_EXISTS;
   static uint32_t Z_NEG_EXISTS;

   static int logInterval;
   static bool includeInnerCellsInFieldLog;
   static Real dx;
   static Box box;
   static Real dV;
   static Real R_object;
   static Real R2_fieldObstacle;
   static Real R2_particleObstacle;
   static Real R2_cellEpObstacle;
#ifdef USE_XMIN_BOUNDARY
   static Real xMinBoundary;
#endif
#ifdef USE_CONIC_INNER_BOUNDARY
   static Real l_conicInnerBoundary;
   static Real e_conicInnerBoundary;
   static Real eta_conicInnerBoundary;
#endif
   static Real M_object;
   static Real maxUe2;
   static Real maxVi2;
   static Real maxVi;
   static Real terminateLimitMaxB;
   static Real minRhoQi;
   static OuterBoundaryZone outerBoundaryZone;
#ifdef USE_ECUT
   static Real Ecut2;
#endif
#ifdef USE_MAXVW
   static Real maxVw;
#endif
#ifdef USE_RESISTIVITY
   static Real resistivityEta;
   static Real resistivityEtaC;
   static Real resistivityR2;
   static Real resistivityGridUnit;
   static Real (*resistivityProfilePtr)(Simulation& sim,SimulationClasses&,const Real x,const Real y,const Real z);
#endif
   static bool useHallElectricField;
   static bool useElectronPressureElectricField;
   static bool useAdiabaticElectronPressure;
   static Real electronTemperature;
   static Real electronPressureCoeff;
   static Real swMacroParticlesCellPerDt;
   static int Efilter;
   static Real EfilterNodeGaussSigma;
   static Real EfilterNodeGaussCoeffs[4];
   static Real IMFBx,IMFBy,IMFBz;
#if defined(USE_B_INITIAL) || defined(USE_B_CONSTANT)
   static Real laminarR2,laminarR3,coeffDip,coeffQuad,dipSurfB,dipSurfR,dipMinR2,dipMomCoeff,xDip,yDip,zDip,thetaDip,phiDip;
   static void (*magneticFieldProfilePtr)(const Real x,const Real y,const Real z,Real B[3]);
#endif
   static unsigned int N_populations;
   static unsigned int N_ionospherePopulations;
   static unsigned int N_exospherePopulations;
   static unsigned int N_outputPopVars;
   static std::vector<particlePopulation> allPops;
   static std::vector<solarWindPopulation> swPops;
   static std::vector<std::string> populationNames;
   static std::vector<std::string> outputPopVarStr;
   static std::vector<int> outputPopVarId;
   static std::vector< std::vector<unsigned int> > outputPopVarIdVector;
   static std::vector<unsigned int> outputPlasmaPopId;
   static std::map<std::string,bool> outputCellParams;

   static std::vector<std::ofstream*> plog;
   static std::ofstream flog;
   static std::vector<Real> particleCounterEscape;
   static std::vector<Real> particleCounterImpact;
   static std::vector<Real> particleCounterInject;
   static std::vector<Real> particleCounterInjectMacroparticles;
   static Real particleCounterTimeStart;

   static bool filterParticlesAfterRestartDone;
   
#ifdef WRITE_POPULATION_AVERAGES
   static pargrid::DataID dataCellAverageBID;
   static std::vector<pargrid::DataID> dataCellAverageDensityID;
   static std::vector<pargrid::DataID> dataCellAverageVelocityID;
   static int averageCounter;
#endif
   static int repartitionCheckIntervalTmp;
};

#endif
