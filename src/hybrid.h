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

#include <cmath>
#include <cstdlib>
#include <utility>
#include <vector>
#include <map>
#include <sstream>
#include <type_traits>
#include <simulationclasses.h>
#include <detectors.h>

template<class T> T sqr(const T& x) { return x*x; }
template<class T> T cube(const T& x) { return x*x*x; }
template<class T> typename std::remove_reference<decltype(std::declval<T>()[0])>::type vecsqr(const T& x) { return sqr(x[0])+sqr(x[1])+sqr(x[2]); }
template<class T> typename std::remove_reference<decltype(std::declval<T>()[0])>::type normvec(const T& x) { using std::sqrt; return sqrt(vecsqr(x)); }

// convert Real to string with given precicion
inline std::string real2str(Real x,unsigned int prec) {
    std::stringstream ss;
    ss.precision(prec);
    ss << x;
    return ss.str();
}

// convert int to string with at least N digits and filled with leading zeros
inline std::string int2str(int x,unsigned int N) {
    std::stringstream ss;
    ss.width(N);
    ss.fill('0');
    ss << x;
    return ss.str();
}

// convert string velocity vector cfg variable to Real vector
inline bool convertConfigFileVariableVelocity(std::string velStr,std::vector<Real>& vel) {
   // check velocity format: (Ux,Uy,Uz)
   bool velStrOk = true;
   if (count(velStr.begin(),velStr.end(),'(') != 1 ||
      count(velStr.begin(),velStr.end(),')') != 1 ||
      count(velStr.begin(),velStr.end(),',') != 2 ||
      velStr.find_first_not_of("(),+-0123456789.e ") != std::string::npos) {
      velStrOk = false;
   }
   // remove non-numeral characters from string
   std::string velStrEdit(velStr);
   replace(velStrEdit.begin(),velStrEdit.end(),'(',' ');
   replace(velStrEdit.begin(),velStrEdit.end(),')',' ');
   replace(velStrEdit.begin(),velStrEdit.end(),',',' ');
   // remove multiple whitespace: not needed
   //velStrEdit.erase(unique(velStrEdit.begin(),velStrEdit.end(),[](char a,char b) { return isspace(a) && isspace(b); } ),velStrEdit.end() );
   // convert string to Reals
   vel.clear();
   Real vtmp;
   std::stringstream ss(velStrEdit);
   while (ss >> vtmp) { vel.push_back(vtmp); }
   // check velocity format: (Ux,Uy,Uz)
   if (vel.size() != 3) { velStrOk = false; }
   return velStrOk;
}

// cross product of 3D vectors
inline void cross(const Real a[3], const Real b[3], Real result[3]) {
   result[0] = a[1]*b[2] - a[2]*b[1];
   result[1] = a[2]*b[0] - a[0]*b[2];
   result[2] = a[0]*b[1] - a[1]*b[0];
}

#ifdef USE_OUTER_BOUNDARY_ZONE
struct OuterBoundaryZone {
   int typeEta=0,typeMinRhoQi=0;
   Real sizeEta=0.0,sizeMinRhoQi=0.0,minRhoQi=0.0,eta=0.0;
   bool constUe = false;
};
#endif

#ifdef USE_NEW_VARIBLE_HANDLING
template<typename T>
struct HybridVariable {
   std::string name = "";
   std::string type = "";
   std::size_t vectorDim = 0;
   pargrid::DataID dataID = pargrid::INVALID_DATAID;
   pargrid::StencilID stencilID = pargrid::INVALID_STENCILID;
   T* ptr = NULL;
   //T* const constPtr = NULL;
};
#endif

struct Hybrid {
#ifdef USE_NEW_VARIBLE_HANDLING
   static std::map< std::string, HybridVariable<Real> > varReal;
   static std::map< std::string, HybridVariable<bool> > varBool;
#endif

   // face data
#ifndef USE_NEW_VARIBLE_HANDLING
   static pargrid::DataID dataFaceBID;
#endif
   static pargrid::DataID dataFaceJID;

   // cell data
   static pargrid::DataID dataCellRhoQiID;
#ifdef USE_BACKGROUND_CHARGE_DENSITY
   static pargrid::DataID dataCellRhoQiBgID;
#endif
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

   // grid constraint counters
#ifdef USE_GRID_CONSTRAINT_COUNTERS
   static pargrid::DataID dataGridCounterCellMaxUeID;
   static pargrid::DataID dataGridCounterCellMaxViID;
   static pargrid::DataID dataGridCounterCellMinRhoQiID;
   static pargrid::DataID dataGridCounterNodeMaxEID;
   static pargrid::DataID dataGridCounterNodeMaxVwID;
#endif

   // stencils
   static pargrid::StencilID accumulationStencilID; // ParGrid Stencil used to exchange accumulation array(s)

   // flags
   static pargrid::DataID dataInnerFlagFieldID;
   static pargrid::DataID dataInnerFlagNodeID;
   static pargrid::DataID dataInnerFlagParticleID;
   static pargrid::DataID dataInnerFlagCellEpID;
#ifdef USE_OUTER_BOUNDARY_ZONE
   static pargrid::DataID dataOuterBoundaryFlagID;
   static pargrid::DataID dataOuterBoundaryFlagNodeID;
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

   static unsigned int logInterval;
   static unsigned int mainLogDiagnosticsInterval;
   static bool writeMainLogDiagnosticsAfterLogStep;
   static bool includeInnerCellsInFieldLog;
   static unsigned int saveReducedStateInterval;
   static unsigned int saveReducedStateNstride;
   static bool saveReducedStateParticles;
   static bool saveParticles;
   static unsigned int saveParticlesNstride;
   static Real dx;
   static Real dV;
   static Real R_object;
   static Real R2_fieldObstacle;
   static Real fieldObstacleUe[3];
   static Real R2_particleObstacle;
   static Real R2_cellEpObstacle;
   static Real upstreamBulkU;
   static Real M_object;
   static Real GMdt;
   static bool initialFlowThrough;
   static Real initialFlowThroughPeriod;
   static Real maxUe2;
   static Real maxVi2;
   static Real maxVi;
   static Real terminateLimitMaxB;
   static Real minRhoQi;
   static Real maxE2;
   static Real maxVw;
#ifdef USE_RESISTIVITY
   static Real resistivityEta;
   static Real resistivityR2;
   static std::vector<Real> resistivitySphericalEta;
   static std::vector<Real> resistivitySphericalR2;
   static Real (*resistivityProfilePtr)(Simulation& sim,SimulationClasses&,const Real x,const Real y,const Real z);
#endif
#ifdef USE_OUTER_BOUNDARY_ZONE
   static OuterBoundaryZone outerBoundaryZone;
#endif
   static bool useHallElectricField;
#ifdef USE_B_CONSTANT
   static bool includeConstantB0InFaradaysLaw;
#endif
   static bool useElectronPressureElectricField;
   static bool useAdiabaticElectronPressure;
   static Real electronTemperature;
   static Real electronPressureCoeff;
   static Real swMacroParticlesCellPerDt;
   static bool useGravity;
   static int Efilter;
   static Real EfilterNodeGaussSigma;
   static Real EfilterNodeGaussCoeffs[4];
   static Real IMFBx,IMFBy,IMFBz;
   static bool IMFBoundaryCellB[6];
   static bool IMFBoundaryFaceB[6];
#if defined(USE_B_INITIAL) || defined(USE_B_CONSTANT)
   static Real dBx,dBy,dBz;
   static Real laminarR2,laminarR3,coeffDip,coeffQuad,dipSurfB,dipSurfR,dipMinR2,dipMomCoeff,xDip,yDip,zDip,thetaDip,phiDip;
   static std::vector<Real> xDipMirror,yDipMirror,zDipMirror;
   static void (*magneticFieldProfilePtr)(Simulation& sim,SimulationClasses& simClasses,const Real x,const Real y,const Real z,Real B[3]);
#endif
   static unsigned int N_populations;
   static unsigned int N_ionospherePopulations;
   static unsigned int N_exospherePopulations;
   static unsigned int N_outputPopVars;
   static std::vector<std::string> populationNames;
   static std::vector<std::string> outputPopVarStr;
   static std::vector<int> outputPopVarId;
   static std::vector< std::vector<unsigned int> > outputPopVarIdVector;
   static std::vector<unsigned int> outputPlasmaPopId;
   static std::map<std::string,bool> outputCellParams;

   // particle population and field logs and their counters
   static std::vector<std::ofstream*> logParticle;
   static std::ofstream logField;
   static std::vector<Real> logCounterParticleEscape;
   static std::vector<Real> logCounterParticleImpact;
   static std::vector<Real> logCounterParticleInject;
   static std::vector<Real> logCounterParticleInjectMacroparticles;
   static std::vector<Real> logCounterParticleEscapeKineticEnergy;
   static std::vector<Real> logCounterParticleImpactKineticEnergy;
   static std::vector<Real> logCounterParticleInjectKineticEnergy;
   static std::vector<Real> logCounterParticleMaxVi;
   static Real logCounterFieldMaxCellUe;
   static Real logCounterFieldMaxNodeUe;
   static Real logCounterFieldMaxVw;
   static Real logCounterFieldMaxE;
   static Real logCounterFieldMinCellRhoQi;
   static Real logCounterFieldMinNodeRhoQi;
   static Real logCounterTimeStart;

   static bool filterParticlesAfterRestartDone;

#ifdef WRITE_GRID_TEMPORAL_AVERAGES
   static pargrid::DataID dataCellAverageBID;
   static std::vector<pargrid::DataID> dataCellAverageDensityID;
   static std::vector<pargrid::DataID> dataCellAverageVelocityID;
   static int gridTemporalAverageCounter;
#endif
   static int repartitionCheckIntervalTmp;
};

#endif
