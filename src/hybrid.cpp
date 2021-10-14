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

#include "hybrid.h"

using namespace std;

// Init static variables:

map< string , HybridVariable<Real> > Hybrid::varReal;
map< string , HybridVariable<bool> > Hybrid::varBool;

// face data
pargrid::DataID Hybrid::dataFaceBID;
pargrid::DataID Hybrid::dataFaceJID;

// cell data
pargrid::DataID Hybrid::dataCellRhoQiID;
pargrid::DataID Hybrid::dataCellRhoQiBgID;
pargrid::DataID Hybrid::dataCellBID;
pargrid::DataID Hybrid::dataCellJID;
pargrid::DataID Hybrid::dataCellUeID;
pargrid::DataID Hybrid::dataCellJiID;
pargrid::DataID Hybrid::dataCellEpID;
pargrid::DataID Hybrid::dataCellIonosphereID;
pargrid::DataID Hybrid::dataCellExosphereID;

// node data
pargrid::DataID Hybrid::dataNodeRhoQiID;
pargrid::DataID Hybrid::dataNodeEID;
pargrid::DataID Hybrid::dataNodeBID;
pargrid::DataID Hybrid::dataNodeJID;
pargrid::DataID Hybrid::dataNodeUeID;
pargrid::DataID Hybrid::dataNodeJiID;

#ifdef USE_RESISTIVITY
pargrid::DataID Hybrid::dataNodeEtaID;
#endif

// counters
pargrid::DataID Hybrid::dataCounterCellMaxUeID;
pargrid::DataID Hybrid::dataCounterCellMaxViID;
pargrid::DataID Hybrid::dataCounterCellMinRhoQiID;
#ifdef USE_ECUT
pargrid::DataID Hybrid::dataCounterNodeEcutID;
#endif
#ifdef USE_MAXVW
pargrid::DataID Hybrid::dataCounterNodeMaxVwID;
#endif

// stencils
pargrid::StencilID Hybrid::accumulationStencilID;

// flags
pargrid::DataID Hybrid::dataInnerFlagFieldID;
pargrid::DataID Hybrid::dataInnerFlagNodeID;
pargrid::DataID Hybrid::dataInnerFlagParticleID;
pargrid::DataID Hybrid::dataInnerFlagCellEpID;
pargrid::DataID Hybrid::dataOuterBoundaryFlagID;
#ifdef USE_XMIN_BOUNDARY
pargrid::DataID Hybrid::dataXminFlagID;
#endif
#ifdef USE_DETECTORS
// detector: particle
pargrid::DataID Hybrid::dataDetectorParticleFlagID;
Real Hybrid::detParticleStartTime;
Real Hybrid::detParticleEndTime;
Real Hybrid::N_detParticleMaxFileLines;
Real Hybrid::detParticleWriteInterval;
Real Hybrid::detParticleTimestepCnt;
Real Hybrid::detParticleFileLineCnt;
bool Hybrid::detParticleRecording = false;
vector<Real> Hybrid::detParticleOutput;
// detector: bulk parameters
pargrid::DataID Hybrid::dataDetectorBulkParamFlagID;
Real Hybrid::detBulkParamStartTime;
Real Hybrid::detBulkParamEndTime;
Real Hybrid::N_detBulkParamMaxFileLines;
Real Hybrid::detBulkParamWriteInterval;
Real Hybrid::detBulkParamTimestepCnt;
Real Hybrid::detBulkParamFileLineCnt;
bool Hybrid::detBulkParamRecording = false;
vector<Real> Hybrid::detBulkParamOutput;
#endif

// bit masks to check the existence of +x, -x, +y, -y, +z, -z neighbour cell
uint32_t Hybrid::X_POS_EXISTS;   
uint32_t Hybrid::X_NEG_EXISTS;   
uint32_t Hybrid::Y_POS_EXISTS;
uint32_t Hybrid::Y_NEG_EXISTS;
uint32_t Hybrid::Z_POS_EXISTS;
uint32_t Hybrid::Z_NEG_EXISTS;

int Hybrid::logInterval;
bool Hybrid::includeInnerCellsInFieldLog;
Real Hybrid::dx;
Box Hybrid::box;
Real Hybrid::dV;
Real Hybrid::R_object;
Real Hybrid::R2_fieldObstacle;
Real Hybrid::R2_particleObstacle;
Real Hybrid::R2_cellEpObstacle;
#ifdef USE_XMIN_BOUNDARY
Real Hybrid::xMinBoundary;
#endif
#ifdef USE_CONIC_INNER_BOUNDARY
Real Hybrid::l_conicInnerBoundary;
Real Hybrid::e_conicInnerBoundary;
Real Hybrid::eta_conicInnerBoundary;
#endif
Real Hybrid::M_object;
Real Hybrid::maxUe2;
Real Hybrid::maxVi2;
Real Hybrid::maxVi;
Real Hybrid::terminateLimitMaxB;
Real Hybrid::minRhoQi;
OuterBoundaryZone Hybrid::outerBoundaryZone;
#ifdef USE_ECUT
Real Hybrid::Ecut2;
#endif
#ifdef USE_MAXVW
Real Hybrid::maxVw;
#endif
bool Hybrid::useHallElectricField;
bool Hybrid::useElectronPressureElectricField;
bool Hybrid::useAdiabaticElectronPressure;
Real Hybrid::electronTemperature;
Real Hybrid::electronPressureCoeff;
Real Hybrid::swMacroParticlesCellPerDt;
int Hybrid::Efilter;
Real Hybrid::EfilterNodeGaussSigma;
Real Hybrid::EfilterNodeGaussCoeffs[4];
#ifdef USE_RESISTIVITY
Real Hybrid::resistivityEta;
Real Hybrid::resistivityEtaC;
Real Hybrid::resistivityR2;
Real Hybrid::resistivityGridUnit;
Real (*Hybrid::resistivityProfilePtr)(Simulation& sim,SimulationClasses&,const Real x,const Real y,const Real z);
#endif
Real Hybrid::IMFBx;
Real Hybrid::IMFBy;
Real Hybrid::IMFBz;
#if defined(USE_B_INITIAL) || defined(USE_B_CONSTANT)
Real Hybrid::laminarR2;
Real Hybrid::laminarR3;
Real Hybrid::coeffDip;
Real Hybrid::coeffQuad;
Real Hybrid::dipSurfB;
Real Hybrid::dipSurfR;
Real Hybrid::dipMinR2;
Real Hybrid::dipMomCoeff;
Real Hybrid::xDip;
Real Hybrid::yDip;
Real Hybrid::zDip;
Real Hybrid::thetaDip;
Real Hybrid::phiDip;
void (*Hybrid::magneticFieldProfilePtr)(const Real x,const Real y,const Real z,Real B[3]);
#endif
// total number of particle populations
unsigned int Hybrid::N_populations;
// number of ionospheric particle populations
unsigned int Hybrid::N_ionospherePopulations;
// number of exospheric particle populations
unsigned int Hybrid::N_exospherePopulations;
// properties of all particle populations
vector<particlePopulation> Hybrid::allPops;
// properties of all solar wind populations
vector<solarWindPopulation> Hybrid::swPops;
// names of particle populations
vector<string> Hybrid::populationNames;
// number of output particle variables
unsigned int Hybrid::N_outputPopVars;
// names of output particle variables
vector<string> Hybrid::outputPopVarStr;
// ids of output particle variables for each population (<0 means no output)
vector<int> Hybrid::outputPopVarId;
// ids (=popid-1) of particle populations included in each output particle variable
vector< vector<unsigned int> > Hybrid::outputPopVarIdVector;
// popids of particle populations included in the total plasma variables
vector<unsigned int> Hybrid::outputPlasmaPopId;
// output cell variables
map<string,bool> Hybrid::outputCellParams;

vector<ofstream*> Hybrid::plog;
ofstream Hybrid::flog;

vector<Real> Hybrid::particleCounterEscape;
vector<Real> Hybrid::particleCounterImpact;
vector<Real> Hybrid::particleCounterInject;
vector<Real> Hybrid::particleCounterInjectMacroparticles;
Real Hybrid::particleCounterTimeStart;

bool Hybrid::filterParticlesAfterRestartDone = true;

#ifdef WRITE_POPULATION_AVERAGES
pargrid::DataID Hybrid::dataCellAverageBID;
vector<pargrid::DataID> Hybrid::dataCellAverageDensityID;
vector<pargrid::DataID> Hybrid::dataCellAverageVelocityID;
int Hybrid::averageCounter;
#endif

int Hybrid::repartitionCheckIntervalTmp = -1;
