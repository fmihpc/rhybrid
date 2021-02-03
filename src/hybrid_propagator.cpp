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

#include <cstdlib>
#include <iostream>
#include <cmath>

#include "hybrid.h"
#include "hybrid_propagator.h"
#include "particle_definition.h"
#ifdef USE_B_CONSTANT
#include "magnetic_field.h"
#endif

using namespace std;

static int totalID = -1;
static int profIntpolID = -1;
static int profPropagFieldID = -1;
static int profBoundCondsID = -1;
static int mpiWaitID = -1;
static int setupGetFieldsID = -1;

static bool saveStepHappened=false;

bool propagateB(Simulation& sim,SimulationClasses& simClasses,vector<ParticleListBase*>& particleLists) {
   bool success = true;
   profile::start("propagateB",totalID);   
   
   // get data array pointers
   //Real* faceB_ = Hybrid::varReal["faceB_"].ptr;
   Real* faceB               = simClasses.pargrid.getUserDataStatic<Real>(Hybrid::dataFaceBID);
   Real* faceJ               = simClasses.pargrid.getUserDataStatic<Real>(Hybrid::dataFaceJID);
   Real* cellRhoQi           = simClasses.pargrid.getUserDataStatic<Real>(Hybrid::dataCellRhoQiID);
   Real* cellB               = simClasses.pargrid.getUserDataStatic<Real>(Hybrid::dataCellBID);
   Real* cellJ               = simClasses.pargrid.getUserDataStatic<Real>(Hybrid::dataCellJID);
   Real* cellUe              = simClasses.pargrid.getUserDataStatic<Real>(Hybrid::dataCellUeID);
   Real* cellJi              = simClasses.pargrid.getUserDataStatic<Real>(Hybrid::dataCellJiID);
   Real* cellEp              = simClasses.pargrid.getUserDataStatic<Real>(Hybrid::dataCellEpID);
   Real* nodeRhoQi           = simClasses.pargrid.getUserDataStatic<Real>(Hybrid::dataNodeRhoQiID);
   Real* nodeE               = simClasses.pargrid.getUserDataStatic<Real>(Hybrid::dataNodeEID);
   Real* nodeB               = simClasses.pargrid.getUserDataStatic<Real>(Hybrid::dataNodeBID);
   Real* nodeJ               = simClasses.pargrid.getUserDataStatic<Real>(Hybrid::dataNodeJID);
   Real* nodeUe              = simClasses.pargrid.getUserDataStatic<Real>(Hybrid::dataNodeUeID);
   Real* nodeJi              = simClasses.pargrid.getUserDataStatic<Real>(Hybrid::dataNodeJiID);
#ifdef USE_RESISTIVITY
   Real* nodeEta             = simClasses.pargrid.getUserDataStatic<Real>(Hybrid::dataNodeEtaID);
#endif
   Real* counterCellMaxUe    = simClasses.pargrid.getUserDataStatic<Real>(Hybrid::dataCounterCellMaxUeID);
   Real* counterCellMaxVi    = simClasses.pargrid.getUserDataStatic<Real>(Hybrid::dataCounterCellMaxViID);
   Real* counterCellMinRhoQi = simClasses.pargrid.getUserDataStatic<Real>(Hybrid::dataCounterCellMinRhoQiID);
#ifdef USE_ECUT
   Real* counterNodeEcut     = simClasses.pargrid.getUserDataStatic<Real>(Hybrid::dataCounterNodeEcutID);
#endif
#ifdef USE_MAXVW
   Real* counterNodeMaxVw    = simClasses.pargrid.getUserDataStatic<Real>(Hybrid::dataCounterNodeMaxVwID);
#endif
   bool* innerFlag           = simClasses.pargrid.getUserDataStatic<bool>(Hybrid::dataInnerFlagFieldID);
   bool* innerFlagNode       = simClasses.pargrid.getUserDataStatic<bool>(Hybrid::dataInnerFlagNodeID);
   bool* innerFlagParticle   = simClasses.pargrid.getUserDataStatic<bool>(Hybrid::dataInnerFlagParticleID);
   bool* innerFlagCellEp     = simClasses.pargrid.getUserDataStatic<bool>(Hybrid::dataInnerFlagCellEpID);
   bool* outerBoundaryFlag   = simClasses.pargrid.getUserDataStatic<bool>(Hybrid::dataOuterBoundaryFlagID);
   
   if(faceB               == NULL) {cerr << "ERROR: obtained NULL faceB array!"        << endl; exit(1);}
   if(faceJ               == NULL) {cerr << "ERROR: obtained NULL faceJ array!"        << endl; exit(1);}
   if(cellRhoQi           == NULL) {cerr << "ERROR: obtained NULL cellRhoQi array!"    << endl; exit(1);}
   if(cellB               == NULL) {cerr << "ERROR: obtained NULL cellB array!"        << endl; exit(1);}
   if(cellJ               == NULL) {cerr << "ERROR: obtained NULL cellJ array!"        << endl; exit(1);}
   if(cellUe              == NULL) {cerr << "ERROR: obtained NULL cellUe array!"       << endl; exit(1);}
   if(cellJi              == NULL) {cerr << "ERROR: obtained NULL cellJi array!"       << endl; exit(1);}
   if(cellEp              == NULL) {cerr << "ERROR: obtained NULL cellEp array!"       << endl; exit(1);}
   if(nodeRhoQi           == NULL) {cerr << "ERROR: obtained NULL nodeRhoQi array!"    << endl; exit(1);}
   if(nodeE               == NULL) {cerr << "ERROR: obtained NULL nodeE array!"        << endl; exit(1);}
   if(nodeB               == NULL) {cerr << "ERROR: obtained NULL nodeB array!"        << endl; exit(1);}
   if(nodeJ               == NULL) {cerr << "ERROR: obtained NULL nodeJ array!"        << endl; exit(1);}
   if(nodeUe              == NULL) {cerr << "ERROR: obtained NULL nodeUe array!"       << endl; exit(1);}
   if(nodeJi              == NULL) {cerr << "ERROR: obtained NULL nodeJi array!"       << endl; exit(1);}
#ifdef USE_RESISTIVITY
   if(nodeEta             == NULL) {cerr << "ERROR: obtained NULL nodeEta array!"      << endl; exit(1);}
#endif
   if(counterCellMaxUe    == NULL) {cerr << "ERROR: obtained NULL counterCellMaxUe array!"    << endl; exit(1);}
   if(counterCellMaxVi    == NULL) {cerr << "ERROR: obtained NULL counterCellMaxVi array!"    << endl; exit(1);}
   if(counterCellMinRhoQi == NULL) {cerr << "ERROR: obtained NULL counterCellMinRhoQi array!" << endl; exit(1);}
#ifdef USE_ECUT
   if(counterNodeEcut     == NULL) {cerr << "ERROR: obtained NULL counterNodeEcut array!"     << endl; exit(1);}
#endif
#ifdef USE_MAXVW
   if(counterNodeMaxVw    == NULL) {cerr << "ERROR: obtained NULL counterNodeMaxVw array!"    << endl; exit(1);}
#endif
   if(innerFlag           == NULL) {cerr << "ERROR: obtained NULL innerFlag array!"    << endl; exit(1);}
   if(innerFlagParticle   == NULL) {cerr << "ERROR: obtained NULL innerFlagParticle array!" << endl; exit(1);}
   if(innerFlagCellEp     == NULL) {cerr << "ERROR: obtained NULL innerFlagCellEp array!" << endl; exit(1);}
   if(innerFlagNode       == NULL) {cerr << "ERROR: obtained NULL innerFlagNode array!"<< endl; exit(1);}
   if(outerBoundaryFlag   == NULL) {cerr << "ERROR: obtained NULL outerBoundaryFlag array!"<< endl; exit(1);}
   
   // get block vectors
   const vector<pargrid::CellID>& innerBlocks = simClasses.pargrid.getInnerCells(pargrid::DEFAULT_STENCIL);
   const vector<pargrid::CellID>& boundaryBlocks = simClasses.pargrid.getBoundaryCells(pargrid::DEFAULT_STENCIL);
   const vector<pargrid::CellID>& exteriorBlocks = simClasses.pargrid.getExteriorCells();

   // zero diagnostic variables
   if(saveStepHappened == true) {
      for(pargrid::CellID b=0;b<simClasses.pargrid.getNumberOfLocalCells();++b) for(int k=0;k<block::WIDTH_Z;++k) for(int j=0;j<block::WIDTH_Y;++j) for(int i=0;i<block::WIDTH_X;++i) {
	 const int n = (b*block::SIZE+block::index(i,j,k));
	 counterCellMaxUe[n]=0.0;
	 counterCellMaxVi[n]=0.0;
	 counterCellMinRhoQi[n]=0.0;
#ifdef USE_ECUT
         counterNodeEcut[n]=0.0;
#endif
#ifdef USE_MAXVW
         counterNodeMaxVw[n]=0.0;
#endif
      }
      saveStepHappened = false;
   }
   
   // face->cell B
   simClasses.pargrid.startNeighbourExchange(pargrid::DEFAULT_STENCIL,Hybrid::dataFaceBID);
   profile::start("intpol",profIntpolID);
   for(pargrid::CellID b=0; b<innerBlocks.size(); ++b) { face2Cell(faceB,cellB,sim,simClasses,innerBlocks[b]); }
   profile::stop();
   profile::start("MPI waits",mpiWaitID);
   simClasses.pargrid.wait(pargrid::DEFAULT_STENCIL,Hybrid::dataFaceBID);
   profile::stop();
   profile::start("intpol",profIntpolID);
   for(pargrid::CellID b=0; b<boundaryBlocks.size(); ++b) { face2Cell(faceB,cellB,sim,simClasses,boundaryBlocks[b]); }
   profile::stop();
   
   // SET FIELD BOUNDARY CONDITIONS
   
   // loop thru ghost cells
   profile::start("BoundaryConds",profBoundCondsID);
   simClasses.pargrid.startNeighbourExchange(pargrid::DEFAULT_STENCIL,Hybrid::dataCellBID);
   simClasses.pargrid.startNeighbourExchange(pargrid::DEFAULT_STENCIL,Hybrid::dataCellJiID);
   simClasses.pargrid.startNeighbourExchange(pargrid::DEFAULT_STENCIL,Hybrid::dataCellRhoQiID);
   profile::start("MPI waits",mpiWaitID);
   simClasses.pargrid.wait(pargrid::DEFAULT_STENCIL,Hybrid::dataCellBID);
   simClasses.pargrid.wait(pargrid::DEFAULT_STENCIL,Hybrid::dataCellJiID);
   simClasses.pargrid.wait(pargrid::DEFAULT_STENCIL,Hybrid::dataCellRhoQiID);
   profile::stop();
   neumannCell(cellB,    sim,simClasses,exteriorBlocks,3);
   neumannCell(cellJi,   sim,simClasses,exteriorBlocks,3);
   neumannCell(cellRhoQi,sim,simClasses,exteriorBlocks,1);
   setIMF(cellB,sim,simClasses,exteriorBlocks);
   profile::stop();
   
#ifdef WRITE_POPULATION_AVERAGES
   // add cellB to cellAverageB and increase average counter
   Real* cellAverageB = simClasses.pargrid.getUserDataStatic<Real>(Hybrid::dataCellAverageBID);
   if(cellAverageB == NULL) { cerr << "ERROR: obtained NULL cellAverageB array!" << endl; exit(1); }
   for(pargrid::CellID b=0;b<simClasses.pargrid.getNumberOfLocalCells();++b) for(int k=0;k<block::WIDTH_Z;++k) for(int j=0;j<block::WIDTH_Y;++j) for(int i=0;i<block::WIDTH_X;++i) {
      const int n = (b*block::SIZE+block::index(i,j,k));
      const int n3 = n*3;
      for(int l=0;l<3;++l) { cellAverageB[n3+l] += cellB[n3+l]; }
   }
   Hybrid::averageCounter++;
#endif
   
   // cell->node B
   profile::start("intpol",profIntpolID);   
   for(pargrid::CellID b=0; b<simClasses.pargrid.getNumberOfLocalCells(); ++b) { cell2Node(cellB,nodeB,sim,simClasses,b); }
   profile::stop();

#ifdef USE_NODE_UE
   // cell -> node RhoQi and Ji
   simClasses.pargrid.startNeighbourExchange(pargrid::DEFAULT_STENCIL,Hybrid::dataCellRhoQiID);
   simClasses.pargrid.startNeighbourExchange(pargrid::DEFAULT_STENCIL,Hybrid::dataCellJiID);
   profile::start("intpol",profIntpolID);
   for(pargrid::CellID b=0; b<innerBlocks.size(); ++b) { cell2Node(cellRhoQi,nodeRhoQi,sim,simClasses,innerBlocks[b],1); }
   for(pargrid::CellID b=0; b<innerBlocks.size(); ++b) { cell2Node(cellJi,nodeJi,sim,simClasses,innerBlocks[b]); }
   profile::stop();
   profile::start("MPI waits",mpiWaitID);   
   simClasses.pargrid.wait(pargrid::DEFAULT_STENCIL,Hybrid::dataCellRhoQiID);
   simClasses.pargrid.wait(pargrid::DEFAULT_STENCIL,Hybrid::dataCellJiID);
   profile::stop();
   profile::start("intpol",profIntpolID);
   for(pargrid::CellID b=0; b<boundaryBlocks.size(); ++b) { cell2Node(cellRhoQi,nodeRhoQi,sim,simClasses,boundaryBlocks[b],1); }
   for(pargrid::CellID b=0; b<boundaryBlocks.size(); ++b) { cell2Node(cellJi,nodeJi,sim,simClasses,boundaryBlocks[b]); }
   profile::stop();
#endif
   
   // calculate J
#ifdef USE_EDGE_J
   neumannFace(faceB,sim,simClasses,exteriorBlocks);
   setIMFFace(faceB,sim,simClasses,exteriorBlocks);
   // nodeJ = avg(edgeJ) = avg(curl(faceB)/mu0)
   simClasses.pargrid.startNeighbourExchange(pargrid::DEFAULT_STENCIL,Hybrid::dataFaceBID);
   profile::start("field propag",profPropagFieldID);
   for(pargrid::CellID b=0; b<innerBlocks.size(); ++b) {
      calcNodeJ(faceB,nodeB,nodeRhoQi,nodeJ,
#ifdef USE_MAXVW
                counterNodeMaxVw,
#endif
                sim,simClasses,innerBlocks[b]);
   }
   profile::stop();
   profile::start("MPI waits",mpiWaitID);
   simClasses.pargrid.wait(pargrid::DEFAULT_STENCIL,Hybrid::dataFaceBID);
   profile::stop();
   profile::start("field propag",profPropagFieldID);
   for(pargrid::CellID b=0; b<boundaryBlocks.size(); ++b) {
      calcNodeJ(faceB,nodeB,nodeRhoQi,nodeJ,
#ifdef USE_MAXVW
                counterNodeMaxVw,
#endif
                sim,simClasses,boundaryBlocks[b]);
   }
   profile::stop(); 
   // node->cell J
   simClasses.pargrid.startNeighbourExchange(pargrid::DEFAULT_STENCIL,Hybrid::dataNodeJID);
   profile::start("intpol",profIntpolID);
   for(pargrid::CellID b=0; b<innerBlocks.size(); ++b) { node2Cell(nodeJ,cellJ,sim,simClasses,innerBlocks[b]); }
   profile::stop();
   profile::start("MPI waits",mpiWaitID);
   simClasses.pargrid.wait(pargrid::DEFAULT_STENCIL,Hybrid::dataNodeJID);
   profile::stop();
   profile::start("intpol",profIntpolID);
   for(pargrid::CellID b=0; b<boundaryBlocks.size(); ++b) { node2Cell(nodeJ,cellJ,sim,simClasses,boundaryBlocks[b]); }
   profile::stop();
#else
   // Ampere: faceJ = curl(nodeB)/mu0
   simClasses.pargrid.startNeighbourExchange(pargrid::DEFAULT_STENCIL,Hybrid::dataNodeBID);
   profile::start("field propag",profPropagFieldID);
   for(pargrid::CellID b=0; b<innerBlocks.size(); ++b) { faceCurl(nodeB,faceJ,false,sim,simClasses,innerBlocks[b]); }
   profile::stop();
   profile::start("MPI waits",mpiWaitID);
   simClasses.pargrid.wait(pargrid::DEFAULT_STENCIL,Hybrid::dataNodeBID);
   profile::stop();
   profile::start("field propag",profPropagFieldID);
   for(pargrid::CellID b=0; b<boundaryBlocks.size(); ++b) { faceCurl(nodeB,faceJ,false,sim,simClasses,boundaryBlocks[b]); }
   profile::stop();
   
   // face->cell J
   simClasses.pargrid.startNeighbourExchange(pargrid::DEFAULT_STENCIL,Hybrid::dataFaceJID);
   profile::start("intpol",profIntpolID);
   for(pargrid::CellID b=0; b<innerBlocks.size(); ++b) { face2Cell(faceJ,cellJ,sim,simClasses,innerBlocks[b]); }
   profile::stop();
   profile::start("MPI waits",mpiWaitID);
   simClasses.pargrid.wait(pargrid::DEFAULT_STENCIL,Hybrid::dataFaceJID);
   profile::stop();
   profile::start("intpol",profIntpolID);
   for(pargrid::CellID b=0; b<boundaryBlocks.size(); ++b) { face2Cell(faceJ,cellJ,sim,simClasses,boundaryBlocks[b]); }
   profile::stop();

   simClasses.pargrid.startNeighbourExchange(pargrid::DEFAULT_STENCIL,Hybrid::dataCellJID);
   simClasses.pargrid.wait(pargrid::DEFAULT_STENCIL,Hybrid::dataCellJID);
   neumannCell(cellJ,sim,simClasses,exteriorBlocks,3);
   
   // cell->node J
   simClasses.pargrid.startNeighbourExchange(pargrid::DEFAULT_STENCIL,Hybrid::dataCellJID);
   profile::start("intpol",profIntpolID);
   for(pargrid::CellID b=0; b<innerBlocks.size(); ++b) { cell2Node(cellJ,nodeJ,sim,simClasses,innerBlocks[b]); }
   profile::stop();
   profile::start("MPI waits",mpiWaitID);
   simClasses.pargrid.wait(pargrid::DEFAULT_STENCIL,Hybrid::dataCellJID);
   profile::stop();
   profile::start("intpol",profIntpolID);
   for(pargrid::CellID b=0; b<boundaryBlocks.size(); ++b) { cell2Node(cellJ,nodeJ,sim,simClasses,boundaryBlocks[b]); }
   profile::stop();
#endif

   // calculate cellUe
   profile::start("field propag",profPropagFieldID);
   for(pargrid::CellID b=0; b<simClasses.pargrid.getNumberOfLocalCells(); ++b) { calcCellUe(cellJ,cellJi,cellRhoQi,cellUe,innerFlag,counterCellMaxUe,sim,simClasses,b); }
   profile::stop();
   
   // loop thru ghost cells
   profile::start("BoundaryConds",profBoundCondsID);
   simClasses.pargrid.startNeighbourExchange(pargrid::DEFAULT_STENCIL,Hybrid::dataCellUeID);
   profile::start("MPI waits",mpiWaitID);
   simClasses.pargrid.wait(pargrid::DEFAULT_STENCIL,Hybrid::dataCellUeID);
   profile::stop();
   neumannCell(cellUe,sim,simClasses,exteriorBlocks,3);
   profile::stop();

   // calculate nodeUe
#ifdef USE_NODE_UE
   // cell -> node RhoQi and Ji
   /*simClasses.pargrid.startNeighbourExchange(pargrid::DEFAULT_STENCIL,Hybrid::dataCellRhoQiID);
   simClasses.pargrid.startNeighbourExchange(pargrid::DEFAULT_STENCIL,Hybrid::dataCellJiID);
   profile::start("intpol",profIntpolID);
   for(pargrid::CellID b=0; b<innerBlocks.size(); ++b) { cell2Node(cellRhoQi,nodeRhoQi,sim,simClasses,innerBlocks[b],1); }
   for(pargrid::CellID b=0; b<innerBlocks.size(); ++b) { cell2Node(cellJi,nodeJi,sim,simClasses,innerBlocks[b]); }
   profile::stop();
   profile::start("MPI waits",mpiWaitID);   
   simClasses.pargrid.wait(pargrid::DEFAULT_STENCIL,Hybrid::dataCellRhoQiID);
   simClasses.pargrid.wait(pargrid::DEFAULT_STENCIL,Hybrid::dataCellJiID);
   profile::stop();
   profile::start("intpol",profIntpolID);
   for(pargrid::CellID b=0; b<boundaryBlocks.size(); ++b) { cell2Node(cellRhoQi,nodeRhoQi,sim,simClasses,boundaryBlocks[b],1); }
   for(pargrid::CellID b=0; b<boundaryBlocks.size(); ++b) { cell2Node(cellJi,nodeJi,sim,simClasses,boundaryBlocks[b]); }
   profile::stop();*/
   profile::start("field propag",profPropagFieldID);
   for(pargrid::CellID b=0; b<simClasses.pargrid.getNumberOfLocalCells(); ++b) { calcNodeUe(nodeRhoQi,nodeJi,nodeJ,nodeUe,innerFlagNode,counterCellMaxUe,sim,simClasses,b); }
   profile::stop();
#else
   // cell->node Ue
   simClasses.pargrid.startNeighbourExchange(pargrid::DEFAULT_STENCIL,Hybrid::dataCellUeID);
   profile::start("intpol",profIntpolID);
   for(pargrid::CellID b=0; b<innerBlocks.size(); ++b) { cell2Node(cellUe,nodeUe,sim,simClasses,innerBlocks[b]); }
   profile::stop();
   profile::start("MPI waits",mpiWaitID);
   simClasses.pargrid.wait(pargrid::DEFAULT_STENCIL,Hybrid::dataCellUeID);
   profile::stop();
   profile::start("intpol",profIntpolID);
   for(pargrid::CellID b=0; b<boundaryBlocks.size(); ++b) { cell2Node(cellUe,nodeUe,sim,simClasses,boundaryBlocks[b]); }
   profile::stop();
#endif
   
   // upwind nodeB
   profile::start("field propag",profPropagFieldID);
   for(pargrid::CellID b=0; b<simClasses.pargrid.getNumberOfLocalCells(); ++b) { upwindNodeB(cellB,nodeUe,nodeB,sim,simClasses,b); }
   profile::stop();
   
   // calculate nodeE
   profile::start("field propag",profPropagFieldID);
   for(pargrid::CellID b=0; b<simClasses.pargrid.getNumberOfLocalCells(); ++b) {
      calcNodeE(nodeUe,nodeB,
#ifdef USE_RESISTIVITY
      nodeEta,
#endif
      nodeJ,nodeE,
#ifdef USE_ECUT
      counterNodeEcut,
#endif
      innerFlagNode,sim,simClasses,b);
   }
   profile::stop();

   // nodeE filter
   for(int i=0;i<Hybrid::Efilter;i++) {
      // node->cell E
      simClasses.pargrid.startNeighbourExchange(pargrid::DEFAULT_STENCIL,Hybrid::dataNodeEID);
      profile::start("intpol",profIntpolID);
      for(pargrid::CellID b=0; b<innerBlocks.size(); ++b) { node2Cell(nodeE,cellJ,sim,simClasses,innerBlocks[b]); }
      profile::stop();
      profile::start("MPI waits",mpiWaitID);
      simClasses.pargrid.wait(pargrid::DEFAULT_STENCIL,Hybrid::dataNodeEID);
      profile::stop();
      profile::start("intpol",profIntpolID);
      for(pargrid::CellID b=0; b<boundaryBlocks.size(); ++b) { node2Cell(nodeE,cellJ,sim,simClasses,boundaryBlocks[b]); }
      profile::stop();
      // Neumann boundary conditions
      simClasses.pargrid.startNeighbourExchange(pargrid::DEFAULT_STENCIL,Hybrid::dataCellJID);
      profile::start("MPI waits",mpiWaitID);
      simClasses.pargrid.wait(pargrid::DEFAULT_STENCIL,Hybrid::dataCellJID);
      profile::stop();
      neumannCell(cellJ,sim,simClasses,exteriorBlocks,3); 
      // cell->node E
      simClasses.pargrid.startNeighbourExchange(pargrid::DEFAULT_STENCIL,Hybrid::dataCellJID);
      profile::start("intpol",profIntpolID);
      for(pargrid::CellID b=0; b<innerBlocks.size(); ++b) { cell2Node(cellJ,nodeE,sim,simClasses,innerBlocks[b]); }
      profile::stop();
      profile::start("MPI waits",mpiWaitID);
      simClasses.pargrid.wait(pargrid::DEFAULT_STENCIL,Hybrid::dataCellJID);
      profile::stop();
      profile::start("intpol",profIntpolID);
      for(pargrid::CellID b=0; b<boundaryBlocks.size(); ++b) { cell2Node(cellJ,nodeE,sim,simClasses,boundaryBlocks[b]); }
      profile::stop();
      // zero cellJ
      for(pargrid::CellID b=0;b<simClasses.pargrid.getNumberOfAllCells();++b) for(int k=0;k<block::WIDTH_Z;++k) for(int j=0;j<block::WIDTH_Y;++j) for(int i=0;i<block::WIDTH_X;++i) {
	 const int n3 = (b*block::SIZE+block::index(i,j,k))*3;
	 for(int l=0;l<3;++l) {
	    cellJ[n3+l]  = 0.0;
	 }
      }
   }
   // nodeE gaussian filter
   if(Hybrid::EfilterNodeGaussSigma > 0) {
      simClasses.pargrid.startNeighbourExchange(pargrid::DEFAULT_STENCIL,Hybrid::dataNodeEID);
      simClasses.pargrid.wait(pargrid::DEFAULT_STENCIL,Hybrid::dataNodeEID);
      const size_t N = simClasses.pargrid.getNumberOfAllCells()*simClasses.pargrid.getUserDataStaticElements(Hybrid::dataNodeEID);
      Real* nodeEOld = new Real[N];
      *nodeEOld = *nodeE;
      for(size_t i=0;i<N;++i) { nodeEOld[i] = nodeE[i]; }
      for(pargrid::CellID b=0; b<innerBlocks.size(); ++b) { nodeAvg(nodeEOld,nodeE,sim,simClasses,innerBlocks[b]); }
      for(pargrid::CellID b=0; b<boundaryBlocks.size(); ++b) { nodeAvg(nodeEOld,nodeE,sim,simClasses,boundaryBlocks[b]); }
      delete [] nodeEOld;
      nodeEOld = NULL;
   }

   // propagate faceB by Faraday's law using nodeE
   simClasses.pargrid.startNeighbourExchange(pargrid::DEFAULT_STENCIL,Hybrid::dataNodeEID);
   profile::start("field propag",profPropagFieldID);
   for(pargrid::CellID b=0; b<innerBlocks.size(); ++b) { faceCurl(nodeE,faceB,true,sim,simClasses,innerBlocks[b]); }
   profile::stop();
   profile::start("MPI waits",mpiWaitID);
   simClasses.pargrid.wait(pargrid::DEFAULT_STENCIL,Hybrid::dataNodeEID);
   profile::stop();
   profile::start("field propag",profPropagFieldID);
   for(pargrid::CellID b=0; b<boundaryBlocks.size(); ++b) { faceCurl(nodeE,faceB,true,sim,simClasses,boundaryBlocks[b]); }
   profile::stop();

   // calculated cellEp
   if(Hybrid::useElectronPressureElectricField == true) {
      simClasses.pargrid.startNeighbourExchange(pargrid::DEFAULT_STENCIL,Hybrid::dataNodeRhoQiID);
      profile::start("intpol",profIntpolID);
      for(pargrid::CellID b=0; b<innerBlocks.size(); ++b) { calcCellEp(nodeRhoQi,cellRhoQi,innerFlagCellEp,cellEp,sim,simClasses,innerBlocks[b]); }
      profile::stop();
      profile::start("MPI waits",mpiWaitID);
      simClasses.pargrid.wait(pargrid::DEFAULT_STENCIL,Hybrid::dataNodeRhoQiID);
      profile::stop();
      profile::start("intpol",profIntpolID);
      for(pargrid::CellID b=0; b<boundaryBlocks.size(); ++b) { calcCellEp(nodeRhoQi,cellRhoQi,innerFlagCellEp,cellEp,sim,simClasses,boundaryBlocks[b]); }
      profile::stop();
   }
   
   if(sim.atDataSaveStep == true) {
      saveStepHappened = true;
   }

   profile::stop();
   return success;
}

// neumann zero-gradient boundary conditions
void neumannCell(Real* cellData,Simulation& sim,SimulationClasses& simClasses,const vector<pargrid::CellID>& exteriorBlocks,const int vectorDim)
{
   const unsigned int tempArraySize = (block::WIDTH_X+2)*(block::WIDTH_Y+2)*(block::WIDTH_Z+2);
   Real* tempArrayCellData = new Real[tempArraySize*vectorDim];
   const std::vector<uint32_t>& neighbourFlags = simClasses.pargrid.getNeighbourFlags();
   
   for(pargrid::CellID eb=0;eb<exteriorBlocks.size();++eb) {
      const pargrid::CellID b = exteriorBlocks[eb];
      fetchData(cellData,tempArrayCellData,simClasses,b,vectorDim);
      int di=0; if(block::WIDTH_X > 1) { di=block::WIDTH_X-1; }
      int dj=0; if(block::WIDTH_Y > 1) { dj=block::WIDTH_Y-1; }
      int dk=0; if(block::WIDTH_Z > 1) { dk=block::WIDTH_Z-1; }
      const uint32_t& nf = neighbourFlags[b];
      // back (-x) wall
      if((nf & Hybrid::X_NEG_EXISTS) == 0) {
	 for(int k=0; k<block::WIDTH_Z; ++k) for(int j=0; j<block::WIDTH_Y; ++j) {
	    const int i = 0+di;
	    const int n = (b*block::SIZE+block::index(i,j,k))*vectorDim;
	    const int m = block::arrayIndex(i+2,j+1,k+1)*vectorDim;
	    for(int l=0;l<vectorDim;++l) { cellData[n+l] = tempArrayCellData[m+l]; }
	 }
      }
      // side (+y) wall
      if((nf & Hybrid::Y_POS_EXISTS) == 0) {
	 for(int k=0; k<block::WIDTH_Z; ++k) for(int i=0; i<block::WIDTH_X; ++i) {
	    const int j = block::WIDTH_Y-1-dj;
	    const int n = (b*block::SIZE+block::index(i,j,k))*vectorDim;
	    const int m = block::arrayIndex(i+1,j+0,k+1)*vectorDim;
	    for(int l=0;l<vectorDim;++l) { cellData[n+l] = tempArrayCellData[m+l]; }
	 }
      }
      // side (-y) wall
      if((nf & Hybrid::Y_NEG_EXISTS) == 0) {
	 for(int k=0; k<block::WIDTH_Z; ++k) for(int i=0; i<block::WIDTH_X; ++i) {
	    const int j = 0+dj;
	    const int n = (b*block::SIZE+block::index(i,j,k))*vectorDim;
	    const int m = block::arrayIndex(i+1,j+2,k+1)*vectorDim;
	    for(int l=0;l<vectorDim;++l) { cellData[n+l] = tempArrayCellData[m+l]; }
	 }
      }
      // side (+z) wall
      if((nf & Hybrid::Z_POS_EXISTS) == 0) {
	 for(int j=0; j<block::WIDTH_Y; ++j) for(int i=0; i<block::WIDTH_X; ++i) {
	    const int k = block::WIDTH_Z-1-dk;
	    const int n = (b*block::SIZE+block::index(i,j,k))*vectorDim;
	    const int m = block::arrayIndex(i+1,j+1,k+0)*vectorDim;
	    for(int l=0;l<vectorDim;++l) { cellData[n+l] = tempArrayCellData[m+l]; }
	 }
      }
      // side (-z) wall
      if((nf & Hybrid::Z_NEG_EXISTS) == 0) {
	 for(int j=0; j<block::WIDTH_Y; ++j) for(int i=0; i<block::WIDTH_X; ++i) {
	    const int k = 0+dk;
	    const int n = (b*block::SIZE+block::index(i,j,k))*vectorDim;
	    const int m = block::arrayIndex(i+1,j+1,k+2)*vectorDim;
	    for(int l=0;l<vectorDim;++l) { cellData[n+l] = tempArrayCellData[m+l]; }
	 }
      }
      // front (+x) wall
      if((nf & Hybrid::X_POS_EXISTS) == 0) {
	 for(int k=0; k<block::WIDTH_Z; ++k) for(int j=0; j<block::WIDTH_Y; ++j) {
	    const int i = 0;
	    const int n = (b*block::SIZE+block::index(i,j,k))*vectorDim;
	    const int m = block::arrayIndex(i+0,j+1,k+1)*vectorDim;
	    for(int l=0;l<vectorDim;++l) { cellData[n+l] = tempArrayCellData[m+l]; }
	 }
      }
      
      // NOTE: block indices not correct for these yet

      // edge (-x,-y)
      if((nf & Hybrid::X_NEG_EXISTS) == 0 && (nf & Hybrid::Y_NEG_EXISTS) == 0) {
	 const int i=0,j=0,k=0;
	 const int n = (b*block::SIZE+block::index(i,j,k))*vectorDim;
	 const int m = block::arrayIndex(i+2,j+2,k+1)*vectorDim;
	 for(int l=0;l<vectorDim;++l) { cellData[n+l] = tempArrayCellData[m+l]; }
      }
      // edge (-x,+y)
      if((nf & Hybrid::X_NEG_EXISTS) == 0 && (nf & Hybrid::Y_POS_EXISTS) == 0) {
	 const int i=0,j=0,k=0;
	 const int n = (b*block::SIZE+block::index(i,j,k))*vectorDim;
	 const int m = block::arrayIndex(i+2,j+0,k+1)*vectorDim;
	 for(int l=0;l<vectorDim;++l) { cellData[n+l] = tempArrayCellData[m+l]; }
      }
      // edge (+x,-y)
      if((nf & Hybrid::X_POS_EXISTS) == 0 && (nf & Hybrid::Y_NEG_EXISTS) == 0) {
	 const int i=0,j=0,k=0;
	 const int n = (b*block::SIZE+block::index(i,j,k))*vectorDim;
	 const int m = block::arrayIndex(i+0,j+2,k+1)*vectorDim;
	 for(int l=0;l<vectorDim;++l) { cellData[n+l] = tempArrayCellData[m+l]; }
      }
      // edge (+x,+y)
      if((nf & Hybrid::X_POS_EXISTS) == 0 && (nf & Hybrid::Y_POS_EXISTS) == 0) {
	 const int i=0,j=0,k=0;
	 const int n = (b*block::SIZE+block::index(i,j,k))*vectorDim;
	 const int m = block::arrayIndex(i+0,j+0,k+1)*vectorDim;
	 for(int l=0;l<vectorDim;++l) { cellData[n+l] = tempArrayCellData[m+l]; }
      }
      // edge (-x,-z)
      if((nf & Hybrid::X_NEG_EXISTS) == 0 && (nf & Hybrid::Z_NEG_EXISTS) == 0) {
	 const int i=0,j=0,k=0;
	 const int n = (b*block::SIZE+block::index(i,j,k))*vectorDim;
	 const int m = block::arrayIndex(i+2,j+1,k+2)*vectorDim;
	 for(int l=0;l<vectorDim;++l) { cellData[n+l] = tempArrayCellData[m+l]; }
      }
      // edge (-x,+z)
      if((nf & Hybrid::X_NEG_EXISTS) == 0 && (nf & Hybrid::Z_POS_EXISTS) == 0) {
	 const int i=0,j=0,k=0;
	 const int n = (b*block::SIZE+block::index(i,j,k))*vectorDim;
	 const int m = block::arrayIndex(i+2,j+1,k+0)*vectorDim;
	 for(int l=0;l<vectorDim;++l) { cellData[n+l] = tempArrayCellData[m+l]; }
      }
      // edge (+x,-z)
      if((nf & Hybrid::X_POS_EXISTS) == 0 && (nf & Hybrid::Z_NEG_EXISTS) == 0) {
	 const int i=0,j=0,k=0;
	 const int n = (b*block::SIZE+block::index(i,j,k))*vectorDim;
	 const int m = block::arrayIndex(i+0,j+1,k+2)*vectorDim;
	 for(int l=0;l<vectorDim;++l) { cellData[n+l] = tempArrayCellData[m+l]; }
      }
      // edge (+x,+z)
      if((nf & Hybrid::X_POS_EXISTS) == 0 && (nf & Hybrid::Z_POS_EXISTS) == 0) {
	 const int i=0,j=0,k=0;
	 const int n = (b*block::SIZE+block::index(i,j,k))*vectorDim;
	 const int m = block::arrayIndex(i+0,j+1,k+0)*vectorDim;
	 for(int l=0;l<vectorDim;++l) { cellData[n+l] = tempArrayCellData[m+l]; }
      }
      // edge (-y,-z)
      if((nf & Hybrid::Y_NEG_EXISTS) == 0 && (nf & Hybrid::Z_NEG_EXISTS) == 0) {
	 const int i=0,j=0,k=0;
	 const int n = (b*block::SIZE+block::index(i,j,k))*vectorDim;
	 const int m = block::arrayIndex(i+1,j+2,k+2)*vectorDim;
	 for(int l=0;l<vectorDim;++l) { cellData[n+l] = tempArrayCellData[m+l]; }
      }
      // edge (-y,+z)
      if((nf & Hybrid::Y_NEG_EXISTS) == 0 && (nf & Hybrid::Z_POS_EXISTS) == 0) {
	 const int i=0,j=0,k=0;
	 const int n = (b*block::SIZE+block::index(i,j,k))*vectorDim;
	 const int m = block::arrayIndex(i+1,j+2,k+0)*vectorDim;
	 for(int l=0;l<vectorDim;++l) { cellData[n+l] = tempArrayCellData[m+l]; }
      }
      // edge (+y,-z)
      if((nf & Hybrid::Y_POS_EXISTS) == 0 && (nf & Hybrid::Z_NEG_EXISTS) == 0) {
	 const int i=0,j=0,k=0;
	 const int n = (b*block::SIZE+block::index(i,j,k))*vectorDim;
	 const int m = block::arrayIndex(i+1,j+0,k+2)*vectorDim;
	 for(int l=0;l<vectorDim;++l) { cellData[n+l] = tempArrayCellData[m+l]; }
      }
      // edge (+y,+z)
      if((nf & Hybrid::Y_POS_EXISTS) == 0 && (nf & Hybrid::Z_POS_EXISTS) == 0) {
	 const int i=0,j=0,k=0;
	 const int n = (b*block::SIZE+block::index(i,j,k))*vectorDim;
	 const int m = block::arrayIndex(i+1,j+0,k+0)*vectorDim;
	 for(int l=0;l<vectorDim;++l) { cellData[n+l] = tempArrayCellData[m+l]; }
      }
      // corner (-x,-y,-z)
      if((nf & Hybrid::X_NEG_EXISTS) == 0 && (nf & Hybrid::Y_NEG_EXISTS) == 0 && (nf & Hybrid::Z_NEG_EXISTS) == 0) {
	 const int i=0,j=0,k=0;
	 const int n = (b*block::SIZE+block::index(i,j,k))*vectorDim;
	 const int m = block::arrayIndex(i+2,j+2,k+2)*vectorDim;
	 for(int l=0;l<vectorDim;++l) { cellData[n+l] = tempArrayCellData[m+l]; }
      }
      // corner (-x,-y,+z)
      if((nf & Hybrid::X_NEG_EXISTS) == 0 && (nf & Hybrid::Y_NEG_EXISTS) == 0 && (nf & Hybrid::Z_POS_EXISTS) == 0) {
	 const int i=0,j=0,k=0;
	 const int n = (b*block::SIZE+block::index(i,j,k))*vectorDim;
	 const int m = block::arrayIndex(i+2,j+2,k+0)*vectorDim;
	 for(int l=0;l<vectorDim;++l) { cellData[n+l] = tempArrayCellData[m+l]; }
      }
      // corner (-x,+y,-z)
      if((nf & Hybrid::X_NEG_EXISTS) == 0 && (nf & Hybrid::Y_POS_EXISTS) == 0 && (nf & Hybrid::Z_NEG_EXISTS) == 0) {
	 const int i=0,j=0,k=0;
	 const int n = (b*block::SIZE+block::index(i,j,k))*vectorDim;
	 const int m = block::arrayIndex(i+2,j+0,k+2)*vectorDim;
	 for(int l=0;l<vectorDim;++l) { cellData[n+l] = tempArrayCellData[m+l]; }
      }
      // corner (-x,+y,+z)
      if((nf & Hybrid::X_NEG_EXISTS) == 0 && (nf & Hybrid::Y_POS_EXISTS) == 0 && (nf & Hybrid::Z_POS_EXISTS) == 0) {
	 const int i=0,j=0,k=0;
	 const int n = (b*block::SIZE+block::index(i,j,k))*vectorDim;
	 const int m = block::arrayIndex(i+2,j+0,k+0)*vectorDim;
	 for(int l=0;l<vectorDim;++l) { cellData[n+l] = tempArrayCellData[m+l]; }
      }
      // corner (+x,-y,-z)
      if((nf & Hybrid::X_POS_EXISTS) == 0 && (nf & Hybrid::Y_NEG_EXISTS) == 0 && (nf & Hybrid::Z_NEG_EXISTS) == 0) {
	 const int i=0,j=0,k=0;
	 const int n = (b*block::SIZE+block::index(i,j,k))*vectorDim;
	 const int m = block::arrayIndex(i+0,j+2,k+2)*vectorDim;
	 for(int l=0;l<vectorDim;++l) { cellData[n+l] = tempArrayCellData[m+l]; }
      }
      // corner (+x,-y,+z)
      if((nf & Hybrid::X_POS_EXISTS) == 0 && (nf & Hybrid::Y_NEG_EXISTS) == 0 && (nf & Hybrid::Z_POS_EXISTS) == 0) {
	 const int i=0,j=0,k=0;
	 const int n = (b*block::SIZE+block::index(i,j,k))*vectorDim;
	 const int m = block::arrayIndex(i+0,j+2,k+0)*vectorDim;
	 for(int l=0;l<vectorDim;++l) { cellData[n+l] = tempArrayCellData[m+l]; }
      }
      // corner (+x,+y,-z)
      if((nf & Hybrid::X_POS_EXISTS) == 0 && (nf & Hybrid::Y_POS_EXISTS) == 0 && (nf & Hybrid::Z_NEG_EXISTS) == 0) {
	 const int i=0,j=0,k=0;
	 const int n = (b*block::SIZE+block::index(i,j,k))*vectorDim;
	 const int m = block::arrayIndex(i+0,j+0,k+2)*vectorDim;
	 for(int l=0;l<vectorDim;++l) { cellData[n+l] = tempArrayCellData[m+l]; }
      }
      // corner (+x,+y,+z)
      if((nf & Hybrid::X_POS_EXISTS) == 0 && (nf & Hybrid::Y_POS_EXISTS) == 0 && (nf & Hybrid::Z_POS_EXISTS) == 0) {
	 const int i=0,j=0,k=0;
	 const int n = (b*block::SIZE+block::index(i,j,k))*vectorDim;
	 const int m = block::arrayIndex(i+0,j+0,k+0)*vectorDim;
	 for(int l=0;l<vectorDim;++l) { cellData[n+l] = tempArrayCellData[m+l]; }
      }
   }
   delete [] tempArrayCellData; tempArrayCellData = NULL;
}


// neumann zero-gradient boundary conditions
void neumannFace(Real* faceData,Simulation& sim,SimulationClasses& simClasses,const vector<pargrid::CellID>& exteriorBlocks)
{
   const int vectorDim = 3;
   const unsigned int tempArraySize = (block::WIDTH_X+2)*(block::WIDTH_Y+2)*(block::WIDTH_Z+2);
   Real* tempArrayFaceData = new Real[tempArraySize*vectorDim];
   const std::vector<uint32_t>& neighbourFlags = simClasses.pargrid.getNeighbourFlags();
   
   for(pargrid::CellID eb=0;eb<exteriorBlocks.size();++eb) {
      const pargrid::CellID b = exteriorBlocks[eb];
      fetchData(faceData,tempArrayFaceData,simClasses,b,vectorDim);
      int di=0; if(block::WIDTH_X > 1) { di=block::WIDTH_X-1; }
      int dj=0; if(block::WIDTH_Y > 1) { dj=block::WIDTH_Y-1; }
      int dk=0; if(block::WIDTH_Z > 1) { dk=block::WIDTH_Z-1; }
      const uint32_t& nf = neighbourFlags[b];
      // back (-x) wall
      if((nf & Hybrid::X_NEG_EXISTS) == 0) {
	 for(int k=0; k<block::WIDTH_Z; ++k) for(int j=0; j<block::WIDTH_Y; ++j) {
	    const int i = 0+di;
	    const int n = (b*block::SIZE+block::index(i,j,k))*vectorDim;
	    const int m = block::arrayIndex(i+2,j+1,k+1)*vectorDim;
	    faceData[n+1] = tempArrayFaceData[m+1];
            faceData[n+2] = tempArrayFaceData[m+2];
	 }
      }
      // side (+y) wall
      if((nf & Hybrid::Y_POS_EXISTS) == 0) {
	 for(int k=0; k<block::WIDTH_Z; ++k) for(int i=0; i<block::WIDTH_X; ++i) {
	    const int j = block::WIDTH_Y-1-dj;
	    const int n = (b*block::SIZE+block::index(i,j,k))*vectorDim;
	    const int m = block::arrayIndex(i+1,j+0,k+1)*vectorDim;
	    faceData[n+0] = tempArrayFaceData[m+0];
            faceData[n+2] = tempArrayFaceData[m+2];
	 }
      }
      // side (-y) wall
      if((nf & Hybrid::Y_NEG_EXISTS) == 0) {
	 for(int k=0; k<block::WIDTH_Z; ++k) for(int i=0; i<block::WIDTH_X; ++i) {
	    const int j = 0+dj;
	    const int n = (b*block::SIZE+block::index(i,j,k))*vectorDim;
	    const int m = block::arrayIndex(i+1,j+2,k+1)*vectorDim;
	    faceData[n+0] = tempArrayFaceData[m+0];
            faceData[n+2] = tempArrayFaceData[m+2];
	 }
      }
      // side (+z) wall
      if((nf & Hybrid::Z_POS_EXISTS) == 0) {
	 for(int j=0; j<block::WIDTH_Y; ++j) for(int i=0; i<block::WIDTH_X; ++i) {
	    const int k = block::WIDTH_Z-1-dk;
	    const int n = (b*block::SIZE+block::index(i,j,k))*vectorDim;
	    const int m = block::arrayIndex(i+1,j+1,k+0)*vectorDim;
	    faceData[n+0] = tempArrayFaceData[m+0];
            faceData[n+1] = tempArrayFaceData[m+1];
	 }
      }
      // side (-z) wall
      if((nf & Hybrid::Z_NEG_EXISTS) == 0) {
	 for(int j=0; j<block::WIDTH_Y; ++j) for(int i=0; i<block::WIDTH_X; ++i) {
	    const int k = 0+dk;
	    const int n = (b*block::SIZE+block::index(i,j,k))*vectorDim;
	    const int m = block::arrayIndex(i+1,j+1,k+2)*vectorDim;
	    faceData[n+0] = tempArrayFaceData[m+0];
            faceData[n+1] = tempArrayFaceData[m+1];
	 }
      }
      // front (+x) wall
      if((nf & Hybrid::X_POS_EXISTS) == 0) {
	 for(int k=0; k<block::WIDTH_Z; ++k) for(int j=0; j<block::WIDTH_Y; ++j) {
	    const int i = 0;
	    const int n = (b*block::SIZE+block::index(i,j,k))*vectorDim;
	    const int m = block::arrayIndex(i+0,j+1,k+1)*vectorDim;
	    faceData[n+1] = tempArrayFaceData[m+1];
            faceData[n+2] = tempArrayFaceData[m+2];
	 }
      }
      // NOTE: block indices not correct for these yet

      // edge (-x,-y)
      if((nf & Hybrid::X_NEG_EXISTS) == 0 && (nf & Hybrid::Y_NEG_EXISTS) == 0) {
	 const int i=0,j=0,k=0;
	 const int n = (b*block::SIZE+block::index(i,j,k))*vectorDim;
	 const int m121 = block::arrayIndex(i+1,j+2,k+1)*vectorDim;
         const int m211 = block::arrayIndex(i+2,j+1,k+1)*vectorDim;
         const int m221 = block::arrayIndex(i+2,j+2,k+1)*vectorDim;
         faceData[n+0] = tempArrayFaceData[m121+0];
         faceData[n+1] = tempArrayFaceData[m211+1];
         faceData[n+2] = tempArrayFaceData[m221+2];
      }
      // edge (-x,+y)
      if((nf & Hybrid::X_NEG_EXISTS) == 0 && (nf & Hybrid::Y_POS_EXISTS) == 0) {
	 const int i=0,j=0,k=0;
	 const int n = (b*block::SIZE+block::index(i,j,k))*vectorDim;
         const int m101 = block::arrayIndex(i+1,j+0,k+1)*vectorDim;
         const int m201 = block::arrayIndex(i+2,j+0,k+1)*vectorDim;
         faceData[n+0] = tempArrayFaceData[m101+0];
         faceData[n+2] = tempArrayFaceData[m201+2];
      }
      // edge (+x,-y)
      if((nf & Hybrid::X_POS_EXISTS) == 0 && (nf & Hybrid::Y_NEG_EXISTS) == 0) {
	 const int i=0,j=0,k=0;
	 const int n = (b*block::SIZE+block::index(i,j,k))*vectorDim;
	 const int m011 = block::arrayIndex(i+0,j+1,k+1)*vectorDim;
         const int m021 = block::arrayIndex(i+0,j+2,k+1)*vectorDim;
         faceData[n+1] = tempArrayFaceData[m011+1];
         faceData[n+2] = tempArrayFaceData[m021+2];
      }
      // edge (+x,+y)
      if((nf & Hybrid::X_POS_EXISTS) == 0 && (nf & Hybrid::Y_POS_EXISTS) == 0) {
	 const int i=0,j=0,k=0;
	 const int n = (b*block::SIZE+block::index(i,j,k))*vectorDim;
	 const int m001 = block::arrayIndex(i+0,j+0,k+1)*vectorDim;
         faceData[n+2] = tempArrayFaceData[m001+2];
      }
      // edge (-x,-z)
      if((nf & Hybrid::X_NEG_EXISTS) == 0 && (nf & Hybrid::Z_NEG_EXISTS) == 0) {
	 const int i=0,j=0,k=0;
	 const int n = (b*block::SIZE+block::index(i,j,k))*vectorDim;
         const int m112 = block::arrayIndex(i+1,j+1,k+2)*vectorDim;
	 const int m212 = block::arrayIndex(i+2,j+1,k+2)*vectorDim;
         const int m211 = block::arrayIndex(i+2,j+1,k+1)*vectorDim;
         faceData[n+0] = tempArrayFaceData[m112+0];
         faceData[n+1] = tempArrayFaceData[m212+1];
         faceData[n+2] = tempArrayFaceData[m211+2];
      }
      // edge (-x,+z)
      if((nf & Hybrid::X_NEG_EXISTS) == 0 && (nf & Hybrid::Z_POS_EXISTS) == 0) {
	 const int i=0,j=0,k=0;
	 const int n = (b*block::SIZE+block::index(i,j,k))*vectorDim;
         const int m110 = block::arrayIndex(i+1,j+1,k+0)*vectorDim;
	 const int m210 = block::arrayIndex(i+2,j+1,k+0)*vectorDim;
         faceData[n+0] = tempArrayFaceData[m110+0];
         faceData[n+1] = tempArrayFaceData[m210+1];
      }
      // edge (+x,-z)
      if((nf & Hybrid::X_POS_EXISTS) == 0 && (nf & Hybrid::Z_NEG_EXISTS) == 0) {
	 const int i=0,j=0,k=0;
	 const int n = (b*block::SIZE+block::index(i,j,k))*vectorDim;
         const int m012 = block::arrayIndex(i+0,j+1,k+2)*vectorDim;
         const int m011 = block::arrayIndex(i+0,j+1,k+1)*vectorDim;
         faceData[n+1] = tempArrayFaceData[m012+1];
         faceData[n+2] = tempArrayFaceData[m011+2];
      }
      // edge (+x,+z)
      if((nf & Hybrid::X_POS_EXISTS) == 0 && (nf & Hybrid::Z_POS_EXISTS) == 0) {
	 const int i=0,j=0,k=0;
	 const int n = (b*block::SIZE+block::index(i,j,k))*vectorDim;
	 const int m010 = block::arrayIndex(i+0,j+1,k+0)*vectorDim;
         faceData[n+1] = tempArrayFaceData[m010+1];
      }
      // edge (-y,-z)
      if((nf & Hybrid::Y_NEG_EXISTS) == 0 && (nf & Hybrid::Z_NEG_EXISTS) == 0) {
	 const int i=0,j=0,k=0;
	 const int n = (b*block::SIZE+block::index(i,j,k))*vectorDim;
	 const int m122 = block::arrayIndex(i+1,j+2,k+2)*vectorDim;
         const int m112 = block::arrayIndex(i+1,j+1,k+2)*vectorDim;
         const int m121 = block::arrayIndex(i+1,j+2,k+1)*vectorDim;
         faceData[n+0] = tempArrayFaceData[m122+0];
         faceData[n+1] = tempArrayFaceData[m112+1];
         faceData[n+2] = tempArrayFaceData[m121+2];
      }
      // edge (-y,+z)
      if((nf & Hybrid::Y_NEG_EXISTS) == 0 && (nf & Hybrid::Z_POS_EXISTS) == 0) {
	 const int i=0,j=0,k=0;
	 const int n = (b*block::SIZE+block::index(i,j,k))*vectorDim;
	 const int m120 = block::arrayIndex(i+1,j+2,k+0)*vectorDim;
         const int m110 = block::arrayIndex(i+1,j+1,k+0)*vectorDim;
         faceData[n+0] = tempArrayFaceData[m120+0];
         faceData[n+1] = tempArrayFaceData[m110+1];
      }
      // edge (+y,-z)
      if((nf & Hybrid::Y_POS_EXISTS) == 0 && (nf & Hybrid::Z_NEG_EXISTS) == 0) {
	 const int i=0,j=0,k=0;
	 const int n = (b*block::SIZE+block::index(i,j,k))*vectorDim;
	 const int m102 = block::arrayIndex(i+1,j+0,k+2)*vectorDim;
         const int m101 = block::arrayIndex(i+1,j+0,k+1)*vectorDim;
         faceData[n+0] = tempArrayFaceData[m102+0];
         faceData[n+2] = tempArrayFaceData[m101+2];
      }
      // edge (+y,+z)
      if((nf & Hybrid::Y_POS_EXISTS) == 0 && (nf & Hybrid::Z_POS_EXISTS) == 0) {
	 const int i=0,j=0,k=0;
	 const int n = (b*block::SIZE+block::index(i,j,k))*vectorDim;
	 const int m100 = block::arrayIndex(i+1,j+0,k+0)*vectorDim;
         faceData[n+0] = tempArrayFaceData[m100+0];
      }
      // corner (-x,-y,-z)
      if((nf & Hybrid::X_NEG_EXISTS) == 0 && (nf & Hybrid::Y_NEG_EXISTS) == 0 && (nf & Hybrid::Z_NEG_EXISTS) == 0) {
	 const int i=0,j=0,k=0;
	 const int n = (b*block::SIZE+block::index(i,j,k))*vectorDim;
	 const int m122 = block::arrayIndex(i+1,j+2,k+2)*vectorDim;
         const int m212 = block::arrayIndex(i+2,j+1,k+2)*vectorDim;
         const int m221 = block::arrayIndex(i+2,j+2,k+1)*vectorDim;
         faceData[n+0] = tempArrayFaceData[m122+0];
         faceData[n+1] = tempArrayFaceData[m212+1];
         faceData[n+2] = tempArrayFaceData[m221+2];
      }
      // corner (-x,-y,+z)
      if((nf & Hybrid::X_NEG_EXISTS) == 0 && (nf & Hybrid::Y_NEG_EXISTS) == 0 && (nf & Hybrid::Z_POS_EXISTS) == 0) {
	 const int i=0,j=0,k=0;
	 const int n = (b*block::SIZE+block::index(i,j,k))*vectorDim;
	 const int m120 = block::arrayIndex(i+1,j+2,k+0)*vectorDim;
         const int m210 = block::arrayIndex(i+2,j+1,k+0)*vectorDim;
         faceData[n+0] = tempArrayFaceData[m120+0];
         faceData[n+1] = tempArrayFaceData[m210+1];
      }
      // corner (-x,+y,-z)
      if((nf & Hybrid::X_NEG_EXISTS) == 0 && (nf & Hybrid::Y_POS_EXISTS) == 0 && (nf & Hybrid::Z_NEG_EXISTS) == 0) {
	 const int i=0,j=0,k=0;
	 const int n = (b*block::SIZE+block::index(i,j,k))*vectorDim;
	 const int m102 = block::arrayIndex(i+1,j+0,k+2)*vectorDim;
         const int m201 = block::arrayIndex(i+2,j+0,k+1)*vectorDim;
         faceData[n+0] = tempArrayFaceData[m102+0];
         faceData[n+2] = tempArrayFaceData[m201+2];
      }
      // corner (-x,+y,+z)
      if((nf & Hybrid::X_NEG_EXISTS) == 0 && (nf & Hybrid::Y_POS_EXISTS) == 0 && (nf & Hybrid::Z_POS_EXISTS) == 0) {
	 const int i=0,j=0,k=0;
	 const int n = (b*block::SIZE+block::index(i,j,k))*vectorDim;
	 const int m100 = block::arrayIndex(i+1,j+0,k+0)*vectorDim;
         faceData[n+0] = tempArrayFaceData[m100+0];
      }
      // corner (+x,-y,-z)
      if((nf & Hybrid::X_POS_EXISTS) == 0 && (nf & Hybrid::Y_NEG_EXISTS) == 0 && (nf & Hybrid::Z_NEG_EXISTS) == 0) {
	 const int i=0,j=0,k=0;
	 const int n = (b*block::SIZE+block::index(i,j,k))*vectorDim;
	 const int m012 = block::arrayIndex(i+0,j+1,k+2)*vectorDim;
         const int m021 = block::arrayIndex(i+0,j+2,k+1)*vectorDim;
         faceData[n+1] = tempArrayFaceData[m012+1];
         faceData[n+2] = tempArrayFaceData[m021+2];
      }
      // corner (+x,-y,+z)
      if((nf & Hybrid::X_POS_EXISTS) == 0 && (nf & Hybrid::Y_NEG_EXISTS) == 0 && (nf & Hybrid::Z_POS_EXISTS) == 0) {
	 const int i=0,j=0,k=0;
	 const int n = (b*block::SIZE+block::index(i,j,k))*vectorDim;
	 const int m010 = block::arrayIndex(i+0,j+1,k+0)*vectorDim;
         faceData[n+1] = tempArrayFaceData[m010+1];
      }
      // corner (+x,+y,-z)
      if((nf & Hybrid::X_POS_EXISTS) == 0 && (nf & Hybrid::Y_POS_EXISTS) == 0 && (nf & Hybrid::Z_NEG_EXISTS) == 0) {
	 const int i=0,j=0,k=0;
	 const int n = (b*block::SIZE+block::index(i,j,k))*vectorDim;
	 const int m001 = block::arrayIndex(i+0,j+0,k+1)*vectorDim;
         faceData[n+2] = tempArrayFaceData[m001+2];
      }
      // corner (+x,+y,+z)
      //if((nf & Hybrid::X_POS_EXISTS) == 0 && (nf & Hybrid::Y_POS_EXISTS) == 0 && (nf & Hybrid::Z_POS_EXISTS) == 0) {
      // all ghost outer faces
      //}
   }
   delete [] tempArrayFaceData; tempArrayFaceData = NULL;
}

void setIMF(Real* cellB,Simulation& sim,SimulationClasses& simClasses,const vector<pargrid::CellID>& exteriorBlocks)
{
   const std::vector<uint32_t>& neighbourFlags = simClasses.pargrid.getNeighbourFlags();
   for(pargrid::CellID eb=0; eb<exteriorBlocks.size(); ++eb) {
      const pargrid::CellID b = exteriorBlocks[eb];
      const uint32_t& nf = neighbourFlags[b];
      // front (+x) wall
      if((nf & Hybrid::X_POS_EXISTS) == 0) {
	 for(int k=0; k<block::WIDTH_Z; ++k) for(int j=0; j<block::WIDTH_Y; ++j) {
	    const int i = 0;
	    const int n = (b*block::SIZE+block::index(i,j,k))*3;
	    cellB[n+0] = Hybrid::IMFBx; // IMF Bx
	    cellB[n+1] = Hybrid::IMFBy; // IMF By
	    cellB[n+2] = Hybrid::IMFBz; // IMF Bz
	 }
      }
   }
}

void setIMFFace(Real* faceB,Simulation& sim,SimulationClasses& simClasses,const vector<pargrid::CellID>& exteriorBlocks)
{
   const std::vector<uint32_t>& neighbourFlags = simClasses.pargrid.getNeighbourFlags();
   for(pargrid::CellID eb=0; eb<exteriorBlocks.size(); ++eb) {
      const pargrid::CellID b = exteriorBlocks[eb];
      const uint32_t& nf = neighbourFlags[b];
      // front (+x) wall
      if((nf & Hybrid::X_POS_EXISTS) == 0) {
	 for(int k=0; k<block::WIDTH_Z; ++k) for(int j=0; j<block::WIDTH_Y; ++j) {
	    const int i = 0;
	    const int n = (b*block::SIZE+block::index(i,j,k))*3;
	    faceB[n+1] = Hybrid::IMFBy; // IMF By
	    faceB[n+2] = Hybrid::IMFBz; // IMF Bz
	 }
      }
   }
}

// interpolation from faces to cells
void face2Cell(Real* faceData,Real* cellData,Simulation& sim,SimulationClasses& simClasses,pargrid::CellID blockID)
{
   if(simClasses.pargrid.getNeighbourFlags(blockID) != pargrid::ALL_NEIGHBOURS_EXIST) return;
   const unsigned int size = (block::WIDTH_X+2)*(block::WIDTH_Y+2)*(block::WIDTH_Z+2);
   Real array[size*3];
   fetchData(faceData,array,simClasses,blockID,3);
   for(int k=0; k<block::WIDTH_Z; ++k) for(int j=0; j<block::WIDTH_Y; ++j) for(int i=0; i<block::WIDTH_X; ++i) {
      const int n = (blockID*block::SIZE+block::index(i,j,k))*3;
      cellData[n+0] = 0.5*(array[(block::arrayIndex(i+1,j+1,k+1))*3+0] + array[(block::arrayIndex(i+0,j+1,k+1))*3+0]);
      cellData[n+1] = 0.5*(array[(block::arrayIndex(i+1,j+1,k+1))*3+1] + array[(block::arrayIndex(i+1,j+0,k+1))*3+1]);
      cellData[n+2] = 0.5*(array[(block::arrayIndex(i+1,j+1,k+1))*3+2] + array[(block::arrayIndex(i+1,j+1,k+0))*3+2]);
   }
}

// inteprolation from cells to nodes
void cell2Node(Real* cellData,Real* nodeData,Simulation& sim,SimulationClasses& simClasses,pargrid::CellID blockID,const int vectorDim)
{
   int di=0;
   int dj=0;
   int dk=0;
   if(simClasses.pargrid.getNeighbourFlags(blockID) != pargrid::ALL_NEIGHBOURS_EXIST) {
      if((simClasses.pargrid.getNeighbourFlags()[blockID] & Hybrid::X_POS_EXISTS) == 0) return;
      if((simClasses.pargrid.getNeighbourFlags()[blockID] & Hybrid::Y_POS_EXISTS) == 0) return;
      if((simClasses.pargrid.getNeighbourFlags()[blockID] & Hybrid::Z_POS_EXISTS) == 0) return;
      if((simClasses.pargrid.getNeighbourFlags()[blockID] & Hybrid::X_NEG_EXISTS) == 0 && block::WIDTH_X > 1) di = block::WIDTH_X-1;
      if((simClasses.pargrid.getNeighbourFlags()[blockID] & Hybrid::Y_NEG_EXISTS) == 0 && block::WIDTH_Y > 1) dj = block::WIDTH_Y-1;
      if((simClasses.pargrid.getNeighbourFlags()[blockID] & Hybrid::Z_NEG_EXISTS) == 0 && block::WIDTH_Z > 1) dk = block::WIDTH_Z-1;
   }
   
   const unsigned int size = (block::WIDTH_X+2)*(block::WIDTH_Y+2)*(block::WIDTH_Z+2);
   Real array[size*vectorDim];
   fetchData(cellData,array,simClasses,blockID,vectorDim);

   for(int k=0+dk; k<block::WIDTH_Z; ++k) for(int j=0+dj; j<block::WIDTH_Y; ++j) for(int i=0+di; i<block::WIDTH_X; ++i) {
      const int n = (blockID*block::SIZE+block::index(i,j,k))*vectorDim;
      for(int l=0;l<vectorDim;++l) {
	 nodeData[n+l] = 0.125*(array[(block::arrayIndex(i+1,j+1,k+1))*vectorDim+l] +
				array[(block::arrayIndex(i+1,j+1,k+2))*vectorDim+l] +
				array[(block::arrayIndex(i+1,j+2,k+1))*vectorDim+l] +
				array[(block::arrayIndex(i+2,j+1,k+1))*vectorDim+l] +
				array[(block::arrayIndex(i+1,j+2,k+2))*vectorDim+l] +
				array[(block::arrayIndex(i+2,j+2,k+1))*vectorDim+l] +
				array[(block::arrayIndex(i+2,j+1,k+2))*vectorDim+l] +
				array[(block::arrayIndex(i+2,j+2,k+2))*vectorDim+l]);
      }
   }
}

// interpolation from nodes to cells
void node2Cell(Real* nodeData,Real* cellData,Simulation& sim,SimulationClasses& simClasses,pargrid::CellID blockID)
{
   if(simClasses.pargrid.getNeighbourFlags(blockID) != pargrid::ALL_NEIGHBOURS_EXIST) return;
   const unsigned int size = (block::WIDTH_X+2)*(block::WIDTH_Y+2)*(block::WIDTH_Z+2);
   Real array[size*3];
   fetchData(nodeData,array,simClasses,blockID,3);
   for(int k=0; k<block::WIDTH_Z; ++k) for(int j=0; j<block::WIDTH_Y; ++j) for(int i=0; i<block::WIDTH_X; ++i) {
      const int n = (blockID*block::SIZE+block::index(i,j,k))*3;
      for(int l=0;l<3;++l) {
	 cellData[n+l] = 0.125*(array[(block::arrayIndex(i+1,j+1,k+1))*3+l] + 
				array[(block::arrayIndex(i+0,j+1,k+1))*3+l] +
				array[(block::arrayIndex(i+0,j+0,k+1))*3+l] +
				array[(block::arrayIndex(i+1,j+0,k+1))*3+l] +
				array[(block::arrayIndex(i+1,j+1,k+0))*3+l] +
				array[(block::arrayIndex(i+1,j+0,k+0))*3+l] +
				array[(block::arrayIndex(i+0,j+1,k+0))*3+l] +
				array[(block::arrayIndex(i+0,j+0,k+0))*3+l]);
      }
   }
}

// node data average from all neighbors
void nodeAvg(Real* nodeDataOld,Real* nodeData,Simulation& sim,SimulationClasses& simClasses,pargrid::CellID blockID,const int vectorDim)
{
   int di=0;
   int dj=0;
   int dk=0;
   if(simClasses.pargrid.getNeighbourFlags(blockID) != pargrid::ALL_NEIGHBOURS_EXIST) {
      if((simClasses.pargrid.getNeighbourFlags()[blockID] & Hybrid::X_POS_EXISTS) == 0) return;
      if((simClasses.pargrid.getNeighbourFlags()[blockID] & Hybrid::Y_POS_EXISTS) == 0) return;
      if((simClasses.pargrid.getNeighbourFlags()[blockID] & Hybrid::Z_POS_EXISTS) == 0) return;
      if((simClasses.pargrid.getNeighbourFlags()[blockID] & Hybrid::X_NEG_EXISTS) == 0 && block::WIDTH_X > 1) di = block::WIDTH_X-1;
      if((simClasses.pargrid.getNeighbourFlags()[blockID] & Hybrid::Y_NEG_EXISTS) == 0 && block::WIDTH_Y > 1) dj = block::WIDTH_Y-1;
      if((simClasses.pargrid.getNeighbourFlags()[blockID] & Hybrid::Z_NEG_EXISTS) == 0 && block::WIDTH_Z > 1) dk = block::WIDTH_Z-1;
   }
   const unsigned int arraySize = (block::WIDTH_X+2)*(block::WIDTH_Y+2)*(block::WIDTH_Z+2)*vectorDim;
   const Real initVal = numeric_limits<Real>::max();
   Real array[arraySize] = {initVal};
   fetchData(nodeDataOld,array,simClasses,blockID,vectorDim);
   // coefficients
   const Real C1 = Hybrid::EfilterNodeGaussCoeffs[0]; // node itself
   const Real C2 = Hybrid::EfilterNodeGaussCoeffs[1]; // direct neighbors (distance = dx)
   const Real C3 = Hybrid::EfilterNodeGaussCoeffs[2]; // diagonal neighbors (distance = sqrt(2)*dx)
   const Real C4 = Hybrid::EfilterNodeGaussCoeffs[3]; // diagonal neighbors (distance = sqrt(3)*dx)
   // loop through all nodes in the simulation domain
   for(int k=0+dk; k<block::WIDTH_Z; ++k) for(int j=0+dj; j<block::WIDTH_Y; ++j) for(int i=0+di; i<block::WIDTH_X; ++i) {
      const int n = (blockID*block::SIZE+block::index(i,j,k))*vectorDim;
      for(int l=0;l<vectorDim;++l) {
         // electric field component l (x=0,y=1,z=2) at 27 nodes
         Real nodeEl111 = array[(block::arrayIndex(i+1,j+1,k+1))*vectorDim+l];
         Real nodeEl110 = array[(block::arrayIndex(i+1,j+1,k+0))*vectorDim+l];
         Real nodeEl011 = array[(block::arrayIndex(i+0,j+1,k+1))*vectorDim+l];
         Real nodeEl121 = array[(block::arrayIndex(i+1,j+2,k+1))*vectorDim+l];
         Real nodeEl211 = array[(block::arrayIndex(i+2,j+1,k+1))*vectorDim+l];
         Real nodeEl101 = array[(block::arrayIndex(i+1,j+0,k+1))*vectorDim+l];
         Real nodeEl112 = array[(block::arrayIndex(i+1,j+1,k+2))*vectorDim+l];
         Real nodeEl010 = array[(block::arrayIndex(i+0,j+1,k+0))*vectorDim+l];
         Real nodeEl100 = array[(block::arrayIndex(i+1,j+0,k+0))*vectorDim+l];
         Real nodeEl120 = array[(block::arrayIndex(i+1,j+2,k+0))*vectorDim+l];
         Real nodeEl210 = array[(block::arrayIndex(i+2,j+1,k+0))*vectorDim+l];
         Real nodeEl001 = array[(block::arrayIndex(i+0,j+0,k+1))*vectorDim+l];
         Real nodeEl021 = array[(block::arrayIndex(i+0,j+2,k+1))*vectorDim+l];
         Real nodeEl221 = array[(block::arrayIndex(i+2,j+2,k+1))*vectorDim+l];
         Real nodeEl201 = array[(block::arrayIndex(i+2,j+0,k+1))*vectorDim+l];
         Real nodeEl012 = array[(block::arrayIndex(i+0,j+1,k+2))*vectorDim+l];
         Real nodeEl122 = array[(block::arrayIndex(i+1,j+2,k+2))*vectorDim+l];
         Real nodeEl212 = array[(block::arrayIndex(i+2,j+1,k+2))*vectorDim+l];
         Real nodeEl102 = array[(block::arrayIndex(i+1,j+0,k+2))*vectorDim+l];
         Real nodeEl002 = array[(block::arrayIndex(i+0,j+0,k+2))*vectorDim+l];
         Real nodeEl022 = array[(block::arrayIndex(i+0,j+2,k+2))*vectorDim+l];
         Real nodeEl202 = array[(block::arrayIndex(i+2,j+0,k+2))*vectorDim+l];
         Real nodeEl222 = array[(block::arrayIndex(i+2,j+2,k+2))*vectorDim+l];
         Real nodeEl000 = array[(block::arrayIndex(i+0,j+0,k+0))*vectorDim+l];
         Real nodeEl020 = array[(block::arrayIndex(i+0,j+2,k+0))*vectorDim+l];
         Real nodeEl200 = array[(block::arrayIndex(i+2,j+0,k+0))*vectorDim+l];
         Real nodeEl220 = array[(block::arrayIndex(i+2,j+2,k+0))*vectorDim+l];
         
         // weight coefficients for each node (if the electric field value is exactly
         // zero, the node is excluedd as it should be a boudary node)
         bool zeroFound = false;
         Real C111 = C1; if(nodeEl111 == initVal) { nodeEl111 = 0.0; C111 = 0.0; zeroFound = true; }
         Real C110 = C2; if(nodeEl110 == initVal) { nodeEl110 = 0.0; C110 = 0.0; zeroFound = true; }
         Real C011 = C2; if(nodeEl011 == initVal) { nodeEl011 = 0.0; C011 = 0.0; zeroFound = true; }
         Real C121 = C2; if(nodeEl121 == initVal) { nodeEl121 = 0.0; C121 = 0.0; zeroFound = true; }
         Real C211 = C2; if(nodeEl211 == initVal) { nodeEl211 = 0.0; C211 = 0.0; zeroFound = true; }
         Real C101 = C2; if(nodeEl101 == initVal) { nodeEl101 = 0.0; C101 = 0.0; zeroFound = true; }
         Real C112 = C2; if(nodeEl112 == initVal) { nodeEl112 = 0.0; C112 = 0.0; zeroFound = true; }
         Real C010 = C3; if(nodeEl010 == initVal) { nodeEl010 = 0.0; C010 = 0.0; zeroFound = true; }
         Real C100 = C3; if(nodeEl100 == initVal) { nodeEl100 = 0.0; C100 = 0.0; zeroFound = true; }
         Real C120 = C3; if(nodeEl120 == initVal) { nodeEl120 = 0.0; C120 = 0.0; zeroFound = true; }
         Real C210 = C3; if(nodeEl210 == initVal) { nodeEl210 = 0.0; C210 = 0.0; zeroFound = true; }
         Real C001 = C3; if(nodeEl001 == initVal) { nodeEl001 = 0.0; C001 = 0.0; zeroFound = true; }
         Real C021 = C3; if(nodeEl021 == initVal) { nodeEl021 = 0.0; C021 = 0.0; zeroFound = true; }
         Real C221 = C3; if(nodeEl221 == initVal) { nodeEl221 = 0.0; C221 = 0.0; zeroFound = true; }
         Real C201 = C3; if(nodeEl201 == initVal) { nodeEl201 = 0.0; C201 = 0.0; zeroFound = true; }
         Real C012 = C3; if(nodeEl012 == initVal) { nodeEl012 = 0.0; C012 = 0.0; zeroFound = true; }
         Real C122 = C3; if(nodeEl122 == initVal) { nodeEl122 = 0.0; C122 = 0.0; zeroFound = true; }
         Real C212 = C3; if(nodeEl212 == initVal) { nodeEl212 = 0.0; C212 = 0.0; zeroFound = true; }
         Real C102 = C3; if(nodeEl102 == initVal) { nodeEl102 = 0.0; C102 = 0.0; zeroFound = true; }
         Real C002 = C4; if(nodeEl002 == initVal) { nodeEl002 = 0.0; C002 = 0.0; zeroFound = true; }
         Real C022 = C4; if(nodeEl022 == initVal) { nodeEl022 = 0.0; C022 = 0.0; zeroFound = true; }
         Real C202 = C4; if(nodeEl202 == initVal) { nodeEl202 = 0.0; C202 = 0.0; zeroFound = true; }
         Real C222 = C4; if(nodeEl222 == initVal) { nodeEl222 = 0.0; C222 = 0.0; zeroFound = true; }
         Real C000 = C4; if(nodeEl000 == initVal) { nodeEl000 = 0.0; C000 = 0.0; zeroFound = true; }
         Real C020 = C4; if(nodeEl020 == initVal) { nodeEl020 = 0.0; C020 = 0.0; zeroFound = true; }
         Real C200 = C4; if(nodeEl200 == initVal) { nodeEl200 = 0.0; C200 = 0.0; zeroFound = true; }
         Real C220 = C4; if(nodeEl220 == initVal) { nodeEl220 = 0.0; C220 = 0.0; zeroFound = true; }
         // renormalize if boundary node found
         if(zeroFound == true) {
            const Real Csum = C111+C110+C011+C121+C211+C101+C112+C010+C100+C120+C210+C001+C021+C221+C201+C012+C122+C212+C102+C002+C022+C202+C222+C000+C020+C200+C220;
            if(Csum > 0.0) {
               C111 /= Csum;
               C110 /= Csum;
               C011 /= Csum;
               C121 /= Csum;
               C211 /= Csum;
               C101 /= Csum;
               C112 /= Csum;
               C010 /= Csum;
               C100 /= Csum;
               C120 /= Csum;
               C210 /= Csum;
               C001 /= Csum;
               C021 /= Csum;
               C221 /= Csum;
               C201 /= Csum;
               C012 /= Csum;
               C122 /= Csum;
               C212 /= Csum;
               C102 /= Csum;
               C002 /= Csum;
               C022 /= Csum;
               C202 /= Csum;
               C222 /= Csum;
               C000 /= Csum;
               C020 /= Csum;
               C200 /= Csum;
               C220 /= Csum;
            }
         }
	 // calculate average
         nodeData[n+l] =
           C111*nodeEl111 +
           C110*nodeEl110 +
           C011*nodeEl011 +
           C121*nodeEl121 +
           C211*nodeEl211 +
           C101*nodeEl101 +
           C112*nodeEl112 +
           C010*nodeEl010 +
           C100*nodeEl100 +
           C120*nodeEl120 +
           C210*nodeEl210 +
           C001*nodeEl001 +
           C021*nodeEl021 +
           C221*nodeEl221 +
           C201*nodeEl201 +
           C012*nodeEl012 +
           C122*nodeEl122 +
           C212*nodeEl212 +
           C102*nodeEl102 +
           C002*nodeEl002 +
           C022*nodeEl022 +
           C202*nodeEl202 +
           C222*nodeEl222 +
           C000*nodeEl000 +
           C020*nodeEl020 +
           C200*nodeEl200 +
           C220*nodeEl220;
      }
   }
}

// upwind nodeB using cellData and nodeUe
void upwindNodeB(Real* cellB,Real* nodeUe,Real* nodeB,Simulation& sim,SimulationClasses& simClasses,pargrid::CellID blockID) {
   int di=0;
   int dj=0;
   int dk=0;
   if(simClasses.pargrid.getNeighbourFlags(blockID) != pargrid::ALL_NEIGHBOURS_EXIST) {
      if((simClasses.pargrid.getNeighbourFlags()[blockID] & Hybrid::X_POS_EXISTS) == 0) return;
      if((simClasses.pargrid.getNeighbourFlags()[blockID] & Hybrid::Y_POS_EXISTS) == 0) return;
      if((simClasses.pargrid.getNeighbourFlags()[blockID] & Hybrid::Z_POS_EXISTS) == 0) return;
      if((simClasses.pargrid.getNeighbourFlags()[blockID] & Hybrid::X_NEG_EXISTS) == 0 && block::WIDTH_X > 1) di = block::WIDTH_X-1;
      if((simClasses.pargrid.getNeighbourFlags()[blockID] & Hybrid::Y_NEG_EXISTS) == 0 && block::WIDTH_Y > 1) dj = block::WIDTH_Y-1;
      if((simClasses.pargrid.getNeighbourFlags()[blockID] & Hybrid::Z_NEG_EXISTS) == 0 && block::WIDTH_Z > 1) dk = block::WIDTH_Z-1;
   }
   const unsigned int tempArraySize = (block::WIDTH_X+2)*(block::WIDTH_Y+2)*(block::WIDTH_Z+2);
   Real tempArray[tempArraySize*3];
   fetchData(cellB,tempArray,simClasses,blockID,3);
   for(int k=0+dk; k<block::WIDTH_Z; ++k) for(int j=0+dj; j<block::WIDTH_Y; ++j) for(int i=0+di; i<block::WIDTH_X; ++i) {
      const int n = (blockID*block::SIZE+block::index(i,j,k))*3;
     
      // upwind displacement vector = (-0.5,0,0) if Ue = 0
      Real xUpwind = -0.5;
      Real yUpwind = 0;
      Real zUpwind = 0;
      
      // upwind displacement vector = -0.5*nodeUe/|nodeUe| (dimensionaless, dx=1)
      const Real Ue = sqrt(sqr(nodeUe[n+0]) + sqr(nodeUe[n+1]) + sqr(nodeUe[n+2]));
      if(Ue > 0) {
	 const Real a = -0.5/Ue;
	 xUpwind = nodeUe[n+0]*a;
	 yUpwind = nodeUe[n+1]*a;
	 zUpwind = nodeUe[n+2]*a;
      }
      
      // cell centroid local dimensionless coordinates
      const Real xCell111 = -0.5; const Real yCell111 = -0.5; const Real zCell111 = -0.5;
      const Real xCell112 = -0.5; const Real yCell112 = -0.5; const Real zCell112 = +0.5;
      const Real xCell121 = -0.5; const Real yCell121 = +0.5; const Real zCell121 = -0.5;
      const Real xCell211 = +0.5; const Real yCell211 = -0.5; const Real zCell211 = -0.5;
      const Real xCell122 = -0.5; const Real yCell122 = +0.5; const Real zCell122 = +0.5;
      const Real xCell221 = +0.5; const Real yCell221 = +0.5; const Real zCell221 = -0.5;
      const Real xCell212 = +0.5; const Real yCell212 = -0.5; const Real zCell212 = +0.5;
      const Real xCell222 = +0.5; const Real yCell222 = +0.5; const Real zCell222 = +0.5;
      
      // weighting factors for the eight cells around the node
      const Real w111 = 1/sqrt(sqr(xCell111-xUpwind) + sqr(yCell111-yUpwind) + sqr(zCell111-zUpwind));
      const Real w112 = 1/sqrt(sqr(xCell112-xUpwind) + sqr(yCell112-yUpwind) + sqr(zCell112-zUpwind));
      const Real w121 = 1/sqrt(sqr(xCell121-xUpwind) + sqr(yCell121-yUpwind) + sqr(zCell121-zUpwind));
      const Real w211 = 1/sqrt(sqr(xCell211-xUpwind) + sqr(yCell211-yUpwind) + sqr(zCell211-zUpwind));
      const Real w122 = 1/sqrt(sqr(xCell122-xUpwind) + sqr(yCell122-yUpwind) + sqr(zCell122-zUpwind));
      const Real w221 = 1/sqrt(sqr(xCell221-xUpwind) + sqr(yCell221-yUpwind) + sqr(zCell221-zUpwind));
      const Real w212 = 1/sqrt(sqr(xCell212-xUpwind) + sqr(yCell212-yUpwind) + sqr(zCell212-zUpwind));
      const Real w222 = 1/sqrt(sqr(xCell222-xUpwind) + sqr(yCell222-yUpwind) + sqr(zCell222-zUpwind));
      const Real wsum = w111+w112+w121+w211+w122+w221+w212+w222;

      if(wsum > 0) {
	 for(int l=0;l<3;++l) {
	    nodeB[n+l] = (w111*tempArray[(block::arrayIndex(i+1,j+1,k+1))*3+l] +
			  w112*tempArray[(block::arrayIndex(i+1,j+1,k+2))*3+l] +
			  w121*tempArray[(block::arrayIndex(i+1,j+2,k+1))*3+l] +
			  w211*tempArray[(block::arrayIndex(i+2,j+1,k+1))*3+l] +
			  w122*tempArray[(block::arrayIndex(i+1,j+2,k+2))*3+l] +
			  w221*tempArray[(block::arrayIndex(i+2,j+2,k+1))*3+l] +
			  w212*tempArray[(block::arrayIndex(i+2,j+1,k+2))*3+l] +
			  w222*tempArray[(block::arrayIndex(i+2,j+2,k+2))*3+l])/wsum;
	 }
      }
      else {
	 nodeB[n+0] = nodeB[n+1] = nodeB[n+2] = 0.0;
      }
   }
}

// calculate cellUe
void calcCellUe(Real* cellJ,Real* cellJi,Real* cellRhoQi,Real* cellUe,bool* innerFlag,Real* counterCellMaxUe,Simulation& sim,SimulationClasses& simClasses,pargrid::CellID blockID)
{
   if(simClasses.pargrid.getNeighbourFlags(blockID) != pargrid::ALL_NEIGHBOURS_EXIST) return;
   for(int k=0; k<block::WIDTH_Z; ++k) for(int j=0; j<block::WIDTH_Y; ++j) for(int i=0; i<block::WIDTH_X; ++i) {
      const int n = (blockID*block::SIZE+block::index(i,j,k));
      const int n3 = n*3;
      // inner boundary condition for Ue
      if(innerFlag[n] == true) {
	 cellUe[n3+0] = cellUe[n3+1] = cellUe[n3+2] = 0.0;
	 continue;
      }
      // calc Ue = (J - Ji)/rhoqi
      for(int l=0;l<3;++l) {
	 if(fabs(cellRhoQi[n]) > 0) {
            if(Hybrid::useHallElectricField == true) {
               cellUe[n3+l] = -(cellJ[n3+l] - cellJi[n3+l])/cellRhoQi[n];
            }
            else {
               cellUe[n3+l] = cellJi[n3+l]/cellRhoQi[n];
            }
	 }
	 else {
	    cellUe[n3+l] = 0.0;
	 }
      }
      
      // check max Ue
      const Real Ue2 = sqr(cellUe[n3+0]) + sqr(cellUe[n3+1]) + sqr(cellUe[n3+2]);
      if(Ue2 > Hybrid::maxUe2) {
	 const Real norm = sqrt(Hybrid::maxUe2/Ue2);
	 cellUe[n3+0] *= norm;
	 cellUe[n3+1] *= norm;
	 cellUe[n3+2] *= norm;
         counterCellMaxUe[n]++;
      }
   }
}
   

// E = -Ue x B + eta*J
void calcNodeE(Real* nodeUe,Real* nodeB,
#ifdef USE_RESISTIVITY
Real* nodeEta,
#endif
Real* nodeJ,Real* nodeE,
#ifdef USE_ECUT
Real* counterNodeEcut,
#endif
bool* innerFlag,Simulation& sim,SimulationClasses& simClasses,pargrid::CellID blockID)
{
   int di=0;
   int dj=0;
   int dk=0;
   if(simClasses.pargrid.getNeighbourFlags(blockID) != pargrid::ALL_NEIGHBOURS_EXIST) {
      if((simClasses.pargrid.getNeighbourFlags()[blockID] & Hybrid::X_POS_EXISTS) == 0) return;
      if((simClasses.pargrid.getNeighbourFlags()[blockID] & Hybrid::Y_POS_EXISTS) == 0) return;
      if((simClasses.pargrid.getNeighbourFlags()[blockID] & Hybrid::Z_POS_EXISTS) == 0) return;
      if((simClasses.pargrid.getNeighbourFlags()[blockID] & Hybrid::X_NEG_EXISTS) == 0 && block::WIDTH_X > 1) di = block::WIDTH_X-1;
      if((simClasses.pargrid.getNeighbourFlags()[blockID] & Hybrid::Y_NEG_EXISTS) == 0 && block::WIDTH_Y > 1) dj = block::WIDTH_Y-1;
      if((simClasses.pargrid.getNeighbourFlags()[blockID] & Hybrid::Z_NEG_EXISTS) == 0 && block::WIDTH_Z > 1) dk = block::WIDTH_Z-1;
   }

#ifdef USE_B_CONSTANT
   const Real* crd = getBlockCoordinateArray(sim,simClasses);
   const size_t b3 = 3*blockID;
#endif
   for(int k=0+dk; k<block::WIDTH_Z; ++k) for(int j=0+dj; j<block::WIDTH_Y; ++j) for(int i=0+di; i<block::WIDTH_X; ++i) {
      const int n = (blockID*block::SIZE+block::index(i,j,k));
      const int n3 = n*3;
      Real Btot[3];
      Btot[0] = nodeB[n3+0];
      Btot[1] = nodeB[n3+1];
      Btot[2] = nodeB[n3+2];
#ifdef USE_B_CONSTANT
      const Real xNode = crd[b3+0] + (i+1.0)*Hybrid::dx;
      const Real yNode = crd[b3+1] + (j+1.0)*Hybrid::dx;
      const Real zNode = crd[b3+2] + (k+1.0)*Hybrid::dx;
      addConstantB(xNode,yNode,zNode,Btot);
#endif
      nodeE[n3+0] = -(nodeUe[n3+1]*Btot[2] - nodeUe[n3+2]*Btot[1]);
      nodeE[n3+1] = -(nodeUe[n3+2]*Btot[0] - nodeUe[n3+0]*Btot[2]);
      nodeE[n3+2] = -(nodeUe[n3+0]*Btot[1] - nodeUe[n3+1]*Btot[0]);
#ifdef USE_RESISTIVITY
      nodeE[n3+0] += nodeEta[n]*nodeJ[n3+0];
      nodeE[n3+1] += nodeEta[n]*nodeJ[n3+1];
      nodeE[n3+2] += nodeEta[n]*nodeJ[n3+2];
#endif
#ifdef USE_ECUT
      if(Hybrid::Ecut2 > 0) {
         const Real E2 = sqr(nodeE[n3+0]) + sqr(nodeE[n3+1]) + sqr(nodeE[n3+2]);
         if (E2 > Hybrid::Ecut2) {
            Real scaling = sqrt(Hybrid::Ecut2/E2);
            nodeE[n3+0] *= scaling;
            nodeE[n3+1] *= scaling;
            nodeE[n3+2] *= scaling;
            counterNodeEcut[n]++;
         }
      }
#endif
   }
}

// nodeJ = nabla x faceB/mu0
void calcNodeJ(Real* faceB,Real* nodeB,Real* nodeRhoQi,Real* nodeJ,
#ifdef USE_MAXVW
Real* counterNodeMaxVw,
#endif
Simulation& sim,SimulationClasses& simClasses,pargrid::CellID blockID)
{
   int di=0;
   int dj=0;
   int dk=0;
   if(simClasses.pargrid.getNeighbourFlags(blockID) != pargrid::ALL_NEIGHBOURS_EXIST) {
      if((simClasses.pargrid.getNeighbourFlags()[blockID] & Hybrid::X_POS_EXISTS) == 0) return;
      if((simClasses.pargrid.getNeighbourFlags()[blockID] & Hybrid::Y_POS_EXISTS) == 0) return;
      if((simClasses.pargrid.getNeighbourFlags()[blockID] & Hybrid::Z_POS_EXISTS) == 0) return;
      if((simClasses.pargrid.getNeighbourFlags()[blockID] & Hybrid::X_NEG_EXISTS) == 0 && block::WIDTH_X > 1) di = block::WIDTH_X-1;
      if((simClasses.pargrid.getNeighbourFlags()[blockID] & Hybrid::Y_NEG_EXISTS) == 0 && block::WIDTH_Y > 1) dj = block::WIDTH_Y-1;
      if((simClasses.pargrid.getNeighbourFlags()[blockID] & Hybrid::Z_NEG_EXISTS) == 0 && block::WIDTH_Z > 1) dk = block::WIDTH_Z-1;
   }
   const unsigned int size = (block::WIDTH_X+2)*(block::WIDTH_Y+2)*(block::WIDTH_Z+2);
   Real Bf[size*3];
   fetchData(faceB,Bf,simClasses,blockID,3);
   for(int k=0+dk; k<block::WIDTH_Z; ++k) for(int j=0+dj; j<block::WIDTH_Y; ++j) for(int i=0+di; i<block::WIDTH_X; ++i) {
      const int n = (blockID*block::SIZE+block::index(i,j,k));
      const int n3 = n*3;
      const Real edgeJx1 =
        - Bf[(block::arrayIndex(i+1,j+1,k+1))*3+2]
        + Bf[(block::arrayIndex(i+1,j+1,k+1))*3+1]
        + Bf[(block::arrayIndex(i+1,j+2,k+1))*3+2]
        - Bf[(block::arrayIndex(i+1,j+1,k+2))*3+1];
      const Real edgeJx2 =
        - Bf[(block::arrayIndex(i+2,j+1,k+1))*3+2]
        + Bf[(block::arrayIndex(i+2,j+1,k+1))*3+1]
        + Bf[(block::arrayIndex(i+2,j+2,k+1))*3+2]
        - Bf[(block::arrayIndex(i+2,j+1,k+2))*3+1];
      const Real edgeJy1 =
        + Bf[(block::arrayIndex(i+1,j+1,k+1))*3+2]
        + Bf[(block::arrayIndex(i+1,j+1,k+2))*3+0]
        - Bf[(block::arrayIndex(i+2,j+1,k+1))*3+2]
        - Bf[(block::arrayIndex(i+1,j+1,k+1))*3+0]; 
      const Real edgeJy2 =
        + Bf[(block::arrayIndex(i+1,j+2,k+1))*3+2]
        + Bf[(block::arrayIndex(i+1,j+2,k+2))*3+0]
        - Bf[(block::arrayIndex(i+2,j+2,k+1))*3+2]
        - Bf[(block::arrayIndex(i+1,j+2,k+1))*3+0];
      const Real edgeJz1 =
        - Bf[(block::arrayIndex(i+1,j+1,k+1))*3+1]
        + Bf[(block::arrayIndex(i+1,j+1,k+1))*3+0]
        + Bf[(block::arrayIndex(i+2,j+1,k+1))*3+1]
        - Bf[(block::arrayIndex(i+1,j+2,k+1))*3+0]; 
      const Real edgeJz2 =
        - Bf[(block::arrayIndex(i+1,j+1,k+2))*3+1]
        + Bf[(block::arrayIndex(i+1,j+1,k+2))*3+0]
        + Bf[(block::arrayIndex(i+2,j+1,k+2))*3+1]
        - Bf[(block::arrayIndex(i+1,j+2,k+2))*3+0];
      Real d = 1.0;
#ifdef USE_MAXVW
      Real vw = 0.0; // fastest whistler signal p. 28 Alho (2016)
      const Real ne = nodeRhoQi[n]/constants::CHARGE_ELEMENTARY;
      const Real Btot = sqrt( sqr(nodeB[n3+0]) + sqr(nodeB[n3+1]) + sqr(nodeB[n3+2]) );
      if(ne > 0.0 && Hybrid::dx > 0.0) {
	vw = 2.0*Btot*M_PI/( constants::PERMEABILITY*ne*constants::CHARGE_ELEMENTARY*Hybrid::dx );
      }
      if (vw > Hybrid::maxVw && Hybrid::maxVw > 0.0) {
	d = vw/Hybrid::maxVw;
	counterNodeMaxVw[n]++;
      }
#endif
      const Real a = 0.5/(Hybrid::dx*constants::PERMEABILITY*d);
      nodeJ[n3+0] = (edgeJx1 + edgeJx2)*a;
      nodeJ[n3+1] = (edgeJy1 + edgeJy2)*a;
      nodeJ[n3+2] = (edgeJz1 + edgeJz2)*a;
   }
}

void calcNodeUe(Real* nodeRhoQi,Real* nodeJi,Real* nodeJ,Real* nodeUe,bool* innerFlag,Real* counterCellMaxUe,Simulation& sim,SimulationClasses& simClasses,pargrid::CellID blockID)
{
   int di=0;
   int dj=0;
   int dk=0;
   if(simClasses.pargrid.getNeighbourFlags(blockID) != pargrid::ALL_NEIGHBOURS_EXIST) {
      if((simClasses.pargrid.getNeighbourFlags()[blockID] & Hybrid::X_POS_EXISTS) == 0) return;
      if((simClasses.pargrid.getNeighbourFlags()[blockID] & Hybrid::Y_POS_EXISTS) == 0) return;
      if((simClasses.pargrid.getNeighbourFlags()[blockID] & Hybrid::Z_POS_EXISTS) == 0) return;
      if((simClasses.pargrid.getNeighbourFlags()[blockID] & Hybrid::X_NEG_EXISTS) == 0 && block::WIDTH_X > 1) di = block::WIDTH_X-1;
      if((simClasses.pargrid.getNeighbourFlags()[blockID] & Hybrid::Y_NEG_EXISTS) == 0 && block::WIDTH_Y > 1) dj = block::WIDTH_Y-1;
      if((simClasses.pargrid.getNeighbourFlags()[blockID] & Hybrid::Z_NEG_EXISTS) == 0 && block::WIDTH_Z > 1) dk = block::WIDTH_Z-1;
   }
   for(int k=0+dk; k<block::WIDTH_Z; ++k) for(int j=0+dj; j<block::WIDTH_Y; ++j) for(int i=0+di; i<block::WIDTH_X; ++i) {
      const int n = (blockID*block::SIZE+block::index(i,j,k));
      const int n3 = n*3;
      // inner boundary condition for Ue
      if(innerFlag[n] == true) {
	 nodeUe[n3+0] = nodeUe[n3+1] = nodeUe[n3+2] = 0.0;
	 continue;
      }
      // check min nodeRhoQi
      if(nodeRhoQi[n] < Hybrid::minRhoQi) { nodeRhoQi[n] = Hybrid::minRhoQi; }
      // calc Ue = (J - Ji)/rhoqi
      for(int l=0;l<3;++l) {
	 if(fabs(nodeRhoQi[n]) > 0) {
            if(Hybrid::useHallElectricField == true) {
               nodeUe[n3+l] = -(nodeJ[n3+l] - nodeJi[n3+l])/nodeRhoQi[n];
            }
            else {
               nodeUe[n3+l] = nodeJi[n3+l]/nodeRhoQi[n];
            }
	 }
	 else {
	    nodeUe[n3+l] = 0.0;
	 }
      }
      
      // check max Ue
      const Real Ue2 = sqr(nodeUe[n3+0]) + sqr(nodeUe[n3+1]) + sqr(nodeUe[n3+2]);
      if(Ue2 > Hybrid::maxUe2) {
	 const Real norm = sqrt(Hybrid::maxUe2/Ue2);
	 nodeUe[n3+0] *= norm;
	 nodeUe[n3+1] *= norm;
	 nodeUe[n3+2] *= norm;
         counterCellMaxUe[n]++; // using cell array here to avoid introducing a new node array
      }
   }
}

// faceData = curl(nodeData), doFaraday: true = Faraday's law, false = Ampere's law
void faceCurl(Real* nodeData,Real* faceData,bool doFaraday,Simulation& sim,SimulationClasses& simClasses,pargrid::CellID blockID)
{
   bool doXFace=true;
   bool doYFace=true;
   bool doZFace=true;
   int di=0;
   int dj=0;
   int dk=0;
   if(simClasses.pargrid.getNeighbourFlags(blockID) != pargrid::ALL_NEIGHBOURS_EXIST) {
      if((simClasses.pargrid.getNeighbourFlags()[blockID] & Hybrid::X_POS_EXISTS) == 0) return;
      if((simClasses.pargrid.getNeighbourFlags()[blockID] & Hybrid::Y_POS_EXISTS) == 0) return;
      if((simClasses.pargrid.getNeighbourFlags()[blockID] & Hybrid::Z_POS_EXISTS) == 0) return;
      if((simClasses.pargrid.getNeighbourFlags()[blockID] & Hybrid::X_NEG_EXISTS) == 0) { doYFace = doZFace = false; }
      if((simClasses.pargrid.getNeighbourFlags()[blockID] & Hybrid::Y_NEG_EXISTS) == 0) { doXFace = doZFace = false; }
      if((simClasses.pargrid.getNeighbourFlags()[blockID] & Hybrid::Z_NEG_EXISTS) == 0) { doXFace = doYFace = false; }
      if((simClasses.pargrid.getNeighbourFlags()[blockID] & Hybrid::X_NEG_EXISTS) == 0 && block::WIDTH_X > 1) di = block::WIDTH_X-1;
      if((simClasses.pargrid.getNeighbourFlags()[blockID] & Hybrid::Y_NEG_EXISTS) == 0 && block::WIDTH_Y > 1) dj = block::WIDTH_Y-1;
      if((simClasses.pargrid.getNeighbourFlags()[blockID] & Hybrid::Z_NEG_EXISTS) == 0 && block::WIDTH_Z > 1) dk = block::WIDTH_Z-1;
   }
   
#ifdef USE_XMIN_BOUNDARY
   bool* xMinFlag = simClasses.pargrid.getUserDataStatic<bool>(Hybrid::dataXminFlagID);
   if(xMinFlag == NULL) {cerr << "ERROR: obtained NULL xMinFlag array!" << endl; exit(1);}
#endif

   const unsigned int size = (block::WIDTH_X+2)*(block::WIDTH_Y+2)*(block::WIDTH_Z+2);
   Real array[size*3];
   fetchData(nodeData,array,simClasses,blockID,3);
   
   for(int k=0+dk; k<block::WIDTH_Z; ++k) for(int j=0+dj; j<block::WIDTH_Y; ++j) for(int i=0+di; i<block::WIDTH_X; ++i) {
      const int n = (blockID*block::SIZE+block::index(i,j,k));
      const int n3 = 3*n;
#ifdef USE_XMIN_BOUNDARY
      // no field propagation at x < xmin
      if(xMinFlag[n] == true && doFaraday == true) { continue; }
#endif

      Real node1x = array[(block::arrayIndex(i+0,j+0,k+1))*3+0];
      Real node1y = array[(block::arrayIndex(i+0,j+0,k+1))*3+1];
      //Real node1z = array[(block::arrayIndex(i+0,j+0,k+1))*3+2];
      
      Real node3x = array[(block::arrayIndex(i+0,j+1,k+0))*3+0];
      //Real node3y = array[(block::arrayIndex(i+0,j+1,k+0))*3+1];
      Real node3z = array[(block::arrayIndex(i+0,j+1,k+0))*3+2];
      
      Real node4x = array[(block::arrayIndex(i+0,j+1,k+1))*3+0];
      Real node4y = array[(block::arrayIndex(i+0,j+1,k+1))*3+1];
      Real node4z = array[(block::arrayIndex(i+0,j+1,k+1))*3+2];      

      Real node5x = array[(block::arrayIndex(i+1,j+0,k+1))*3+0];
      Real node5y = array[(block::arrayIndex(i+1,j+0,k+1))*3+1];
      Real node5z = array[(block::arrayIndex(i+1,j+0,k+1))*3+2];

      //Real node6x = array[(block::arrayIndex(i+1,j+0,k+0))*3+0];
      Real node6y = array[(block::arrayIndex(i+1,j+0,k+0))*3+1];
      Real node6z = array[(block::arrayIndex(i+1,j+0,k+0))*3+2];
      
      Real node7x = array[(block::arrayIndex(i+1,j+1,k+0))*3+0];
      Real node7y = array[(block::arrayIndex(i+1,j+1,k+0))*3+1];
      Real node7z = array[(block::arrayIndex(i+1,j+1,k+0))*3+2];
      
      Real node8x = array[(block::arrayIndex(i+1,j+1,k+1))*3+0];
      Real node8y = array[(block::arrayIndex(i+1,j+1,k+1))*3+1];
      Real node8z = array[(block::arrayIndex(i+1,j+1,k+1))*3+2];
      
      if(doXFace == true) {
	 Real curlX = 0.5*(+node5z+node6z-node6y-node7y-node7z-node8z+node8y+node5y)/Hybrid::dx;
	 if(doFaraday == true) {
	    faceData[n3+0] += sim.dt*curlX; // Faraday's law
	 }
	 else {
	    faceData[n3+0] = -curlX/constants::PERMEABILITY; // Ampere's law
	 }
      }      
      if(doYFace == true) {
	 Real curlY = 0.5*(-node4x-node8x+node8z+node7z+node7x+node3x-node3z-node4z)/Hybrid::dx;
	 if(doFaraday == true) {
	    faceData[n3+1] += sim.dt*curlY;
	 }
	 else {
	    faceData[n3+1] = -curlY/constants::PERMEABILITY;
	 }
      }
      if(doZFace == true) {
	 Real curlZ = 0.5*(-node1x-node5x-node5y-node8y+node8x+node4x+node4y+node1y)/Hybrid::dx;
	 if(doFaraday == true) {
	    faceData[n3+2] += sim.dt*curlZ;
	 }
	 else {
	    faceData[n3+2] = -curlZ/constants::PERMEABILITY;
	 }
      }
   }
}

// calculate the electron pressure gradient term of the electric field
void calcCellEp(Real* nodeRhoQi,Real* cellRhoQi,bool* innerFlagCellEp,Real* cellEp,Simulation& sim,SimulationClasses& simClasses,pargrid::CellID blockID)
{
   if(simClasses.pargrid.getNeighbourFlags(blockID) != pargrid::ALL_NEIGHBOURS_EXIST) return;
   const unsigned int size = (block::WIDTH_X+2)*(block::WIDTH_Y+2)*(block::WIDTH_Z+2);
   Real array[size*3];
   fetchData(nodeRhoQi,array,simClasses,blockID,1);
   for(int k=0; k<block::WIDTH_Z; ++k) for(int j=0; j<block::WIDTH_Y; ++j) for(int i=0; i<block::WIDTH_X; ++i) {
      const int n = (blockID*block::SIZE+block::index(i,j,k));
      const int n3 = n*3;
      cellEp[n3+0] = 0.0;
      cellEp[n3+1] = 0.0;
      cellEp[n3+2] = 0.0;
      if(innerFlagCellEp[n] == true) { continue; }
      // node2face interpolation: nodeRhoQi -> faceRhoQi
      const Real xPosFaceRhoQi =
        0.25*(array[(block::arrayIndex(i+1,j+1,k+1))] +
              array[(block::arrayIndex(i+1,j+1,k+0))] +
              array[(block::arrayIndex(i+1,j+0,k+0))] +
              array[(block::arrayIndex(i+1,j+0,k+1))]);
      const Real xNegFaceRhoQi =
        0.25*(array[(block::arrayIndex(i+0,j+1,k+1))] +
              array[(block::arrayIndex(i+0,j+1,k+0))] +
              array[(block::arrayIndex(i+0,j+0,k+0))] +
              array[(block::arrayIndex(i+0,j+0,k+1))]);
      const Real yPosFaceRhoQi =
        0.25*(array[(block::arrayIndex(i+1,j+1,k+1))] +
              array[(block::arrayIndex(i+1,j+1,k+0))] +
              array[(block::arrayIndex(i+0,j+1,k+0))] +
              array[(block::arrayIndex(i+0,j+1,k+1))]);
      const Real yNegFaceRhoQi =
        0.25*(array[(block::arrayIndex(i+1,j+0,k+1))] +
              array[(block::arrayIndex(i+1,j+0,k+0))] +
              array[(block::arrayIndex(i+0,j+0,k+0))] +
              array[(block::arrayIndex(i+0,j+0,k+1))]);
      const Real zPosFaceRhoQi =
        0.25*(array[(block::arrayIndex(i+1,j+1,k+1))] +
              array[(block::arrayIndex(i+1,j+0,k+1))] +
              array[(block::arrayIndex(i+0,j+0,k+1))] +
              array[(block::arrayIndex(i+0,j+1,k+1))]);
      const Real zNegFaceRhoQi =
        0.25*(array[(block::arrayIndex(i+1,j+1,k+0))] +
              array[(block::arrayIndex(i+1,j+0,k+0))] +
              array[(block::arrayIndex(i+0,j+0,k+0))] +
              array[(block::arrayIndex(i+0,j+1,k+0))]);
      // calculate linear gradient of faceRhoQi as a cell volume average
      const Real xCellGradRhoQi = (xPosFaceRhoQi - xNegFaceRhoQi)/Hybrid::dx;
      const Real yCellGradRhoQi = (yPosFaceRhoQi - yNegFaceRhoQi)/Hybrid::dx;
      const Real zCellGradRhoQi = (zPosFaceRhoQi - zNegFaceRhoQi)/Hybrid::dx;
      // calculate electron pressure term of the electric field
      // adiabatic electrons
      if(Hybrid::useAdiabaticElectronPressure == true) {
         cellEp[n3+0] = -Hybrid::electronPressureCoeff*xCellGradRhoQi;
         cellEp[n3+1] = -Hybrid::electronPressureCoeff*yCellGradRhoQi;
         cellEp[n3+2] = -Hybrid::electronPressureCoeff*zCellGradRhoQi;
      }
      else {
         // isothermal electrons
         if(fabs(cellRhoQi[n]) > 0) {
            cellEp[n3+0] = -Hybrid::electronPressureCoeff*xCellGradRhoQi/cellRhoQi[n];
            cellEp[n3+1] = -Hybrid::electronPressureCoeff*yCellGradRhoQi/cellRhoQi[n];
            cellEp[n3+2] = -Hybrid::electronPressureCoeff*zCellGradRhoQi/cellRhoQi[n];
         }
      }
   }
}

// interpolation from faces to arbitrary point r (in block's local coordinates)
void face2r(Real* r,Real* faceData,Simulation& sim,SimulationClasses& simClasses,pargrid::CellID blockID,Real* result)
{
   if(simClasses.pargrid.getNeighbourFlags(blockID) != pargrid::ALL_NEIGHBOURS_EXIST) return;
   const unsigned int size = (block::WIDTH_X+2)*(block::WIDTH_Y+2)*(block::WIDTH_Z+2);
   Real array[size*3];
   fetchData(faceData,array,simClasses,blockID,3);
   
   // dimensionless local coordinates
   const Real x=r[0]/Hybrid::dx;
   const Real y=r[1]/Hybrid::dx;
   const Real z=r[2]/Hybrid::dx;

   // local indices of the cell the point is in
   const int i = static_cast<int>(floor(x));
   const int j = static_cast<int>(floor(y));
   const int k = static_cast<int>(floor(z));
      
   // face center local dimensionless coordinates / indices
   const Real xFaceNeg = i+0;
   const Real xFacePos = i+1;
   const Real yFaceNeg = j+0;
   const Real yFacePos = j+1;
   const Real zFaceNeg = k+0;
   const Real zFacePos = k+1;

   // weight factors
   const Real w111x = x - xFaceNeg;
   const Real w111y = y - yFaceNeg;
   const Real w111z = z - zFaceNeg;
   const Real w011x = xFacePos - x;
   const Real w101y = yFacePos - y;
   const Real w110z = zFacePos - z;

   // linear interpolation between faces
   result[0] = w111x*array[(block::arrayIndex(i+1,j+1,k+1))*3+0] + w011x*array[(block::arrayIndex(i+0,j+1,k+1))*3+0];
   result[1] = w111y*array[(block::arrayIndex(i+1,j+1,k+1))*3+1] + w101y*array[(block::arrayIndex(i+1,j+0,k+1))*3+1];
   result[2] = w111z*array[(block::arrayIndex(i+1,j+1,k+1))*3+2] + w110z*array[(block::arrayIndex(i+1,j+1,k+0))*3+2];
}

// interpolation from nodes to arbitrary point r (in block's local coordinates)
void node2r(Real* r,Real* nodeData,Simulation& sim,SimulationClasses& simClasses,pargrid::CellID blockID,Real* result)
{
   if(simClasses.pargrid.getNeighbourFlags(blockID) != pargrid::ALL_NEIGHBOURS_EXIST) return;
   const unsigned int size = (block::WIDTH_X+2)*(block::WIDTH_Y+2)*(block::WIDTH_Z+2);
   Real array[size*3];
   fetchData(nodeData,array,simClasses,blockID,3);
   
   // dimensionless local coordinates
   const Real x=r[0]/Hybrid::dx;
   const Real y=r[1]/Hybrid::dx;
   const Real z=r[2]/Hybrid::dx;

   // local indices of the cell the point is in
   const int i = static_cast<int>(floor(x));
   const int j = static_cast<int>(floor(y));
   const int k = static_cast<int>(floor(z));
   
   // node local dimensionless coordinates / indices
   const Real xNode111 = i+1; const Real yNode111 = j+1; const Real zNode111 = k+1;
   const Real xNode011 = i+0; const Real yNode011 = j+1; const Real zNode011 = k+1;
   const Real xNode101 = i+1; const Real yNode101 = j+0; const Real zNode101 = k+1;
   const Real xNode001 = i+0; const Real yNode001 = j+0; const Real zNode001 = k+1;
   const Real xNode110 = i+1; const Real yNode110 = j+1; const Real zNode110 = k+0;
   const Real xNode010 = i+0; const Real yNode010 = j+1; const Real zNode010 = k+0;
   const Real xNode100 = i+1; const Real yNode100 = j+0; const Real zNode100 = k+0;
   const Real xNode000 = i+0; const Real yNode000 = j+0; const Real zNode000 = k+0;
   
   // weight factors
   const Real w111 = 1/sqrt(sqr(xNode111-x) + sqr(yNode111-y) + sqr(zNode111-z));
   const Real w011 = 1/sqrt(sqr(xNode011-x) + sqr(yNode011-y) + sqr(zNode011-z));
   const Real w101 = 1/sqrt(sqr(xNode101-x) + sqr(yNode101-y) + sqr(zNode101-z));
   const Real w001 = 1/sqrt(sqr(xNode001-x) + sqr(yNode001-y) + sqr(zNode001-z));
   const Real w110 = 1/sqrt(sqr(xNode110-x) + sqr(yNode110-y) + sqr(zNode110-z));
   const Real w010 = 1/sqrt(sqr(xNode010-x) + sqr(yNode010-y) + sqr(zNode010-z));
   const Real w100 = 1/sqrt(sqr(xNode100-x) + sqr(yNode100-y) + sqr(zNode100-z));
   const Real w000 = 1/sqrt(sqr(xNode000-x) + sqr(yNode000-y) + sqr(zNode000-z));
   const Real wsum = w111+w011+w101+w001+w110+w010+w100+w000;

   if(wsum > 0) {
      for(int l=0;l<3;++l) {
	 result[l] = (w111*array[(block::arrayIndex(i+1,j+1,k+1))*3+l] +
		      w011*array[(block::arrayIndex(i+0,j+1,k+1))*3+l] +
		      w101*array[(block::arrayIndex(i+1,j+0,k+1))*3+l] +
		      w001*array[(block::arrayIndex(i+0,j+0,k+1))*3+l] +
		      w110*array[(block::arrayIndex(i+1,j+1,k+0))*3+l] +
		      w010*array[(block::arrayIndex(i+0,j+1,k+0))*3+l] +
		      w100*array[(block::arrayIndex(i+1,j+0,k+0))*3+l] +
		      w000*array[(block::arrayIndex(i+0,j+0,k+0))*3+l])/wsum;
      }
   }
   else {
      result[0] = result[1] = result[2] = 0.0;
   }
}

// zero order interpolation from cells to arbitrary point r (in block's local coordinates)
void cell2r(Real* r,Real* cellData,Simulation& sim,SimulationClasses& simClasses,pargrid::CellID blockID,Real* result) {
   if(simClasses.pargrid.getNeighbourFlags(blockID) != pargrid::ALL_NEIGHBOURS_EXIST) return;
   const unsigned int size = (block::WIDTH_X+2)*(block::WIDTH_Y+2)*(block::WIDTH_Z+2);
   Real array[size*3];
   fetchData(cellData,array,simClasses,blockID,3);
   
   // dimensionless local coordinates
   const Real x=r[0]/Hybrid::dx;
   const Real y=r[1]/Hybrid::dx;
   const Real z=r[2]/Hybrid::dx;

   // local indices of the cell the point is in
   const int i = static_cast<int>(floor(x));
   const int j = static_cast<int>(floor(y));
   const int k = static_cast<int>(floor(z));
   
   for(int l=0;l<3;++l) {
      result[l] = array[(block::arrayIndex(i+1,j+1,k+1))*3+l];
   }
}

// exchange faceB and nodeE
void setupGetFields(Simulation& sim,SimulationClasses& simClasses) {
   profile::start("setupGetFields",setupGetFieldsID);
   Real* faceB = simClasses.pargrid.getUserDataStatic<Real>(Hybrid::dataFaceBID);
   Real* nodeE = simClasses.pargrid.getUserDataStatic<Real>(Hybrid::dataNodeEID);
   if(faceB == NULL) {cerr << "ERROR: obtained NULL faceB array!" << endl; exit(1);}
   if(nodeE == NULL) {cerr << "ERROR: obtained NULL nodeE array!" << endl; exit(1);}
   simClasses.pargrid.startNeighbourExchange(pargrid::DEFAULT_STENCIL,Hybrid::dataFaceBID);
   simClasses.pargrid.startNeighbourExchange(pargrid::DEFAULT_STENCIL,Hybrid::dataNodeEID);
   profile::start("MPI waits",mpiWaitID);
   simClasses.pargrid.wait(pargrid::DEFAULT_STENCIL,Hybrid::dataFaceBID);
   simClasses.pargrid.wait(pargrid::DEFAULT_STENCIL,Hybrid::dataNodeEID);
   profile::stop();
   profile::stop();
}

// get E, B and Ue fields at arbitrary point r
void getFields(Real* r,Real* B,Real* Ue,Real* Ep,Simulation& sim,SimulationClasses& simClasses,pargrid::CellID blockID) {
   Real* faceB  = reinterpret_cast<Real*>(simClasses.pargrid.getUserData(Hybrid::dataFaceBID));
   //Real* nodeE  = reinterpret_cast<Real*>(simClasses.pargrid.getUserData(Hybrid::dataNodeEID));
   Real* cellUe = reinterpret_cast<Real*>(simClasses.pargrid.getUserData(Hybrid::dataCellUeID));
   Real* cellEp = reinterpret_cast<Real*>(simClasses.pargrid.getUserData(Hybrid::dataCellEpID));
   face2r(r,faceB,sim,simClasses,blockID,B);
   //node2r(r,nodeE,sim,simClasses,blockID,E);
   cell2r(r,cellUe,sim,simClasses,blockID,Ue);
   if(Hybrid::useElectronPressureElectricField == true) {
      cell2r(r,cellEp,sim,simClasses,blockID,Ep);
   }
}

// fetch ALL neighours
void fetchData(Real* data,Real* array,SimulationClasses& simClasses,pargrid::CellID blockID,int vectorDim) {
   pargrid::CellID nbrLID = simClasses.pargrid.invalid();
   const pargrid::CellID* const nbrs = simClasses.pargrid.getCellNeighbourIDs(blockID);   
   
   // THIS BLOCK
      
   for(int k=0;k<block::WIDTH_Z;++k) for(int j=0;j<block::WIDTH_Y;++j) for(int i=0;i<block::WIDTH_X;++i) for(int l=0;l<vectorDim;++l) {
      array[block::arrayIndex(i+1,j+1,k+1)*vectorDim+l] = data[(blockID*block::SIZE+block::index(i,j,k))*vectorDim+l];
   }
   
   // FACE NEIGHBOURS

   // -x
   nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(-1,+0,+0)];
   if(nbrLID != simClasses.pargrid.invalid()) {
      for(int k=0;k<block::WIDTH_Z;++k) for(int j=0;j<block::WIDTH_Y;++j) for(int l=0;l<vectorDim;++l) {
	 array[block::arrayIndex(0,j+1,k+1)*vectorDim+l] = data[(nbrLID*block::SIZE+block::index(block::WIDTH_X-1,j,k))*vectorDim+l];
      }
   }
   // -y
   nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+0,-1,+0)];
   if(nbrLID != simClasses.pargrid.invalid()) {
      for(int k=0;k<block::WIDTH_Z;++k) for(int i=0;i<block::WIDTH_X;++i) for(int l=0;l<vectorDim;++l) {
	 array[block::arrayIndex(i+1,0,k+1)*vectorDim+l] = data[(nbrLID*block::SIZE+block::index(i,block::WIDTH_Y-1,k))*vectorDim+l];
      }
   }
   // -z
   nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+0,+0,-1)];
   if(nbrLID != simClasses.pargrid.invalid()) {
      for(int j=0;j<block::WIDTH_Y;++j) for(int i=0;i<block::WIDTH_X;++i) for(int l=0;l<vectorDim;++l) {
	 array[block::arrayIndex(i+1,j+1,0)*vectorDim+l] = data[(nbrLID*block::SIZE+block::index(i,j,block::WIDTH_Z-1))*vectorDim+l];
      }
   }
   // +x
   nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+1,+0,+0)];
   if(nbrLID != simClasses.pargrid.invalid()) {
      for(int k=0;k<block::WIDTH_Z;++k) for(int j=0;j<block::WIDTH_Y;++j) for(int l=0;l<vectorDim;++l) {
	 array[block::arrayIndex(block::WIDTH_X+1,j+1,k+1)*vectorDim+l] = data[(nbrLID*block::SIZE+block::index(0,j,k))*vectorDim+l];
      }
   }
   // +y
   nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+0,+1,+0)];
   if(nbrLID != simClasses.pargrid.invalid()) {
      for(int k=0;k<block::WIDTH_Z;++k) for(int i=0;i<block::WIDTH_X;++i) for(int l=0;l<vectorDim;++l) {
	 array[block::arrayIndex(i+1,block::WIDTH_Y+1,k+1)*vectorDim+l] = data[(nbrLID*block::SIZE+block::index(i,0,k))*vectorDim+l];
      }
   }
   // +z
   nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+0,+0,+1)];
   if(nbrLID != simClasses.pargrid.invalid()) {
      for(int j=0;j<block::WIDTH_Y;++j) for(int i=0;i<block::WIDTH_X;++i) for(int l=0;l<vectorDim;++l) {
	 array[block::arrayIndex(i+1,j+1,block::WIDTH_Z+1)*vectorDim+l] = data[(nbrLID*block::SIZE+block::index(i,j,0))*vectorDim+l];
      }
   }
      
   // EDGE NEIGHBOURS

   // -x,-y
   nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(-1,-1,+0)];
   if(nbrLID != simClasses.pargrid.invalid()) {
      for(int k=0;k<block::WIDTH_Z;++k) for(int l=0;l<vectorDim;++l) {
	 array[block::arrayIndex(0,0,k+1)*vectorDim+l] = data[(nbrLID*block::SIZE+block::index(block::WIDTH_X-1,block::WIDTH_Y-1,k))*vectorDim+l];
      }
   }
   // +x,-y
   nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+1,-1,+0)];
   if(nbrLID != simClasses.pargrid.invalid()) {
      for(int k=0;k<block::WIDTH_Z;++k) for(int l=0;l<vectorDim;++l) {
	 array[block::arrayIndex(block::WIDTH_X+1,0,k+1)*vectorDim+l] = data[(nbrLID*block::SIZE+block::index(0,block::WIDTH_Y-1,k))*vectorDim+l];
      }
   }
   // -x,+y
   nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(-1,+1,+0)];
   if(nbrLID != simClasses.pargrid.invalid()) {
      for(int k=0;k<block::WIDTH_Z;++k) for(int l=0;l<vectorDim;++l) {
	 array[block::arrayIndex(0,block::WIDTH_Y+1,k+1)*vectorDim+l] = data[(nbrLID*block::SIZE+block::index(block::WIDTH_X-1,0,k))*vectorDim+l];
      }
   }
   // +x,+y
   nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+1,+1,+0)];
   if(nbrLID != simClasses.pargrid.invalid()) {
      for(int k=0;k<block::WIDTH_Z;++k) for(int l=0;l<vectorDim;++l) {
	 array[block::arrayIndex(block::WIDTH_X+1,block::WIDTH_Y+1,k+1)*vectorDim+l] = data[(nbrLID*block::SIZE+block::index(0,0,k))*vectorDim+l];
      }
   }
   // -y,-z
   nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+0,-1,-1)];
   if(nbrLID != simClasses.pargrid.invalid()) {
      for(int i=0;i<block::WIDTH_X;++i) for(int l=0;l<vectorDim;++l) {
	 array[block::arrayIndex(i+1,0,0)*vectorDim+l] = data[(nbrLID*block::SIZE+block::index(i,block::WIDTH_Y-1,block::WIDTH_Z-1))*vectorDim+l];
      }
   }
   // +y,-z
   nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+0,+1,-1)];
   if(nbrLID != simClasses.pargrid.invalid()) {
      for(int i=0;i<block::WIDTH_X;++i) for(int l=0;l<vectorDim;++l) {
	 array[block::arrayIndex(i+1,block::WIDTH_Y+1,0)*vectorDim+l] = data[(nbrLID*block::SIZE+block::index(i,0,block::WIDTH_Z-1))*vectorDim+l];
      }
   }
   // -y,+z
   nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+0,-1,+1)];
   if(nbrLID != simClasses.pargrid.invalid()) {
      for(int i=0;i<block::WIDTH_X;++i) for(int l=0;l<vectorDim;++l) {
	 array[block::arrayIndex(i+1,0,block::WIDTH_Z+1)*vectorDim+l] = data[(nbrLID*block::SIZE+block::index(i,block::WIDTH_Y-1,0))*vectorDim+l];
      }
   }
   // +y,+z
   nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+0,+1,+1)];
   if(nbrLID != simClasses.pargrid.invalid()) {
      for(int i=0;i<block::WIDTH_X;++i) for(int l=0;l<vectorDim;++l) {
	 array[block::arrayIndex(i+1,block::WIDTH_Y+1,block::WIDTH_Z+1)*vectorDim+l] = data[(nbrLID*block::SIZE+block::index(i,0,0))*vectorDim+l];
      }
   }
   // -x,-z
   nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(-1,+0,-1)];
   if(nbrLID != simClasses.pargrid.invalid()) {
      for(int j=0;j<block::WIDTH_Y;++j) for(int l=0;l<vectorDim;++l) {
	 array[block::arrayIndex(0,j+1,0)*vectorDim+l] = data[(nbrLID*block::SIZE+block::index(block::WIDTH_X-1,j,block::WIDTH_Z-1))*vectorDim+l];
      }
   }
   // +x,-z
   nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+1,+0,-1)];
   if(nbrLID != simClasses.pargrid.invalid()) {
      for(int j=0;j<block::WIDTH_Y;++j) for(int l=0;l<vectorDim;++l) {
	 array[block::arrayIndex(block::WIDTH_X+1,j+1,0)*vectorDim+l] = data[(nbrLID*block::SIZE+block::index(0,j,block::WIDTH_Z-1))*vectorDim+l];
      }
   }
   // -x,+z
   nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(-1,+0,+1)];
   if(nbrLID != simClasses.pargrid.invalid()) {
      for(int j=0;j<block::WIDTH_Y;++j) for(int l=0;l<vectorDim;++l) {
	 array[block::arrayIndex(0,j+1,block::WIDTH_Z+1)*vectorDim+l] = data[(nbrLID*block::SIZE+block::index(block::WIDTH_X-1,j,0))*vectorDim+l];
      }
   }
   // +x,+z
   nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+1,+0,+1)];
   if(nbrLID != simClasses.pargrid.invalid()) {
      for(int j=0;j<block::WIDTH_Y;++j) for(int l=0;l<vectorDim;++l) {
	 array[block::arrayIndex(block::WIDTH_X+1,j+1,block::WIDTH_Z+1)*vectorDim+l] = data[(nbrLID*block::SIZE+block::index(0,j,0))*vectorDim+l];
      }
   }
   
   // CORNER NEIGHBOURS
    
   // -x,-y,-z
   nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(-1,-1,-1)];
   if(nbrLID != simClasses.pargrid.invalid()) {
      for(int l=0;l<vectorDim;++l) {
	 array[block::arrayIndex(0,0,0)*vectorDim+l] = data[(nbrLID*block::SIZE+block::index(block::WIDTH_X-1,block::WIDTH_Y-1,block::WIDTH_Z-1))*vectorDim+l];
      }
   }
   // +x,-y,-z
   nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+1,-1,-1)];
   if(nbrLID != simClasses.pargrid.invalid()) {
      for(int l=0;l<vectorDim;++l) {
	 array[block::arrayIndex(block::WIDTH_X+1,0,0)*vectorDim+l] = data[(nbrLID*block::SIZE+block::index(0,block::WIDTH_Y-1,block::WIDTH_Z-1))*vectorDim+l];
      }
   }
   // -x,+y,-z
   nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(-1,+1,-1)];
   if(nbrLID != simClasses.pargrid.invalid()) {
      for(int l=0;l<vectorDim;++l) {
	 array[block::arrayIndex(0,block::WIDTH_Y+1,0)*vectorDim+l] = data[(nbrLID*block::SIZE+block::index(block::WIDTH_X-1,0,block::WIDTH_Z-1))*vectorDim+l];
      }
   }
   // +x,+y,-z
   nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+1,+1,-1)];
   if(nbrLID != simClasses.pargrid.invalid()) {
      for(int l=0;l<vectorDim;++l) {
	 array[block::arrayIndex(block::WIDTH_X+1,block::WIDTH_Y+1,0)*vectorDim+l] = data[(nbrLID*block::SIZE+block::index(0,0,block::WIDTH_Z-1))*vectorDim+l];
      }
   }
   // -x,-y,+z
   nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(-1,-1,+1)];
   if(nbrLID != simClasses.pargrid.invalid()) {
      for(int l=0;l<vectorDim;++l) {
	 array[block::arrayIndex(0,0,block::WIDTH_Z+1)*vectorDim+l] = data[(nbrLID*block::SIZE+block::index(block::WIDTH_X-1,block::WIDTH_Y-1,0))*vectorDim+l];
      }
   }
   // +x,-y,+z
   nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+1,-1,+1)];
   if(nbrLID != simClasses.pargrid.invalid()) {
      for(int l=0;l<vectorDim;++l) {
	 array[block::arrayIndex(block::WIDTH_X+1,0,block::WIDTH_Z+1)*vectorDim+l] = data[(nbrLID*block::SIZE+block::index(0,block::WIDTH_Y-1,0))*vectorDim+l];
      }
   }
   // -x,+y,+z
   nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(-1,+1,+1)];
   if(nbrLID != simClasses.pargrid.invalid()) {
      for(int l=0;l<vectorDim;++l) {
	 array[block::arrayIndex(0,block::WIDTH_Y+1,block::WIDTH_Z+1)*vectorDim+l] = data[(nbrLID*block::SIZE+block::index(block::WIDTH_X-1,0,0))*vectorDim+l];
      }
   }
   // +x,+y,+z
   nbrLID = nbrs[simClasses.pargrid.calcNeighbourTypeID(+1,+1,+1)];
   if(nbrLID != simClasses.pargrid.invalid()) {
      for(int l=0;l<vectorDim;++l) {
	 array[block::arrayIndex(block::WIDTH_X+1,block::WIDTH_Y+1,block::WIDTH_Z+1)*vectorDim+l] = data[(nbrLID*block::SIZE+block::index(0,0,0))*vectorDim+l];
      }
   }
}
