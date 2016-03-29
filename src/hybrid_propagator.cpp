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
   Real* faceB        = simClasses.pargrid.getUserDataStatic<Real>(Hybrid::dataFaceBID);
   Real* faceJ        = simClasses.pargrid.getUserDataStatic<Real>(Hybrid::dataFaceJID);
   Real* cellRhoQi    = simClasses.pargrid.getUserDataStatic<Real>(Hybrid::dataCellRhoQiID);
   Real* cellB        = simClasses.pargrid.getUserDataStatic<Real>(Hybrid::dataCellBID);
   Real* cellJ        = simClasses.pargrid.getUserDataStatic<Real>(Hybrid::dataCellJID);
   Real* cellUe       = simClasses.pargrid.getUserDataStatic<Real>(Hybrid::dataCellUeID);
   Real* cellJi       = simClasses.pargrid.getUserDataStatic<Real>(Hybrid::dataCellJiID);
   Real* cellMaxUe    = simClasses.pargrid.getUserDataStatic<Real>(Hybrid::dataCellMaxUeID);
   Real* cellMaxVi    = simClasses.pargrid.getUserDataStatic<Real>(Hybrid::dataCellMaxViID);
   Real* cellMinRhoQi = simClasses.pargrid.getUserDataStatic<Real>(Hybrid::dataCellMinRhoQiID);
   Real* nodeRhoQi    = simClasses.pargrid.getUserDataStatic<Real>(Hybrid::dataNodeRhoQiID);
   Real* nodeE        = simClasses.pargrid.getUserDataStatic<Real>(Hybrid::dataNodeEID);
   Real* nodeB        = simClasses.pargrid.getUserDataStatic<Real>(Hybrid::dataNodeBID);
   Real* nodeJ        = simClasses.pargrid.getUserDataStatic<Real>(Hybrid::dataNodeJID);
   Real* nodeUe       = simClasses.pargrid.getUserDataStatic<Real>(Hybrid::dataNodeUeID);
   Real* nodeJi       = simClasses.pargrid.getUserDataStatic<Real>(Hybrid::dataNodeJiID);
   bool* innerFlag    = simClasses.pargrid.getUserDataStatic<bool>(Hybrid::dataInnerFlagFieldID);
   bool* innerFlagNode= simClasses.pargrid.getUserDataStatic<bool>(Hybrid::dataInnerFlagNodeID);

   if(faceB        == NULL) {cerr << "ERROR: obtained NULL faceB array!"        << endl; exit(1);}
   if(faceJ        == NULL) {cerr << "ERROR: obtained NULL faceJ array!"        << endl; exit(1);}
   if(cellRhoQi    == NULL) {cerr << "ERROR: obtained NULL cellRhoQi array!"    << endl; exit(1);}
   if(cellB        == NULL) {cerr << "ERROR: obtained NULL cellB array!"        << endl; exit(1);}
   if(cellJ        == NULL) {cerr << "ERROR: obtained NULL cellJ array!"        << endl; exit(1);}
   if(cellUe       == NULL) {cerr << "ERROR: obtained NULL cellUe array!"       << endl; exit(1);}
   if(cellJi       == NULL) {cerr << "ERROR: obtained NULL cellJi array!"       << endl; exit(1);}
   if(cellMaxUe    == NULL) {cerr << "ERROR: obtained NULL cellMaxUe array!"    << endl; exit(1);}
   if(cellMaxVi    == NULL) {cerr << "ERROR: obtained NULL cellMaxVi array!"    << endl; exit(1);}
   if(cellMinRhoQi == NULL) {cerr << "ERROR: obtained NULL cellMinRhoQi array!" << endl; exit(1);}
   if(nodeRhoQi    == NULL) {cerr << "ERROR: obtained NULL nodeRhoQi array!"    << endl; exit(1);}
   if(nodeE        == NULL) {cerr << "ERROR: obtained NULL nodeE array!"        << endl; exit(1);}
   if(nodeB        == NULL) {cerr << "ERROR: obtained NULL nodeB array!"        << endl; exit(1);}
   if(nodeJ        == NULL) {cerr << "ERROR: obtained NULL nodeJ array!"        << endl; exit(1);}
   if(nodeUe       == NULL) {cerr << "ERROR: obtained NULL nodeUe array!"       << endl; exit(1);}
   if(nodeJi       == NULL) {cerr << "ERROR: obtained NULL nodeJi array!"       << endl; exit(1);}
   if(innerFlag    == NULL) {cerr << "ERROR: obtained NULL innerFlag array!"    << endl; exit(1);}
   if(innerFlagNode== NULL) {cerr << "ERROR: obtained NULL innerFlagNode array!"<< endl; exit(1);}
   
   // get block vectors
   const vector<pargrid::CellID>& innerBlocks = simClasses.pargrid.getInnerCells(pargrid::DEFAULT_STENCIL);
   const vector<pargrid::CellID>& boundaryBlocks = simClasses.pargrid.getBoundaryCells(pargrid::DEFAULT_STENCIL);
   const vector<pargrid::CellID>& exteriorBlocks = simClasses.pargrid.getExteriorCells();

   // zero diagnostic variables
   if(saveStepHappened == true) {
      for(pargrid::CellID b=0;b<simClasses.pargrid.getNumberOfLocalCells();++b) for(int k=0;k<block::WIDTH_Z;++k) for(int j=0;j<block::WIDTH_Y;++j) for(int i=0;i<block::WIDTH_X;++i) {
	 const int n = (b*block::SIZE+block::index(i,j,k));
	 cellMaxUe[n]=0.0;
	 cellMaxVi[n]=0.0;
	 cellMinRhoQi[n]=0.0;
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

   // calculate J
#ifdef USE_EDGE_J
   neumannFace(faceB,sim,simClasses,exteriorBlocks);
   // nodeJ = avg(edgeJ) = avg(curl(faceB)/mu0)
   simClasses.pargrid.startNeighbourExchange(pargrid::DEFAULT_STENCIL,Hybrid::dataFaceBID);
   profile::start("field propag",profPropagFieldID);
   for(pargrid::CellID b=0; b<innerBlocks.size(); ++b) { calcNodeJ(faceB,nodeJ,sim,simClasses,innerBlocks[b]); }
   profile::stop();
   profile::start("MPI waits",mpiWaitID);
   simClasses.pargrid.wait(pargrid::DEFAULT_STENCIL,Hybrid::dataFaceBID);
   profile::stop();
   profile::start("field propag",profPropagFieldID);
   for(pargrid::CellID b=0; b<boundaryBlocks.size(); ++b) { calcNodeJ(faceB,nodeJ,sim,simClasses,boundaryBlocks[b]); }
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
   for(pargrid::CellID b=0; b<simClasses.pargrid.getNumberOfLocalCells(); ++b) { calcCellUe(cellJ,cellJi,cellRhoQi,cellUe,innerFlag,cellMaxUe,sim,simClasses,b); }
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
   profile::start("field propag",profPropagFieldID);
   for(pargrid::CellID b=0; b<simClasses.pargrid.getNumberOfLocalCells(); ++b) { calcNodeUe(nodeRhoQi,nodeJi,nodeJ,nodeUe,innerFlagNode,cellMaxUe,sim,simClasses,b); }
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
   for(pargrid::CellID b=0; b<simClasses.pargrid.getNumberOfLocalCells(); ++b) { calcNodeE(nodeUe,nodeB,nodeJ,nodeE,innerFlagNode,sim,simClasses,b); }
   profile::stop();

   // E filtering
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
      // edge (-x+y)
      if((nf & Hybrid::X_NEG_EXISTS) == 0 && (nf & Hybrid::Y_POS_EXISTS) == 0) {
	 const int i=0,j=0,k=0;
	 const int n = (b*block::SIZE+block::index(i,j,k))*vectorDim;
	 const int m = block::arrayIndex(i+2,j+0,k+1)*vectorDim;
	 for(int l=0;l<vectorDim;++l) { cellData[n+l] = tempArrayCellData[m+l]; }
      }
      // edge (-x-y)
      if((nf & Hybrid::X_NEG_EXISTS) == 0 && (nf & Hybrid::Y_NEG_EXISTS) == 0) {
	 const int i=0,j=0,k=0;
	 const int n = (b*block::SIZE+block::index(i,j,k))*vectorDim;
	 const int m = block::arrayIndex(i+2,j+2,k+1)*vectorDim;
	 for(int l=0;l<vectorDim;++l) { cellData[n+l] = tempArrayCellData[m+l]; }
      }
      // edge (-x+z)
      if((nf & Hybrid::X_NEG_EXISTS) == 0 && (nf & Hybrid::Z_POS_EXISTS) == 0) {
	 const int i=0,j=0,k=0;
	 const int n = (b*block::SIZE+block::index(i,j,k))*vectorDim;
	 const int m = block::arrayIndex(i+2,j+1,k+0)*vectorDim;
	 for(int l=0;l<vectorDim;++l) { cellData[n+l] = tempArrayCellData[m+l]; }
      }
      // edge (-x-z)
      if((nf & Hybrid::X_NEG_EXISTS) == 0 && (nf & Hybrid::Z_NEG_EXISTS) == 0) {
	 const int i=0,j=0,k=0;
	 const int n = (b*block::SIZE+block::index(i,j,k))*vectorDim;
	 const int m = block::arrayIndex(i+2,j+1,k+2)*vectorDim;
	 for(int l=0;l<vectorDim;++l) { cellData[n+l] = tempArrayCellData[m+l]; }
      }
      // edge (+y+z)
      if((nf & Hybrid::Y_POS_EXISTS) == 0 && (nf & Hybrid::Z_POS_EXISTS) == 0) {
	 const int i=0,j=0,k=0;
	 const int n = (b*block::SIZE+block::index(i,j,k))*vectorDim;
	 const int m = block::arrayIndex(i+1,j+0,k+0)*vectorDim;
	 for(int l=0;l<vectorDim;++l) { cellData[n+l] = tempArrayCellData[m+l]; }
      }
      // edge (+y-z)
      if((nf & Hybrid::Y_POS_EXISTS) == 0 && (nf & Hybrid::Z_NEG_EXISTS) == 0) {
	 const int i=0,j=0,k=0;
	 const int n = (b*block::SIZE+block::index(i,j,k))*vectorDim;
	 const int m = block::arrayIndex(i+1,j+0,k+2)*vectorDim;
	 for(int l=0;l<vectorDim;++l) { cellData[n+l] = tempArrayCellData[m+l]; }
      }
      // edge (-y+z)
      if((nf & Hybrid::Y_NEG_EXISTS) == 0 && (nf & Hybrid::Z_POS_EXISTS) == 0) {
	 const int i=0,j=0,k=0;
	 const int n = (b*block::SIZE+block::index(i,j,k))*vectorDim;
	 const int m = block::arrayIndex(i+1,j+2,k+0)*vectorDim;
	 for(int l=0;l<vectorDim;++l) { cellData[n+l] = tempArrayCellData[m+l]; }
      }
      // edge (-y-z)
      if((nf & Hybrid::Y_NEG_EXISTS) == 0 && (nf & Hybrid::Z_NEG_EXISTS) == 0) {
	 const int i=0,j=0,k=0;
	 const int n = (b*block::SIZE+block::index(i,j,k))*vectorDim;
	 const int m = block::arrayIndex(i+1,j+2,k+2)*vectorDim;
	 for(int l=0;l<vectorDim;++l) { cellData[n+l] = tempArrayCellData[m+l]; }
      }
      // edge (+x+y)
      if((nf & Hybrid::X_POS_EXISTS) == 0 && (nf & Hybrid::Y_POS_EXISTS) == 0) {
	 const int i=0,j=0,k=0;
	 const int n = (b*block::SIZE+block::index(i,j,k))*vectorDim;
	 const int m = block::arrayIndex(i+0,j+0,k+1)*vectorDim;
	 for(int l=0;l<vectorDim;++l) { cellData[n+l] = tempArrayCellData[m+l]; }
      }
      // edge (+x-y)
      if((nf & Hybrid::X_POS_EXISTS) == 0 && (nf & Hybrid::Y_NEG_EXISTS) == 0) {
	 const int i=0,j=0,k=0;
	 const int n = (b*block::SIZE+block::index(i,j,k))*vectorDim;
	 const int m = block::arrayIndex(i+0,j+2,k+1)*vectorDim;
	 for(int l=0;l<vectorDim;++l) { cellData[n+l] = tempArrayCellData[m+l]; }
      }
      // edge (+x+z)
      if((nf & Hybrid::X_POS_EXISTS) == 0 && (nf & Hybrid::Z_POS_EXISTS) == 0) {
	 const int i=0,j=0,k=0;
	 const int n = (b*block::SIZE+block::index(i,j,k))*vectorDim;
	 const int m = block::arrayIndex(i+0,j+1,k+0)*vectorDim;
	 for(int l=0;l<vectorDim;++l) { cellData[n+l] = tempArrayCellData[m+l]; }
      }
      // edge (+x-z)
      if((nf & Hybrid::X_POS_EXISTS) == 0 && (nf & Hybrid::Z_NEG_EXISTS) == 0) {
	 const int i=0,j=0,k=0;
	 const int n = (b*block::SIZE+block::index(i,j,k))*vectorDim;
	 const int m = block::arrayIndex(i+0,j+1,k+2)*vectorDim;
	 for(int l=0;l<vectorDim;++l) { cellData[n+l] = tempArrayCellData[m+l]; }
      }
      // corner (-x-y-z)
      if((nf & Hybrid::X_NEG_EXISTS) == 0 && (nf & Hybrid::Y_NEG_EXISTS) == 0 && (nf & Hybrid::Z_NEG_EXISTS) == 0) {
	 const int i=0,j=0,k=0;
	 const int n = (b*block::SIZE+block::index(i,j,k))*vectorDim;
	 const int m = block::arrayIndex(i+2,j+2,k+2)*vectorDim;
	 for(int l=0;l<vectorDim;++l) { cellData[n+l] = tempArrayCellData[m+l]; }
      }
      // corner (-x-y+z)
      if((nf & Hybrid::X_NEG_EXISTS) == 0 && (nf & Hybrid::Y_NEG_EXISTS) == 0 && (nf & Hybrid::Z_POS_EXISTS) == 0) {
	 const int i=0,j=0,k=0;
	 const int n = (b*block::SIZE+block::index(i,j,k))*vectorDim;
	 const int m = block::arrayIndex(i+2,j+2,k+0)*vectorDim;
	 for(int l=0;l<vectorDim;++l) { cellData[n+l] = tempArrayCellData[m+l]; }
      }
      // corner (-x+y+z)
      if((nf & Hybrid::X_NEG_EXISTS) == 0 && (nf & Hybrid::Y_POS_EXISTS) == 0 && (nf & Hybrid::Z_POS_EXISTS) == 0) {
	 const int i=0,j=0,k=0;
	 const int n = (b*block::SIZE+block::index(i,j,k))*vectorDim;
	 const int m = block::arrayIndex(i+2,j+0,k+0)*vectorDim;
	 for(int l=0;l<vectorDim;++l) { cellData[n+l] = tempArrayCellData[m+l]; }
      }
      // corner (-x+y-z)
      if((nf & Hybrid::X_NEG_EXISTS) == 0 && (nf & Hybrid::Y_POS_EXISTS) == 0 && (nf & Hybrid::Z_NEG_EXISTS) == 0) {
	 const int i=0,j=0,k=0;
	 const int n = (b*block::SIZE+block::index(i,j,k))*vectorDim;
	 const int m = block::arrayIndex(i+2,j+0,k+2)*vectorDim;
	 for(int l=0;l<vectorDim;++l) { cellData[n+l] = tempArrayCellData[m+l]; }
      }
      // corner (+x-y-z)
      if((nf & Hybrid::X_POS_EXISTS) == 0 && (nf & Hybrid::Y_NEG_EXISTS) == 0 && (nf & Hybrid::Z_NEG_EXISTS) == 0) {
	 const int i=0,j=0,k=0;
	 const int n = (b*block::SIZE+block::index(i,j,k))*vectorDim;
	 const int m = block::arrayIndex(i+0,j+2,k+2)*vectorDim;
	 for(int l=0;l<vectorDim;++l) { cellData[n+l] = tempArrayCellData[m+l]; }
      }
      // corner (+x+y-z)
      if((nf & Hybrid::X_POS_EXISTS) == 0 && (nf & Hybrid::Y_POS_EXISTS) == 0 && (nf & Hybrid::Z_NEG_EXISTS) == 0) {
	 const int i=0,j=0,k=0;
	 const int n = (b*block::SIZE+block::index(i,j,k))*vectorDim;
	 const int m = block::arrayIndex(i+0,j+0,k+2)*vectorDim;
	 for(int l=0;l<vectorDim;++l) { cellData[n+l] = tempArrayCellData[m+l]; }
      }
      // corner (+x-y+z)
      if((nf & Hybrid::X_POS_EXISTS) == 0 && (nf & Hybrid::Y_NEG_EXISTS) == 0 && (nf & Hybrid::Z_POS_EXISTS) == 0) {
	 const int i=0,j=0,k=0;
	 const int n = (b*block::SIZE+block::index(i,j,k))*vectorDim;
	 const int m = block::arrayIndex(i+0,j+2,k+0)*vectorDim;
	 for(int l=0;l<vectorDim;++l) { cellData[n+l] = tempArrayCellData[m+l]; }
      }
      // corner (+x+y+z)
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
   Real* tempArrayCellData = new Real[tempArraySize*vectorDim];
   const std::vector<uint32_t>& neighbourFlags = simClasses.pargrid.getNeighbourFlags();
   
   for(pargrid::CellID eb=0;eb<exteriorBlocks.size();++eb) {
      const pargrid::CellID b = exteriorBlocks[eb];
      fetchData(faceData,tempArrayCellData,simClasses,b,vectorDim);
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
	    faceData[n+1] = tempArrayCellData[m+1];
            faceData[n+2] = tempArrayCellData[m+2];
	 }
      }
      // side (+y) wall
      if((nf & Hybrid::Y_POS_EXISTS) == 0) {
	 for(int k=0; k<block::WIDTH_Z; ++k) for(int i=0; i<block::WIDTH_X; ++i) {
	    const int j = block::WIDTH_Y-1-dj;
	    const int n = (b*block::SIZE+block::index(i,j,k))*vectorDim;
	    const int m = block::arrayIndex(i+1,j+0,k+1)*vectorDim;
	    faceData[n+0] = tempArrayCellData[m+0];
            faceData[n+2] = tempArrayCellData[m+2];
	 }
      }
      // side (-y) wall
      if((nf & Hybrid::Y_NEG_EXISTS) == 0) {
	 for(int k=0; k<block::WIDTH_Z; ++k) for(int i=0; i<block::WIDTH_X; ++i) {
	    const int j = 0+dj;
	    const int n = (b*block::SIZE+block::index(i,j,k))*vectorDim;
	    const int m = block::arrayIndex(i+1,j+2,k+1)*vectorDim;
	    faceData[n+0] = tempArrayCellData[m+0];
            faceData[n+2] = tempArrayCellData[m+2];
	 }
      }
      // side (+z) wall
      if((nf & Hybrid::Z_POS_EXISTS) == 0) {
	 for(int j=0; j<block::WIDTH_Y; ++j) for(int i=0; i<block::WIDTH_X; ++i) {
	    const int k = block::WIDTH_Z-1-dk;
	    const int n = (b*block::SIZE+block::index(i,j,k))*vectorDim;
	    const int m = block::arrayIndex(i+1,j+1,k+0)*vectorDim;
	    faceData[n+0] = tempArrayCellData[m+0];
            faceData[n+1] = tempArrayCellData[m+1];
	 }
      }
      // side (-z) wall
      if((nf & Hybrid::Z_NEG_EXISTS) == 0) {
	 for(int j=0; j<block::WIDTH_Y; ++j) for(int i=0; i<block::WIDTH_X; ++i) {
	    const int k = 0+dk;
	    const int n = (b*block::SIZE+block::index(i,j,k))*vectorDim;
	    const int m = block::arrayIndex(i+1,j+1,k+2)*vectorDim;
	    faceData[n+0] = tempArrayCellData[m+0];
            faceData[n+1] = tempArrayCellData[m+1];
	 }
      }
      // front (+x) wall
      if((nf & Hybrid::X_POS_EXISTS) == 0) {
	 for(int k=0; k<block::WIDTH_Z; ++k) for(int j=0; j<block::WIDTH_Y; ++j) {
	    const int i = 0;
	    const int n = (b*block::SIZE+block::index(i,j,k))*vectorDim;
	    const int m = block::arrayIndex(i+0,j+1,k+1)*vectorDim;
	    faceData[n+1] = tempArrayCellData[m+1];
            faceData[n+2] = tempArrayCellData[m+2];
	 }
      }
   }
   delete [] tempArrayCellData; tempArrayCellData = NULL;
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
void calcCellUe(Real* cellJ,Real* cellJi,Real* cellRhoQi,Real* cellUe,bool* innerFlag,Real* cellMaxUe,Simulation& sim,SimulationClasses& simClasses,pargrid::CellID blockID)
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
         cellMaxUe[n]++;
      }
   }
}
   

// E = -Ue x B + eta*J
void calcNodeE(Real* nodeUe,Real* nodeB,Real* nodeJ,Real* nodeE,bool* innerFlag,Simulation& sim,SimulationClasses& simClasses,pargrid::CellID blockID)
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
      if(innerFlag[n] == false) {
	 nodeE[n3+0] += Hybrid::eta*nodeJ[n3+0];
	 nodeE[n3+1] += Hybrid::eta*nodeJ[n3+1];
	 nodeE[n3+2] += Hybrid::eta*nodeJ[n3+2];
      }
   }
}

// nodeJ = nabla x faceB/mu0
void calcNodeJ(Real* faceB,Real* nodeJ,Simulation& sim,SimulationClasses& simClasses,pargrid::CellID blockID)
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
      const Real a = 0.5/(Hybrid::dx*constants::PERMEABILITY);
      nodeJ[n3+0] = (edgeJx1 + edgeJx2)*a;
      nodeJ[n3+1] = (edgeJy1 + edgeJy2)*a;
      nodeJ[n3+2] = (edgeJz1 + edgeJz2)*a;
   }
}

void calcNodeUe(Real* nodeRhoQi,Real* nodeJi,Real* nodeJ,Real* nodeUe,bool* innerFlag,Real* cellMaxUe,Simulation& sim,SimulationClasses& simClasses,pargrid::CellID blockID)
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
         cellMaxUe[n]++; // using cell array here to avoid introducing a new node array
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
      if(xMinFlag[n] == true) { continue; }
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
void getFields(Real* r,Real* B,Real* Ue,Simulation& sim,SimulationClasses& simClasses,pargrid::CellID blockID) {
   Real* faceB  = reinterpret_cast<Real*>(simClasses.pargrid.getUserData(Hybrid::dataFaceBID));
   //Real* nodeE  = reinterpret_cast<Real*>(simClasses.pargrid.getUserData(Hybrid::dataNodeEID));
   Real* cellUe = reinterpret_cast<Real*>(simClasses.pargrid.getUserData(Hybrid::dataCellUeID));
   face2r(r,faceB,sim,simClasses,blockID,B);
   //node2r(r,nodeE,sim,simClasses,blockID,E);
   cell2r(r,cellUe,sim,simClasses,blockID,Ue);
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
