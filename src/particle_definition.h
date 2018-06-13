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

#ifndef PARTICLE_DEFINITION_H
#define PARTICLE_DEFINITION_H

#include <mpi.h>

namespace particle {
//#ifdef USE_DETECTORS
//   enum STATE {X,Y,Z,VX,VY,VZ,WEIGHT,INI_CELLID,INI_TIME,INI_X,INI_Y,INI_Z,INI_VX,INI_VY,INI_VZ,SIZE};
//#else
   enum STATE {X,Y,Z,VX,VY,VZ,WEIGHT,SIZE};
//#endif
}

template<typename REAL>
struct Particle {
   REAL state[particle::SIZE];
   
   static void getDatatype(MPI_Datatype& datatype);
};

template<typename REAL> inline
void Particle<REAL>::getDatatype(MPI_Datatype& datatype) {
   MPI_Type_contiguous(particle::SIZE,MPI_Type<REAL>(),&datatype);
}

#endif
