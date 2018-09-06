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

#ifndef DETECTORS_H
#define DETECTORS_H
#ifdef USE_DETECTORS

#include <particle_list_base.h>
#include <simulation.h>
#include <simulationclasses.h>

bool writeDetectorParticle(Simulation& sim,SimulationClasses& simClasses);
bool recordDetectorBulkParam(Simulation& sim,SimulationClasses& simClasses);
bool writeDetectorBulkParam(Simulation& sim,SimulationClasses& simClasses);

#endif
#endif
