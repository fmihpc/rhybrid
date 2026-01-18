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

#ifndef NEUTRAL_PROFILES_H
#define NEUTRAL_PROFILES_H

#include <string>
#include <vector>
#include <simulation.h>
#include <simulationclasses.h>
#include "hybrid.h"

struct NeutralProfileArgs{
   Real m,r0,R_exobase,R_shadow;
   std::vector<Real> n0,H0,T0,k0;
};

Real neutralDensityChamberlainH(SimulationClasses& simClasses,Real x,Real y,Real z,Real r0,std::vector<Real> n0,std::vector<Real> H0,Real R_exobase,Real R_shadow) {
   const Real r = sqrt(sqr(x) + sqr(y) + sqr(z));
   if (r < R_exobase) { return 0.0; }
   if (x < 0 && sqr(y) + sqr(z) < sqr(R_shadow)) { return 0.0; }
   if (n0.size() < 1 || n0.size() != H0.size()) {
      simClasses.logger << "(neutralDensityChamberlainH) ERROR: Give same number of n0 and H0 parameters" << std::endl << write;
      return -1.0;
   }
   Real n = 0.0;
   for (size_t i=0;i<n0.size();i++) { n += n0[i]*exp(-H0[i]*(1/r0 - 1/r)); }
   return n;
}

Real neutralDensityChamberlainT(SimulationClasses& simClasses,Real x,Real y,Real z,Real m,Real r0,std::vector<Real> n0,std::vector<Real> T0,Real R_exobase,Real R_shadow) {
   const Real r = sqrt(sqr(x) + sqr(y) + sqr(z));
   if (r < R_exobase) { return 0.0; }
   if (x < 0 && sqr(y) + sqr(z) < sqr(R_shadow)) { return 0.0; }
   if (n0.size() < 1 || n0.size() != T0.size()) {
      simClasses.logger << "(neutralDensityChamberlainT) ERROR: Give same number of n0 and T0 parameters" << std::endl << write;
      return -1.0;
   }
   const Real a = constants::GRAVITY * Hybrid::M_object * m / constants::BOLTZMANN;
   Real n = 0.0;
   for (size_t i=0;i<n0.size();i++) { n += n0[i]*exp(-a/T0[i]*(1/r0 - 1/r)); }
   return n;
}

Real neutralDensityExponential(SimulationClasses& simClasses,Real x,Real y,Real z,Real r0,std::vector<Real> n0,std::vector<Real> H0,Real R_exobase,Real R_shadow) {
   const Real r = sqrt(sqr(x) + sqr(y) + sqr(z));
   if (r < R_exobase) { return 0.0; }
   if (x < 0 && sqr(y) + sqr(z) < sqr(R_shadow)) { return 0.0; }
   if (n0.size() < 1 || n0.size() != H0.size()) {
      simClasses.logger << "(neutralDensityExponential) ERROR: Give same number of n0 and H0 parameters" << std::endl << write;
      return -1.0;
   }
   Real n = 0.0;
   for (size_t i=0;i<n0.size();i++) { n += n0[i]*exp(-(r-r0)/H0[i]); }
   return n;
}

Real neutralDensityPowerLaw(SimulationClasses& simClasses,Real x,Real y,Real z,Real r0,std::vector<Real> n0,std::vector<Real> k0,Real R_exobase,Real R_shadow) {
   const Real r = sqrt(sqr(x) + sqr(y) + sqr(z));
   if (r < R_exobase) { return 0.0; }
   if (x < 0 && sqr(y) + sqr(z) < sqr(R_shadow)) { return 0.0; }
   if (n0.size() < 1 || n0.size() != k0.size()) {
      simClasses.logger << "(neutralDensityPowerLaw) ERROR: Give same number of n0 and k0 parameters" << std::endl << write;
      return -1.0;
   }
   Real n = 0.0;
   for (size_t i=0;i<n0.size();i++) { n += n0[i]*pow(r0/r,k0[i]); }
   return n;
}

Real neutralDensityVenusHydrogen(SimulationClasses& simClasses,Real x,Real y,Real z,Real R_exobase,Real R_shadow) {
   const Real r = sqrt(sqr(x) + sqr(y) + sqr(z));
   if (r < R_exobase) { return 0.0; }
   if (x < 0 && sqr(y) + sqr(z) < sqr(R_shadow)) { return 0.0; }
   const Real sza = acos(x/r);
   const Real exoR = constants::DIST_VENUS_RADIUS + 170e3;
   // thermal hydrogen
   const Real n_0_zip1 = 1.32e11;
   const Real beta_zip1 = constants::GRAVITY * constants::MASS_VENUS * constants::MASS_PROTON/(constants::BOLTZMANN*285);
   const Real n_0_zip5 = 2.59e12;
   const Real beta_zip5 = constants::GRAVITY * constants::MASS_VENUS * constants::MASS_PROTON/(constants::BOLTZMANN*110);
   const Real n_thermal_day = n_0_zip1*exp(-beta_zip1*(1/exoR - 1.0/r));
   const Real n_thermal_night = n_0_zip5*exp(-beta_zip5*(1.0/exoR - 1.0/r));
   const Real n_thermal = (1.0 - (sza/M_PI))*n_thermal_day + (sza/M_PI)*n_thermal_night;
   // hot hydrogen
   const Real a_1_noon = -6.2625e-5/1e3;
   const Real a_2_noon = 15.4817;
   const Real a_3_noon = 3.6414e4*1e3;
   const Real a_1_term = -8.4607e-5/1e3;
   const Real a_2_term = 15.9944;
   const Real a_3_term = 2.9743e4*1e3;
   const Real a_1_midn = -6.2309e-5/1e3;
   const Real a_2_midn = 15.2723;
   const Real a_3_midn = 4.3781e4*1e3;
   const Real n_hot_noon = exp(a_1_noon*r + a_2_noon + a_3_noon/r);
   const Real n_hot_term = exp(a_1_term*r + a_2_term + a_3_term/r);
   const Real n_hot_midn = exp(a_1_midn*r + a_2_midn + a_3_midn/r);
   Real n_hot;
   if (sza <= M_PI/2) {
      n_hot = (1.0 - (sza/(M_PI/2)))*n_hot_noon + (sza/(M_PI/2))*n_hot_term;
   }
   else {
      n_hot = (1.0 - (sza/(M_PI/2) - 1))*n_hot_term + (sza/(M_PI/2) - 1)*n_hot_midn;
   }
   return n_thermal + n_hot;
}

Real neutralDensityVenusOxygen(SimulationClasses& simClasses,Real x,Real y,Real z,Real R_exobase,Real R_shadow) {
   const Real r = sqrt(sqr(x) + sqr(y) + sqr(z));
   if (r < R_exobase) { return 0.0; }
   if (x < 0 && sqr(y) + sqr(z) < sqr(R_shadow)) { return 0.0; }
   const Real sza = acos(x/r);
   const Real r_0_day = constants::DIST_VENUS_RADIUS + 200e3;
   const Real beta_day = constants::GRAVITY * constants::MASS_VENUS * constants::MASS_OXYGEN/(constants::BOLTZMANN*6400);
   const Real n_0_day = 7.5e10;
   const Real r_0_night = constants::DIST_VENUS_RADIUS + 300e3;
   const Real beta_night = constants::GRAVITY * constants::MASS_VENUS * constants::MASS_OXYGEN/(constants::BOLTZMANN*4847);
   const Real n_0_night = 2e9;
   Real n_day;
   Real n_night;
   if (r < r_0_day) {
      n_day = 0;
   }
   else {
      n_day = n_0_day*exp(-beta_day*(1.0/r_0_day - 1.0/r));
   }
   if (r < r_0_night) {
      n_night = 0;
   }
   else {
      n_night = n_0_night*exp(-beta_night*(1.0/r_0_night - 1.0/r));
   }
   return (1.0 - (sza/M_PI))*n_day + (sza/M_PI)*n_night;
}

// for neutralDensityMercurySodiumExner20
Real expexp(Real x,Real y,Real z,Real r,Real Rp,Real n0,Real Hr,Real Htheta,Real lat_nmax,Real lon_nmax) {
   const Real lat = atan2(z,sqrt(sqr(x) + sqr(y))); // latitude in radians
   const Real lon = atan2(y,x); // longitude in radians
   const Real dtheta = acos(sin(lat)*sin(lat_nmax) + cos(lat)*cos(lat_nmax)*cos(lon_nmax - lon));
   const Real res = n0 * exp(-(r-Rp)/Hr) * exp(-sqr(dtheta/Htheta));
   return res;
}

// neutral profile from: Exner+ (2020), Inﬂuence of Mercury's exosphere on the structure of the magnetosphere, J. Geophys. Res., 125, e2019JA027691, doi:10.1029/2019JA027691
Real neutralDensityMercurySodiumExner20(SimulationClasses& simClasses,Real x,Real y,Real z) {
   const Real Rp = constants::DIST_MERCURY_RADIUS;
   const Real r = sqrt(sqr(x) + sqr(y) + sqr(z));
   if (r < Rp) { return 0.0; }
   // surface density [m-3]
   const Real n0_TD  = 8.86e9;
   const Real n0_MIV = 7.84e6;
   const Real n0_PSD = 4.06e10;
   const Real n0_SP  = 5.67e6;

   // ccale height [m]
   const Real Hr_TD  = 100e3;
   const Real Hr_MIV = 431e3;
   const Real Hr_PSD = 232e3;
   const Real Hr_SP  = 748e3;

   // angular width
   const Real Htheta_TD = 15*M_PI/180.0;
   const Real Htheta_PSD = 20*M_PI/180.0;
   const Real Htheta_SP_day = 15*M_PI/180.0;
   const Real Htheta_SP_nigth = 10*M_PI/180.0;

   // lat and lon of maximum density
   const Real lat_nmax_TD = 0*M_PI/180.0;
   const Real lon_nmax_TD = 0*M_PI/180.0;

   const Real lat_nmax_PSD1 = 50*M_PI/180.0;
   const Real lon_nmax_PSD1 = 0*M_PI/180.0;

   const Real lat_nmax_PSD2 = -50*M_PI/180.0;
   const Real lon_nmax_PSD2 = 0*M_PI/180.0;

   const Real lat_nmax_SP_day1 = 80*M_PI/180.0;
   const Real lon_nmax_SP_day1 = 0*M_PI/180.0;

   const Real lat_nmax_SP_day2 = -80*M_PI/180.0;
   const Real lon_nmax_SP_day2 = 0*M_PI/180.0;

   const Real lat_nmax_SP_nigth1 = 15*M_PI/180.0;
   const Real lon_nmax_SP_nigth1 = 180*M_PI/180.0;

   const Real lat_nmax_SP_nigth2 = -15*M_PI/180.0;
   const Real lon_nmax_SP_nigth2 = 180*M_PI/180.0;

   // MIV
   const Real n_miv = n0_MIV * exp(-(r-Rp)/Hr_MIV);

   // TD
   const Real n_td = expexp(x,y,z,r,Rp,n0_TD,Hr_TD,Htheta_TD,lat_nmax_TD,lon_nmax_TD);

   // SP day 1 & 2
   const Real n_sp_day1 = expexp(x,y,z,r,Rp,n0_SP,Hr_SP,Htheta_SP_day,lat_nmax_SP_day1,lon_nmax_SP_day1);
   const Real n_sp_day2 = expexp(x,y,z,r,Rp,n0_SP,Hr_SP,Htheta_SP_day,lat_nmax_SP_day2,lon_nmax_SP_day2);

   // SP nigth 1 & 2
   const Real n_sp_nigth1 = expexp(x,y,z,r,Rp,n0_SP,Hr_SP,Htheta_SP_nigth,lat_nmax_SP_nigth1,lon_nmax_SP_nigth1);
   const Real n_sp_nigth2 = expexp(x,y,z,r,Rp,n0_SP,Hr_SP,Htheta_SP_nigth,lat_nmax_SP_nigth2,lon_nmax_SP_nigth2);

   // PSD 1 & 2
   const Real n_psd1 = expexp(x,y,z,r,Rp,n0_PSD,Hr_PSD,Htheta_PSD,lat_nmax_PSD1,lon_nmax_PSD1);
   const Real n_psd2 = expexp(x,y,z,r,Rp,n0_PSD,Hr_PSD,Htheta_PSD,lat_nmax_PSD2,lon_nmax_PSD2);

   Real n = sqr(Rp/r)*(n_miv + n_td + n_sp_day1 + n_sp_day2 + n_sp_nigth1 + n_sp_nigth2 + n_psd1 + n_psd2);
   // ionization is only 1/5 in the shadow
   if ( x < 0 && ( (sqr(y) + sqr(z)) < sqr(Rp) ) ) { n *= 0.2; }
   return n;
}

Real getNeutralDensity(SimulationClasses& simClasses,std::string name,Real x,Real y,Real z,NeutralProfileArgs a) {
   if (name.compare("ChamberlainH") == 0) {
      return neutralDensityChamberlainH(simClasses,x,y,z,a.r0,a.n0,a.H0,a.R_exobase,a.R_shadow);
   }
   else if (name.compare("ChamberlainT") == 0) {
      return neutralDensityChamberlainT(simClasses,x,y,z,a.m,a.r0,a.n0,a.T0,a.R_exobase,a.R_shadow);
   }
   else if (name.compare("Exponential") == 0) {
      return neutralDensityExponential(simClasses,x,y,z,a.r0,a.n0,a.H0,a.R_exobase,a.R_shadow);
   }
   else if (name.compare("PowerLaw") == 0) {
      return neutralDensityPowerLaw(simClasses,x,y,z,a.r0,a.n0,a.k0,a.R_exobase,a.R_shadow);
   }
   else if (name.compare("VenusHydrogen") == 0) {
      return neutralDensityVenusHydrogen(simClasses,x,y,z,a.R_exobase,a.R_shadow);
   }
   else if (name.compare("VenusOxygen") == 0) {
      return neutralDensityVenusOxygen(simClasses,x,y,z,a.R_exobase,a.R_shadow);
   }
   else if (name.compare("MercurySodiumExner20") == 0) {
      return neutralDensityMercurySodiumExner20(simClasses,x,y,z);
   }
   else {
      simClasses.logger << "(getNeutralDensity) ERROR: unknown name of an exospheric neutral density profile (" << name << ")" << std::endl << write;
      return -1.0;
   }
   return -1.0;
}

#endif
