#ifndef PARAMETERS_H
#define PARAMETERS_H

#include<string>
#include<sstream>
#include "H5Cpp.h"

using namespace std;

#define SYMMETRIC_PT_PTS 		0		// chooses whether or not to use gaussian spacing or symmetric spacing for pt points
#define UNIFORM_SPACING			false		// specifies uniform or non-uniform grid spacing for interpolation
#define ASSUME_ETA_SYMMETRIC 		1		// 1 means integrate only over eta_s = 0..eta_s_max, and multiply by 2 or 0 to speed up calculations
							// 0 means just integrate over given range of eta_s without worrying about symmetry
#define GROUPING_PARTICLES 		0		// set to 1 to perform calculations for similar particles together
#define PARTICLE_DIFF_TOLERANCE 	0.00		// particles with mass and chemical potential (for each FO-cell) difference less than this value
							// will be considered to be identical (b/c Cooper-Frye)
#define TRUNCATE_COOPER_FRYE		false		// ignore contributions to CF integral which are extremely small --> speeds up code by factor of 3-4
#define VERBOSE 			1		// specifies level of output - 0 is lowest (no output)
#define DEBUG				false		// flag for output of debugging statements
#define SPACETIME_MOMENTS_ONLY		true		// duh
#define DO_ALL_DECAY_CHANNELS		false		// duh
#define USE_HDF5			false		// utilizes HDF5 software to store large arrays
#define USE_LAMBDA			false		// fit correlation function with adjustable intercept parameter

#ifndef H5_NO_NAMESPACE
    using namespace H5;
#endif

//HDF information
const int RANK = 1;
const int ntrig = 2;			// for cos or sin

const double hbarC=0.197327053;		//GeV*fm
const double hbarC3=0.00768351405;
const double hbarCm1=5.067728853;
const double twopi = 2.*M_PI;
const double MeVToGeV = 0.001;

//particle information
const int Maxparticle=400;            //size of array for storage of the particles
const int Maxdecaychannel=13;
const int Maxdecaypart=5;

//spatial rapidity information
const int eta_s_npts = 15;
const double eta_s_i = 0.0;
const double eta_s_f = 4.0;

//relative momentum information
const int qonpts = 51;
const int qsnpts = 51;
const int qlnpts = 51;
//const int qnpts = 11;
//const double delta_q = 0.02;
//const double init_q = -5.0*delta_q;
//const int qnpts = 7;
//const double delta_q = 0.02;
//const double init_q = -3.0*delta_q;
//const int qnpts = 6;
//const double delta_q = 0.05;
//const double init_q = -2.5*delta_q;
//const int qnpts = 5;
//const double delta_q = 0.03;
//const double init_q = -2.0 * delta_q;
const int qnpts = 3;
const double delta_q = 0.06;
const double init_q = -delta_q;
//const int qnpts = 2;
//const double delta_q = 0.05;
//const double init_q = -0.5*delta_q;
//const int qnpts = 1;
//const double delta_q = 0.005;
//const double init_q = 0.0;

//all direction-specific q points information here
const int qtnpts = qnpts;
const int qxnpts = qnpts;
const int qynpts = qnpts;
const int qznpts = qnpts;
const double delta_qt = delta_q;
const double delta_qx = delta_q;
const double delta_qy = delta_q;
const double delta_qz = delta_q;
const double init_qt = init_q;
const double init_qx = init_q;
const double init_qy = init_q;
const double init_qz = init_q;
/*const int qtnpts = 7;
const int qxnpts = 13;
const int qynpts = 13;
const int qznpts = 11;
const double delta_qt = -0.02;
const double delta_qx = -0.02;
const double delta_qy = -0.02;
const double delta_qz = -0.02;
const double init_qt = -3.0*delta_qt;
const double init_qx = -6.0*delta_qt;
const double init_qy = -6.0*delta_qt;
const double init_qz = -5.0*delta_qt;*/

//single particle spectra info
const int n_SP_pT = 15;
const int n_SP_pphi = 48;
const double SP_pT_min = 0.0;
const double SP_pT_max = 3.0;

//parameters for interpolation grid
//  - polar
const int n_interp_pT_pts = 15;
const int n_interp_pphi_pts = 48;
const double interp_pT_min = 0.01;
const double interp_pphi_min = 0.0;
const double interp_pT_max = 4.05;
const double interp_pphi_max = 2.*M_PI;
const double Del2_pT = (interp_pT_max - interp_pT_min) / (double)(n_interp_pT_pts-1);
const double Del2_pphi = (interp_pphi_max - interp_pphi_min) / (double)(n_interp_pphi_pts-1);

//correlation function info
const int corrfuncdim = 1;
const bool lambdaflag = true;

//pair momentum info
const int n_localp_T = 14;
const double localp_T_min = 0.05;
const double localp_T_max = 0.7;
const int n_localp_phi = 51;
const double localp_phi_min = 0.0;
const double localp_phi_max = 2*M_PI;

const int n_order = 1;

const double tol = 0.0;		//tolerance
const int flagneg = 1;		//neglect all points that are negative
				//choose flagneg == 0 to agree with iS.e
				//choose flagneg == 1 for the real world

const size_t fit_max_iterations = 1000;  // stop at this point if not converged 
const double fit_tolerance = 1e-6;

//misc. resonance info
const double max_lifetime = 100.;	// fm/c

namespace patch
{
    template < typename T > std::string to_string( const T& n )
    {
        std::ostringstream stm ;
        stm << n ;
        return stm.str() ;
    }
}

#endif
