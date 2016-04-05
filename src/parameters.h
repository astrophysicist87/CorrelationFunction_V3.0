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
#define VERBOSE 			1		// specifies level of output - 0 is lowest (no output)
#define DEBUG				false		// flag for output of debugging statements
#define SPACETIME_MOMENTS_ONLY		false		// duh
#define DO_ALL_DECAY_CHANNELS		false		// duh
#define UNZIP_HDF5			false		// utilizes HDF5 software to store large arrays
#define USE_LAMBDA			true		// fit correlation function with adjustable intercept parameter
#define USE_EXTRAPOLATION		true		// extrapolates results of CF integrals instead of competing them, false just calculates full integrals (slower)
#define EXTRAPOLATION_METHOD		0		// 0 - GSL polynomial fit
							// 1 - direct calculation of rational function fit using ratint in Arsenal.* files (numerator and denominator
							// orders chosen automatically to be n+m+1==number of percentage markers)
							// 2 - best fit rational function of given orders in numerator and denominator
#define PC_MARKER_SPACING		1		// 0 - automatic
							// 1 - use usr_def_pc_markers
							// 2 - use usr_def_pc_markers_thinned
#define COMPUTE_RESONANCE_ARRAYS	true		// alternative is to read them in from a file
#define COMPUTE_RESONANCE_DECAYS	true		// alternative is to read them in from a file
#define IGNORE_LONG_LIVED_RESONANCES	true		// particularly, whether to include eta or eta' in spectra calculations
							// true means C(q=0) ~ 1 + \lambda
#define QT_POINTS_SPACING		1		// 0 - uniform from -qmax to +qmax
							// 1 - Chebyshev nodes from -qmax to +qmax
#define QX_POINTS_SPACING		0
#define QY_POINTS_SPACING		0
#define QZ_POINTS_SPACING		0
//#define VARY_ALPHA			false		// (not yet implemented) feature to treat power in exponential as a fit variable (alpha == 2 <==> traditional Gaussian)
#define Q_AXES_AND_RAYS_ONLY		false		// true - only do points along q-axes (only works for odd points right now)
							// false - do full grid
#define FIT_WITH_PROJECTED_CFVALS	true		// as opposed to unprojected CFvals...
#define FLESH_OUT_CF			false		// refines grid via interpolation before fitting
#define REGULATE_CF			true		// true (false) means (don't) try to catch spurious values of projected
							// or regular CF and replace them with median value in that window

#ifndef H5_NO_NAMESPACE
    using namespace H5;
#endif

//HDF information
const int RANK = 1;
const int ntrig = 2;			// for cos or sin

const double hbarC = 0.197327053;		//GeV*fm
const double hbarC3 = 0.00768351405;
const double hbarCm1 = 5.067728853;
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

//extrapolation information
const int polynomial_fit_order = 4;
const int rational_function_numerator_order = 3;
const int rational_function_denominator_order = 4;
const int UDPMsize = 15;
const int UDPMTsize = 10;
static double usr_def_pc_markers[UDPMsize] = {
					0.00, 0.85, 0.86, 0.87, 0.88, 0.89, 0.90,
					0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98
				};
static double usr_def_pc_markers_thinned[UDPMTsize] = { 0.00, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90 };

//relative momentum information
const int qonpts = 11;
const int qsnpts = 11;
const int qlnpts = 11;
const int qnpts = 1;
const double delta_q = 0.005;
const double init_q = 0.0;

const int new_nqpts = 51;

//all direction-specific q points information here
const int qtnpts = 3;
const int qxnpts = 3;
const int qynpts = 3;
const int qznpts = 3;
const double delta_qt = 0.02;
const double delta_qx = 0.0016;
const double delta_qy = 0.02;
const double delta_qz = 0.02;
const double init_qt = -0.5*double(qtnpts-1)*delta_qt;
const double init_qx = -0.5*double(qxnpts-1)*delta_qx;
const double init_qy = -0.5*double(qynpts-1)*delta_qy;
const double init_qz = -0.5*double(qznpts-1)*delta_qz;

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
const double correlator_minus_one_cutoff = 0.0;		//zero means all calculations happen as usual

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
