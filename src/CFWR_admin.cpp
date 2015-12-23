#include<iostream>
#include<sstream>
#include<string>
#include<fstream>
#include<cmath>
#include<iomanip>
#include<vector>
#include<stdio.h>
#include<time.h>
#include<algorithm>
#include <set>

#include<gsl/gsl_sf_bessel.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>            // gsl random number generators
#include <gsl/gsl_randist.h>        // gsl random number distributions
#include <gsl/gsl_vector.h>         // gsl vector and matrix definitions
#include <gsl/gsl_blas.h>           // gsl linear algebra stuff
#include <gsl/gsl_multifit_nlin.h>  // gsl multidimensional fitting

#include "CFWR.h"
#include "Arsenal.h"
#include "gauss_quadrature.h"

using namespace std;

const double PTCHANGE = 1.0;

template <typename T> int sgn(T val)
{
    return (T(0) < val) - (val < T(0));
}

CorrelationFunction::CorrelationFunction(particle_info* particle, particle_info* all_particles_in, int Nparticle_in,
					FO_surf* FOsurf_ptr, vector<int> chosen_resonances_in, int particle_idx, ofstream& myout)
{
	//set ofstream for output file
	global_out_stream_ptr = &myout;
	
	//particle information (both final-state particle used in HBT and all decay decay_channels)
	particle_name = particle->name;
	particle_mass = particle->mass;
	particle_sign = particle->sign;
	particle_gspin = particle->gspin;
	particle_id = particle_idx;
	target_particle_id = particle_id;
	particle_monval = particle->monval;
	S_prefactor = 1.0/(8.0*(M_PI*M_PI*M_PI))/hbarC/hbarC/hbarC;
	all_particles = all_particles_in;
	for (int icr = 0; icr < (int)chosen_resonances_in.size(); icr++)
		chosen_resonances.push_back(chosen_resonances_in[icr]);
	thermal_pions_only = false;
	Nparticle = Nparticle_in;
	read_in_all_dN_dypTdpTdphi = false;
	output_all_dN_dypTdpTdphi = true;
	currentfolderindex = -1;
	current_level_of_output = 0;
	qspace_cs_slice_length = qnpts*qnpts*qnpts*qnpts*2;		//factor of 2 for sin or cos

	Set_q_points();

	n_zeta_pts = 12;
	n_v_pts = 12;
	n_s_pts = 12;
	v_min = -1.;
	v_max = 1.;
	zeta_min = 0.;
	zeta_max = M_PI;
	
	//default: use delta_f in calculations
	use_delta_f = true;
	no_df_stem = "";
	current_FOsurf_ptr = FOsurf_ptr;

//****************************************************************************************************
//OLD CODE FOR READING IN SELECTED decay_channels...
	current_resonance_mass = 0.0;
	current_resonance_mu = 0.0;
	current_resonance_Gamma = 0.0;
	current_resonance_total_br = 0.0;
	current_resonance_decay_masses = new double [2];
	current_resonance_decay_masses[0] = 0.0;
	current_resonance_decay_masses[1] = 0.0;
	previous_resonance_particle_id = -1;
	previous_decay_channel_idx = -1;				//different for each decay channel
	previous_resonance_mass = 0.0;
	previous_resonance_Gamma = 0.0;
	previous_resonance_total_br = 0.0;
	if (chosen_resonances.size() == 0)
	{
		n_decay_channels = 1;
		n_resonance = 0;
		thermal_pions_only = true;
		if (VERBOSE > 0) *global_out_stream_ptr << "Thermal pion(+) only!" << endl;
		decay_channels = new decay_info [n_decay_channels];
		decay_channels[0].resonance_decay_masses = new double [Maxdecaypart];	// Maxdecaypart == 5
debugger(__LINE__, __FILE__);
	}
	else
	{
		//n_decay_channels is actually total number of decay channels which can generate pions
		//from chosen decay_channels
		n_decay_channels = get_number_of_decay_channels(chosen_resonances, all_particles);
		n_resonance = (int)chosen_resonances.size();
		if (VERBOSE > 0) *global_out_stream_ptr << "Computed n_decay_channels = " << n_decay_channels << endl
							<< "Computed n_resonance = " << n_resonance << endl;
		decay_channels = new decay_info [n_decay_channels];
		int temp_idx = 0;
		for (int icr = 0; icr < n_resonance; icr++)
		{
			particle_info particle_temp = all_particles[chosen_resonances[icr]];
			if (VERBOSE > 0) *global_out_stream_ptr << "Loading resonance: " << particle_temp.name
					<< ", chosen_resonances[" << icr << "] = " << chosen_resonances[icr] << endl;
			for (int idecay = 0; idecay < particle_temp.decays; idecay++)
			{
				if (VERBOSE > 0) *global_out_stream_ptr << "Current temp_idx = " << temp_idx << endl;
				if (temp_idx == n_decay_channels)	//i.e., all contributing decay channels have been loaded
					break;
				decay_channels[temp_idx].resonance_name = particle_temp.name;		// set name of resonance

				//check if effective branching is too small for inclusion in source variances
				bool effective_br_is_too_small = false;
				if (particle_temp.decays_effective_branchratio[idecay] <= 1.e-12)
					effective_br_is_too_small = true;

				decay_channels[temp_idx].resonance_particle_id = chosen_resonances[icr];	// set index of resonance in all_particles
				decay_channels[temp_idx].resonance_idx = icr;					// set index of resonance in chosen_resonances
				decay_channels[temp_idx].resonance_decay_masses = new double [Maxdecaypart];	// Maxdecaypart == 5
				decay_channels[temp_idx].resonance_decay_monvals = new int [Maxdecaypart];	// Maxdecaypart == 5
				decay_channels[temp_idx].resonance_decay_Gammas = new double [Maxdecaypart];	// Maxdecaypart == 5



				//*** SETTING RESONANCE DECAY MASSES DIFFERENTLY FOR NEW ANALYZE SF
				for (int ii = 0; ii < Maxdecaypart; ii++)
				{
					decay_channels[temp_idx].resonance_decay_monvals[ii] = particle_temp.decays_part[idecay][ii];
					if (particle_temp.decays_part[idecay][ii] == 0)
					{
						decay_channels[temp_idx].resonance_decay_masses[ii] = 0.0;
						decay_channels[temp_idx].resonance_decay_Gammas[ii] = 0.0;

					}
					else
					{
						int tempID = lookup_particle_id_from_monval(all_particles, Nparticle, particle_temp.decays_part[idecay][ii]);
						decay_channels[temp_idx].resonance_decay_masses[ii] = all_particles[tempID].mass;
						decay_channels[temp_idx].resonance_decay_Gammas[ii] = all_particles[tempID].width;

					}
				}
				decay_channels[temp_idx].resonance_mu = particle_temp.mu;
				decay_channels[temp_idx].resonance_gspin = particle_temp.gspin;
				decay_channels[temp_idx].resonance_sign = particle_temp.sign;
				decay_channels[temp_idx].resonance_mass = particle_temp.mass;
				decay_channels[temp_idx].nbody = abs(particle_temp.decays_Npart[idecay]);
				decay_channels[temp_idx].resonance_Gamma = particle_temp.width;
				decay_channels[temp_idx].resonance_total_br = particle_temp.decays_effective_branchratio[idecay];
				decay_channels[temp_idx].resonance_direct_br = particle_temp.decays_branchratio[idecay];

				
				//check if particle lifetime is too long for inclusion in source variances
				bool lifetime_is_too_long = false;
				//if (decay_channels[temp_idx].resonance_Gamma < hbarC / max_lifetime)
				//	lifetime_is_too_long = true;		//i.e., for lifetimes longer than 100 fm/c, skip decay channel

				if (VERBOSE > 0) *global_out_stream_ptr << "Resonance = " << decay_channels[temp_idx].resonance_name << ", decay channel " << idecay + 1
						<< ": mu=" << decay_channels[temp_idx].resonance_mu
						<< ", gs=" << decay_channels[temp_idx].resonance_gspin << ", sign=" << decay_channels[temp_idx].resonance_sign
						<< ", M=" << decay_channels[temp_idx].resonance_mass << ", nbody=" << decay_channels[temp_idx].nbody
						<< ", Gamma=" << decay_channels[temp_idx].resonance_Gamma << ", total br=" << decay_channels[temp_idx].resonance_total_br
						<< ", direct br=" << decay_channels[temp_idx].resonance_direct_br << endl;

				if (VERBOSE > 0) *global_out_stream_ptr << "Resonance = " << decay_channels[temp_idx].resonance_name << ": ";
				for (int decay_part_idx = 0; decay_part_idx < decay_channels[temp_idx].nbody; decay_part_idx++)
					if (VERBOSE > 0) *global_out_stream_ptr << "m[" << decay_part_idx << "] = "
						<< decay_channels[temp_idx].resonance_decay_masses[decay_part_idx] << "   "
						<< decay_channels[temp_idx].resonance_decay_monvals[decay_part_idx] << "   ";
				if (VERBOSE > 0) *global_out_stream_ptr << endl << endl;

				// if decay channel parent resonance is not too long-lived
				// and decay channel contains at least one target daughter particle,
				// include channel
				if (DO_ALL_DECAY_CHANNELS)
					decay_channels[temp_idx].include_channel = true;
				else
					decay_channels[temp_idx].include_channel = (!lifetime_is_too_long
																	&& !effective_br_is_too_small);

				temp_idx++;
			}
		}
	}

	current_dN_dypTdpTdphi_moments = new double ****** [n_interp_pT_pts];
	current_ln_dN_dypTdpTdphi_moments = new double ****** [n_interp_pT_pts];
	current_sign_of_dN_dypTdpTdphi_moments = new double ****** [n_interp_pT_pts];
	for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
	{
		current_dN_dypTdpTdphi_moments[ipt] = new double ***** [n_interp_pphi_pts];
		current_ln_dN_dypTdpTdphi_moments[ipt] = new double ***** [n_interp_pphi_pts];
		current_sign_of_dN_dypTdpTdphi_moments[ipt] = new double ***** [n_interp_pphi_pts];
		for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
		{
			current_dN_dypTdpTdphi_moments[ipt][ipphi] = new double **** [qnpts];
			current_ln_dN_dypTdpTdphi_moments[ipt][ipphi] = new double **** [qnpts];
			current_sign_of_dN_dypTdpTdphi_moments[ipt][ipphi] = new double **** [qnpts];
			for (int iqt = 0; iqt < qnpts; ++iqt)
			{
				current_dN_dypTdpTdphi_moments[ipt][ipphi][iqt] = new double *** [qnpts];
				current_ln_dN_dypTdpTdphi_moments[ipt][ipphi][iqt] = new double *** [qnpts];
				current_sign_of_dN_dypTdpTdphi_moments[ipt][ipphi][iqt] = new double *** [qnpts];
				for (int iqx = 0; iqx < qnpts; ++iqx)
				{
					current_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx] = new double ** [qnpts];
					current_ln_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx] = new double ** [qnpts];
					current_sign_of_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx] = new double ** [qnpts];
					for (int iqy = 0; iqy < qnpts; ++iqy)
					{
						current_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx][iqy] = new double * [qnpts];
						current_ln_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx][iqy] = new double * [qnpts];
						current_sign_of_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx][iqy] = new double * [qnpts];
						for (int iqz = 0; iqz < qnpts; ++iqz)
						{
							current_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx][iqy][iqz] = new double [2];
							current_ln_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx][iqy][iqz] = new double [2];
							current_sign_of_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx][iqy][iqz] = new double [2];
							for (int itrig = 0; itrig < 2; ++itrig)
							{
								current_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx][iqy][iqz][itrig] = 0.0;
								current_ln_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx][iqy][iqz][itrig] = 0.0;
								current_sign_of_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx][iqy][iqz][itrig] = 0.0;
							}
						}
					}
				}
			}
		}
	}

	res_log_info = new double ** [n_interp_pT_pts];
	res_sign_info = new double ** [n_interp_pT_pts];
	res_moments_info = new double ** [n_interp_pT_pts];
	for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
	{
		res_log_info[ipt] = new double * [n_interp_pphi_pts];
		res_sign_info[ipt] = new double * [n_interp_pphi_pts];
		res_moments_info[ipt] = new double * [n_interp_pphi_pts];
		for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
		{
			res_log_info[ipt][ipphi] = new double [qspace_cs_slice_length];
			res_sign_info[ipt][ipphi] = new double [qspace_cs_slice_length];
			res_moments_info[ipt][ipphi] = new double [qspace_cs_slice_length];
		}
	}

	// initialize spectra and correlation function arrays
	spectra = new double ** [Nparticle];
	abs_spectra = new double ** [Nparticle];
	for (int ir = 0; ir < Nparticle; ++ir)
	{
		spectra[ir] = new double * [n_interp_pT_pts];
		abs_spectra[ir] = new double * [n_interp_pT_pts];
		for (int ipT = 0; ipT < n_interp_pT_pts; ++ipT)
		{
			spectra[ir][ipT] = new double [n_interp_pphi_pts];
			abs_spectra[ir][ipT] = new double [n_interp_pphi_pts];
			for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
			{
				spectra[ir][ipT][ipphi] = 0.0;
				abs_spectra[ir][ipT][ipphi] = 0.0;
			}
		}
	}

	CFvals = new double **** [n_interp_pT_pts];
	for (int ipT = 0; ipT < n_interp_pT_pts; ++ipT)
	{
		CFvals[ipT] = new double *** [n_interp_pphi_pts];
		for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
		{
			CFvals[ipT][ipphi] = new double ** [qonpts];
			for (int iqo = 0; iqo < qonpts; ++iqo)
			{
				CFvals[ipT][ipphi][iqo] = new double * [qsnpts];
				for (int iqs = 0; iqs < qsnpts; ++iqs)
				{
					CFvals[ipT][ipphi][iqo][iqs] = new double [qlnpts];
					for (int iql = 0; iql < qlnpts; ++iql)
						CFvals[ipT][ipphi][iqo][iqs][iql] = 0.0;
				}
			}
		}
	}

	// used for keeping track of how many FO cells are important for given pT, pphi
	number_of_FOcells_above_cutoff_array = new int * [n_interp_pT_pts];
	for (int ipT = 0; ipT < n_interp_pT_pts; ++ipT)
	{
		number_of_FOcells_above_cutoff_array[ipT] = new int [n_interp_pphi_pts];
		for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
			number_of_FOcells_above_cutoff_array[ipT][ipphi] = 0;
	}


	// set-up integration points for resonance integrals
	s_pts = new double * [n_decay_channels];
	s_wts = new double * [n_decay_channels];
	v_pts = new double [n_v_pts];
	v_wts = new double [n_v_pts];
	zeta_pts = new double [n_zeta_pts];
	zeta_wts = new double [n_zeta_pts];
	for (int idc=0; idc<n_decay_channels; idc++)
	{
		s_pts[idc] = new double [n_s_pts];
		s_wts[idc] = new double [n_s_pts];
	}
	//initialize all gaussian points for resonance integrals
	//syntax: int gauss_quadrature(int order, int kind, double alpha, double beta, double a, double b, double x[], double w[])
	gauss_quadrature(n_zeta_pts, 1, 0.0, 0.0, zeta_min, zeta_max, zeta_pts, zeta_wts);
	gauss_quadrature(n_v_pts, 1, 0.0, 0.0, v_min, v_max, v_pts, v_wts);
	for (int idc = 0; idc < n_decay_channels; idc++)
	{
		//cerr << "working on resonance #" << idc << "..." << endl;
		double Gamma_temp = decay_channels[idc].resonance_Gamma;
		double m2_temp = decay_channels[idc].resonance_decay_masses[0];
		double m3_temp = decay_channels[idc].resonance_decay_masses[1];
		double M_temp = decay_channels[idc].resonance_mass;
		double s_min_temp = (m2_temp + m3_temp)*(m2_temp + m3_temp);
		double s_max_temp = (M_temp - particle_mass)*(M_temp - particle_mass);
		// N.B. - this is only really necessary for 3-body decays,
		//			but doesn't cause any problems for 2-body and is easier/simpler to code...
		gauss_quadrature(n_s_pts, 1, 0.0, 0.0, s_min_temp, s_max_temp, s_pts[idc], s_wts[idc]);
	}	
   //single particle spectra for plane angle determination
   SP_pT = new double [n_SP_pT];
   SP_pT_weight = new double [n_SP_pT];
   gauss_quadrature(n_SP_pT, 1, 0.0, 0.0, SP_pT_min, SP_pT_max, SP_pT, SP_pT_weight);
   SP_pphi = new double [n_SP_pphi];
   SP_pphi_weight = new double [n_SP_pphi];
   gauss_quadrature(n_SP_pphi, 1, 0.0, 0.0, 0.0, 2.*M_PI, SP_pphi, SP_pphi_weight);
   SP_p_y = 0.0e0;
//initialize and set evenly spaced grid of px-py points in transverse plane,
//and corresponding p0 and pz points
	SPinterp_pT = new double [n_interp_pT_pts];
	SPinterp_pT_public = new double [n_interp_pT_pts];
	SPinterp_pphi = new double [n_interp_pphi_pts];
	sin_SPinterp_pphi = new double [n_interp_pphi_pts];
	cos_SPinterp_pphi = new double [n_interp_pphi_pts];
	SPinterp_p0 = new double* [n_interp_pT_pts];
	SPinterp_pz = new double* [n_interp_pT_pts];
	for(int ipt=0; ipt<n_interp_pT_pts; ipt++)
	{
		SPinterp_p0[ipt] = new double [eta_s_npts];
		SPinterp_pz[ipt] = new double [eta_s_npts];
	}
	SPinterp_p0_Transpose = new double * [eta_s_npts];
	SPinterp_pz_Transpose = new double * [eta_s_npts];
	for (int ieta = 0; ieta < eta_s_npts; ieta++)
	{
		SPinterp_p0_Transpose[ieta] = new double [n_interp_pT_pts];
		SPinterp_pz_Transpose[ieta] = new double [n_interp_pT_pts];
	}
	double * dummywts3 = new double [n_interp_pT_pts];
	double * dummywts4 = new double [n_interp_pphi_pts];

	if (UNIFORM_SPACING)
	{
		//use uniformly spaced points in transverse momentum to make
		//interpolation simpler/faster
		for(int ipt=0; ipt<n_interp_pT_pts; ipt++)
		{
			SPinterp_pT[ipt] = interp_pT_min + (double)ipt*Del2_pT;
			SPinterp_pT_public[ipt] = SPinterp_pT[ipt];
		}
		for(int ipphi=0; ipphi<n_interp_pphi_pts; ipphi++)
		{
			SPinterp_pphi[ipphi] = interp_pphi_min + (double)ipphi*Del2_pphi;
			sin_SPinterp_pphi[ipphi] = sin(SPinterp_pphi[ipphi]);
			cos_SPinterp_pphi[ipphi] = cos(SPinterp_pphi[ipphi]);
		}
	}
	else
	{
		//try just gaussian points...
		//syntax:
		//int gauss_quadrature(int order, int kind, double alpha, double beta, double a, double b, double x[], double w[]) 
		gauss_quadrature(n_interp_pT_pts, 5, 0.0, 0.0, 0.0, 13.0, SPinterp_pT, dummywts3);	//use this one to agree with iS.e
		gauss_quadrature(n_interp_pT_pts, 5, 0.0, 0.0, 0.0, 13.0, SPinterp_pT_public, dummywts3);	//use this one to agree with iS.e
		gauss_quadrature(n_interp_pphi_pts, 1, 0.0, 0.0, interp_pphi_min, interp_pphi_max, SPinterp_pphi, dummywts4);
		for(int ipphi=0; ipphi<n_interp_pphi_pts; ipphi++)
		{
			//SPinterp_pphi[ipphi] = interp_pphi_min + (double)ipphi*Del2_pphi;
			sin_SPinterp_pphi[ipphi] = sin(SPinterp_pphi[ipphi]);
			cos_SPinterp_pphi[ipphi] = cos(SPinterp_pphi[ipphi]);
		}
	}

debugger(__LINE__, __FILE__);
	dN_dypTdpTdphi = new double* [n_SP_pT];
	cosine_iorder = new double* [n_SP_pT];
	sine_iorder = new double* [n_SP_pT];
	for(int i=0; i<n_SP_pT; i++)
	{
		dN_dypTdpTdphi[i] = new double [n_SP_pphi];
		cosine_iorder[i] = new double [n_order];
		sine_iorder[i] = new double [n_order];
	}
	dN_dydphi = new double [n_SP_pphi];
	dN_dypTdpT = new double [n_SP_pT];
	pTdN_dydphi = new double [n_SP_pphi];
	for(int i=0; i<n_SP_pphi; i++)
	{
		dN_dydphi[i] = 0.0e0;
		pTdN_dydphi[i] = 0.0e0;
		for(int j=0; j<n_SP_pT; j++)
			dN_dypTdpTdphi[j][i] = 0.0e0;
	}
	for(int i=0; i<n_SP_pT; i++)
	{
		dN_dypTdpT[i] = 0.0e0;
		for(int j=0; j<n_order; j++)
		{
			cosine_iorder[i][j] = 0.0e0;
			sine_iorder[i][j] = 0.0e0;
		}
	}
	plane_angle = new double [n_order];

	//pair momentum
	K_T = new double [n_localp_T];
	double dK_T = (localp_T_max - localp_T_min)/(n_localp_T - 1 + 1e-100);
	for(int i=0; i<n_localp_T; i++) K_T[i] = localp_T_min + i*dK_T;
	//K_y = p_y;
	K_y = 0.;
	ch_K_y = cosh(K_y);
	sh_K_y = sinh(K_y);
	beta_l = sh_K_y/ch_K_y;
	K_phi = new double [n_localp_phi];
	K_phi_weight = new double [n_localp_phi];
	gauss_quadrature(n_localp_phi, 1, 0.0, 0.0, localp_phi_min, localp_phi_max, K_phi, K_phi_weight);

	//spatial rapidity grid
	eta_s = new double [eta_s_npts];
	eta_s_weight = new double [eta_s_npts];
	gauss_quadrature(eta_s_npts, 1, 0.0, 0.0, eta_s_i, eta_s_f, eta_s, eta_s_weight);
	ch_eta_s = new double [eta_s_npts];
	sh_eta_s = new double [eta_s_npts];
	for (int ieta = 0; ieta < eta_s_npts; ieta++)
	{
		ch_eta_s[ieta] = cosh(eta_s[ieta]);
		sh_eta_s[ieta] = sinh(eta_s[ieta]);
	}
debugger(__LINE__, __FILE__);
	R2_side = new double* [n_localp_T];
	R2_side_C = new double* [n_localp_T];
	R2_side_S = new double* [n_localp_T];
	R2_out = new double* [n_localp_T];
	R2_out_C = new double* [n_localp_T];
	R2_out_S = new double* [n_localp_T];
	R2_long = new double* [n_localp_T];
	R2_long_C = new double* [n_localp_T];
	R2_long_S = new double* [n_localp_T];
	R2_outside = new double* [n_localp_T];
	R2_outside_C = new double* [n_localp_T];
	R2_outside_S = new double* [n_localp_T];
	R2_sidelong = new double* [n_localp_T];
	R2_sidelong_C = new double* [n_localp_T];
	R2_sidelong_S = new double* [n_localp_T];
	R2_outlong = new double* [n_localp_T];
	R2_outlong_C = new double* [n_localp_T];
	R2_outlong_S = new double* [n_localp_T];

	R2_side_err = new double* [n_localp_T];
	R2_out_err = new double* [n_localp_T];
	R2_long_err = new double* [n_localp_T];
	R2_outside_err = new double* [n_localp_T];
	R2_sidelong_err = new double* [n_localp_T];
	R2_outlong_err = new double* [n_localp_T];

	lambda_Correl = new double * [n_localp_T];
	lambda_Correl_err = new double * [n_localp_T];

	for(int i=0; i<n_localp_T; i++)
	{
		R2_side[i] = new double [n_localp_phi];
		R2_side_C[i] = new double [n_order];
		R2_side_S[i] = new double [n_order];
		R2_out[i] = new double [n_localp_phi];
		R2_out_C[i] = new double [n_order];
		R2_out_S[i] = new double [n_order];
		R2_outside[i] = new double [n_localp_phi];
		R2_outside_C[i] = new double [n_order];
		R2_outside_S[i] = new double [n_order];
		R2_long[i] = new double [n_localp_phi];
		R2_long_C[i] = new double [n_order];
		R2_long_S[i] = new double [n_order];
		R2_sidelong[i] = new double [n_localp_phi];
		R2_sidelong_C[i] = new double [n_order];
		R2_sidelong_S[i] = new double [n_order];
		R2_outlong[i] = new double [n_localp_phi];
		R2_outlong_C[i] = new double [n_order];
		R2_outlong_S[i] = new double [n_order];

		R2_side_err[i] = new double [n_localp_phi];
		R2_out_err[i] = new double [n_localp_phi];
		R2_long_err[i] = new double [n_localp_phi];
		R2_outside_err[i] = new double [n_localp_phi];
		R2_sidelong_err[i] = new double [n_localp_phi];
		R2_outlong_err[i] = new double [n_localp_phi];

		lambda_Correl[i] = new double [n_localp_phi];
		lambda_Correl_err[i] = new double [n_localp_phi];
	}
debugger(__LINE__, __FILE__);
	//initialize all source variances and HBT radii/coeffs
	for(int i=0; i<n_localp_T; i++)
	{
		for(int j=0; j<n_localp_phi; j++)
		{
			R2_side[i][j] = 0.;
			R2_out[i][j] = 0.;
			R2_long[i][j] = 0.;
			R2_outside[i][j] = 0.;
			R2_sidelong[i][j] = 0.;
			R2_outlong[i][j] = 0.;

			R2_side_err[i][j] = 0.;
			R2_out_err[i][j] = 0.;
			R2_long_err[i][j] = 0.;
			R2_outside_err[i][j] = 0.;
			R2_sidelong_err[i][j] = 0.;
			R2_outlong_err[i][j] = 0.;

			lambda_Correl[i][j] = 0.0;
			lambda_Correl_err[i][j] = 0.0;
		}
		for(int j=0; j<n_order; j++)
		{
			R2_side_C[i][j] = 0.;
			R2_side_S[i][j] = 0.;
			R2_out_C[i][j] = 0.;
			R2_out_S[i][j] = 0.;
			R2_outside_C[i][j] = 0.;
			R2_outside_S[i][j] = 0.;
			R2_long_C[i][j] = 0.;
			R2_long_S[i][j] = 0.;
			R2_sidelong_C[i][j] = 0.;
			R2_sidelong_S[i][j] = 0.;
			R2_outlong_C[i][j] = 0.;
			R2_outlong_S[i][j] = 0.;
		}
	}

   return;
}

void CorrelationFunction::Update_sourcefunction(particle_info* particle, int FOarray_length, int particle_idx)
{
	full_FO_length = FOarray_length * eta_s_npts;

	osc0 = new double *** [FOarray_length];					//to hold cos/sin(q0 t)
	osc1 = new double ** [FOarray_length];					//to hold cos/sin(qx x)
	osc2 = new double ** [FOarray_length];					//to hold cos/sin(qy y)
	osc3 = new double *** [FOarray_length];				//to hold cos/sin(+/- qz z)

	for (int isurf = 0; isurf < FOarray_length; ++isurf)
	{
		FO_surf * surf = &current_FOsurf_ptr[isurf];

		double tau = surf->tau;
		double xpt = surf->xpt;
		double ypt = surf->ypt;

		osc0[isurf] = new double ** [eta_s_npts];
		osc1[isurf] = new double * [qnpts];
		osc2[isurf] = new double * [qnpts];
		osc3[isurf] = new double ** [eta_s_npts];

		for (int iq = 0; iq < qnpts; ++iq)
		{
			osc1[isurf][iq] = new double [2];
			osc2[isurf][iq] = new double [2];

			osc1[isurf][iq][0] = cos(hbarCm1*qx_pts[iq]*xpt);
			osc2[isurf][iq][0] = cos(hbarCm1*qy_pts[iq]*ypt);
			osc1[isurf][iq][1] = sin(hbarCm1*qx_pts[iq]*xpt);
			osc2[isurf][iq][1] = sin(hbarCm1*qy_pts[iq]*ypt);
		}

		for (int ieta = 0; ieta < eta_s_npts; ++ieta)
		{
			double tpt = tau*ch_eta_s[ieta];
			double zpt = tau*sh_eta_s[ieta];

			osc0[isurf][ieta] = new double * [qnpts];
			osc3[isurf][ieta] = new double * [qnpts];

			for (int iq = 0; iq < qnpts; ++iq)
			{
				osc0[isurf][ieta][iq] = new double [2];
				osc3[isurf][ieta][iq] = new double [2];

				osc0[isurf][ieta][iq][0] = cos(hbarCm1*qt_pts[iq]*tpt);
				osc3[isurf][ieta][iq][0] = cos(hbarCm1*qz_pts[iq]*zpt);
				osc0[isurf][ieta][iq][1] = sin(hbarCm1*qt_pts[iq]*tpt);
				osc3[isurf][ieta][iq][1] = sin(hbarCm1*qz_pts[iq]*zpt);
			}
		}
	}

	S_p_withweight_array = new double *** [n_interp_pT_pts];
	zero_FOcell_flag = new bool *** [n_interp_pT_pts];
	for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
	{
		S_p_withweight_array[ipt] = new double ** [n_interp_pphi_pts];
		zero_FOcell_flag[ipt] = new bool ** [n_interp_pphi_pts];
		for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
		{
			//S_p_withweight_array[ipt][ipphi] = new double * [FOarray_length];
			zero_FOcell_flag[ipt][ipphi] = new bool * [FOarray_length];
			for (int isurf = 0; isurf < FOarray_length; ++isurf)
			{
				//S_p_withweight_array[ipt][ipphi][isurf] = new double [eta_s_npts];
				zero_FOcell_flag[ipt][ipphi][isurf] = new bool [eta_s_npts];
				for (int ieta = 0; ieta < eta_s_npts; ++ieta)
				{
					//S_p_withweight_array[ipt][ipphi][isurf][ieta] = 0.0;
					zero_FOcell_flag[ipt][ipphi][isurf][ieta] = false;	//assume by default that all FOcells need to be included, i.e., cannot be zeroed
				}
			}
		}
	}

   //particle information
   particle_name = particle->name;
   particle_mass = particle->mass;
   particle_sign = particle->sign;
   particle_gspin = particle->gspin;
   particle_id = particle_idx;

	*global_out_stream_ptr << "Inside Update_sourcefunction(...): using fraction_of_resonances = " << fraction_of_resonances << endl;

   //erase contents of single - and two-particle spectra
   for(int i=0; i<n_SP_pphi; i++)
   {
      dN_dydphi[i] = 0.0e0;
      pTdN_dydphi[i] = 0.0e0;
      for(int j=0; j<n_SP_pT; j++) dN_dypTdpTdphi[j][i] = 0.0e0;
   }
   //erase anisotropic flows
   for(int i=0; i<n_SP_pT; i++)
   for(int j=0; j<n_order; j++)
   {
      cosine_iorder[i][j] = 0.0e0;
      sine_iorder[i][j] = 0.0e0;
   }

   FO_length = FOarray_length;

	// set the rest later
	most_important_FOcells = new size_t ** [n_interp_pT_pts];
	for(int ipt=0; ipt<n_interp_pT_pts; ipt++)
	{
		most_important_FOcells[ipt] = new size_t * [n_interp_pphi_pts];
		for(int ipphi=0; ipphi<n_interp_pphi_pts; ipphi++)
			most_important_FOcells[ipt][ipphi] = new size_t [FO_length];
	}


//reset only EBE source variances and EBE HBT radii/coeffs
for(int i=0; i<n_localp_T; i++)
{
	for(int j=0; j<n_order; j++)
	{
		R2_side_C[i][j] = 0.;
		R2_side_S[i][j] = 0.;
		R2_out_C[i][j] = 0.;
		R2_out_S[i][j] = 0.;
		R2_outside_C[i][j] = 0.;
		R2_outside_S[i][j] = 0.;
		R2_long_C[i][j] = 0.;
		R2_long_S[i][j] = 0.;
		R2_sidelong_C[i][j] = 0.;
		R2_sidelong_S[i][j] = 0.;
		R2_outlong_C[i][j] = 0.;
		R2_outlong_S[i][j] = 0.;
	}
}

   return;
}

CorrelationFunction::~CorrelationFunction()
{
   delete[] SP_pT;
   delete[] SP_pT_weight;
   delete[] SP_pphi;
   delete[] SP_pphi_weight;
   delete[] dN_dydphi;
   delete[] dN_dypTdpT;
   delete[] pTdN_dydphi;
   for(int i=0; i<n_SP_pT; i++)
   {
      delete[] dN_dypTdpTdphi[i];
      delete[] cosine_iorder[i];
      delete[] sine_iorder[i];
   }
   delete[] dN_dypTdpTdphi;
   delete[] cosine_iorder;
   delete[] sine_iorder;
   delete[] plane_angle;

   delete[] K_T;
   delete[] K_phi;
   delete[] K_phi_weight;
   delete[] eta_s;
   delete[] eta_s_weight;

   for(int i=0; i<n_localp_T; i++)
   {
      delete[] R2_side[i];
   }

   delete[] R2_side;

   return;
}

void CorrelationFunction::Delete_S_p_withweight_array()
{
	for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
	{
		for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
		{
			for (int isurf = 0; isurf < FO_length; ++isurf)
				delete [] zero_FOcell_flag[ipt][ipphi][isurf];
			delete [] zero_FOcell_flag[ipt][ipphi];
		}
		delete [] S_p_withweight_array[ipt];
		delete [] zero_FOcell_flag[ipt];
	}
	delete [] S_p_withweight_array;
	delete [] zero_FOcell_flag;

	return;
}

void CorrelationFunction::Reset_zero_FOcell_flag_array()
{
	for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
	for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
	for (int isurf = 0; isurf < FO_length; ++isurf)
	for (int ieta = 0; ieta < eta_s_npts; ++ieta)
		zero_FOcell_flag[ipt][ipphi][isurf][ieta] = false;	//assume by default that all FOcells need to be included, i.e., cannot be zeroed

	return;
}

void CorrelationFunction::Allocate_decay_channel_info()
{
	if (VERBOSE > 2) *global_out_stream_ptr << "Reallocating memory for decay channel information..." << endl;
	VEC_n2_v_factor = new double [n_v_pts];
	VEC_n2_zeta_factor = new double * [n_v_pts];
	VEC_n2_P_Y = new double [n_v_pts];
	VEC_n2_MTbar = new double [n_v_pts];
	VEC_n2_DeltaMT = new double [n_v_pts];
	VEC_n2_MTp = new double [n_v_pts];
	VEC_n2_MTm = new double [n_v_pts];
	VEC_n2_MT = new double * [n_v_pts];
	VEC_n2_PPhi_tilde = new double * [n_v_pts];
	VEC_n2_PPhi_tildeFLIP = new double * [n_v_pts];
	VEC_n2_PT = new double * [n_v_pts];
	for(int iv = 0; iv < n_v_pts; ++iv)
	{
		VEC_n2_MT[iv] = new double [n_zeta_pts];
		VEC_n2_PPhi_tilde[iv] = new double [n_zeta_pts];
		VEC_n2_PPhi_tildeFLIP[iv] = new double [n_zeta_pts];
		VEC_n2_PT[iv] = new double [n_zeta_pts];
		VEC_n2_zeta_factor[iv] = new double [n_zeta_pts];
	}
	NEW_s_pts = new double [n_s_pts];
	NEW_s_wts = new double [n_s_pts];
	VEC_pstar = new double [n_s_pts];
	VEC_Estar = new double [n_s_pts];
	VEC_DeltaY = new double [n_s_pts];
	VEC_g_s = new double [n_s_pts];
	VEC_s_factor = new double [n_s_pts];
	VEC_v_factor = new double * [n_s_pts];
	VEC_zeta_factor = new double ** [n_s_pts];
	VEC_Yp = new double [n_s_pts];
	VEC_Ym = new double [n_s_pts];
	VEC_P_Y = new double * [n_s_pts];
	VEC_MTbar = new double * [n_s_pts];
	VEC_DeltaMT = new double * [n_s_pts];
	VEC_MTp = new double * [n_s_pts];
	VEC_MTm = new double * [n_s_pts];
	VEC_MT = new double ** [n_s_pts];
	VEC_PPhi_tilde = new double ** [n_s_pts];
	VEC_PPhi_tildeFLIP = new double ** [n_s_pts];
	VEC_PT = new double ** [n_s_pts];
	for(int is = 0; is < n_s_pts; ++is)
	{
		VEC_v_factor[is] = new double [n_v_pts];
		VEC_zeta_factor[is] = new double * [n_v_pts];
		VEC_P_Y[is] = new double [n_v_pts];
		VEC_MTbar[is] = new double [n_v_pts];
		VEC_DeltaMT[is] = new double [n_v_pts];
		VEC_MTp[is] = new double [n_v_pts];
		VEC_MTm[is] = new double [n_v_pts];
		VEC_MT[is] = new double * [n_v_pts];
		VEC_PPhi_tilde[is] = new double * [n_v_pts];
		VEC_PPhi_tildeFLIP[is] = new double * [n_v_pts];
		VEC_PT[is] = new double * [n_v_pts];
		for(int iv = 0; iv < n_v_pts; ++iv)
		{
			VEC_MT[is][iv] = new double [n_zeta_pts];
			VEC_PPhi_tilde[is][iv] = new double [n_zeta_pts];
			VEC_PPhi_tildeFLIP[is][iv] = new double [n_zeta_pts];
			VEC_PT[is][iv] = new double [n_zeta_pts];
			VEC_zeta_factor[is][iv] = new double [n_zeta_pts];
		}
	}
	if (VERBOSE > 2) *global_out_stream_ptr << "Reallocated memory for decay channel information." << endl;

	return;
}

void CorrelationFunction::Delete_decay_channel_info()
{
	if (VERBOSE > 2) *global_out_stream_ptr << "Deleting memory for decay channel information..." << endl;
	for(int iv = 0; iv < n_v_pts; ++iv)
	{
		delete [] VEC_n2_MT[iv];
		delete [] VEC_n2_PPhi_tilde[iv];
		delete [] VEC_n2_PPhi_tildeFLIP[iv];
		delete [] VEC_n2_PT[iv];
		delete [] VEC_n2_zeta_factor[iv];
	}
	delete [] VEC_n2_v_factor;
	delete [] VEC_n2_zeta_factor;
	delete [] VEC_n2_P_Y;
	delete [] VEC_n2_MTbar ;
	delete [] VEC_n2_DeltaMT;
	delete [] VEC_n2_MTp;
	delete [] VEC_n2_MTm;
	delete [] VEC_n2_MT;
	delete [] VEC_n2_PPhi_tilde;
	delete [] VEC_n2_PPhi_tildeFLIP;
	delete [] VEC_n2_PT;

	for(int is = 0; is < n_s_pts; ++is)
	{
		for(int iv = 0; iv < n_v_pts; ++iv)
		{
			delete [] VEC_MT[is][iv];
			delete [] VEC_PPhi_tilde[is][iv];
			delete [] VEC_PPhi_tildeFLIP[is][iv];
			delete [] VEC_PT[is][iv];
			delete [] VEC_zeta_factor[is][iv];
		}
		delete [] VEC_v_factor[is];
		delete [] VEC_zeta_factor[is];
		delete [] VEC_P_Y[is];
		delete [] VEC_MTbar[is];
		delete [] VEC_DeltaMT[is];
		delete [] VEC_MTp[is];
		delete [] VEC_MTm[is];
		delete [] VEC_MT[is];
		delete [] VEC_PPhi_tilde[is];
		delete [] VEC_PPhi_tildeFLIP[is];
		delete [] VEC_PT[is];
	}
	delete [] NEW_s_pts;
	delete [] NEW_s_wts;
	delete [] VEC_pstar;
	delete [] VEC_Estar;
	delete [] VEC_DeltaY;
	delete [] VEC_g_s;
	delete [] VEC_s_factor;
	delete [] VEC_v_factor;
	delete [] VEC_zeta_factor;
	delete [] VEC_Yp;
	delete [] VEC_Ym;
	delete [] VEC_P_Y;
	delete [] VEC_MTbar;
	delete [] VEC_DeltaMT;
	delete [] VEC_MTp;
	delete [] VEC_MTm;
	delete [] VEC_MT;
	delete [] VEC_PPhi_tilde;
	delete [] VEC_PPhi_tildeFLIP;
	delete [] VEC_PT;
	if (VERBOSE > 2) *global_out_stream_ptr << "Deleted memory for decay channel information." << endl;

	return;
}

// sets points in q-space for computing weighted spectra grid
void CorrelationFunction::Set_q_points()
{
	q_pts = new double [qnpts];
	for (int iq = 0; iq < qnpts; ++iq)
		q_pts[iq] = init_q + (double)iq * delta_q;
	q_axes = new double [3];

	qt_pts = new double [qnpts];
	qx_pts = new double [qnpts];
	qy_pts = new double [qnpts];
	qz_pts = new double [qnpts];
	for (int iq = 0; iq < qnpts; ++iq)
	{
		qt_pts[iq] = init_q + (double)iq * delta_q;
		//qt_pts[iq] = 0.0 + 0.5 * (double)iq * delta_q;
		qx_pts[iq] = init_q + (double)iq * delta_q;
		qy_pts[iq] = init_q + (double)iq * delta_q;
		qz_pts[iq] = init_q + (double)iq * delta_q;
		if (abs(qt_pts[iq]) < 1.e-15)
			iqt0 = iq;
		if (abs(qx_pts[iq]) < 1.e-15)
			iqx0 = iq;
		if (abs(qy_pts[iq]) < 1.e-15)
			iqy0 = iq;
		if (abs(qz_pts[iq]) < 1.e-15)
			iqz0 = iq;
	}

	cerr << "Output iq*0 = " << iqt0 << "   " << iqx0 << "   " << iqy0 << "   " << iqz0 << endl;

if (1) exit;

	return;
}

void CorrelationFunction::Set_correlation_function_q_pts()
{
	// initialize q-points for actual correlation function
	qo_pts = new double [qonpts];
	qs_pts = new double [qsnpts];
	ql_pts = new double [qlnpts];
	for (int iq = 0; iq < qonpts; ++iq)
		qo_pts[iq] = init_q + (double)iq * delta_q;
	for (int iq = 0; iq < qsnpts; ++iq)
		qs_pts[iq] = init_q + (double)iq * delta_q;
	for (int iq = 0; iq < qlnpts; ++iq)
		ql_pts[iq] = init_q + (double)iq * delta_q;

	// initialize error matrix
	Correl_3D_err = new double [qonpts];
	for (int iqo = 0; iqo < qonpts; ++iqo)
	{
		Correl_3D_err[iqo] = new double [qsnpts];
		for (int iqs = 0; iqs < qsnpts; ++iqs)
		{
			Correl_3D_err[iqo][iqs] = new double [qlnpts];
			for (int iql = 0; iql < qlnpts; ++iql)
				Correl_3D_err[iqo][iqs][iql] = 0.0;
		}
	}

	return;
}


// returns points in q-space for computing weighted spectra grid corresponding to to given q and K choices
// weighted spectra grid thus needs to be interpolated at point returned in qgridpts
void CorrelationFunction::Get_q_points(double qo, double qs, double ql, double pT, double pphi, double * qgridpts)
{
	double mtarget = all_particles[target_particle_id].mass;
	double xi2 = 0.25*mtarget*mtarget + pT*pT + qo*qo + qs*qs + ql*ql;
	double ckp = cos(pphi), skp = sin(pphi);

	// set qpts at which to interpolate spectra
	qgridpts[0] = sqrt(xi2 + qo*pT) - sqrt(xi2 - qo*pT);	//set qt component
	qgridpts[1] = qo*ckp - qs*skp;							//set qx component
	qgridpts[2] = qo*skp + qs*ckp;							//set qy component
	qgridpts[3] = ql;										//set qz component, since qz = ql

	return;
}

bool CorrelationFunction::fexists(const char *filename)
{
  ifstream ifile(filename);
  return ifile;
}

//print output to output filestream, one line at a time
void CorrelationFunction::Set_ofstream(ofstream& myout)
{
	global_out_stream_ptr = &myout;

	return;
}

//print output to output filestream, one line at a time
void CorrelationFunction::Set_path(string localpath)
{
	global_path = localpath;

	return;
}

void CorrelationFunction::Set_runfolder(string localrunfolder)
{
	global_runfolder = localrunfolder;

	return;
}

void CorrelationFunction::Set_resultsfolder_stem(string usrdef_resultsfolder_stem)
{
	global_resultsfolder_stem = usrdef_resultsfolder_stem;

	return;
}

void CorrelationFunction::Set_use_delta_f(bool usrdef_usedeltaf)
{
	use_delta_f = usrdef_usedeltaf;
	if (!use_delta_f)
		no_df_stem = "_no_df";
	return;
}

void CorrelationFunction::Set_particle_mass(double usrdef_particle_mass)
{
	particle_mass = usrdef_particle_mass;
	return;
}

void CorrelationFunction::Set_current_FOsurf_ptr(FO_surf* FOsurf_ptr)
{
	current_FOsurf_ptr = FOsurf_ptr;
	return;
}

void CorrelationFunction::Get_current_decay_string(int dc_idx, string * decay_string)
{
	// N.B. - dc_idx == 0 is thermal pion(+)s in calling loop, dc_idx > 0 gives resonance decays
	//      ==> need to use dc_idx - 1 here
	*decay_string = decay_channels[dc_idx - 1].resonance_name + " --->> ";
	int temp_monval, tempID;
	for (int decay_part_idx = 0; decay_part_idx < decay_channels[dc_idx - 1].nbody; decay_part_idx++)
	{
		temp_monval = decay_channels[dc_idx - 1].resonance_decay_monvals[decay_part_idx];
		//if (VERBOSE > 0) *global_out_stream_ptr << "Get_current_decay_string(): temp_monval = " << temp_monval << endl;
		if (temp_monval == 0)
			continue;
		else
		{
			tempID = lookup_particle_id_from_monval(all_particles, Nparticle, temp_monval);
			*decay_string += all_particles[tempID].name;
			if (decay_part_idx < decay_channels[dc_idx - 1].nbody - 1) *decay_string += " + ";
		}
	}
	return;
}

int CorrelationFunction::Set_daughter_list(int parent_resonance_index)
{
	// reset list
	daughter_resonance_indices.clear();
	
	// then re-populate it
	particle_info parent = all_particles[parent_resonance_index];
	if (parent.stable == 1 && parent.decays_Npart[0] == 1)
		return (0);									// no daughters to worry about if parent resonance is actually stable
	int number_of_decays = parent.decays;
	for (int k = 0; k < number_of_decays; k++)		// loop through decays for parent resonance
	{
		int nb = abs(parent.decays_Npart[k]);		// for each decay, nb is the number of daughter particles
		for (int l = 0; l < nb; l++)				// loop through each daughter particle
		{
			int pid = lookup_particle_id_from_monval(all_particles, Nparticle, parent.decays_part[k][l]);
			daughter_resonance_indices.insert(pid);		// using a <set> object will automatically remove duplicates and keep pid's in a fixed order
		}
	}

	// return value is total number of daughters found
	return (daughter_resonance_indices.size());
}

int CorrelationFunction::lookup_resonance_idx_from_particle_id(int pid)
{
	// pid - particle index in all_particles array
	// looks up location in chosen_resonances of given value particle_id
	int result = -1;

	for (int ii = 0; ii < (int)chosen_resonances.size(); ii++)
	{
		if (chosen_resonances[ii] == pid)
		{
			result = ii;
			break;
		}
	}

	// if pid is not one of the chosen_resonances, is not the target daughter (pion(+)), is not stable and has a non-zero effective branching ratio
	if (result < 0 && pid != particle_id && all_particles[pid].stable == 0 && all_particles[pid].effective_branchratio >= 1.e-12)
	{
		*global_out_stream_ptr << " *** lookup_resonance_idx_from_particle_id(): Particle_id = " << pid
					<< " (" << all_particles[pid].name <<") not found in chosen_resonances!" << endl
					<< " *** br = " << all_particles[pid].effective_branchratio << endl;
	}
	return (result);
}

void CorrelationFunction::Setup_temp_arrays(double ***** local_temp_moments, double ******* temp_moments_array)
{
	//set first argument
	local_temp_moments = new double **** [qnpts];
	for (int iqt = 0; iqt < qnpts; ++iqt)
	{
		local_temp_moments[iqt] = new double *** [qnpts];
		for (int iqx = 0; iqx < qnpts; ++iqx)
		{
			local_temp_moments[iqt][iqx] = new double ** [qnpts];
			for (int iqy = 0; iqy < qnpts; ++iqy)
			{
				local_temp_moments[iqt][iqx][iqy] = new double * [qnpts];
				for (int iqz = 0; iqz < qnpts; ++iqz)
				{
					local_temp_moments[iqt][iqx][iqy][iqz] = new double [2];
					for (int itrig = 0; itrig < 2; ++itrig)
						local_temp_moments[iqt][iqx][iqy][iqz][itrig] = 0.0;
				}
			}
		}
	}
	
	// set second argument
	temp_moments_array = new double ****** [n_interp_pT_pts];
	for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
	{
		temp_moments_array[ipt] = new double ***** [n_interp_pphi_pts];
		for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
		{
			temp_moments_array[ipt][ipphi] = new double **** [qnpts];
			for (int iqt = 0; iqt < qnpts; ++iqt)
			{
				temp_moments_array[ipt][ipphi][iqt] = new double *** [qnpts];
				for (int iqx = 0; iqx < qnpts; ++iqx)
				{
					temp_moments_array[ipt][ipphi][iqt][iqx] = new double ** [qnpts];
					for (int iqy = 0; iqy < qnpts; ++iqy)
					{
						temp_moments_array[ipt][ipphi][iqt][iqx][iqy] = new double * [qnpts];
						for (int iqz = 0; iqz < qnpts; ++iqz)
						{
							temp_moments_array[ipt][ipphi][iqt][iqx][iqy][iqz] = new double [2];
							for (int itrig = 0; itrig < 2; ++itrig)
								temp_moments_array[ipt][ipphi][iqt][iqx][iqy][iqz][itrig] = 0.0;
						}
					}
				}
			}
		}
	}

	return;
}

void CorrelationFunction::Teardown_temp_arrays(double ***** local_temp_moments, double ******* temp_moments_array)
{
	//clean up first argument
	for (int iqt = 0; iqt < qnpts; ++iqt)
	{
		for (int iqx = 0; iqx < qnpts; ++iqx)
		{
			for (int iqy = 0; iqy < qnpts; ++iqy)
			{
				for (int iqz = 0; iqz < qnpts; ++iqz)
					delete [] local_temp_moments[iqt][iqx][iqy][iqz];
				delete [] local_temp_moments[iqt][iqx][iqy];
			}
			delete [] local_temp_moments[iqt][iqx];
		}
		delete [] local_temp_moments[iqt];
	}
	delete [] local_temp_moments;

	//clean up second argument
	for(int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
	{
		for(int iphi = 0; iphi < n_interp_pphi_pts; ++iphi)
		{
			for (int iqt = 0; iqt < qnpts; ++iqt)
			{
				for (int iqx = 0; iqx < qnpts; ++iqx)
				{
					for (int iqy = 0; iqy < qnpts; ++iqy)
					{
						for (int iqz = 0; iqz < qnpts; ++iqz)
							delete [] temp_moments_array[ipt][iphi][iqt][iqx][iqy][iqz];
						delete [] temp_moments_array[ipt][iphi][iqt][iqx][iqy];
					}
					delete [] temp_moments_array[ipt][iphi][iqt][iqx];
				}
				delete [] temp_moments_array[ipt][iphi][iqt];
			}
			delete [] temp_moments_array[ipt][iphi];
		}
		delete [] temp_moments_array[ipt];
	}
	delete [] temp_moments_array;

	return;
}

void CorrelationFunction::Allocate_resonance_running_sum_vectors()
{
    ssum_vec = new double [qspace_cs_slice_length];
    vsum_vec = new double [qspace_cs_slice_length];
    zetasum_vec = new double [qspace_cs_slice_length];
    Csum_vec = new double [qspace_cs_slice_length];
}

void CorrelationFunction::Delete_resonance_running_sum_vectors()
{
    delete [] ssum_vec;
    delete [] vsum_vec;
    delete [] zetasum_vec;
    delete [] Csum_vec;
}

void CorrelationFunction::Zero_resonance_running_sum_vector(double * vec)
{
	for (int tmp = 0; tmp < qspace_cs_slice_length; ++tmp)
		vec[tmp] = 0.0;
}

void CorrelationFunction::Setup_current_daughters_dN_dypTdpTdphi_moments(int n_daughter)
{
	current_daughters_dN_dypTdpTdphi_moments = new double ******* [n_daughter];
	current_daughters_ln_dN_dypTdpTdphi_moments = new double ******* [n_daughter];
	current_daughters_sign_of_dN_dypTdpTdphi_moments = new double ******* [n_daughter];
	for (int id = 0; id < n_daughter; ++id)
	{
		current_daughters_dN_dypTdpTdphi_moments[id] = new double ****** [n_interp_pT_pts];
		current_daughters_ln_dN_dypTdpTdphi_moments[id] = new double ****** [n_interp_pT_pts];
		current_daughters_sign_of_dN_dypTdpTdphi_moments[id] = new double ****** [n_interp_pT_pts];
		for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
		{
			current_daughters_dN_dypTdpTdphi_moments[id][ipt] = new double ***** [n_interp_pphi_pts];
			current_daughters_ln_dN_dypTdpTdphi_moments[id][ipt] = new double ***** [n_interp_pphi_pts];
			current_daughters_sign_of_dN_dypTdpTdphi_moments[id][ipt] = new double ***** [n_interp_pphi_pts];
			for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
			{
				current_daughters_dN_dypTdpTdphi_moments[id][ipt][ipphi] = new double **** [qnpts];
				current_daughters_ln_dN_dypTdpTdphi_moments[id][ipt][ipphi] = new double **** [qnpts];
				current_daughters_sign_of_dN_dypTdpTdphi_moments[id][ipt][ipphi] = new double **** [qnpts];
				for (int iqt = 0; iqt < qnpts; ++iqt)
				{
					current_daughters_dN_dypTdpTdphi_moments[id][ipt][ipphi][iqt] = new double *** [qnpts];
					current_daughters_ln_dN_dypTdpTdphi_moments[id][ipt][ipphi][iqt] = new double *** [qnpts];
					current_daughters_sign_of_dN_dypTdpTdphi_moments[id][ipt][ipphi][iqt] = new double *** [qnpts];
					for (int iqx = 0; iqx < qnpts; ++iqx)
					{
						current_daughters_dN_dypTdpTdphi_moments[id][ipt][ipphi][iqt][iqx] = new double ** [qnpts];
						current_daughters_ln_dN_dypTdpTdphi_moments[id][ipt][ipphi][iqt][iqx] = new double ** [qnpts];
						current_daughters_sign_of_dN_dypTdpTdphi_moments[id][ipt][ipphi][iqt][iqx] = new double ** [qnpts];
						for (int iqy = 0; iqy < qnpts; ++iqy)
						{
							current_daughters_dN_dypTdpTdphi_moments[id][ipt][ipphi][iqt][iqx][iqy] = new double * [qnpts];
							current_daughters_ln_dN_dypTdpTdphi_moments[id][ipt][ipphi][iqt][iqx][iqy] = new double * [qnpts];
							current_daughters_sign_of_dN_dypTdpTdphi_moments[id][ipt][ipphi][iqt][iqx][iqy] = new double * [qnpts];
							for (int iqz = 0; iqz < qnpts; ++iqz)
							{
								current_daughters_dN_dypTdpTdphi_moments[id][ipt][ipphi][iqt][iqx][iqy][iqz] = new double [2];
								current_daughters_ln_dN_dypTdpTdphi_moments[id][ipt][ipphi][iqt][iqx][iqy][iqz] = new double [2];
								current_daughters_sign_of_dN_dypTdpTdphi_moments[id][ipt][ipphi][iqt][iqx][iqy][iqz] = new double [2];
								for (int itrig = 0; itrig < 2; ++itrig)
								{
									current_daughters_dN_dypTdpTdphi_moments[id][ipt][ipphi][iqt][iqx][iqy][iqz][itrig] = 0.0;
									current_daughters_ln_dN_dypTdpTdphi_moments[id][ipt][ipphi][iqt][iqx][iqy][iqz][itrig] = 0.0;
									current_daughters_sign_of_dN_dypTdpTdphi_moments[id][ipt][ipphi][iqt][iqx][iqy][iqz][itrig] = 0.0;
								}
							}
						}
					}
				}
			}
		}
	}
	return;
}

void CorrelationFunction::Cleanup_current_daughters_dN_dypTdpTdphi_moments(int n_daughter)
{
	for (int id = 0; id < n_daughter; ++id)
	{
		for(int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
		{
			for(int iphi = 0; iphi < n_interp_pphi_pts; ++iphi)
			{
				for (int iqt = 0; iqt < qnpts; ++iqt)
				{
					for (int iqx = 0; iqx < qnpts; ++iqx)
					{
						for (int iqy = 0; iqy < qnpts; ++iqy)
						{
							for (int iqz = 0; iqz < qnpts; ++iqz)
							{
								delete [] current_daughters_dN_dypTdpTdphi_moments[id][ipt][iphi][iqt][iqx][iqy][iqz];
								delete [] current_daughters_ln_dN_dypTdpTdphi_moments[id][ipt][iphi][iqt][iqx][iqy][iqz];
								delete [] current_daughters_sign_of_dN_dypTdpTdphi_moments[id][ipt][iphi][iqt][iqx][iqy][iqz];
							}
							delete [] current_daughters_dN_dypTdpTdphi_moments[id][ipt][iphi][iqt][iqx][iqy];
							delete [] current_daughters_ln_dN_dypTdpTdphi_moments[id][ipt][iphi][iqt][iqx][iqy];
							delete [] current_daughters_sign_of_dN_dypTdpTdphi_moments[id][ipt][iphi][iqt][iqx][iqy];
						}
						delete [] current_daughters_dN_dypTdpTdphi_moments[id][ipt][iphi][iqt][iqx];
						delete [] current_daughters_ln_dN_dypTdpTdphi_moments[id][ipt][iphi][iqt][iqx];
						delete [] current_daughters_sign_of_dN_dypTdpTdphi_moments[id][ipt][iphi][iqt][iqx];
					}
					delete [] current_daughters_dN_dypTdpTdphi_moments[id][ipt][iphi][iqt];
					delete [] current_daughters_ln_dN_dypTdpTdphi_moments[id][ipt][iphi][iqt];
					delete [] current_daughters_sign_of_dN_dypTdpTdphi_moments[id][ipt][iphi][iqt];
				}
				delete [] current_daughters_dN_dypTdpTdphi_moments[id][ipt][iphi];
				delete [] current_daughters_ln_dN_dypTdpTdphi_moments[id][ipt][iphi];
				delete [] current_daughters_sign_of_dN_dypTdpTdphi_moments[id][ipt][iphi];
			}
			delete [] current_daughters_dN_dypTdpTdphi_moments[id][ipt];
			delete [] current_daughters_ln_dN_dypTdpTdphi_moments[id][ipt];
			delete [] current_daughters_sign_of_dN_dypTdpTdphi_moments[id][ipt];
		}
		delete [] current_daughters_dN_dypTdpTdphi_moments[id];
		delete [] current_daughters_ln_dN_dypTdpTdphi_moments[id];
		delete [] current_daughters_sign_of_dN_dypTdpTdphi_moments[id];
	}
	delete [] current_daughters_dN_dypTdpTdphi_moments;
	delete [] current_daughters_ln_dN_dypTdpTdphi_moments;
	delete [] current_daughters_sign_of_dN_dypTdpTdphi_moments;

	return;
}

inline double CorrelationFunction::lin_int(double x_m_x1, double one_by_x2_m_x1, double f1, double f2)
{
	return ( f1 + (f2 - f1) * x_m_x1 * one_by_x2_m_x1 );
}

void CorrelationFunction::Edndp3(double ptr, double phir, double * results)
{
	double phi0, phi1;
	double f1, f2;

	int npphi_max = n_interp_pphi_pts - 1;
	int npT_max = n_interp_pT_pts - 1;

	// locate pT interval
	int npt = 1;
	while ((ptr > SPinterp_pT[npt]) &&
			(npt < npT_max)) ++npt;
	double pT0 = SPinterp_pT[npt-1];
	double pT1 = SPinterp_pT[npt];

	// locate pphi interval
	int nphi = 1, nphim1 = 0;
	if(phir < SPinterp_pphi[0])			//if angle is less than minimum angle grid point
	{
		phi0 = SPinterp_pphi[npphi_max] - 2. * M_PI;
		phi1 = SPinterp_pphi[0];
		nphi = 0;
		nphim1 = npphi_max;
	}
	else if(phir > SPinterp_pphi[npphi_max])	//if angle is greater than maximum angle grid point
	{
		phi0 = SPinterp_pphi[npphi_max];
		phi1 = SPinterp_pphi[0] + 2. * M_PI;
		nphi = 0;
		nphim1 = npphi_max;
	}
	else						//if angle is within grid range
	{
		while ((phir > SPinterp_pphi[nphi]) &&
				(nphi < npphi_max)) ++nphi;
		nphim1 = nphi - 1;
		phi0 = SPinterp_pphi[nphim1];
		phi1 = SPinterp_pphi[nphi];
	}

	if (pT0==pT1 || phi0==phi1)
	{
		cerr << "ERROR in Edndp3(): pT and/or pphi values equal!" << endl;
		exit(1);
	}

	double one_by_pTdiff = 1./(pT1 - pT0), one_by_pphidiff = 1./(phi1 - phi0);

	// choose pt-pphi slice of resonance info arrays
	double * sign_of_f11_arr = res_sign_info[npt-1][nphim1];
	double * sign_of_f12_arr = res_sign_info[npt-1][nphi];
	double * sign_of_f21_arr = res_sign_info[npt][nphim1];
	double * sign_of_f22_arr = res_sign_info[npt][nphi];

	double * log_f11_arr = res_log_info[npt-1][nphim1];
	double * log_f12_arr = res_log_info[npt-1][nphi];
	double * log_f21_arr = res_log_info[npt][nphim1];
	double * log_f22_arr = res_log_info[npt][nphi];

	double * f11_arr = res_moments_info[npt-1][nphim1];
	double * f12_arr = res_moments_info[npt-1][nphi];
	double * f21_arr = res_moments_info[npt][nphim1];
	double * f22_arr = res_moments_info[npt][nphi];

	// set index for looping
	int qpt_cs_idx = 0;

	for (int iqt = 0; iqt < qnpts; ++iqt)
	for (int iqx = 0; iqx < qnpts; ++iqx)
	for (int iqy = 0; iqy < qnpts; ++iqy)
	for (int iqz = 0; iqz < qnpts; ++iqz)
	for (int itrig = 0; itrig < 2; ++itrig)
	{
		// interpolate over pT values first
		if(ptr > PTCHANGE)				// if pT interpolation point is larger than PTCHANGE (currently 1.0 GeV)
		{
			double sign_of_f11 = sign_of_f11_arr[qpt_cs_idx];
			double sign_of_f12 = sign_of_f12_arr[qpt_cs_idx];
			double sign_of_f21 = sign_of_f21_arr[qpt_cs_idx];
			double sign_of_f22 = sign_of_f22_arr[qpt_cs_idx];
			
			//*******************************************************************************************************************
			// set f1 first
			//*******************************************************************************************************************
			// if using extrapolation and spectra at pT1 has larger magnitude than at pT0, just return zero
			if (ptr > pT1 && log_f21_arr[qpt_cs_idx] > log_f11_arr[qpt_cs_idx])
				f1 = 0.0;
			else if (sign_of_f11 * sign_of_f21 > 0)	// if the two points have the same sign in the pT direction, interpolate logs
				f1 = sign_of_f11 * exp( lin_int(ptr-pT0, one_by_pTdiff, log_f11_arr[qpt_cs_idx], log_f21_arr[qpt_cs_idx]) );
			else					// otherwise, just interpolate original vals
				f1 = lin_int(ptr-pT0, one_by_pTdiff, f11_arr[qpt_cs_idx], f21_arr[qpt_cs_idx]);
				
			//*******************************************************************************************************************
			// set f2 next
			//*******************************************************************************************************************
			if (ptr > pT1 && log_f22_arr[qpt_cs_idx] > log_f12_arr[qpt_cs_idx])
				f2 = 0.0;
			else if (sign_of_f12 * sign_of_f22 > 0)	// if the two points have the same sign in the pT direction, interpolate logs
				f2 = sign_of_f12 * exp( lin_int(ptr-pT0, one_by_pTdiff, log_f12_arr[qpt_cs_idx], log_f22_arr[qpt_cs_idx]) );
			else					// otherwise, just interpolate original vals
				f2 = lin_int(ptr-pT0, one_by_pTdiff, f12_arr[qpt_cs_idx], f22_arr[qpt_cs_idx]);
			//*******************************************************************************************************************
		}
		else						// if pT is smaller than PTCHANGE, just use linear interpolation, no matter what
		{
			f1 = lin_int(ptr-pT0, one_by_pTdiff, f11_arr[qpt_cs_idx], f21_arr[qpt_cs_idx]);
			f2 = lin_int(ptr-pT0, one_by_pTdiff, f12_arr[qpt_cs_idx], f22_arr[qpt_cs_idx]);
		}
					
		// now, interpolate f1 and f2 over the pphi direction
		results[qpt_cs_idx] += lin_int(phir-phi0, one_by_pphidiff, f1, f2);
					
		if ( isnan( results[qpt_cs_idx] ) )
		{
			*global_out_stream_ptr << "ERROR in Edndp3(double, double, double*): problems encountered!" << endl
				<< "results(" << iqt << "," << iqx << "," << iqy << "," << iqz << "," << itrig << ") = "
				<< setw(8) << setprecision(15) << results[qpt_cs_idx] << endl
				<< "  --> ptr = " << ptr << endl
				<< "  --> pt0 = " << pT0 << endl
				<< "  --> pt1 = " << pT1 << endl
				<< "  --> phir = " << phir << endl
				<< "  --> phi0 = " << phi0 << endl
				<< "  --> phi1 = " << phi1 << endl
				<< "  --> f11 = " << f11_arr[qpt_cs_idx] << endl
				<< "  --> f12 = " << f12_arr[qpt_cs_idx] << endl
				<< "  --> f21 = " << f21_arr[qpt_cs_idx] << endl
				<< "  --> f22 = " << f22_arr[qpt_cs_idx] << endl
				<< "  --> f1 = " << f1 << endl
				<< "  --> f2 = " << f2 << endl;
							exit(1);
		}
		++qpt_cs_idx;	// step to next cell in results array
	}	// ending all loops at once in linearized version

	return;
}

//End of file
