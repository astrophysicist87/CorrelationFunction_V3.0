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
	current_total_resonance_percentage = 0.0;
	all_particles = all_particles_in;
	for (int icr = 0; icr < (int)chosen_resonances_in.size(); icr++)
		chosen_resonances.push_back(chosen_resonances_in[icr]);
	thermal_pions_only = false;
	Nparticle = Nparticle_in;
	read_in_all_dN_dypTdpTdphi = false;
	output_all_dN_dypTdpTdphi = true;
	currentfolderindex = -1;
	current_level_of_output = 0;
	//qspace_cs_slice_length = qnpts*qnpts*qnpts*qnpts*2;		//factor of 2 for sin or cos
	qspace_cs_slice_length = qtnpts*qxnpts*qynpts*qznpts*2;		//factor of 2 for sin or cos

	gsl_set_error_handler_off();

	switch(PC_MARKER_SPACING)
	{
		case 0:
			number_of_percentage_markers = 101;
			break;
		case 1:
			number_of_percentage_markers = UDPMsize;
			break;
		case 2:
			number_of_percentage_markers = UDPMTsize;
			break;
		default:
			break;
	}

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
//debugger(__LINE__, __FILE__);
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

	thermal_target_dN_dypTdpTdphi_moments = new double ****** [n_interp_pT_pts];
	current_dN_dypTdpTdphi_moments = new double ****** [n_interp_pT_pts];
	current_ln_dN_dypTdpTdphi_moments = new double ****** [n_interp_pT_pts];
	current_sign_of_dN_dypTdpTdphi_moments = new double ****** [n_interp_pT_pts];
	for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
	{
		thermal_target_dN_dypTdpTdphi_moments[ipt] = new double ***** [n_interp_pphi_pts];
		current_dN_dypTdpTdphi_moments[ipt] = new double ***** [n_interp_pphi_pts];
		current_ln_dN_dypTdpTdphi_moments[ipt] = new double ***** [n_interp_pphi_pts];
		current_sign_of_dN_dypTdpTdphi_moments[ipt] = new double ***** [n_interp_pphi_pts];
		for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
		{
			thermal_target_dN_dypTdpTdphi_moments[ipt][ipphi] = new double **** [qtnpts];
			current_dN_dypTdpTdphi_moments[ipt][ipphi] = new double **** [qtnpts];
			current_ln_dN_dypTdpTdphi_moments[ipt][ipphi] = new double **** [qtnpts];
			current_sign_of_dN_dypTdpTdphi_moments[ipt][ipphi] = new double **** [qtnpts];
			for (int iqt = 0; iqt < qtnpts; ++iqt)
			{
				thermal_target_dN_dypTdpTdphi_moments[ipt][ipphi][iqt] = new double *** [qxnpts];
				current_dN_dypTdpTdphi_moments[ipt][ipphi][iqt] = new double *** [qxnpts];
				current_ln_dN_dypTdpTdphi_moments[ipt][ipphi][iqt] = new double *** [qxnpts];
				current_sign_of_dN_dypTdpTdphi_moments[ipt][ipphi][iqt] = new double *** [qxnpts];
				for (int iqx = 0; iqx < qxnpts; ++iqx)
				{
					thermal_target_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx] = new double ** [qynpts];
					current_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx] = new double ** [qynpts];
					current_ln_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx] = new double ** [qynpts];
					current_sign_of_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx] = new double ** [qynpts];
					for (int iqy = 0; iqy < qynpts; ++iqy)
					{
						thermal_target_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx][iqy] = new double * [qznpts];
						current_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx][iqy] = new double * [qznpts];
						current_ln_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx][iqy] = new double * [qznpts];
						current_sign_of_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx][iqy] = new double * [qznpts];
						for (int iqz = 0; iqz < qznpts; ++iqz)
						{
							thermal_target_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx][iqy][iqz] = new double [2];
							current_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx][iqy][iqz] = new double [2];
							current_ln_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx][iqy][iqz] = new double [2];
							current_sign_of_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx][iqy][iqz] = new double [2];
							for (int itrig = 0; itrig < 2; ++itrig)
							{
								thermal_target_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx][iqy][iqz][itrig] = 0.0;
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

	qlist = new double ** [n_interp_pT_pts];
	for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
	{
		int qidx = 0;
		qlist[ipt] = new double * [qtnpts*qxnpts*qynpts*qznpts];
		for (int iqt = 0; iqt < qtnpts; ++iqt)
		for (int iqx = 0; iqx < qxnpts; ++iqx)
		for (int iqy = 0; iqy < qynpts; ++iqy)
		for (int iqz = 0; iqz < qznpts; ++iqz)
			qlist[ipt][qidx++] = new double [4];
	}
	int qidx = 0;
	current_qlist_slice = new double * [qtnpts*qxnpts*qynpts*qznpts];
	for (int iqt = 0; iqt < qtnpts; ++iqt)
	for (int iqx = 0; iqx < qxnpts; ++iqx)
	for (int iqy = 0; iqy < qynpts; ++iqy)
	for (int iqz = 0; iqz < qznpts; ++iqz)
	{
		current_qlist_slice[qidx] = new double [4];
		for (int i = 0; i < 4; ++i)
			current_qlist_slice[qidx][i] = 0.0;
		qidx++;
	}

/*
	spectra_to_subtract = new double * [n_interp_pT_pts];
	for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
	{
		spectra_to_subtract[ipt] = new double [n_interp_pphi_pts];
		for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
			spectra_to_subtract[ipt][ipphi] = 0.0;
	}
*/

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
	thermal_spectra = new double ** [Nparticle];
	log_spectra = new double ** [Nparticle];
	sign_spectra = new double ** [Nparticle];
	for (int ir = 0; ir < Nparticle; ++ir)
	{
		spectra[ir] = new double * [n_interp_pT_pts];
		abs_spectra[ir] = new double * [n_interp_pT_pts];
		thermal_spectra[ir] = new double * [n_interp_pT_pts];
		log_spectra[ir] = new double * [n_interp_pT_pts];
		sign_spectra[ir] = new double * [n_interp_pT_pts];
		for (int ipT = 0; ipT < n_interp_pT_pts; ++ipT)
		{
			spectra[ir][ipT] = new double [n_interp_pphi_pts];
			abs_spectra[ir][ipT] = new double [n_interp_pphi_pts];
			thermal_spectra[ir][ipT] = new double [n_interp_pphi_pts];
			log_spectra[ir][ipT] = new double [n_interp_pphi_pts];
			sign_spectra[ir][ipT] = new double [n_interp_pphi_pts];
			for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
			{
				spectra[ir][ipT][ipphi] = 0.0;
				abs_spectra[ir][ipT][ipphi] = 0.0;
				thermal_spectra[ir][ipT][ipphi] = 0.0;
				log_spectra[ir][ipT][ipphi] = 0.0;
				sign_spectra[ir][ipT][ipphi] = 0.0;
			}
		}
	}

	// used for keeping track of how many FO cells are important for given pT, pphi
	// also set up q-space cutoffs array
	number_of_FOcells_above_cutoff_array = new int * [n_interp_pT_pts];
	current_q_space_cutoff = new double * [n_interp_pT_pts];
	correlator_minus_one_cutoff_norms = new int ** [n_interp_pT_pts];
	for (int ipT = 0; ipT < n_interp_pT_pts; ++ipT)
	{
		number_of_FOcells_above_cutoff_array[ipT] = new int [n_interp_pphi_pts];
		current_q_space_cutoff[ipT] = new double [n_interp_pphi_pts];
		correlator_minus_one_cutoff_norms[ipT] = new int * [n_interp_pphi_pts];
		for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
		{
			correlator_minus_one_cutoff_norms[ipT][ipphi] = new int [4];
			for (int ii = 0; ii < 4; ++ii)
				correlator_minus_one_cutoff_norms[ipT][ipphi][ii] = qtnpts*qtnpts+qxnpts*qxnpts+qynpts*qynpts+qznpts*qznpts;
				// i.e., something larger than maximum q-array ranges, only made smaller if correlator cutoff threshhold reached, making large q-values redundant
			number_of_FOcells_above_cutoff_array[ipT][ipphi] = 0;
			current_q_space_cutoff[ipT][ipphi] = 0.0;
		}
	}

	// set-up integration points for resonance integrals
	v_pts = new double [n_v_pts];
	v_wts = new double [n_v_pts];
	zeta_pts = new double [n_zeta_pts];
	zeta_wts = new double [n_zeta_pts];
	//initialize all gaussian points for resonance integrals
	//syntax: int gauss_quadrature(int order, int kind, double alpha, double beta, double a, double b, double x[], double w[])
	gauss_quadrature(n_zeta_pts, 1, 0.0, 0.0, zeta_min, zeta_max, zeta_pts, zeta_wts);
	gauss_quadrature(n_v_pts, 1, 0.0, 0.0, v_min, v_max, v_pts, v_wts);
   //single particle spectra for plane angle determination
   /*SP_pT = new double [n_SP_pT];
   SP_pT_weight = new double [n_SP_pT];
   gauss_quadrature(n_SP_pT, 1, 0.0, 0.0, SP_pT_min, SP_pT_max, SP_pT, SP_pT_weight);
   SP_pphi = new double [n_SP_pphi];
   SP_pphi_weight = new double [n_SP_pphi];
   gauss_quadrature(n_SP_pphi, 1, 0.0, 0.0, 0.0, 2.*M_PI, SP_pphi, SP_pphi_weight);*/
   SP_p_y = 0.0e0;
//initialize and set evenly spaced grid of px-py points in transverse plane,
//and corresponding p0 and pz points
	SPinterp_pT = new double [n_interp_pT_pts];
	SPinterp_pT_public = new double [n_interp_pT_pts];
	SPinterp_pphi = new double [n_interp_pphi_pts];
	sin_SPinterp_pphi = new double [n_interp_pphi_pts];
	cos_SPinterp_pphi = new double [n_interp_pphi_pts];
	SPinterp_p0 = new double * [n_interp_pT_pts];
	SPinterp_pz = new double * [n_interp_pT_pts];
	for(int ipt=0; ipt<n_interp_pT_pts; ipt++)
	{
		SPinterp_p0[ipt] = new double [eta_s_npts];
		SPinterp_pz[ipt] = new double [eta_s_npts];
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
		//for(int ipt=0; ipt<n_interp_pT_pts; ipt++)
		//	cout << "PointCheck, pT: " << scientific << setprecision(17) << setw(20) << SPinterp_pT[ipt] << "   " << dummywts3[ipt] << endl;
		gauss_quadrature(n_interp_pphi_pts, 1, 0.0, 0.0, interp_pphi_min, interp_pphi_max, SPinterp_pphi, dummywts4);
		//for(int ipphi=0; ipphi<n_interp_pphi_pts; ipphi++)
		//	cout << "PointCheck, pphi: " << scientific << setprecision(17) << setw(20) << SPinterp_pphi[ipphi] << "   " << dummywts4[ipphi] << endl;
		for(int ipphi=0; ipphi<n_interp_pphi_pts; ipphi++)
		{
			//SPinterp_pphi[ipphi] = interp_pphi_min + (double)ipphi*Del2_pphi;
			sin_SPinterp_pphi[ipphi] = sin(SPinterp_pphi[ipphi]);
			cos_SPinterp_pphi[ipphi] = cos(SPinterp_pphi[ipphi]);
		}
	}

//debugger(__LINE__, __FILE__);

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

	//for (int ieta = 0; ieta < eta_s_npts; ieta++)
	//	cout << "PointCheck, eta_s: " << scientific << setprecision(17) << setw(20) << eta_s[ieta] << "   " << 2.*eta_s_weight[ieta] << endl;


//debugger(__LINE__, __FILE__);
	R2_side = new double * [n_interp_pT_pts];
	R2_out = new double * [n_interp_pT_pts];
	R2_long = new double * [n_interp_pT_pts];
	R2_outside = new double * [n_interp_pT_pts];
	R2_sidelong = new double * [n_interp_pT_pts];
	R2_outlong = new double * [n_interp_pT_pts];

	R2_side_err = new double * [n_interp_pT_pts];
	R2_out_err = new double * [n_interp_pT_pts];
	R2_long_err = new double * [n_interp_pT_pts];
	R2_outside_err = new double * [n_interp_pT_pts];
	R2_sidelong_err = new double * [n_interp_pT_pts];
	R2_outlong_err = new double * [n_interp_pT_pts];

	lambda_Correl = new double * [n_interp_pT_pts];
	lambda_Correl_err = new double * [n_interp_pT_pts];

	for(int ipt=0; ipt<n_interp_pT_pts; ipt++)
	{
		R2_side[ipt] = new double [n_interp_pphi_pts];
		R2_out[ipt] = new double [n_interp_pphi_pts];
		R2_outside[ipt] = new double [n_interp_pphi_pts];
		R2_long[ipt] = new double [n_interp_pphi_pts];
		R2_sidelong[ipt] = new double [n_interp_pphi_pts];
		R2_outlong[ipt] = new double [n_interp_pphi_pts];

		R2_side_err[ipt] = new double [n_interp_pphi_pts];
		R2_out_err[ipt] = new double [n_interp_pphi_pts];
		R2_long_err[ipt] = new double [n_interp_pphi_pts];
		R2_outside_err[ipt] = new double [n_interp_pphi_pts];
		R2_sidelong_err[ipt] = new double [n_interp_pphi_pts];
		R2_outlong_err[ipt] = new double [n_interp_pphi_pts];

		lambda_Correl[ipt] = new double [n_interp_pphi_pts];
		lambda_Correl_err[ipt] = new double [n_interp_pphi_pts];
	}
//debugger(__LINE__, __FILE__);
	//initialize all source variances and HBT radii/coeffs
	for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
	{
		for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
		{
			R2_side[ipt][ipphi] = 0.;
			R2_out[ipt][ipphi] = 0.;
			R2_long[ipt][ipphi] = 0.;
			R2_outside[ipt][ipphi] = 0.;
			R2_sidelong[ipt][ipphi] = 0.;
			R2_outlong[ipt][ipphi] = 0.;

			R2_side_err[ipt][ipphi] = 0.;
			R2_out_err[ipt][ipphi] = 0.;
			R2_long_err[ipt][ipphi] = 0.;
			R2_outside_err[ipt][ipphi] = 0.;
			R2_sidelong_err[ipt][ipphi] = 0.;
			R2_outlong_err[ipt][ipphi] = 0.;

			lambda_Correl[ipt][ipphi] = 0.0;
			lambda_Correl_err[ipt][ipphi] = 0.0;
		}
	}

   return;
}

void CorrelationFunction::Update_sourcefunction(particle_info* particle, int FOarray_length, int particle_idx)
{
	full_FO_length = FOarray_length * eta_s_npts;

	eiqtt = new double * [qtnpts];					//to hold cos/sin(q0 t)
	eiqxx = new double * [qxnpts];					//to hold cos/sin(qx x)
	eiqyy = new double * [qynpts];					//to hold cos/sin(qy y)
	eiqzz = new double * [qznpts];					//to hold cos/sin(+/- qz z)

	for (int iqt = 0; iqt < qtnpts; ++iqt)
		eiqtt[iqt] = new double [FOarray_length * eta_s_npts * 2];
	for (int iqx = 0; iqx < qxnpts; ++iqx)
		eiqxx[iqx] = new double [FOarray_length * 2];
	for (int iqy = 0; iqy < qynpts; ++iqy)
		eiqyy[iqy] = new double [FOarray_length * 2];
	for (int iqz = 0; iqz < qznpts; ++iqz)
		eiqzz[iqz] = new double [FOarray_length * eta_s_npts * 2];

	qt_PTdep_pts = new double * [n_interp_pT_pts];
	qx_PTdep_pts = new double * [n_interp_pT_pts];
	qy_PTdep_pts = new double * [n_interp_pT_pts];
	qz_PTdep_pts = new double * [n_interp_pT_pts];

	S_p_withweight_array = new double ** [n_interp_pT_pts];
	for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
		S_p_withweight_array[ipt] = new double * [n_interp_pphi_pts];

   //particle information
   particle_name = particle->name;
   particle_mass = particle->mass;
   particle_sign = particle->sign;
   particle_gspin = particle->gspin;
   particle_id = particle_idx;

	*global_out_stream_ptr << "Inside Update_sourcefunction(...): using fraction_of_resonances = " << fraction_of_resonances << endl;

   FO_length = FOarray_length;

	// set the rest later
	most_important_FOcells = new size_t ** [n_interp_pT_pts];
	for(int ipt=0; ipt<n_interp_pT_pts; ipt++)
	{
		most_important_FOcells[ipt] = new size_t * [n_interp_pphi_pts];
		//for(int ipphi=0; ipphi<n_interp_pphi_pts; ipphi++)
		//	most_important_FOcells[ipt][ipphi] = new size_t [FO_length * eta_s_npts];
	}

   return;
}

void CorrelationFunction::Fill_out_pts(double * pointsarray, int numpoints, double max_val, int spacing_type)
{
	// spacing_type:		0 - uniform spacing
	//						1 - Chebyshev-node based spacing
	if (numpoints == 1)
		pointsarray[0] = 0.0;
	else
	{
		// if I want q-points equally spaced...
		if (spacing_type == 0)
		{
			for (int iqd = 0; iqd < numpoints; ++iqd)
				pointsarray[iqd] = -max_val + (double)iqd * 2.*max_val / double(numpoints - 1+1e-100);
		}
		// else, use Chebyshev nodes instead...
		else if (spacing_type == 1)
		{
			//double local_scale = max_val / cos(M_PI / (2.*qtnpts));
			double local_scale = -max_val;
			for (int iqd = 0; iqd < numpoints; ++iqd)
				pointsarray[iqd] = local_scale * cos( M_PI*(2.*(iqd+1.) - 1.) / (2.*numpoints) );
		}
	}
	return;
}

void CorrelationFunction::Set_q_pTdep_pts(int ipt, double qxw, double qyw, double qzw)
{
	double pT_local = SPinterp_pT[ipt];

	qt_PTdep_pts[ipt] = new double [qtnpts];
	qx_PTdep_pts[ipt] = new double [qxnpts];
	qy_PTdep_pts[ipt] = new double [qynpts];
	qz_PTdep_pts[ipt] = new double [qznpts];

	double mpion = all_particles[target_particle_id].mass;
	double eps = 2.e-1;									//specifies approximate CF value at which to truncate calculation
														// (used for computing q(i)max)
	double ln_one_by_eps = hbarC*sqrt(log(1./eps));

////////////////////////////////////////////////////////////////////////
	double qxmax = ln_one_by_eps / sqrt(qxw*qxw+qyw*qyw);
	double qymax = ln_one_by_eps / sqrt(qxw*qxw+qyw*qyw);
	double qzmax = ln_one_by_eps / qzw;
//double qzmax = qxmax;
	double xi2 = mpion*mpion + pT_local*pT_local + 2.0*0.25*qxmax*qxmax;	//pretend that Kphi == 0, qx == qo and qs == ql == 0, to maximize qtmax
	double qtmax = sqrt(xi2 + sqrt(2.0)*pT_local*qxmax) - sqrt(xi2 - sqrt(2.0)*pT_local*qxmax) + 1.e-10;

//cout << "pt and qmax list: " << pT_local << "   " << qtmax << "   " << qxmax << "   " << qymax << "   " << qzmax << endl;
	//double qtmax = 2.0 * qxmax;
////////////////////////////////////////////////////////////////////////

	//finally, set q(i)_PTdep_pts
	Fill_out_pts(qt_PTdep_pts[ipt], qtnpts, qtmax, QT_POINTS_SPACING);
	Fill_out_pts(qx_PTdep_pts[ipt], qxnpts, qxmax, QX_POINTS_SPACING);
	Fill_out_pts(qy_PTdep_pts[ipt], qynpts, qymax, QY_POINTS_SPACING);
	Fill_out_pts(qz_PTdep_pts[ipt], qznpts, qzmax, QZ_POINTS_SPACING);

	int qidx = 0;
	//qlist[ipt] = new double * [qtnpts*qxnpts*qynpts*qznpts];
	for (int iqt = 0; iqt < qtnpts; ++iqt)
	for (int iqx = 0; iqx < qxnpts; ++iqx)
	for (int iqy = 0; iqy < qynpts; ++iqy)
	for (int iqz = 0; iqz < qznpts; ++iqz)
	//for (int itrig = 0; itrig < 2; ++itrig)
	{
		qlist[ipt][qidx][0] = qt_PTdep_pts[ipt][iqt];
		qlist[ipt][qidx][1] = qx_PTdep_pts[ipt][iqx];
		qlist[ipt][qidx][2] = qy_PTdep_pts[ipt][iqy];
		qlist[ipt][qidx][3] = qz_PTdep_pts[ipt][iqz];
		qidx++;
	}

	return;
}

void CorrelationFunction::Set_eiqx_with_q_pTdep_pts(int ipt)
{
	int iFO = 0;
	double * current_qt_slice = qt_PTdep_pts[ipt];
	double * current_qx_slice = qx_PTdep_pts[ipt];
	double * current_qy_slice = qy_PTdep_pts[ipt];
	double * current_qz_slice = qz_PTdep_pts[ipt];

	for (int isurf = 0; isurf < FO_length; ++isurf)
	{
		FO_surf * surf = &current_FOsurf_ptr[isurf];
		double tau = surf->tau;

		for (int ieta = 0; ieta < eta_s_npts; ++ieta)
		{
			double tpt = tau*ch_eta_s[ieta];
			double zpt = tau*sh_eta_s[ieta];

			for (int iqt = 0; iqt < qtnpts; ++iqt)
			{
				eiqtt[iqt][iFO] = cos(hbarCm1*current_qt_slice[iqt]*tpt);
				eiqtt[iqt][iFO+1] = sin(hbarCm1*current_qt_slice[iqt]*tpt);
			}
			for (int iqz = 0; iqz < qznpts; ++iqz)
			{
				eiqzz[iqz][iFO] = cos(hbarCm1*current_qz_slice[iqz]*zpt);
				eiqzz[iqz][iFO+1] = sin(hbarCm1*current_qz_slice[iqz]*zpt);
			}

			iFO += 2;
		}

		double xpt = surf->xpt;
		double ypt = surf->ypt;
		for (int iqx = 0; iqx < qxnpts; ++iqx)
		{
			eiqxx[iqx][2*isurf] = cos(hbarCm1*current_qx_slice[iqx]*xpt);
			eiqxx[iqx][2*isurf+1] = sin(hbarCm1*current_qx_slice[iqx]*xpt);
		}
		for (int iqy = 0; iqy < qynpts; ++iqy)
		{
			eiqyy[iqy][2*isurf] = cos(hbarCm1*current_qy_slice[iqy]*ypt);
			eiqyy[iqy][2*isurf+1] = sin(hbarCm1*current_qy_slice[iqy]*ypt);
		}
	}

	//dump results to binary files, which is way faster than recalculating everything each time
	Dump_phases_to_binary('t', ipt, eiqtt, qtnpts, 2*FO_length*eta_s_npts);
	Dump_phases_to_binary('x', ipt, eiqxx, qxnpts, 2*FO_length);
	Dump_phases_to_binary('y', ipt, eiqyy, qynpts, 2*FO_length);
	Dump_phases_to_binary('z', ipt, eiqzz, qznpts, 2*FO_length*eta_s_npts);

	return;
}

void CorrelationFunction::Load_eiqx_with_q_pTdep_pts(int ipt)
{
	Load_phases_from_binary('t', ipt, eiqtt, qtnpts, 2*FO_length*eta_s_npts);
	Load_phases_from_binary('x', ipt, eiqxx, qxnpts, 2*FO_length);
	Load_phases_from_binary('y', ipt, eiqyy, qynpts, 2*FO_length);
	Load_phases_from_binary('z', ipt, eiqzz, qznpts, 2*FO_length*eta_s_npts);

	return;
}

CorrelationFunction::~CorrelationFunction()
{
   delete [] K_T;
   delete [] K_phi;
   delete [] K_phi_weight;
   delete [] eta_s;
   delete [] eta_s_weight;

	for(int ipt=0; ipt<n_interp_pT_pts; ipt++)
	{
		delete [] lambda_Correl[ipt];
		delete [] R2_side[ipt];
		delete [] R2_out[ipt];
		delete [] R2_long[ipt];
		delete [] R2_outside[ipt];
		delete [] R2_sidelong[ipt];
		delete [] R2_outlong[ipt];

		delete [] lambda_Correl_err[ipt];
		delete [] R2_side_err[ipt];
		delete [] R2_out_err[ipt];
		delete [] R2_long_err[ipt];
		delete [] R2_outside_err[ipt];
		delete [] R2_sidelong_err[ipt];
		delete [] R2_outlong_err[ipt];
	}

	delete [] R2_side;
	delete [] R2_out;
	delete [] R2_long;
	delete [] R2_outside;
	delete [] R2_sidelong;
	delete [] R2_outlong;

	delete [] lambda_Correl_err;
	delete [] R2_side_err;
	delete [] R2_out_err;
	delete [] R2_long_err;
	delete [] R2_outside_err;
	delete [] R2_sidelong_err;
	delete [] R2_outlong_err;

   return;
}

void CorrelationFunction::Allocate_CFvals()
{
	CFvals = new double **** [n_interp_pT_pts];
	for (int ipT = 0; ipT < n_interp_pT_pts; ++ipT)
	{
		CFvals[ipT] = new double *** [n_interp_pphi_pts];
		for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
		{
			CFvals[ipT][ipphi] = new double ** [qxnpts];
			for (int iqx = 0; iqx < qxnpts; ++iqx)
			{
				CFvals[ipT][ipphi][iqx] = new double * [qynpts];
				for (int iqy = 0; iqy < qynpts; ++iqy)
				{
					CFvals[ipT][ipphi][iqx][iqy] = new double [qznpts];
					for (int iqz = 0; iqz < qznpts; ++iqz)
						CFvals[ipT][ipphi][iqx][iqy][iqz] = 0.0;
				}
			}
		}
	}

	return;
}

void CorrelationFunction::Delete_S_p_withweight_array()
{
	for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
		delete [] S_p_withweight_array[ipt];
	delete [] S_p_withweight_array;

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
	VEC_n2_Ppm = new double *** [n_v_pts];
	for (int iv = 0; iv < n_v_pts; ++iv)
	{
		VEC_n2_MT[iv] = new double [n_zeta_pts];
		VEC_n2_PPhi_tilde[iv] = new double [n_zeta_pts];
		VEC_n2_PPhi_tildeFLIP[iv] = new double [n_zeta_pts];
		VEC_n2_PT[iv] = new double [n_zeta_pts];
		VEC_n2_zeta_factor[iv] = new double [n_zeta_pts];
		VEC_n2_Ppm[iv] = new double ** [n_zeta_pts];
		for (int izeta = 0; izeta < n_zeta_pts; ++izeta)
		{
			VEC_n2_Ppm[iv][izeta] = new double * [2];		//two corresponds to +/-
			for (int i = 0; i < 2; ++i)
				VEC_n2_Ppm[iv][izeta][i] = new double [4];	//four corresponds to space-time components
		}
	}
	s_pts = new double [n_s_pts];
	s_wts = new double [n_s_pts];
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
	VEC_Ppm = new double **** [n_s_pts];
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
		VEC_Ppm[is] = new double *** [n_v_pts];
		for(int iv = 0; iv < n_v_pts; ++iv)
		{
			VEC_MT[is][iv] = new double [n_zeta_pts];
			VEC_PPhi_tilde[is][iv] = new double [n_zeta_pts];
			VEC_PPhi_tildeFLIP[is][iv] = new double [n_zeta_pts];
			VEC_PT[is][iv] = new double [n_zeta_pts];
			VEC_zeta_factor[is][iv] = new double [n_zeta_pts];
			VEC_Ppm[is][iv] = new double ** [n_zeta_pts];
			for (int izeta = 0; izeta < n_zeta_pts; ++izeta)
			{
				VEC_Ppm[is][iv][izeta] = new double * [2];		//two corresponds to +/-
				for (int i = 0; i < 2; ++i)
					VEC_Ppm[is][iv][izeta][i] = new double [4];	//four corresponds to space-time components
			}
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
		for (int izeta = 0; izeta < n_zeta_pts; ++izeta)
		{
			for (int i = 0; i < 2; ++i)
				delete [] VEC_n2_Ppm[iv][izeta][i];
			delete [] VEC_n2_Ppm[iv][izeta];
		}
		delete [] VEC_n2_Ppm[iv];
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
	delete [] VEC_n2_Ppm;

	for(int is = 0; is < n_s_pts; ++is)
	{
		for(int iv = 0; iv < n_v_pts; ++iv)
		{
			delete [] VEC_MT[is][iv];
			delete [] VEC_PPhi_tilde[is][iv];
			delete [] VEC_PPhi_tildeFLIP[is][iv];
			delete [] VEC_PT[is][iv];
			delete [] VEC_zeta_factor[is][iv];
			for (int izeta = 0; izeta < n_zeta_pts; ++izeta)
			{
				for (int i = 0; i < 2; ++i)
					delete [] VEC_Ppm[is][iv][izeta][i];
				delete [] VEC_Ppm[is][iv][izeta];
			}
			delete [] VEC_Ppm[is][iv];
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
		delete [] VEC_Ppm[is];
	}
	delete [] s_pts;
	delete [] s_wts;
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
	delete [] VEC_Ppm;
	if (VERBOSE > 2) *global_out_stream_ptr << "Deleted memory for decay channel information." << endl;

	return;
}

void CorrelationFunction::Set_correlation_function_q_pts()
{
	q1npts = qxnpts;
	q2npts = qynpts;
	q3npts = qznpts;

	qx_pts = new double [qxnpts];
	qy_pts = new double [qynpts];
	qz_pts = new double [qznpts];
	
	// initialize error matrix
	Correl_3D_err = new double ** [qxnpts];
	for (int iqx = 0; iqx < qxnpts; ++iqx)
	{
		Correl_3D_err[iqx] = new double * [qynpts];
		for (int iqy = 0; iqy < qynpts; ++iqy)
		{
			Correl_3D_err[iqx][iqy] = new double [qznpts];
			for (int iqz = 0; iqz < qznpts; ++iqz)
				Correl_3D_err[iqx][iqy][iqz] = 1e-3;	//naive choice for now
		}
	}

	return;
}


// returns points in q-space for computing weighted spectra grid corresponding to to given q and K choices
// weighted spectra grid thus needs to be interpolated at point returned in qgridpts
void CorrelationFunction::Get_q_points(double q1, double q2, double q3, double pT, double pphi, double * qgridpts)
{
	double mtarget = all_particles[target_particle_id].mass;
	double xi2 = mtarget*mtarget + pT*pT + 0.25*(q1*q1 + q2*q2 + q3*q3);
	double ckp = cos(pphi), skp = sin(pphi);

	double qo = ckp * q1 + skp * q2;
		
	// set qpts at which to interpolate spectra
	qgridpts[0] = sqrt(xi2 + qo*pT) - sqrt(xi2 - qo*pT);	//set qt component
	qgridpts[1] = q1;										//set qx component
	qgridpts[2] = q2;										//set qy component
	qgridpts[3] = q3;										//set qz component, since qz = ql

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
				current_daughters_dN_dypTdpTdphi_moments[id][ipt][ipphi] = new double **** [qtnpts];
				current_daughters_ln_dN_dypTdpTdphi_moments[id][ipt][ipphi] = new double **** [qtnpts];
				current_daughters_sign_of_dN_dypTdpTdphi_moments[id][ipt][ipphi] = new double **** [qtnpts];
				for (int iqt = 0; iqt < qtnpts; ++iqt)
				{
					current_daughters_dN_dypTdpTdphi_moments[id][ipt][ipphi][iqt] = new double *** [qxnpts];
					current_daughters_ln_dN_dypTdpTdphi_moments[id][ipt][ipphi][iqt] = new double *** [qxnpts];
					current_daughters_sign_of_dN_dypTdpTdphi_moments[id][ipt][ipphi][iqt] = new double *** [qxnpts];
					for (int iqx = 0; iqx < qxnpts; ++iqx)
					{
						current_daughters_dN_dypTdpTdphi_moments[id][ipt][ipphi][iqt][iqx] = new double ** [qynpts];
						current_daughters_ln_dN_dypTdpTdphi_moments[id][ipt][ipphi][iqt][iqx] = new double ** [qynpts];
						current_daughters_sign_of_dN_dypTdpTdphi_moments[id][ipt][ipphi][iqt][iqx] = new double ** [qynpts];
						for (int iqy = 0; iqy < qynpts; ++iqy)
						{
							current_daughters_dN_dypTdpTdphi_moments[id][ipt][ipphi][iqt][iqx][iqy] = new double * [qznpts];
							current_daughters_ln_dN_dypTdpTdphi_moments[id][ipt][ipphi][iqt][iqx][iqy] = new double * [qznpts];
							current_daughters_sign_of_dN_dypTdpTdphi_moments[id][ipt][ipphi][iqt][iqx][iqy] = new double * [qznpts];
							for (int iqz = 0; iqz < qznpts; ++iqz)
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
			for(int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
			{
				for (int iqt = 0; iqt < qtnpts; ++iqt)
				{
					for (int iqx = 0; iqx < qxnpts; ++iqx)
					{
						for (int iqy = 0; iqy < qynpts; ++iqy)
						{
							for (int iqz = 0; iqz < qznpts; ++iqz)
							{
								delete [] current_daughters_dN_dypTdpTdphi_moments[id][ipt][ipphi][iqt][iqx][iqy][iqz];
								delete [] current_daughters_ln_dN_dypTdpTdphi_moments[id][ipt][ipphi][iqt][iqx][iqy][iqz];
								delete [] current_daughters_sign_of_dN_dypTdpTdphi_moments[id][ipt][ipphi][iqt][iqx][iqy][iqz];
							}
							delete [] current_daughters_dN_dypTdpTdphi_moments[id][ipt][ipphi][iqt][iqx][iqy];
							delete [] current_daughters_ln_dN_dypTdpTdphi_moments[id][ipt][ipphi][iqt][iqx][iqy];
							delete [] current_daughters_sign_of_dN_dypTdpTdphi_moments[id][ipt][ipphi][iqt][iqx][iqy];
						}
						delete [] current_daughters_dN_dypTdpTdphi_moments[id][ipt][ipphi][iqt][iqx];
						delete [] current_daughters_ln_dN_dypTdpTdphi_moments[id][ipt][ipphi][iqt][iqx];
						delete [] current_daughters_sign_of_dN_dypTdpTdphi_moments[id][ipt][ipphi][iqt][iqx];
					}
					delete [] current_daughters_dN_dypTdpTdphi_moments[id][ipt][ipphi][iqt];
					delete [] current_daughters_ln_dN_dypTdpTdphi_moments[id][ipt][ipphi][iqt];
					delete [] current_daughters_sign_of_dN_dypTdpTdphi_moments[id][ipt][ipphi][iqt];
				}
				delete [] current_daughters_dN_dypTdpTdphi_moments[id][ipt][ipphi];
				delete [] current_daughters_ln_dN_dypTdpTdphi_moments[id][ipt][ipphi];
				delete [] current_daughters_sign_of_dN_dypTdpTdphi_moments[id][ipt][ipphi];
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

void CorrelationFunction::Set_target_pphiavgd_CFs()
{
	double N = (double)n_interp_pphi_pts;
	double factor1 = 1./N;
	double factor2 = 2./(N*(N-1.));

	target_pphiavgd_CFs = new double **** [n_interp_pT_pts];
	target_pphivar_CFs = new double **** [n_interp_pT_pts];
	for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
	{
		target_pphiavgd_CFs[ipt] = new double *** [qtnpts];
		target_pphivar_CFs[ipt] = new double *** [qtnpts];
		for (int iqt = 0; iqt < qtnpts; ++iqt)
		{
			target_pphiavgd_CFs[ipt][iqt] = new double ** [qxnpts];
			target_pphivar_CFs[ipt][iqt] = new double ** [qxnpts];
			for (int iqx = 0; iqx < qxnpts; ++iqx)
			{
				target_pphiavgd_CFs[ipt][iqt][iqx] = new double * [qynpts];
				target_pphivar_CFs[ipt][iqt][iqx] = new double * [qynpts];
				for (int iqy = 0; iqy < qynpts; ++iqy)
				{
					target_pphiavgd_CFs[ipt][iqt][iqx][iqy] = new double [qznpts];
					target_pphivar_CFs[ipt][iqt][iqx][iqy] = new double [qznpts];
					for (int iqz = 0; iqz < qznpts; ++iqz)
					{
						target_pphiavgd_CFs[ipt][iqt][iqx][iqy][iqz] = 0.0;
						target_pphivar_CFs[ipt][iqt][iqx][iqy][iqz] = 0.0;
						double tmpsum = 0.0, tmp2sum = 0.0;
						for(int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
						{
							double tmpC = current_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx][iqy][iqz][0];
							double tmpS = current_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx][iqy][iqz][1];
							double tmpspec = spectra[target_particle_id][ipt][ipphi];
							double tmp = 1.0 + (tmpC*tmpC+tmpS*tmpS) / (tmpspec*tmpspec);
							target_pphiavgd_CFs[ipt][iqt][iqx][iqy][iqz] += tmp;
							target_pphivar_CFs[ipt][iqt][iqx][iqy][iqz] -= tmp*tmpsum;
							tmpsum += tmp;
							tmp2sum += tmp*tmp;
						}
						target_pphiavgd_CFs[ipt][iqt][iqx][iqy][iqz] *= factor1;			//set pphi-averaged CF
						target_pphivar_CFs[ipt][iqt][iqx][iqy][iqz] *= factor2;
						target_pphivar_CFs[ipt][iqt][iqx][iqy][iqz] += factor1*tmp2sum;		//set "pphi-varianced" CF
					}
				}
			}
		}
	}

	return;
}


inline double CorrelationFunction::lin_int(double x_m_x1, double one_by_x2_m_x1, double f1, double f2)
{
	return ( f1 + (f2 - f1) * x_m_x1 * one_by_x2_m_x1 );
}

//dots 4-vectors together
inline double CorrelationFunction::dot_four_vectors(double * a, double * b)
{
	double sum = a[0]*b[0];
	for (size_t i=1; i<4; ++i) sum -= a[i]*b[i];
	return (sum);
}

void CorrelationFunction::Edndp3(double ptr, double phir, double * result, int loc_verb /*==0*/)
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
	double log_f11 = spec_log_info[npt-1][nphim1];
	double log_f12 = spec_log_info[npt-1][nphi];
	double log_f21 = spec_log_info[npt][nphim1];
	double log_f22 = spec_log_info[npt][nphi];

	double f11 = spec_vals_info[npt-1][nphim1];
	double f12 = spec_vals_info[npt-1][nphi];
	double f21 = spec_vals_info[npt][nphim1];
	double f22 = spec_vals_info[npt][nphi];

	//double sign_of_f11 = spec_sign_info[npt-1][nphim1];
	//double sign_of_f12 = spec_sign_info[npt-1][nphi];
	//double sign_of_f21 = spec_sign_info[npt][nphim1];
	//double sign_of_f22 = spec_sign_info[npt][nphi];

	/////////////////////////////////////////////////////////////////
	// interpolate over pT values first
	/////////////////////////////////////////////////////////////////
	if(ptr > PTCHANGE)				// if pT interpolation point is larger than PTCHANGE (currently 1.0 GeV)
	{
		double sign_of_f11 = spec_sign_info[npt-1][nphim1];
		double sign_of_f12 = spec_sign_info[npt-1][nphi];
		double sign_of_f21 = spec_sign_info[npt][nphim1];
		double sign_of_f22 = spec_sign_info[npt][nphi];
		
		//*******************************************************************************************************************
		// set f1 first
		//*******************************************************************************************************************
		// if using extrapolation and spectra at pT1 has larger magnitude than at pT0 (or the signs are different), just return zero
		if ( ptr > pT1 && ( log_f21 > log_f11 || sign_of_f11 * sign_of_f21 < 0 ) )
		{
			f1 = 0.0;
			//if ( loc_verb ) *global_out_stream_ptr << "Chose branch 1Aa!" << endl;
		}
		else if (sign_of_f11 * sign_of_f21 > 0)	// if the two points have the same sign in the pT direction, interpolate logs
		{
			f1 = sign_of_f11 * exp( lin_int(ptr-pT0, one_by_pTdiff, log_f11, log_f21) );
			//if ( loc_verb ) *global_out_stream_ptr << "Chose branch 1Ba!" << endl;
		}
		else					// otherwise, just interpolate original vals
		{
			f1 = lin_int(ptr-pT0, one_by_pTdiff, f11, f21);
			//if ( loc_verb ) *global_out_stream_ptr << "Chose branch 1Ca!" << endl;
		}
			
		//*******************************************************************************************************************
		// set f2 next
		//*******************************************************************************************************************
		if ( ptr > pT1 && ( log_f22 > log_f12 || sign_of_f12 * sign_of_f22 < 0 ) )
		{
			f2 = 0.0;
			//if ( loc_verb ) *global_out_stream_ptr << "Chose branch 1Ab!" << endl;
		}
		else if (sign_of_f12 * sign_of_f22 > 0)	// if the two points have the same sign in the pT direction, interpolate logs
		{
			f2 = sign_of_f12 * exp( lin_int(ptr-pT0, one_by_pTdiff, log_f12, log_f22) );
			//if ( loc_verb ) *global_out_stream_ptr << "Chose branch 1Bb!" << endl;
		}
		else					// otherwise, just interpolate original vals
		{
			f2 = lin_int(ptr-pT0, one_by_pTdiff, f12, f22);
			//if ( loc_verb ) *global_out_stream_ptr << "Chose branch 1Cb!" << endl;
		}
		//*******************************************************************************************************************
	}
	else						// if pT is smaller than PTCHANGE, just use linear interpolation, no matter what
	{
		f1 = lin_int(ptr-pT0, one_by_pTdiff, f11, f21);
		f2 = lin_int(ptr-pT0, one_by_pTdiff, f12, f22);
		//if ( loc_verb ) *global_out_stream_ptr << "Chose branch 2!" << endl;
	}
				
	// now, interpolate f1 and f2 over the pphi direction
	*result += lin_int(phir-phi0, one_by_pphidiff, f1, f2);
	//double tmp = lin_int(phir-phi0, one_by_pphidiff, f1, f2);
	//*result += tmp;

/*if ( loc_verb )
		{
			*global_out_stream_ptr << "INFODUMP in Edndp3(double, double, double*):" << endl
				<< setw(25) << setprecision(20) 
				<< "  --> result = " << *result << endl
				<< "  --> ptr = " << ptr << endl
				<< "  --> pt0 = " << pT0 << endl
				<< "  --> pt1 = " << pT1 << endl
				<< "  --> phir = " << phir << endl
				<< "  --> phi0 = " << phi0 << endl
				<< "  --> phi1 = " << phi1 << endl
				<< "  --> f11 = " << f11 << endl
				<< "  --> f12 = " << f12 << endl
				<< "  --> f21 = " << f21 << endl
				<< "  --> f22 = " << f22 << endl
				<< "  --> lf11 = " << log_f11 << endl
				<< "  --> lf12 = " << log_f12 << endl
				<< "  --> lf21 = " << log_f21 << endl
				<< "  --> lf22 = " << log_f22 << endl
				<< "  --> sf11 = " << sign_of_f11 << endl
				<< "  --> sf12 = " << sign_of_f12 << endl
				<< "  --> sf21 = " << sign_of_f21 << endl
				<< "  --> sf22 = " << sign_of_f22 << endl
				<< "  --> f1 = " << f1 << endl
				<< "  --> f2 = " << f2 << endl
				<< "  --> tmp = " << tmp << endl;
		}*/

	//return (result);
	return;
	//return ( lin_int(phir-phi0, one_by_pphidiff, f1, f2) );
}


void CorrelationFunction::eiqxEdndp3(double ptr, double phir, double * results, int loc_verb /*==0*/)
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
		cerr << "ERROR in eiqxEdndp3(): pT and/or pphi values equal!" << endl;
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
	int qlist_idx = 0;

	for (int iqt = 0; iqt < qtnpts; ++iqt)
	for (int iqx = 0; iqx < qxnpts; ++iqx)
	for (int iqy = 0; iqy < qynpts; ++iqy)
	for (int iqz = 0; iqz < qznpts; ++iqz)
	//for (int itrig = 0; itrig < 2; ++itrig)
	{
		double arg = one_by_Gamma_Mres * dot_four_vectors(current_qlist_slice[qlist_idx], currentPpm);
		double akr = 1./(1.+arg*arg);
		double aki = arg/(1.+arg*arg);

		/////////////////////////////////////////////////////////////////
		// DO COSINE PART FIRST
		/////////////////////////////////////////////////////////////////
		// interpolate over pT values first
		/////////////////////////////////////////////////////////////////
		if(ptr > PTCHANGE)				// if pT interpolation point is larger than PTCHANGE (currently 1.0 GeV)
		{
			double sign_of_f11 = sign_of_f11_arr[qpt_cs_idx];
			double sign_of_f12 = sign_of_f12_arr[qpt_cs_idx];
			double sign_of_f21 = sign_of_f21_arr[qpt_cs_idx];
			double sign_of_f22 = sign_of_f22_arr[qpt_cs_idx];
			
			//*******************************************************************************************************************
			// set f1 first
			//*******************************************************************************************************************
			// if using extrapolation and spectra at pT1 has larger magnitude than at pT0 (or the signs are different), just return zero
			if ( ptr > pT1 && ( log_f21_arr[qpt_cs_idx] > log_f11_arr[qpt_cs_idx] || sign_of_f11 * sign_of_f21 < 0 ) )
				f1 = 0.0;
			else if (sign_of_f11 * sign_of_f21 > 0)	// if the two points have the same sign in the pT direction, interpolate logs
				f1 = sign_of_f11 * exp( lin_int(ptr-pT0, one_by_pTdiff, log_f11_arr[qpt_cs_idx], log_f21_arr[qpt_cs_idx]) );
			else					// otherwise, just interpolate original vals
				f1 = lin_int(ptr-pT0, one_by_pTdiff, f11_arr[qpt_cs_idx], f21_arr[qpt_cs_idx]);
				
			//*******************************************************************************************************************
			// set f2 next
			//*******************************************************************************************************************
			if ( ptr > pT1 && ( log_f22_arr[qpt_cs_idx] > log_f12_arr[qpt_cs_idx] || sign_of_f12 * sign_of_f22 < 0 ) )
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
		// note "+=", since must sum over "+" and "-" roots in Eq. (A19) in Wiedemann & Heinz (1997)
		//results[qpt_cs_idx] += lin_int(phir-phi0, one_by_pphidiff, f1, f2);
		double Zkr = lin_int(phir-phi0, one_by_pphidiff, f1, f2);

		double f1s = f1;
		double f2s = f2;

		/////////////////////////////////////////////////////////////////
		// DO SINE PART NEXT
		/////////////////////////////////////////////////////////////////
		// interpolate over pT values first
		/////////////////////////////////////////////////////////////////
		if(ptr > PTCHANGE)				// if pT interpolation point is larger than PTCHANGE (currently 1.0 GeV)
		{
			double sign_of_f11 = sign_of_f11_arr[qpt_cs_idx+1];
			double sign_of_f12 = sign_of_f12_arr[qpt_cs_idx+1];
			double sign_of_f21 = sign_of_f21_arr[qpt_cs_idx+1];
			double sign_of_f22 = sign_of_f22_arr[qpt_cs_idx+1];
			
			//*******************************************************************************************************************
			// set f1 first
			//*******************************************************************************************************************
			// if using extrapolation and spectra at pT1 has larger magnitude than at pT0 (or the signs are different), just return zero
			if ( ptr > pT1 && ( log_f21_arr[qpt_cs_idx+1] > log_f11_arr[qpt_cs_idx+1] || sign_of_f11 * sign_of_f21 < 0 ) )
				f1 = 0.0;
			else if (sign_of_f11 * sign_of_f21 > 0)	// if the two points have the same sign in the pT direction, interpolate logs
				f1 = sign_of_f11 * exp( lin_int(ptr-pT0, one_by_pTdiff, log_f11_arr[qpt_cs_idx+1], log_f21_arr[qpt_cs_idx+1]) );
			else					// otherwise, just interpolate original vals
				f1 = lin_int(ptr-pT0, one_by_pTdiff, f11_arr[qpt_cs_idx+1], f21_arr[qpt_cs_idx+1]);
				
			//*******************************************************************************************************************
			// set f2 next
			//*******************************************************************************************************************
			if ( ptr > pT1 && ( log_f22_arr[qpt_cs_idx+1] > log_f12_arr[qpt_cs_idx+1] || sign_of_f12 * sign_of_f22 < 0 ) )
				f2 = 0.0;
			else if (sign_of_f12 * sign_of_f22 > 0)	// if the two points have the same sign in the pT direction, interpolate logs
				f2 = sign_of_f12 * exp( lin_int(ptr-pT0, one_by_pTdiff, log_f12_arr[qpt_cs_idx+1], log_f22_arr[qpt_cs_idx+1]) );
			else					// otherwise, just interpolate original vals
				f2 = lin_int(ptr-pT0, one_by_pTdiff, f12_arr[qpt_cs_idx+1], f22_arr[qpt_cs_idx+1]);
			//*******************************************************************************************************************
		}
		else						// if pT is smaller than PTCHANGE, just use linear interpolation, no matter what
		{
			f1 = lin_int(ptr-pT0, one_by_pTdiff, f11_arr[qpt_cs_idx+1], f21_arr[qpt_cs_idx+1]);
			f2 = lin_int(ptr-pT0, one_by_pTdiff, f12_arr[qpt_cs_idx+1], f22_arr[qpt_cs_idx+1]);
		}
					
		// now, interpolate f1 and f2 over the pphi direction
		// note "+=", since must sum over "+" and "-" roots in Eq. (A19) in Wiedemann & Heinz (1997)
		//results[qpt_cs_idx+1] += lin_int(phir-phi0, one_by_pphidiff, f1, f2);
		double Zki = lin_int(phir-phi0, one_by_pphidiff, f1, f2);
					
		//Finally, update results vectors appropriately
		//--> update the real part of weighted daughter spectra
		results[qpt_cs_idx] += akr*Zkr-aki*Zki;
		//--> update the imaginary part of weighted daughter spectra
		results[qpt_cs_idx+1] += akr*Zki+aki*Zkr;


		/*if ( loc_verb || isinf( results[qpt_cs_idx] ) || isnan( results[qpt_cs_idx] ) || isinf( results[qpt_cs_idx+1] ) || isnan( results[qpt_cs_idx+1] ) )
		{
			*global_out_stream_ptr << "ERROR in eiqxEdndp3(double, double, double*): problems encountered!" << endl
				<< "results(" << iqt << "," << iqx << "," << iqy << "," << iqz << ") = "
				<< setw(25) << setprecision(20) << results[qpt_cs_idx] << ",   " << results[qpt_cs_idx+1] << endl
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
				<< "  --> f1 = " << f1s << endl
				<< "  --> f2 = " << f2s << endl
				<< "  --> akr = " << akr << endl
				<< "  --> aki = " << aki << endl
				<< "  --> Zkr = " << Zkr << endl
				<< "  --> Zki = " << Zki << endl
				<< "  --> akr*Zkr-aki*Zki = " << akr*Zkr-aki*Zki << endl
				<< "  --> akr*Zki+aki*Zkr = " << akr*Zki+aki*Zkr << endl;
							//exit(1);
		}*/

		//++qpt_cs_idx;	// step to next cell in results array
		qpt_cs_idx += 2;
		qlist_idx++;
	}	// ending all loops at once in linearized version

	return;
}

//End of file
