#include<iostream>
#include<sstream>
#include<string>
#include<fstream>
#include<cmath>
#include<iomanip>
#include<vector>
#include<set>
#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<queue>
#include<map>
#include<numeric>

#include "CFWR.h"
#include "Arsenal.h"
#include "Stopwatch.h"
#include "CPStopwatch.h"
#include "gauss_quadrature.h"

using namespace std;

// only need to calculated interpolation grid of spacetime moments for each resonance, NOT each decay channel!
bool recycle_previous_moments = false;
bool recycle_similar_moments = false;
int reso_particle_id_of_moments_to_recycle = -1;
string reso_name_of_moments_to_recycle = "NULL";
string current_decay_channel_string = "NULL";

template < typename T >
void check_for_NaNs(string variable_name, const T variable_value, ofstream& localout)
{
	if (isnan(variable_value))
		localout << "ERROR: " << variable_name << " = " << variable_value << endl;
	return;
}

double CorrelationFunction::place_in_range(double phi, double min, double max)
{
	while (phi < min || phi > max)
	{
		if (phi < min) phi += twopi;
		else phi -= twopi;
	}

	return (phi);
}

// ************************************************************
// Compute correlation function at all specified q points for all resonances here
// ************************************************************
void CorrelationFunction::Compute_correlation_function(FO_surf* FOsurf_ptr)
{
	Stopwatch BIGsw;
	int decay_channel_loop_cutoff = n_decay_channels;			//loop over direct pions and decay_channels

	int HDFInitializationSuccess = Initialize_resonance_HDF_array();
	if (HDFInitializationSuccess < 0)
	{
		cerr << "Failed to initialize HDF array of resonances!  Exiting..." << endl;
		exit;
	}

	*global_out_stream_ptr << "Setting spacetime moments grid..." << endl;
	BIGsw.tic();
	// ************************************************************
	// loop over decay_channels (idc == 0 corresponds to thermal pions)
	// ************************************************************
	for (int idc = 0; idc <= decay_channel_loop_cutoff; ++idc)				//this is inefficient, but will do the job for now
	{
		// ************************************************************
		// check whether to do this decay channel
		// ************************************************************
		if (idc > 0 && thermal_pions_only)
			break;
		else if (!Do_this_decay_channel(idc))
			continue;

		// ************************************************************
		// if so, set decay channel info
		// ************************************************************
		Set_current_particle_info(idc);

		// ************************************************************
		// decide whether to recycle old moments or calculate new moments
		// ************************************************************
		Get_spacetime_moments(FOsurf_ptr, idc);
	}	//computing all resonances' spacetime moments here first
		//THEN do phase-space integrals

	// once all spacetime moments have been computed, get rid of weighted S_p array to save space
	Delete_S_p_withweight_array();

	if (VERBOSE > 0) *global_out_stream_ptr << endl << "************************************************************"
											<< endl << "* Computed all (thermal) space-time moments!" << endl
											<< "************************************************************" << endl << endl;
	BIGsw.toc();
	*global_out_stream_ptr << "\t ...finished all (thermal) space-time moments in " << BIGsw.takeTime() << " seconds." << endl;
	
	if (SPACETIME_MOMENTS_ONLY)
		return;
	else if (thermal_pions_only)
		goto correlation_function_calculation;

	*global_out_stream_ptr << "Computing all phase-space integrals..." << endl;
	BIGsw.tic();
	// ************************************************************
	// Compute feeddown with heaviest resonances first
	// ************************************************************
	for (int idc = 1; idc <= decay_channel_loop_cutoff; ++idc)
	{
		// ************************************************************
		// check whether to do this decay channel
		// ************************************************************
		if (decay_channels[idc-1].resonance_particle_id == target_particle_id || thermal_pions_only)
			break;
		else if (!Do_this_decay_channel(idc))
			continue;

		// ************************************************************
		// if so, set decay channel info
		// ************************************************************
		Set_current_particle_info(idc);
		Load_resonance_and_daughter_spectra(decay_channels[idc-1].resonance_particle_id);

		// ************************************************************
		// begin source variances calculations here...
		// ************************************************************
		Allocate_decay_channel_info();				// allocate needed memory

		for (int idc_DI = 0; idc_DI < current_reso_nbody; ++idc_DI)
		{
			int daughter_resonance_particle_id = -1;
			if (!Do_this_daughter_particle(idc, idc_DI, &daughter_resonance_particle_id))
				continue;
			Set_current_daughter_info(idc, idc_DI);
			Do_resonance_integrals(current_resonance_particle_id, daughter_resonance_particle_id, idc);
		}

		Update_daughter_spectra();
		Delete_decay_channel_info();				// free up memory
	}											// END of decay channel loop

correlation_function_calculation:
	// Now, with all resonance contributions to correlation function computed, do the actual calculation
	Cal_correlationfunction();

   return;
}

bool CorrelationFunction::Do_this_decay_channel(int dc_idx)
{
	string local_name = "Thermal pion(+)";
	if (dc_idx == 0)
	{
		if (VERBOSE > 0) *global_out_stream_ptr << endl << local_name << ": doing this one." << endl;
		return true;
	}
	else
	{
		local_name = decay_channels[dc_idx-1].resonance_name;
		Get_current_decay_string(dc_idx, &current_decay_channel_string);
	}
	if (VERBOSE > 0) *global_out_stream_ptr << endl << local_name << ": skipping decay " << current_decay_channel_string << "." << endl;

	return (decay_channels[dc_idx-1].include_channel);
}

// ************************************************************
// Checks whether to do daughter particle for any given decay channel
// ************************************************************
bool CorrelationFunction::Do_this_daughter_particle(int dc_idx, int daughter_idx, int * daughter_resonance_pid)
{
	// assume dc_idx > 0
	string local_name = decay_channels[dc_idx-1].resonance_name;

	// look up daughter particle info
	int temp_monval = decay_channels[dc_idx-1].resonance_decay_monvals[daughter_idx];

	if (temp_monval == 0)
		return false;

	int temp_ID = lookup_particle_id_from_monval(all_particles, Nparticle, temp_monval);
	//*daughter_resonance_idx = lookup_resonance_idx_from_particle_id(temp_ID) + 1;
	*daughter_resonance_pid = temp_ID;
	// if daughter was found in chosen_resonances or is pion(+), this is correct
	particle_info temp_daughter = all_particles[temp_ID];

	if (*daughter_resonance_pid < 0 && temp_daughter.monval != particle_monval && temp_daughter.effective_branchratio >= 1.e-12)
		*global_out_stream_ptr << "Couldn't find " << temp_daughter.name << " in chosen_resonances!  Results are probably not reliable..." << endl;

	//bool daughter_does_not_contribute = ( (temp_daughter.stable == 1 || temp_daughter.effective_branchratio < 1.e-12) && temp_daughter.monval != particle_monval );
	bool daughter_does_not_contribute = ( (temp_daughter.decays_Npart[0] == 1 || temp_daughter.effective_branchratio < 1.e-12) && temp_daughter.monval != particle_monval );

	// if daughter particle gives no contribution to final pion spectra
	if (daughter_does_not_contribute)
	{
		if (VERBOSE > 0) *global_out_stream_ptr << "\t * " << local_name << ": in decay " << current_decay_channel_string << ", skipping " << temp_daughter.name
												<< " (daughter_resonance_pid = " << *daughter_resonance_pid << ")." << endl;
		return false;
	}
	else
	{
		if (VERBOSE > 0) *global_out_stream_ptr << "\t * " << local_name << ": in decay " << current_decay_channel_string << ", doing " << temp_daughter.name
												<< " (daughter_resonance_pid = " << *daughter_resonance_pid << ")." << endl;
		return true;
	}
}

void CorrelationFunction::Set_current_particle_info(int dc_idx)
{
	if (dc_idx == 0)
	{
		muRES = particle_mu;
		signRES = particle_sign;
		gRES = particle_gspin;
		current_resonance_particle_id = target_particle_id;
		
		return;
	}
	else
	{
		// assume dc_idx > 0
		string local_name = decay_channels[dc_idx-1].resonance_name;

		if (VERBOSE > 0) *global_out_stream_ptr << local_name << ": doing decay " << current_decay_channel_string << "." << endl
			<< "\t * " << local_name << ": setting information for this decay channel..." << endl;

		if (dc_idx > 1)
		{
			//cerr << "Setting previous decay channel information for dc_idx = " << dc_idx << endl;
			previous_resonance_particle_id = current_resonance_particle_id;		//for look-up in all_particles
			previous_decay_channel_idx = current_decay_channel_idx;			//different for each decay channel
			previous_resonance_idx = current_resonance_idx;				//different for each decay channel
			previous_resonance_mass = current_resonance_mass;
			previous_resonance_Gamma = current_resonance_Gamma;
			previous_resonance_total_br = current_resonance_total_br;
			previous_resonance_direct_br = current_resonance_direct_br;
			previous_reso_nbody = current_reso_nbody;
		}
		//cerr << "Setting current decay channel information for dc_idx = " << dc_idx << endl;
		current_decay_channel_idx = dc_idx;
		current_resonance_idx = decay_channels[dc_idx-1].resonance_idx;
		current_resonance_particle_id = decay_channels[dc_idx-1].resonance_particle_id;
		current_resonance_mass = decay_channels[dc_idx-1].resonance_mass;
		current_resonance_Gamma = decay_channels[dc_idx-1].resonance_Gamma;
		current_resonance_total_br = decay_channels[dc_idx-1].resonance_total_br;
		current_resonance_direct_br = decay_channels[dc_idx-1].resonance_direct_br;
		current_reso_nbody = decay_channels[dc_idx-1].nbody;
		
		// might want to rename these for notational consistency...
		muRES = decay_channels[dc_idx-1].resonance_mu;
		signRES = decay_channels[dc_idx-1].resonance_sign;
		gRES = decay_channels[dc_idx-1].resonance_gspin;
		
		if (dc_idx > 1)
		{
			int similar_particle_idx = -1;
			int temp_reso_idx = decay_channels[dc_idx-1].resonance_idx;
			
			if ( current_resonance_particle_id == previous_resonance_particle_id )
			{
				//previous resonance is the same as this one...
				recycle_previous_moments = true;
				recycle_similar_moments = false;
			}
			else if ( Search_for_similar_particle( temp_reso_idx, &similar_particle_idx ) )
			{
				//previous resonance is NOT the same as this one BUT this one is sufficiently similar to some preceding one...
				recycle_previous_moments = false;
				recycle_similar_moments = true;
				reso_particle_id_of_moments_to_recycle = chosen_resonances[similar_particle_idx];
			}
			else
			{
				recycle_previous_moments = false;
				recycle_similar_moments = false;
				reso_particle_id_of_moments_to_recycle = -1;	//guarantees it won't be used spuriously
			}
		}
	}
	
	return;
}

void CorrelationFunction::Set_current_daughter_info(int dc_idx, int daughter_idx)
{
	if (dc_idx > 1)
	{
		previous_resonance_particle_id = current_resonance_particle_id;		//for look-up in all_particles
		previous_decay_channel_idx = current_decay_channel_idx;			//different for each decay channel
		previous_resonance_idx = current_resonance_idx;
		previous_resonance_mass = current_resonance_mass;
		previous_resonance_Gamma = current_resonance_Gamma;
		previous_m2_Gamma = current_m2_Gamma;
		previous_m3_Gamma = current_m3_Gamma;
		previous_resonance_total_br = current_resonance_total_br;
		previous_resonance_direct_br = current_resonance_direct_br;
		previous_reso_nbody = current_reso_nbody;
		previous_daughter_mass = current_daughter_mass;
		previous_daughter_Gamma = current_daughter_Gamma;
	}
	current_decay_channel_idx = dc_idx;
	current_resonance_idx = decay_channels[dc_idx-1].resonance_idx;
	current_resonance_particle_id = decay_channels[dc_idx-1].resonance_particle_id;
	current_resonance_mass = decay_channels[dc_idx-1].resonance_mass;
	current_resonance_Gamma = decay_channels[dc_idx-1].resonance_Gamma;
	current_resonance_total_br = decay_channels[dc_idx-1].resonance_total_br;
	current_resonance_direct_br = decay_channels[dc_idx-1].resonance_direct_br;
	current_reso_nbody = decay_channels[dc_idx-1].nbody;
	current_daughter_mass = decay_channels[dc_idx-1].resonance_decay_masses[daughter_idx];
	current_daughter_Gamma = decay_channels[dc_idx-1].resonance_decay_Gammas[daughter_idx];

	// might want to rename these for notational consistency...
	muRES = decay_channels[dc_idx-1].resonance_mu;
	signRES = decay_channels[dc_idx-1].resonance_sign;
	gRES = decay_channels[dc_idx-1].resonance_gspin;

	// set non-daughter decay masses for computing contributions to spectra of daughter
	double m2ex = 0.0, m3ex = 0.0, m4ex = 0.0;
	switch(current_reso_nbody)
	{
		case 1:
			break;
		case 2:
			current_resonance_decay_masses[1] = 0.0;
			current_m3_Gamma = decay_channels[dc_idx-1].resonance_decay_Gammas[0];
			if (daughter_idx == 0)
			{
				current_resonance_decay_masses[0] = decay_channels[dc_idx-1].resonance_decay_masses[1];
				current_m2_Gamma = decay_channels[dc_idx-1].resonance_decay_Gammas[1];
			}
			else
			{
				current_resonance_decay_masses[0] = decay_channels[dc_idx-1].resonance_decay_masses[0];
				current_m2_Gamma = decay_channels[dc_idx-1].resonance_decay_Gammas[0];
			}
			break;
		case 3:
			if (daughter_idx == 0)
			{
				current_resonance_decay_masses[0] = decay_channels[dc_idx-1].resonance_decay_masses[1];
				current_resonance_decay_masses[1] = decay_channels[dc_idx-1].resonance_decay_masses[2];
				current_m2_Gamma = decay_channels[dc_idx-1].resonance_decay_Gammas[1];
				current_m3_Gamma = decay_channels[dc_idx-1].resonance_decay_Gammas[2];
			}
			else if (daughter_idx == 1)
			{
				current_resonance_decay_masses[0] = decay_channels[dc_idx-1].resonance_decay_masses[0];
				current_resonance_decay_masses[1] = decay_channels[dc_idx-1].resonance_decay_masses[2];
				current_m2_Gamma = decay_channels[dc_idx-1].resonance_decay_Gammas[0];
				current_m3_Gamma = decay_channels[dc_idx-1].resonance_decay_Gammas[2];
			}
			else
			{
				current_resonance_decay_masses[0] = decay_channels[dc_idx-1].resonance_decay_masses[0];
				current_resonance_decay_masses[1] = decay_channels[dc_idx-1].resonance_decay_masses[1];
				current_m2_Gamma = decay_channels[dc_idx-1].resonance_decay_Gammas[0];
				current_m3_Gamma = decay_channels[dc_idx-1].resonance_decay_Gammas[1];
			}
			break;
		case 4:
			if (daughter_idx == 0)
			{
				m2ex = decay_channels[dc_idx-1].resonance_decay_masses[1];
				m3ex = decay_channels[dc_idx-1].resonance_decay_masses[2];
				m4ex = decay_channels[dc_idx-1].resonance_decay_masses[3];
			}
			else if (daughter_idx == 1)
			{
				m2ex = decay_channels[dc_idx-1].resonance_decay_masses[0];
				m3ex = decay_channels[dc_idx-1].resonance_decay_masses[2];
				m4ex = decay_channels[dc_idx-1].resonance_decay_masses[3];
			}
			else if (daughter_idx == 2)
			{
				m2ex = decay_channels[dc_idx-1].resonance_decay_masses[0];
				m3ex = decay_channels[dc_idx-1].resonance_decay_masses[1];
				m4ex = decay_channels[dc_idx-1].resonance_decay_masses[3];
			}
			else
			{
				m2ex = decay_channels[dc_idx-1].resonance_decay_masses[0];
				m3ex = decay_channels[dc_idx-1].resonance_decay_masses[1];
				m4ex = decay_channels[dc_idx-1].resonance_decay_masses[2];
			}
			current_resonance_decay_masses[0] = m2ex;
			current_resonance_decay_masses[1] = 0.5 * (m3ex + m4ex + current_resonance_mass - current_daughter_mass - m2ex);
			// approximation obtained from earlier resonances code
			//*global_out_stream_ptr << "Current decay " << current_decay_channel_string << ", br = " << current_resonance_direct_br
			//						<< ": {m2ex, m3ex, m4ex, m3eff} = {"
			//						<< m2ex << ", " << m3ex << ", " << m4ex << ", " << current_resonance_decay_masses[1] << "}" << endl;
			break;
		default:
			cerr << "Set_current_daughter_info(): shouldn't have ended up here, bad value of current_reso_nbody!" << endl;
			exit(1);
	}
}

bool CorrelationFunction::Search_for_similar_particle(int reso_idx, int * result)
{
	// for the timebeing, just search from beginning of decay_channels until similar particle is found;
	// should be more careful, since could lead to small numerical discrepancies if similar particle was
	// already recycled by some other (dissimilar) particle, but ignore this complication for now...
	*result = -1;
	
	for (int local_ir = 0; local_ir < reso_idx; ++local_ir)
	{// only need to search decay_channels that have already been calculated
		if (particles_are_the_same(local_ir, reso_idx))
		{
			*result = local_ir;
			break;
		}
	}
	
	return (*result >= 0);
}

//**********************************************************************************************
bool CorrelationFunction::particles_are_the_same(int reso_idx1, int reso_idx2)
{
	int icr1 = chosen_resonances[reso_idx1];
	int icr2 = chosen_resonances[reso_idx2];
	//int icr1 = reso_idx1;
	//int icr2 = reso_idx2;
	if (all_particles[icr1].sign != all_particles[icr2].sign)
		return false;
	if (abs(all_particles[icr1].mass-all_particles[icr2].mass) / (all_particles[icr2].mass+1.e-30) > PARTICLE_DIFF_TOLERANCE)
		return false;
	//assume chemical potential mu is constant over entire FO surface
	double chem1 = all_particles[icr1].mu, chem2 = all_particles[icr2].mu;
	if (2.*abs(chem1 - chem2)/(chem1 + chem2 + 1.e-30) > PARTICLE_DIFF_TOLERANCE)
		return false;

	return true;
}

void CorrelationFunction::Recycle_spacetime_moments()
{
	int HDFcopyChunkSuccess = Copy_chunk(current_resonance_particle_id, reso_particle_id_of_moments_to_recycle);

	return;
}

//**************************************************************
//**************************************************************
void CorrelationFunction::Load_resonance_and_daughter_spectra(int local_pid)
{
	// get parent resonance spectra, set logs and signs arrays that are needed for interpolation
	int getHDFresonanceSpectra = Get_resonance_from_HDF_array(local_pid, current_dN_dypTdpTdphi_moments);
	Set_current_resonance_logs_and_signs();

	// get spectra for all daughters, set all of the logs and signs arrays that are needed for interpolation
	int n_daughters = Set_daughter_list(local_pid);

	if (n_daughters > 0)
	{
		int d_idx = 0;

		Setup_current_daughters_dN_dypTdpTdphi_moments(n_daughters);

		for (set<int>::iterator it = daughter_resonance_indices.begin(); it != daughter_resonance_indices.end(); ++it)
		{
			int daughter_pid = *it;		//daughter pid is pointed to by iterator
			getHDFresonanceSpectra = Get_resonance_from_HDF_array(daughter_pid, current_daughters_dN_dypTdpTdphi_moments[d_idx]);
			++d_idx;
		}
	
		Set_current_daughters_resonance_logs_and_signs(n_daughters);
	}
	else
	{
		cerr << "Particle is stable, shouldn't have ended up here!  Something went wrong..." << endl;
		exit;
	}

	return;
}

void CorrelationFunction::Update_daughter_spectra()
{
	int d_idx = 0;
	for (set<int>::iterator it = daughter_resonance_indices.begin(); it != daughter_resonance_indices.end(); ++it)
	{
		int daughter_pid = *it;		//daughter pid is pointed to by iterator
		int setHDFresonanceSpectra = Set_resonance_in_HDF_array(daughter_pid, current_daughters_dN_dypTdpTdphi_moments[d_idx]);

		for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
		for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
			spectra[daughter_pid][ipt][ipphi] = current_daughters_dN_dypTdpTdphi_moments[d_idx][ipt][ipphi][iqt0][iqx0][iqy0][iqz0][0];

		++d_idx;
	}

	// cleanup previous iteration and setup the new one
	Cleanup_current_daughters_dN_dypTdpTdphi_moments(daughter_resonance_indices.size());

	return;
}

void CorrelationFunction::Set_current_resonance_logs_and_signs()
{
	for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
	for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
	for (int iqt = 0; iqt < qtnpts; ++iqt)
	for (int iqx = 0; iqx < qxnpts; ++iqx)
	for (int iqy = 0; iqy < qynpts; ++iqy)
	for (int iqz = 0; iqz < qznpts; ++iqz)
	for (int itrig = 0; itrig < 2; ++itrig)
	{
		double temp = current_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx][iqy][iqz][itrig];
		current_ln_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx][iqy][iqz][itrig] = log(abs(temp)+1.e-100);
		current_sign_of_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx][iqy][iqz][itrig] = sgn(temp);
	}

	return;
}

void CorrelationFunction::Set_current_daughters_resonance_logs_and_signs(int n_daughters)
{
	for (int idaughter = 0; idaughter < n_daughters; ++idaughter)
	for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
	for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
	for (int iqt = 0; iqt < qtnpts; ++iqt)
	for (int iqx = 0; iqx < qxnpts; ++iqx)
	for (int iqy = 0; iqy < qynpts; ++iqy)
	for (int iqz = 0; iqz < qznpts; ++iqz)
	for (int itrig = 0; itrig < 2; ++itrig)
	{
		double temp = current_daughters_dN_dypTdpTdphi_moments[idaughter][ipt][ipphi][iqt][iqx][iqy][iqz][itrig];
		current_daughters_ln_dN_dypTdpTdphi_moments[idaughter][ipt][ipphi][iqt][iqx][iqy][iqz][itrig] = log(abs(temp)+1.e-100);
		current_daughters_sign_of_dN_dypTdpTdphi_moments[idaughter][ipt][ipphi][iqt][iqx][iqy][iqz][itrig] = sgn(temp);
	}

	return;
}
//**************************************************************
//**************************************************************

void CorrelationFunction::Get_spacetime_moments(FO_surf* FOsurf_ptr, int dc_idx)
{
//**************************************************************
//Set resonance name
//**************************************************************
	string local_name = "Thermal pion(+)";
	if (dc_idx > 0)
		local_name = decay_channels[dc_idx-1].resonance_name;
//**************************************************************
//Decide what to do with this resonance / decay channel
//**************************************************************
	if (recycle_previous_moments && dc_idx > 1)	// similar (but different) earlier resonance
	{
		if (VERBOSE > 0) *global_out_stream_ptr << local_name
			<< ": new parent resonance (" << decay_channels[current_decay_channel_idx-1].resonance_name << ", dc_idx = " << current_decay_channel_idx
			<< " of " << n_decay_channels << ") same as preceding parent resonance \n\t\t--> reusing old dN_dypTdpTdphi_moments!" << endl;
	}
	else if (recycle_similar_moments && dc_idx > 1)	// similar (but different) earlier resonance
	{
		if (VERBOSE > 0) *global_out_stream_ptr << local_name
			<< ": new parent resonance (" << decay_channels[current_decay_channel_idx-1].resonance_name << ", dc_idx = " << current_decay_channel_idx
			<< " of " << n_decay_channels << ") sufficiently close to preceding parent resonance (" << all_particles[reso_particle_id_of_moments_to_recycle].name
			<< ", reso_particle_id = " << reso_particle_id_of_moments_to_recycle << ") \n\t\t--> reusing old dN_dypTdpTdphi_moments!" << endl;
		Recycle_spacetime_moments();
	}
	else
	{
		if (dc_idx == 0)	//if it's thermal pions
		{
			if (VERBOSE > 0) *global_out_stream_ptr << "  --> Computing dN_dypTdpTdphi_moments for thermal pion(+)!" << endl;
		}
		else if (dc_idx == 1)	//if it's the first resonance
		{
			if (VERBOSE > 0) *global_out_stream_ptr << "  --> Computing dN_dypTdpTdphi_moments for " << local_name << endl;
		}
		else			//if it's a later resonance
		{
			if (VERBOSE > 0 && !recycle_previous_moments && !recycle_similar_moments) *global_out_stream_ptr << local_name
				<< ": new parent resonance (" << decay_channels[current_decay_channel_idx-1].resonance_name << ", dc_idx = " << current_decay_channel_idx
				<< " of " << n_decay_channels << ") dissimilar from all preceding decay_channels \n\t\t--> calculating new dN_dypTdpTdphi_moments!" << endl;
			else
			{
				cerr << "You shouldn't have ended up here!" << endl;
				exit(1);
			}
		}
		Set_dN_dypTdpTdphi_moments(FOsurf_ptr, current_resonance_particle_id);
	}
//**************************************************************
//Spacetime moments now set
//**************************************************************
	return;
}

void CorrelationFunction::Set_dN_dypTdpTdphi_moments(FO_surf* FOsurf_ptr, int local_pid)
{
	double localmass = all_particles[local_pid].mass;
	string local_name = "Thermal pion(+)";
	if (local_pid != target_particle_id)
		local_name = all_particles[local_pid].name;

	for(int i=0; i<eta_s_npts; ++i)
	{
		double local_eta_s = eta_s[i];
		double local_cosh = cosh(SP_p_y - local_eta_s);
		double local_sinh = sinh(SP_p_y - local_eta_s);

		for(int ipt=0; ipt<n_interp_pT_pts; ++ipt)
		{
			double mT = sqrt(localmass*localmass + SPinterp_pT[ipt]*SPinterp_pT[ipt]);
			SPinterp_p0[ipt][i] = mT*local_cosh;
			SPinterp_pz[ipt][i] = mT*local_sinh;
			SPinterp_p0_Transpose[i][ipt] = mT*local_cosh;
			SPinterp_pz_Transpose[i][ipt] = mT*local_sinh;
		}
	}
	
	// get spectra at each fluid cell, sort by importance
	*global_out_stream_ptr << "Computing spectra..." << endl;
	CPStopwatch sw;
	sw.Start();
	Cal_dN_dypTdpTdphi_heap(FOsurf_ptr, local_pid, 0.9);
	sw.Stop();
	*global_out_stream_ptr << "CP#1: Took " << sw.printTime() << " seconds." << endl;

	// get weighted spectra with only most important fluid cells, up to given threshhold
	*global_out_stream_ptr << "Computing weighted spectra..." << endl;
	sw.Reset();
	Cal_dN_dypTdpTdphi_with_weights(FOsurf_ptr, local_pid);
	sw.Stop();
	*global_out_stream_ptr << "CP#2: Took " << sw.printTime() << " seconds." << endl;

	// store in HDF5 file
	int setHDFresonanceSpectra = Set_resonance_in_HDF_array(local_pid, current_dN_dypTdpTdphi_moments);
	if (setHDFresonanceSpectra < 0)
	{
		cerr << "Failed to initialize HDF array of resonances!  Exiting..." << endl;
		exit;
	}


	return;
}

void CorrelationFunction::form_trig_sign_z(int isurf, int ieta, int iqt, int iqx, int iqy, int iqz, int ii, double * results)
{
	double zfactor = 1.0 - 2.0 * double(ii);		// related to symmetry in z-direction
	double cosA0 = osc0[isurf][ieta][iqt][0], cosA1 = osc1[isurf][iqx][0], cosA2 = osc2[isurf][iqy][0], cosA3 = osc3[isurf][ieta][iqz][0];
	double sinA0 = osc0[isurf][ieta][iqt][1], sinA1 = osc1[isurf][iqx][1], sinA2 = osc2[isurf][iqy][1], sinA3 = osc3[isurf][ieta][iqz][1];
	results[0] = cosA0*cosA1*cosA2*cosA3
				+ cosA2*cosA3*sinA0*sinA1
				+ cosA1*cosA3*sinA0*sinA2
				- cosA0*cosA3*sinA1*sinA2
				+ zfactor*cosA1*cosA2*sinA0*sinA3
				- zfactor*cosA0*cosA2*sinA1*sinA3
				- zfactor*cosA0*cosA1*sinA2*sinA3
				- zfactor*sinA0*sinA1*sinA2*sinA3;
	results[1] = cosA1*cosA2*cosA3*sinA0
				- cosA0*cosA2*cosA3*sinA1
				- cosA0*cosA1*cosA3*sinA2
				- cosA3*sinA0*sinA1*sinA2
				- zfactor*cosA0*cosA1*cosA2*sinA3
				- zfactor*cosA2*sinA0*sinA1*sinA3
				- zfactor*cosA1*sinA0*sinA2*sinA3
				+ zfactor*cosA0*sinA1*sinA2*sinA3;

	return;
}

inline void CorrelationFunction::addElementToQueue(priority_queue<pair<double, size_t> >& p, pair<double, size_t> elem, size_t max_size)
{
	if ( ( max_size <= p.size() ) && ( elem >= p.top() ) )
		return; // nothing to do.
	p.push(elem);
	//while ( max_size < p.size() )
	if( max_size < p.size() )
		p.pop();
	return;
};

void CorrelationFunction::Cal_dN_dypTdpTdphi_heap(FO_surf* FOsurf_ptr, int local_pid, double cutoff)
{
	double ** temp_moments_array = new double * [n_interp_pT_pts];
	double ** abs_temp_moments_array = new double * [n_interp_pT_pts];
	for (int ipt = 0; ipt < n_interp_pT_pts; ipt++)
	{
		temp_moments_array[ipt] = new double [n_interp_pphi_pts];
		abs_temp_moments_array[ipt] = new double [n_interp_pphi_pts];
		for (int ipphi = 0; ipphi < n_interp_pphi_pts; ipphi++)
		{
			temp_moments_array[ipt][ipphi] = 0.0;
			abs_temp_moments_array[ipt][ipphi] = 0.0;
		}
	}

	// set particle information
	double sign = all_particles[local_pid].sign;
	double degen = all_particles[local_pid].gspin;
	double localmass = all_particles[local_pid].mass;
	double mu = all_particles[local_pid].mu;

	// set some freeze-out surface information that's constant the whole time
	double prefactor = 1.0*degen/(8.0*M_PI*M_PI*M_PI)/(hbarC*hbarC*hbarC);
	double Tdec = (&FOsurf_ptr[0])->Tdec;
	double Pdec = (&FOsurf_ptr[0])->Pdec;
	double Edec = (&FOsurf_ptr[0])->Edec;
	double one_by_Tdec = 1./Tdec;
	double deltaf_prefactor = 0.;
	if (use_delta_f)
		deltaf_prefactor = 1./(2.0*Tdec*Tdec*(Edec+Pdec));

	double eta_s_symmetry_factor = 2.0;

CPStopwatch sw, sw2, sw3;

sw2.Start();
cutoff_FOcells.resize( n_interp_pT_pts * n_interp_pphi_pts );
cutoff_FOcell_vals_C.resize( n_interp_pT_pts * n_interp_pphi_pts );
cutoff_FOcell_vals_S.resize( n_interp_pT_pts * n_interp_pphi_pts );
pc_cutoff_vals.resize( number_of_percentage_markers );
	for(int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
	{
		//if (ipt > 0) continue;
		double pT = SPinterp_pT[ipt];
		double * p0_pTslice = SPinterp_p0[ipt];
		double * pz_pTslice = SPinterp_pz[ipt];

		for(int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
		{
			//if (ipphi > 0) continue;
			int ptphi_index = ipt * n_interp_pphi_pts + ipphi;
			double sin_pphi = sin_SPinterp_pphi[ipphi];
			double cos_pphi = cos_SPinterp_pphi[ipphi];
			//*global_out_stream_ptr << "\t\t--> Doing pT = " << pT << ", pphi = " << SPinterp_pphi[ipphi] << "..." << endl;
			double px = pT*cos_pphi;
			double py = pT*sin_pphi;
			number_of_FOcells_above_cutoff_array[ipt][ipphi] = FO_length * eta_s_npts;
			int current_num_of_FOcells_ACA = number_of_FOcells_above_cutoff_array[ipt][ipphi];

			double tempsum = 0.0, tempabssum = 0.0;
			double * tmp_S_p_withweight_array = new double [FO_length * eta_s_npts];
			int iFOcell = 0;
			priority_queue<pair<double, size_t> > FOcells_PQ;

			for(int isurf = 0; isurf < FO_length; ++isurf)
			{
				FO_surf * surf = &FOsurf_ptr[isurf];

				double tau = surf->tau;

				double vx = surf->vx;
				double vy = surf->vy;
				double gammaT = surf->gammaT;

				double da0 = surf->da0;
				double da1 = surf->da1;
				double da2 = surf->da2;

				double pi00 = surf->pi00;
				double pi01 = surf->pi01;
				double pi02 = surf->pi02;
				double pi11 = surf->pi11;
				double pi12 = surf->pi12;
				double pi22 = surf->pi22;
				double pi33 = surf->pi33;

				double temp_running_sum = 0.0;

				for(int ieta = 0; ieta < eta_s_npts; ++ieta)
				{
					double p0 = p0_pTslice[ieta];
					double pz = pz_pTslice[ieta];

					//now get distribution function, emission function, etc.
					double f0 = 1./(exp( (gammaT*(p0*1. - px*vx - py*vy) - mu)*one_by_Tdec ) + sign);	//thermal equilibrium distributions
	
					//viscous corrections
					double deltaf = 0.;
					if (use_delta_f)
						deltaf = deltaf_prefactor * (1. - sign*f0)
								* (p0*p0*pi00 - 2.0*p0*px*pi01 - 2.0*p0*py*pi02 + px*px*pi11 + 2.0*px*py*pi12 + py*py*pi22 + pz*pz*pi33);

					//p^mu d^3sigma_mu factor: The plus sign is due to the fact that the DA# variables are for the covariant surface integration
					double S_p = prefactor*(p0*da0 + px*da1 + py*da2)*f0*(1.+deltaf);

					// store values to recycle later
					double S_p_withweight = S_p*tau*eta_s_weight[ieta];		//don't include eta_s_symmetry_factor here, for consistency with later calculations...

					//ignore points where delta f is large or emission function goes negative from pdsigma
					if ((1. + deltaf < 0.0) || (flagneg == 1 && S_p < tol))
					{
						S_p_withweight = 0.0;
						tmp_S_p_withweight_array[iFOcell] = 0.0;
						++iFOcell;			// N.B. - iFOcell == isurf * eta_s_npts + ieta
						--current_num_of_FOcells_ACA;	// might save some time, since forces FOcells_PQ to be smaller?
						continue;
					}

					double absS_p_withweight = abs(S_p_withweight);
					tmp_S_p_withweight_array[iFOcell] = S_p_withweight;
					tempsum += S_p_withweight;
					tempabssum += absS_p_withweight;
					addElementToQueue(FOcells_PQ, pair<double, size_t>(-absS_p_withweight, iFOcell), current_num_of_FOcells_ACA);
					++iFOcell;			// N.B. - iFOcell == isurf * eta_s_npts + ieta

				}	//end of ieta loop
			}		//end of isurf loop

			temp_moments_array[ipt][ipphi] = tempsum;
			abs_temp_moments_array[ipt][ipphi] = tempabssum;

			// start timing this part
			sw3.Start();
			
			// figure out just how many cells I need to reach cutoff
			vector<size_t> most_impt_FOcells_vec;
			vector<double> most_impt_FOcells_vals_vec;
			size_t FOcells_PQ_size = FOcells_PQ.size();

			// read out full priority_queue into vector first
			most_impt_FOcells_vec.reserve(FOcells_PQ_size);
			most_impt_FOcells_vals_vec.reserve(FOcells_PQ_size);
			for (int ii = 0; ii < FOcells_PQ_size; ii++)
			{
				pair<double, size_t> tmp = FOcells_PQ.top();
				most_impt_FOcells_vals_vec.push_back(-tmp.first);	//recall: tmp.first <= 0
				most_impt_FOcells_vec.push_back(tmp.second);
				FOcells_PQ.pop();
			}

			// get vectors in the right order
			reverse( most_impt_FOcells_vec.begin(), most_impt_FOcells_vec.end() );
			reverse( most_impt_FOcells_vals_vec.begin(), most_impt_FOcells_vals_vec.end() );

			// loop through vectors from largest FOcell contribution to smallest, terminate when cutoff is reached
			double running_sum = 0.0;
			int breaker = FOcells_PQ_size;
			double abs_cutoff = cutoff * abs_temp_moments_array[ipt][ipphi];

			int current_iPC = 0;
			cutoff_FOcells[ptphi_index].reserve( number_of_percentage_markers );
			cutoff_FOcell_vals_C[ptphi_index].reserve( number_of_percentage_markers );
			cutoff_FOcell_vals_S[ptphi_index].reserve( number_of_percentage_markers );
			cutoff_FOcells[ptphi_index].push_back(0);		//always use sum over zero FO cells as trivial point

			//set the cutoff %-age values here...
			//just makes them equally spaced from 0% to cutoff
			linspace(pc_cutoff_vals, 0.0, cutoff);

			for (int ii = 0; ii < FOcells_PQ_size; ii++)
			{
				running_sum += most_impt_FOcells_vals_vec[ii];	//starts with largest first...
				if (running_sum >= abs_cutoff)	//cutoff-dependence only enters here
				{
					breaker = ii + 1;	//marks where the final cutoff was reached
					break;
				}
				else if (running_sum > pc_cutoff_vals[current_iPC + 1] * tempabssum)
				{
					++current_iPC;
					cutoff_FOcells[ptphi_index].push_back(ii+1);
				}
			}
			cutoff_FOcells[ptphi_index].push_back(breaker);

			// finally, chop off whatever is left
			most_impt_FOcells_vec.erase ( most_impt_FOcells_vec.begin() + breaker, most_impt_FOcells_vec.end() );
			most_impt_FOcells_vals_vec.erase ( most_impt_FOcells_vals_vec.begin() + breaker, most_impt_FOcells_vals_vec.end() );

			//!!! - RESIZING THE CUTOFF COUNTER TO HOLD EXACTLY ONLY AS MANY CELLS AS NECESSARY
			number_of_FOcells_above_cutoff_array[ipt][ipphi] = most_impt_FOcells_vec.size();
			FOcells_PQ_size = number_of_FOcells_above_cutoff_array[ipt][ipphi];

			// define this array so that it's smaller at run-time...
			S_p_withweight_array[ipt][ipphi] = new double [FOcells_PQ_size];
			most_important_FOcells[ipt][ipphi] = new size_t [FOcells_PQ_size];

			// copy over sorted results
			for (int ii = 0; ii < FOcells_PQ_size; ++ii)
			{
				size_t topFOcell = most_impt_FOcells_vec[ii];
				most_important_FOcells[ipt][ipphi][ii] = topFOcell;
				S_p_withweight_array[ipt][ipphi][ii] = tmp_S_p_withweight_array[topFOcell];
			}

			//end of the section to time
			sw3.Stop();

			delete [] tmp_S_p_withweight_array;
		}	// end of pphi loop
	}		// end of pt loop
*global_out_stream_ptr << "\t\t\t*** Took total of " << sw3.printTime() << " seconds on ordering and copying." << endl;

	cerr << "Working with the following %-age cutoffs:" << endl;
	for (int ii = 0; ii < number_of_percentage_markers; ++ii)
		cerr << ii << "   " << 100.0 * pc_cutoff_vals[ii] << "   " << cutoff_FOcells[0][ii] << " of " << number_of_FOcells_above_cutoff_array[0][0] << endl;

sw2.Stop();
*global_out_stream_ptr << "\t\t\t*** Took " << sw2.printTime() << " seconds for whole function." << endl;

	for(int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
	for(int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
	{
		spectra[local_pid][ipt][ipphi] = eta_s_symmetry_factor * temp_moments_array[ipt][ipphi];
		abs_spectra[local_pid][ipt][ipphi] = abs_temp_moments_array[ipt][ipphi];
	}

	for(int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
	{
		delete [] temp_moments_array[ipt];
		delete [] abs_temp_moments_array[ipt];
	}
	delete [] temp_moments_array;
	delete [] abs_temp_moments_array;

//if(1)
//	exit(1);

	return;
}


// this function sets the FT phase factor which depends on all q pts. and all x pts., but no K or p pts.
void CorrelationFunction::Set_giant_arrays(int iqt, int iqx, int iqy, int iqz)
{
	giant_array_C = new double [FO_length * eta_s_npts];
	giant_array_S = new double [FO_length * eta_s_npts];

	double * tmp_results_ii0 = new double [ntrig];
	double * tmp_results_ii1 = new double [ntrig];

	int iFO = 0;
	for (int isurf = 0; isurf < FO_length; ++isurf)
	for (int ieta = 0; ieta < eta_s_npts; ++ieta)
	{
		form_trig_sign_z(isurf, ieta, iqt, iqx, iqy, iqz, 0, tmp_results_ii0);
		form_trig_sign_z(isurf, ieta, iqt, iqx, iqy, iqz, 1, tmp_results_ii1);
		giant_array_C[iFO] = tmp_results_ii0[0] + tmp_results_ii1[0];
		giant_array_S[iFO] = tmp_results_ii0[1] + tmp_results_ii1[1];
		++iFO;		//N.B. - iFO == isurf * eta_s_npts + ieta
	}

	delete [] tmp_results_ii0;
	delete [] tmp_results_ii1;

	return;
}

void CorrelationFunction::Cal_dN_dypTdpTdphi_with_weights(FO_surf* FOsurf_ptr, int local_pid)
{
	CPStopwatch debug_sw, debug_sw2, sw_set_giant_array_slice;
	debug_sw.Start();

	int temp_moms_lin_arr_length = n_interp_pT_pts * n_interp_pphi_pts * ntrig;

	int lin_TMLAL_idx = 0;

	double * temp_moms_linear_array = new double [temp_moms_lin_arr_length];

	for (int iqt = 0; iqt < (qtnpts / 2) + 1; ++iqt)	//assumes each central q point is zero!!!
	for (int iqx = 0; iqx < qxnpts; ++iqx)
	for (int iqy = 0; iqy < qynpts; ++iqy)
	for (int iqz = 0; iqz < qznpts; ++iqz)
	{	
		*global_out_stream_ptr << "Working on (iqt, iqx, iqy, iqz) = (" << iqt << ", " << iqx << ", " << iqy << ", " << iqz << ") DIRECTLY..." << endl;
		debug_sw2.Reset();
	
		sw_set_giant_array_slice.Start();
		Set_giant_arrays(iqt, iqx, iqy, iqz);
		sw_set_giant_array_slice.Stop();

		for (int i = 0; i < temp_moms_lin_arr_length; ++i)
			temp_moms_linear_array[i] = 0.0;

		for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
		{
			//if (ipt > 0) continue;
			double ** slice1 = S_p_withweight_array[ipt];
			size_t ** mif1 = most_important_FOcells[ipt];
			for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
			{
				//if (ipphi > 0) continue;
				double * slice2 = slice1[ipphi];
				size_t * most_important_FOcells_for_current_pt_and_pphi = mif1[ipphi];

				int ptphi_index = ipt * n_interp_pphi_pts + ipphi;
				int maxFOnum = number_of_FOcells_above_cutoff_array[ipt][ipphi];
				double running_sum = 0.0;
				double tmla_C = 0.0, tmla_S = 0.0;
				double current_abs_spectra = abs_spectra[local_pid][ipt][ipphi];
				vector<double> runsumvals;
				vector<int> tmpvec = cutoff_FOcells[ptphi_index];
	
				for (int ipc = 0; ipc < tmpvec.size() - 1; ++ipc)
				{
					cutoff_FOcell_vals_C[ptphi_index].push_back(tmla_C);
					cutoff_FOcell_vals_S[ptphi_index].push_back(tmla_S);
					runsumvals.push_back(running_sum);

					for (int iFO = tmpvec[ipc]; iFO < tmpvec[ipc+1]; ++iFO)
					{
						size_t next_most_important_FOindex = most_important_FOcells_for_current_pt_and_pphi[iFO];
						double S_p_withweight = slice2[iFO];
						//if (iqt==0 && iqx==0 && iqy==0 && iqz==0)
						//{
						//	cerr << "iFO = " << iFO << endl;
						//	cerr << "next_most_important_FOindex = " << next_most_important_FOindex << endl;
						//}

						tmla_C += giant_array_C[next_most_important_FOindex] * S_p_withweight;
						tmla_S += giant_array_S[next_most_important_FOindex] * S_p_withweight;
	
						running_sum += abs(S_p_withweight);
					}
				}
				cutoff_FOcell_vals_C[ptphi_index].push_back(tmla_C);
				cutoff_FOcell_vals_S[ptphi_index].push_back(tmla_S);
				runsumvals.push_back(running_sum);

				//calculate ***PROJECTED*** tmla_C and tmla_S
				double chisqC = 0.0, chisqS = 0.0;
				double proj_tmla_C = gsl_polynomial_fit(pc_cutoff_vals, cutoff_FOcell_vals_C[ptphi_index], 4, chisqC);
				double proj_tmla_S = gsl_polynomial_fit(pc_cutoff_vals, cutoff_FOcell_vals_S[ptphi_index], 4, chisqS);
				//double proj_tmla_C = 0.0;
				//ratint(pc_cutoff_vals, cutoff_FOcell_vals_C[ptphi_index], 1.0, &proj_tmla_C, &chisqC);
				//double proj_tmla_S = 0.0;
				//ratint(pc_cutoff_vals, cutoff_FOcell_vals_S[ptphi_index], 1.0, &proj_tmla_S, &chisqS);

				temp_moms_linear_array[ntrig * ptphi_index + 0] = proj_tmla_C;
				temp_moms_linear_array[ntrig * ptphi_index + 1] = proj_tmla_S;

				//cutoff_FOcells[ptphi_index].clear();
				cutoff_FOcell_vals_C[ptphi_index].clear();
				cutoff_FOcell_vals_S[ptphi_index].clear();
if (ipt==0 && ipphi == 0 && iqt==0 && iqx==0 && iqy==0 && iqz==0)
	for (int ii = 0; ii < runsumvals.size(); ++ii)
		cerr << "Now in Cal_dN_dypTdpTdphi_with_weights(): " << ii << "   " << runsumvals[ii]/current_abs_spectra << endl;
			}	//end of pphi-loop
		}		//end of pt-loop
	
		int iidx = 0;
		for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
		for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
		for (int itrig = 0; itrig < ntrig; ++itrig)
		{
			current_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx][iqy][iqz][itrig] = temp_moms_linear_array[iidx];
			++iidx;
		}

		// Clean up
		delete [] giant_array_C;
		delete [] giant_array_S;

		debug_sw2.Stop();
		*global_out_stream_ptr << "   --> Finished (iqt, iqx, iqy, iqz) = ("
								<< iqt << ", " << iqx << ", " << iqy << ", " << iqz
								<< ") in " << debug_sw2.printTime() << " seconds." << endl;
	}		//end of first set of q-loops

	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// !!!!!  USE SYMMETRY IN Q-SPACE TO GET SPECTRA FOR ALL +VE QT POINTS FROM -VE QT POINTS  !!!!!
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
	for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
	for (int iqt = (qtnpts / 2) + 1; iqt < qtnpts; ++iqt)	//assumes each central q point is zero!!!
	for (int iqx = 0; iqx < qxnpts; ++iqx)
	for (int iqy = 0; iqy < qynpts; ++iqy)
	for (int iqz = 0; iqz < qznpts; ++iqz)
	for (int itrig = 0; itrig < ntrig; ++itrig)
	{
		//*global_out_stream_ptr << "Working on (iqt, iqx, iqy, iqz) = (" << iqt << ", " << iqx << ", " << iqy << ", " << iqz << ") INDIRECTLY..." << endl;
		current_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx][iqy][iqz][itrig]
			= (1.0 - 2.0 * itrig) * current_dN_dypTdpTdphi_moments[ipt][ipphi][qtnpts - iqt - 1][qxnpts - iqx - 1][qynpts - iqy - 1][qznpts - iqz - 1][itrig];
	}		//end of second set of q-loops

	*global_out_stream_ptr << "Spent " << sw_set_giant_array_slice.printTime() << " seconds setting giant_array_slice." << endl;

	for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
	for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
	{
		delete [] S_p_withweight_array[ipt][ipphi];
		delete [] most_important_FOcells[ipt][ipphi];
		int ptphi_index = ipt * n_interp_pphi_pts + ipphi;
		cutoff_FOcells[ptphi_index].clear();
	}

	delete [] temp_moms_linear_array;

	debug_sw.Stop();
	*global_out_stream_ptr << "Total function call took " << debug_sw.printTime() << " seconds." << endl;

	return;
}

void CorrelationFunction::Load_decay_channel_info(int dc_idx, double K_T_local, double K_phi_local)
{
	Mres = current_resonance_mass;
	Gamma = current_resonance_Gamma;
	double one_by_Gamma_Mres = hbarC/(Gamma*Mres);
	mass = current_daughter_mass;
	br = current_resonance_direct_br;	//doesn't depend on target daughter particle, just parent resonance and decay channel
	m2 = current_resonance_decay_masses[0];
	m3 = current_resonance_decay_masses[1];
	pT = K_T_local;
	current_K_phi = K_phi_local;
	n_body = current_reso_nbody;
	if (n_body == 2)
	{
		// some particles may decay to particles with more total mass than originally
		// --> broaden with resonance widths
		while ((mass + m2) > Mres)
		{
			Mres += 0.25 * current_resonance_Gamma;
			mass -= 0.5 * current_daughter_Gamma;
			m2 -= 0.5 * current_m2_Gamma;
		}

		mT = sqrt(mass*mass + pT*pT);

		//set up vectors of points to speed-up integrals...
		double s_loc = m2*m2;
		VEC_n2_spt = s_loc;
		double pstar_loc = sqrt( ((Mres+mass)*(Mres+mass) - s_loc)*((Mres-mass)*(Mres-mass) - s_loc) )/(2.0*Mres);
		VEC_n2_pstar = pstar_loc;
		check_for_NaNs("pstar_loc", pstar_loc, *global_out_stream_ptr);
		double g_s_loc = g(s_loc);	//for n_body == 2, doesn't actually use s_loc since result is just a factor * delta(...); just returns factor
		VEC_n2_s_factor = br/(4.*M_PI*VEC_n2_pstar);	//==g_s_loc
		double Estar_loc = sqrt(mass*mass + pstar_loc*pstar_loc);
		VEC_n2_Estar = Estar_loc;
		double psBmT = pstar_loc / mT;
		VEC_n2_psBmT = psBmT;
		double DeltaY_loc = log(psBmT + sqrt(1.+psBmT*psBmT));
		VEC_n2_DeltaY = DeltaY_loc;
		p_y = 0.0;
		VEC_n2_Yp = p_y + DeltaY_loc;
		VEC_n2_Ym = p_y - DeltaY_loc;
		check_for_NaNs("DeltaY_loc", DeltaY_loc, *global_out_stream_ptr);
		for(int iv = 0; iv < n_v_pts; ++iv)
		{
			double v_loc = v_pts[iv];
			double P_Y_loc = p_y + v_loc*DeltaY_loc;
			VEC_n2_P_Y[iv] = P_Y_loc;
			double mT_ch_P_Y_p_y = mT*cosh(v_loc*DeltaY_loc);
			double x2 = mT_ch_P_Y_p_y*mT_ch_P_Y_p_y - pT*pT;
			VEC_n2_v_factor[iv] = v_wts[iv]*DeltaY_loc/sqrt(x2);
			double MTbar_loc = Estar_loc*Mres*mT_ch_P_Y_p_y/x2;
			VEC_n2_MTbar[iv] = MTbar_loc;
			double DeltaMT_loc = Mres*pT*sqrt(Estar_loc*Estar_loc - x2)/x2;
			VEC_n2_DeltaMT[iv] = DeltaMT_loc;
			VEC_n2_MTp[iv] = MTbar_loc + DeltaMT_loc;
			VEC_n2_MTm[iv] = MTbar_loc - DeltaMT_loc;
			check_for_NaNs("MTbar_loc", MTbar_loc, *global_out_stream_ptr);
			check_for_NaNs("DeltaMT_loc", DeltaMT_loc, *global_out_stream_ptr);

			for(int izeta = 0; izeta < n_zeta_pts; ++izeta)
			{
				double zeta_loc = zeta_pts[izeta];
				double MT_loc = MTbar_loc + cos(zeta_loc)*DeltaMT_loc;
				VEC_n2_MT[iv][izeta] = MT_loc;
				VEC_n2_zeta_factor[iv][izeta] = zeta_wts[izeta]*MT_loc;
				double PT_loc = sqrt(MT_loc*MT_loc - Mres*Mres);
				double temp_cos_PPhi_tilde_loc = (mT*MT_loc*cosh(P_Y_loc-p_y) - Estar_loc*Mres)/(pT*PT_loc);
				//assume that PPhi_tilde is +ve in next step...
				double temp_sin_PPhi_tilde_loc = sqrt(1. - temp_cos_PPhi_tilde_loc*temp_cos_PPhi_tilde_loc);
				double PPhi_tilde_loc = place_in_range( atan2(temp_sin_PPhi_tilde_loc, temp_cos_PPhi_tilde_loc), interp_pphi_min, interp_pphi_max);
				VEC_n2_PPhi_tilde[iv][izeta] = place_in_range( K_phi_local + PPhi_tilde_loc, interp_pphi_min, interp_pphi_max);
				VEC_n2_PPhi_tildeFLIP[iv][izeta] = place_in_range( K_phi_local - PPhi_tilde_loc, interp_pphi_min, interp_pphi_max);
				VEC_n2_PT[iv][izeta] = PT_loc;
				check_for_NaNs("PT_loc", PT_loc, *global_out_stream_ptr);
				check_for_NaNs("PPhi_tilde_loc", PPhi_tilde_loc, *global_out_stream_ptr);
			}
		}
	}
	else
	{
		mT = sqrt(mass*mass + pT*pT);
		double s_min_temp = (m2 + m3)*(m2 + m3);
		double s_max_temp = (Mres - mass)*(Mres - mass);
		gauss_quadrature(n_s_pts, 1, 0.0, 0.0, s_min_temp, s_max_temp, NEW_s_pts, NEW_s_wts);
		Qfunc = get_Q();
		for (int is = 0; is < n_s_pts; ++is)
		{
			double s_loc = NEW_s_pts[is];
			double g_s_loc = g(s_loc);
			VEC_g_s[is] = g_s_loc;
			VEC_s_factor[is] = NEW_s_wts[is]*g_s_loc;
			double pstar_loc = sqrt(((Mres+mass)*(Mres+mass) - s_loc)*((Mres-mass)*(Mres-mass) - s_loc))/(2.0*Mres);
			VEC_pstar[is] = pstar_loc;
			double Estar_loc = sqrt(mass*mass + pstar_loc*pstar_loc);
			VEC_Estar[is] = Estar_loc;
			double psBmT = pstar_loc / mT;
			double DeltaY_loc = log(psBmT + sqrt(1.+psBmT*psBmT));
			VEC_DeltaY[is] = DeltaY_loc;
			p_y = 0.0;
			VEC_Yp[is] = p_y + DeltaY_loc;
			VEC_Ym[is] = p_y - DeltaY_loc;
			for(int iv = 0; iv < n_v_pts; ++iv)
			{
				double v_loc = v_pts[iv];
				double P_Y_loc = p_y + v_loc*DeltaY_loc;
				VEC_P_Y[is][iv] = P_Y_loc;
				double mT_ch_P_Y_p_y = mT*cosh(v_loc*DeltaY_loc);
				double x2 = mT_ch_P_Y_p_y*mT_ch_P_Y_p_y - pT*pT;
				VEC_v_factor[is][iv] = v_wts[iv]*DeltaY_loc/sqrt(x2);
				double MTbar_loc = Estar_loc*Mres*mT_ch_P_Y_p_y/x2;
				VEC_MTbar[is][iv] = MTbar_loc;
				double DeltaMT_loc = Mres*pT*sqrt(Estar_loc*Estar_loc - x2)/x2;
				VEC_DeltaMT[is][iv] = DeltaMT_loc;
				VEC_MTp[is][iv] = MTbar_loc + DeltaMT_loc;
				VEC_MTm[is][iv] = MTbar_loc - DeltaMT_loc;
				for(int izeta = 0; izeta < n_zeta_pts; ++izeta)
				{
					double zeta_loc = zeta_pts[izeta];
					double MT_loc = MTbar_loc + cos(zeta_loc)*DeltaMT_loc;
					VEC_MT[is][iv][izeta] = MT_loc;
					VEC_zeta_factor[is][iv][izeta] = zeta_wts[izeta]*MT_loc;
					double PT_loc = sqrt(MT_loc*MT_loc - Mres*Mres);
					double temp_cos_PPhi_tilde_loc = (mT*MT_loc*cosh(P_Y_loc-p_y) - Estar_loc*Mres)/(pT*PT_loc);
					//assume that PPhi_tilde is +ve in next step...
					double temp_sin_PPhi_tilde_loc = sqrt(1. - temp_cos_PPhi_tilde_loc*temp_cos_PPhi_tilde_loc);
					double PPhi_tilde_loc = place_in_range( atan2(temp_sin_PPhi_tilde_loc, temp_cos_PPhi_tilde_loc), interp_pphi_min, interp_pphi_max);
					VEC_PPhi_tilde[is][iv][izeta] = place_in_range( K_phi_local + PPhi_tilde_loc, interp_pphi_min, interp_pphi_max);
					VEC_PPhi_tildeFLIP[is][iv][izeta] = place_in_range( K_phi_local - PPhi_tilde_loc, interp_pphi_min, interp_pphi_max);
					VEC_PT[is][iv][izeta] = PT_loc;
				}
			}
		}
	}

	return;
}

//***************************************************************************************************

double CorrelationFunction::Cal_dN_dypTdpTdphi_function(FO_surf* FOsurf_ptr, int local_pid, double pT, double pphi)
{
	// set particle information
	double sign = all_particles[local_pid].sign;
	double degen = all_particles[local_pid].gspin;
	double localmass = all_particles[local_pid].mass;
	double mu = all_particles[local_pid].mu;

	// set some freeze-out surface information that's constant the whole time
	double prefactor = 1.0*degen/(8.0*M_PI*M_PI*M_PI)/(hbarC*hbarC*hbarC);
	double Tdec = (&FOsurf_ptr[0])->Tdec;
	double Pdec = (&FOsurf_ptr[0])->Pdec;
	double Edec = (&FOsurf_ptr[0])->Edec;
	double one_by_Tdec = 1./Tdec;
	double deltaf_prefactor = 0.;
	if (use_delta_f)
		deltaf_prefactor = 1./(2.0*Tdec*Tdec*(Edec+Pdec));

	// set the rapidity-integration symmetry factor
	double eta_odd_factor = 1.0, eta_even_factor = 1.0;
	if (ASSUME_ETA_SYMMETRIC)
	{
		eta_odd_factor = 0.0;
		eta_even_factor = 2.0;
	}

	double sin_pphi = sin(pphi);
	double cos_pphi = cos(pphi);
	double px = pT*cos_pphi;
	double py = pT*sin_pphi;

	double dN_dypTdpTdphi = 0.0;

	for(int isurf=0; isurf<FO_length; ++isurf)
	{
		FO_surf*surf = &FOsurf_ptr[isurf];

		double tau = surf->tau;
		double r = surf->r;
		double sin_temp_phi = surf->sin_phi;
		double cos_temp_phi = surf->cos_phi;

		double vx = surf->vx;
		double vy = surf->vy;
		double gammaT = surf->gammaT;

		double da0 = surf->da0;
		double da1 = surf->da1;
		double da2 = surf->da2;

		double pi00 = surf->pi00;
		double pi01 = surf->pi01;
		double pi02 = surf->pi02;
		double pi11 = surf->pi11;
		double pi12 = surf->pi12;
		double pi22 = surf->pi22;
		double pi33 = surf->pi33;

		for(int ieta=0; ieta < eta_s_npts; ++ieta)
		{
			double p0 = sqrt(pT*pT+localmass*localmass)*cosh(SP_p_y - eta_s[ieta]);
			double pz = sqrt(pT*pT+localmass*localmass)*sinh(SP_p_y - eta_s[ieta]);
			double expon, f0;
	
			//now get distribution function, emission function, etc.
			if (TRUNCATE_COOPER_FRYE)
			{
				expon = (gammaT*(p0*1. - px*vx - py*vy) - mu)*one_by_Tdec;
				if (expon > 20.) continue;
				f0 = 1./(exp(expon)+sign);	//thermal equilibrium distributions
			}
			else
				f0 = 1./(exp( one_by_Tdec*(gammaT*(p0*1. - px*vx - py*vy) - mu) )+sign);	//thermal equilibrium distributions
	
			//viscous corrections
			double deltaf = 0.;
			if (use_delta_f)
				deltaf = deltaf_prefactor * (1. - sign*f0)
							* (p0*p0*pi00 - 2.0*p0*px*pi01 - 2.0*p0*py*pi02 + px*px*pi11 + 2.0*px*py*pi12 + py*py*pi22 + pz*pz*pi33);

			//p^mu d^3sigma_mu factor: The plus sign is due to the fact that the DA# variables are for the covariant surface integration
			double S_p = prefactor*(p0*da0 + px*da1 + py*da2)*f0*(1.+deltaf);

			//ignore points where delta f is large or emission function goes negative from pdsigma
			if ((1. + deltaf < 0.0) || (flagneg == 1 && S_p < tol))
				S_p = 0.0;

			dN_dypTdpTdphi += eta_even_factor*S_p*tau*eta_s_weight[ieta];
		}
	}

	return dN_dypTdpTdphi;
}

//performs extrapolation of running_sum of FO integrals to unity (1) by polynomial fit
double CorrelationFunction::gsl_polynomial_fit(const vector<double> &data_x, const vector<double> &data_y, const int order, double & chisq, bool verbose /* == false*/)
{
	const int n = data_x.size();
	double * in_data_x = new double [n];
	double * in_data_y = new double [n];
	gsl_vector *y, *c;
	gsl_matrix *X, *cov;
	y = gsl_vector_alloc (n);
	c = gsl_vector_alloc (order+1);
	X   = gsl_matrix_alloc (n, order+1);
	cov = gsl_matrix_alloc (order+1, order+1);

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < order+1; j++)
		{
			in_data_x[i] = data_x[i];
			gsl_matrix_set (X, i, j, pow(in_data_x[i],j));
		}
		in_data_y[i] = data_y[i];
		gsl_vector_set (y, i, in_data_y[i]);
	}

	gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc (n, order+1);
	gsl_multifit_linear (X, y, c, cov, &chisq, work);
	gsl_multifit_linear_free (work);

	vector<double> vc;
	for (int i = 0; i < order+1; i++)
	{
		vc.push_back(gsl_vector_get(c,i));
		if (verbose) cerr << "In gsl_polynomial_fit(): vc[" << i << "] = " << vc[i] << endl;
	}

	gsl_vector_free (y);
	gsl_vector_free (c);
	gsl_matrix_free (X);
	gsl_matrix_free (cov);

	return ( accumulate(vc.begin(), vc.end(), 0.0) );
}

//End of file
