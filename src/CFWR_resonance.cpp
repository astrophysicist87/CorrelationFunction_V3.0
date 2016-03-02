#include<iostream>
#include<sstream>
#include<string>
#include<fstream>
#include<cmath>
#include<iomanip>
#include<vector>
#include<stdio.h>
#include<time.h>

#include "CFWR.h"
#include "Arsenal.h"
#include "gauss_quadrature.h"

using namespace std;

double CorrelationFunction::get_Q()
{
	double smin = (m2+m3)*(m2+m3);
	double smax = (Mres-mass)*(Mres-mass);
	double sum = 0.;
	
	for (int is = 0; is < n_s_pts; ++is)
	{
		double sp = NEW_s_pts[is];
		double f1 = (Mres+mass)*(Mres+mass) - sp;
		double f2 = smax - sp;
		double f3 = smin - sp;
		double f4 = (m2-m3)*(m2-m3) - sp;
		sum += NEW_s_wts[is]*sqrt(f1*f2*f3*f4)/(sp+1.e-15);
	}

	return sum;
}

double CorrelationFunction::g(double s)
{
	double gs_pstar_loc = sqrt( ((Mres+mass)*(Mres+mass) - s)*((Mres-mass)*(Mres-mass) - s) )/(2.0*Mres);
	double g_res = br/(4.*M_PI*gs_pstar_loc);
	if (n_body == 3 || n_body == 4)		//both set up to work the same way
	{
		double pre_f = (Mres*br)/(2.*M_PI*s);
		double num = sqrt( (s - (m2+m3)*(m2+m3)) * (s - (m2-m3)*(m2-m3)) );
		double den = Qfunc;
		g_res = pre_f * num / den;
	}

	return g_res;
}

void CorrelationFunction::Do_resonance_integrals(int parent_resonance_particle_id, int daughter_particle_id, int decay_channel)
{
	time_t rawtime;
  	struct tm * timeinfo;

	int daughter_lookup_idx = distance(daughter_resonance_indices.begin(), daughter_resonance_indices.find(daughter_particle_id));

	Allocate_resonance_running_sum_vectors();

	Flatten_dN_dypTdpTdphi_moments();

	//set these for quick look-up in EdNd3p() routine
	spec_vals_info = spectra[parent_resonance_particle_id];
	spec_log_info = log_spectra[parent_resonance_particle_id];
	spec_sign_info = sign_spectra[parent_resonance_particle_id];

	int tmp_parent_monval = all_particles[parent_resonance_particle_id].monval;
	n_body = current_reso_nbody;
	//double tmp_spec_save = current_daughters_dN_dypTdpTdphi_moments[daughter_lookup_idx][0][0][0][0][0][0][0];

	//double ssumarray[n_s_pts], vsumarray[n_v_pts], zetasumarray[n_zeta_pts], Csumarray[];

	if (n_body == 2)
	{
		for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
		for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
		{
			double local_pT = SPinterp_pT[ipt];
			double local_pphi = SPinterp_pphi[ipphi];
			current_ipt = ipt;
			current_ipphi = ipphi;
			current_qlist_slice = qlist[ipt];
			Zero_resonance_running_sum_vector(ssum_vec);
			Zero_resonance_running_sum_vector(vsum_vec);
			Zero_resonance_running_sum_vector(zetasum_vec);
			Zero_resonance_running_sum_vector(Csum_vec);
			Load_decay_channel_info(decay_channel, local_pT, local_pphi);	// set decay channel information

			//then g(s) is delta-function, skip s-integration entirely
			//double s_loc = m2*m2;
			double vsum = 0.0;
			for (int iv = 0; iv < n_v_pts; ++iv)
			{
				time (&rawtime);
				timeinfo = localtime (&rawtime);
				Zero_resonance_running_sum_vector(zetasum_vec);
				double zetasum = 0.0;
				for (int izeta = 0; izeta < n_zeta_pts; ++izeta)
				{
					Zero_resonance_running_sum_vector(Csum_vec);
					double Csum = 0.0;
					double PKT = VEC_n2_PT[iv][izeta];
					double PKY = VEC_n2_P_Y[iv];
					double PKphi = VEC_n2_PPhi_tilde[iv][izeta];
					for (int tempidx = 0; tempidx <= 1; ++tempidx)
					{
						if (tempidx != 0)
							PKphi = VEC_n2_PPhi_tildeFLIP[iv][izeta];		//also takes Pp --> Pm
						currentPpm = VEC_n2_Ppm[iv][izeta][tempidx];
						Csum += Edndp3(PKT, PKphi);							//now set spectra and weights separately
						eiqxEdndp3(PKT, PKphi, Csum_vec);
					}												// end of tempidx sum
					for (int qpt_cs_idx = 0; qpt_cs_idx < qspace_cs_slice_length; ++qpt_cs_idx)
						zetasum_vec[qpt_cs_idx] += VEC_n2_zeta_factor[iv][izeta]*Csum_vec[qpt_cs_idx];
					zetasum += VEC_n2_zeta_factor[iv][izeta]*Csum;
				}													// end of zeta sum
				for (int qpt_cs_idx = 0; qpt_cs_idx < qspace_cs_slice_length; ++qpt_cs_idx)
					vsum_vec[qpt_cs_idx] += VEC_n2_v_factor[iv]*zetasum_vec[qpt_cs_idx];
				vsum += VEC_n2_v_factor[iv]*zetasum;
			}														// end of v sum
			for (int qpt_cs_idx = 0; qpt_cs_idx < qspace_cs_slice_length; ++qpt_cs_idx)
				ssum_vec[qpt_cs_idx] += Mres*VEC_n2_s_factor*vsum_vec[qpt_cs_idx];
			double ssum = Mres*VEC_n2_s_factor*vsum;

			//update all gridpoints for all daughter moments
			int qpt_cs_idx = 0;
			for (int iqt = 0; iqt < qtnpts; ++iqt)
			for (int iqx = 0; iqx < qxnpts; ++iqx)
			for (int iqy = 0; iqy < qynpts; ++iqy)
			for (int iqz = 0; iqz < qznpts; ++iqz)
			for (int itrig = 0; itrig < 2; ++itrig)
			{
				//current_daughters_dN_dypTdpTdphi_moments[daughter_lookup_idx][ipt][ipphi][iqt][iqx][iqy][iqz][itrig] += ssum_vec[qpt_cs_idx] / fraction_of_resonances;
				current_daughters_dN_dypTdpTdphi_moments[daughter_lookup_idx][ipt][ipphi][iqt][iqx][iqy][iqz][itrig] += ssum_vec[qpt_cs_idx];
				//double temp = current_daughters_dN_dypTdpTdphi_moments[daughter_lookup_idx][ipt][ipphi][iqt][iqx][iqy][iqz][itrig];
				++qpt_cs_idx;
			}

			//update daughter spectra separately
			spectra[daughter_particle_id][ipt][ipphi] += ssum;
			log_spectra[daughter_particle_id][ipt][ipphi] = log(abs(spectra[daughter_particle_id][ipt][ipphi])+1.e-100);
			sign_spectra[daughter_particle_id][ipt][ipphi] = sgn(spectra[daughter_particle_id][ipt][ipphi]);

/*
			//now, if ignoring long-lived resonances, take them out of the correlator numerator, but keep them in the denominator (AKA, spectra)
			if (IGNORE_LONG_LIVED_RESONANCES && Gamma < 1.e-3)		//i.e., the eta and eta'(958) mesons
			{
				*global_out_stream_ptr << "RESONANCES: Leaving out long-lived resonance " << all_particles[parent_resonance_particle_id].name << " from q=0!" << endl;
				int tmp_qpt_cs_idx = ( ( ( iqt0 * qxnpts + iqx0 ) * qynpts + iqy0 ) * qznpts + iqz0 ) * 2;
				current_daughters_dN_dypTdpTdphi_moments[daughter_lookup_idx][ipt][ipphi][iqt0][iqx0][iqy0][iqz0][0] -= ssum_vec[tmp_qpt_cs_idx];
				current_daughters_dN_dypTdpTdphi_moments[daughter_lookup_idx][ipt][ipphi][iqt0][iqx0][iqy0][iqz0][1] -= ssum_vec[tmp_qpt_cs_idx+1];
			}
*/
	
			if (isnan(current_daughters_dN_dypTdpTdphi_moments[daughter_lookup_idx][ipt][ipphi][0][0][0][0][0]
					+ current_daughters_dN_dypTdpTdphi_moments[daughter_lookup_idx][ipt][ipphi][0][0][0][0][1]))
			{
				*global_out_stream_ptr << "ERROR: NaNs encountered!" << endl
										<< "current_daughters_dN_dypTdpTdphi_moments[" << daughter_lookup_idx << "][" << ipt << "][" << ipphi << "][0][0][0][0][0] = "
										<< setw(8) << setprecision(15)
										<< current_daughters_dN_dypTdpTdphi_moments[daughter_lookup_idx][ipt][ipphi][0][0][0][0][0] << endl
										<< "current_daughters_dN_dypTdpTdphi_moments[" << daughter_lookup_idx << "][" << ipt << "][" << ipphi << "][0][0][0][0][1] = "
										<< setw(8) << setprecision(15)
										<< current_daughters_dN_dypTdpTdphi_moments[daughter_lookup_idx][ipt][ipphi][0][0][0][0][1] << endl
										<< "  --> pt = " << local_pT << std::endl
										<< "  --> pphi = " << local_pphi << std::endl
										<< "daughter_particle_id = " << daughter_particle_id << endl
										<< "parent_resonance_particle_id = " << parent_resonance_particle_id << endl
										<< "  --> Qfunc = " << Qfunc << endl
										<< "  --> n_body = " << n_body << endl
										<< "  --> gRES = " << gRES << endl
										<< "  --> Mres = " << Mres << endl
										<< "  --> mass = " << mass << endl
										<< "  --> Gamma = " << Gamma << endl
										<< "  --> br = " << br << endl
										<< "  --> m2 = " << m2 << endl
										<< "  --> m3 = " << m3 << endl << endl;
				exit(1);
			}
		}											// end of pT, pphi loops
	}												// end of nbody == 2
	else
	{
		for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
		for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
		{
			double local_pT = SPinterp_pT[ipt];
			double local_pphi = SPinterp_pphi[ipphi];
			current_ipt = ipt;
			current_ipphi = ipphi;
			current_qlist_slice = qlist[ipt];
			Zero_resonance_running_sum_vector(ssum_vec);
			Zero_resonance_running_sum_vector(vsum_vec);
			Zero_resonance_running_sum_vector(zetasum_vec);
			Zero_resonance_running_sum_vector(Csum_vec);
			Load_decay_channel_info(decay_channel, local_pT, local_pphi);	// set decay channel information

			double ssum = 0.0;
			for (int is = 0; is < n_s_pts; ++is)
			{
				double vsum = 0.0;
 		  		Zero_resonance_running_sum_vector(vsum_vec);
				for (int iv = 0; iv < n_v_pts; ++iv)
				{
					double zetasum = 0.0;
					Zero_resonance_running_sum_vector(zetasum_vec);
					for (int izeta = 0; izeta < n_zeta_pts; ++izeta)
					{
						double Csum = 0.0;
						Zero_resonance_running_sum_vector(Csum_vec);
						double PKT = VEC_PT[is][iv][izeta];
						double PKY = VEC_P_Y[is][iv];
						double PKphi = VEC_PPhi_tilde[is][iv][izeta];
						for (int tempidx = 0; tempidx <= 1; ++tempidx)
						{
							if (tempidx != 0)
								PKphi = VEC_PPhi_tildeFLIP[is][iv][izeta];		//also takes Pp --> Pm
							currentPpm = VEC_Ppm[is][iv][izeta][tempidx];
							Csum += Edndp3(PKT, PKphi);							//now set spectra and weights separately
							eiqxEdndp3(PKT, PKphi, Csum_vec);
						}										// end of tempidx sum
						for (int qpt_cs_idx = 0; qpt_cs_idx < qspace_cs_slice_length; ++qpt_cs_idx)
							zetasum_vec[qpt_cs_idx] += VEC_zeta_factor[is][iv][izeta]*Csum_vec[qpt_cs_idx];
						zetasum += VEC_zeta_factor[is][iv][izeta]*Csum;
					}											// end of zeta sum
					for (int qpt_cs_idx = 0; qpt_cs_idx < qspace_cs_slice_length; ++qpt_cs_idx)
					    vsum_vec[qpt_cs_idx] += VEC_v_factor[is][iv]*zetasum_vec[qpt_cs_idx];
					vsum += VEC_v_factor[is][iv]*zetasum;
				}												// end of v sum
				for (int qpt_cs_idx = 0; qpt_cs_idx < qspace_cs_slice_length; ++qpt_cs_idx)
					ssum_vec[qpt_cs_idx] += Mres*VEC_s_factor[is]*vsum_vec[qpt_cs_idx];
				ssum += Mres*VEC_s_factor[is]*vsum;
			}													// end of s sum
			//update all gridpoints for daughter moments
			int qpt_cs_idx = 0;
			for (int iqt = 0; iqt < qtnpts; ++iqt)
			for (int iqx = 0; iqx < qxnpts; ++iqx)
			for (int iqy = 0; iqy < qynpts; ++iqy)
			for (int iqz = 0; iqz < qznpts; ++iqz)
			for (int itrig = 0; itrig < 2; ++itrig)
			{
				//current_daughters_dN_dypTdpTdphi_moments[daughter_lookup_idx][ipt][ipphi][iqt][iqx][iqy][iqz][itrig] += ssum_vec[qpt_cs_idx] / fraction_of_resonances;
				current_daughters_dN_dypTdpTdphi_moments[daughter_lookup_idx][ipt][ipphi][iqt][iqx][iqy][iqz][itrig] += ssum_vec[qpt_cs_idx];
				//double temp = current_daughters_dN_dypTdpTdphi_moments[daughter_lookup_idx][ipt][ipphi][iqt][iqx][iqy][iqz][itrig];
				++qpt_cs_idx;
			}

			//update daughter spectra separately
			spectra[daughter_particle_id][ipt][ipphi] += ssum;
			log_spectra[daughter_particle_id][ipt][ipphi] = log(abs(spectra[daughter_particle_id][ipt][ipphi])+1.e-100);
			sign_spectra[daughter_particle_id][ipt][ipphi] = sgn(spectra[daughter_particle_id][ipt][ipphi]);

/*
			//now, if ignoring long-lived resonances, take them out of the correlator numerator, but keep them in the denominator (AKA, spectra)
			if (IGNORE_LONG_LIVED_RESONANCES && Gamma < 1.e-3)		//i.e., the eta and eta'(958) mesons
			{
				*global_out_stream_ptr << "RESONANCES: Leaving out long-lived resonance " << all_particles[parent_resonance_particle_id].name << " from q=0!" << endl;
				int tmp_qpt_cs_idx = ( ( ( iqt0 * qxnpts + iqx0 ) * qynpts + iqy0 ) * qznpts + iqz0 ) * 2;
				current_daughters_dN_dypTdpTdphi_moments[daughter_lookup_idx][ipt][ipphi][iqt0][iqx0][iqy0][iqz0][0] -= ssum_vec[tmp_qpt_cs_idx];
				current_daughters_dN_dypTdpTdphi_moments[daughter_lookup_idx][ipt][ipphi][iqt0][iqx0][iqy0][iqz0][1] -= ssum_vec[tmp_qpt_cs_idx+1];
			}
*/

			if (isnan(current_daughters_dN_dypTdpTdphi_moments[daughter_lookup_idx][ipt][ipphi][0][0][0][0][0]
					+ current_daughters_dN_dypTdpTdphi_moments[daughter_lookup_idx][ipt][ipphi][0][0][0][0][1]))
			{
				*global_out_stream_ptr << "ERROR: NaNs encountered!" << endl
										<< "current_daughters_dN_dypTdpTdphi_moments[" << daughter_lookup_idx << "][" << ipt << "][" << ipphi << "][0][0][0][0][0] = "
										<< setw(8) << setprecision(15)
										<< current_daughters_dN_dypTdpTdphi_moments[daughter_lookup_idx][ipt][ipphi][0][0][0][0][0] << endl
										<< "current_daughters_dN_dypTdpTdphi_moments[" << daughter_lookup_idx << "][" << ipt << "][" << ipphi << "][0][0][0][0][1] = "
										<< setw(8) << setprecision(15)
										<< current_daughters_dN_dypTdpTdphi_moments[daughter_lookup_idx][ipt][ipphi][0][0][0][0][1] << endl
										<< "  --> pt = " << local_pT << endl
										<< "  --> pphi = " << local_pphi << endl
										<< "daughter_particle_id = " << daughter_particle_id << endl
										<< "parent_resonance_particle_id = " << parent_resonance_particle_id << endl
										<< "  --> Qfunc = " << Qfunc << endl
										<< "  --> n_body = " << n_body << endl
										<< "  --> gRES = " << gRES << endl
										<< "  --> Mres = " << Mres << endl
										<< "  --> mass = " << mass << endl
										<< "  --> Gamma = " << Gamma << endl
										<< "  --> br = " << br << endl
										<< "  --> m2 = " << m2 << endl
										<< "  --> m3 = " << m3 << endl << endl;
				exit(1);
			}
		}								// end of pT, pphi loops
	}										// end of nbody == 3

	/*cout << "Resonance integrals (n_body = " << n_body << "): added " << scientific << setprecision(17) << setw(20)
			<< current_daughters_dN_dypTdpTdphi_moments[daughter_lookup_idx][0][0][0][0][0][0][0] - tmp_spec_save << endl
			<< "  --> pt = " << SPinterp_pT[0] << endl
			<< "  --> pphi = " << SPinterp_pphi[0] << endl
			<< "daughter_particle_id = " << daughter_particle_id << endl
			<< "parent_resonance_particle_id = " << parent_resonance_particle_id << endl
			<< "  --> Qfunc = " << Qfunc << endl
			<< "  --> n_body = " << n_body << endl
			<< "  --> gRES = " << gRES << endl
			<< "  --> Mres = " << Mres << endl
			<< "  --> mass = " << mass << endl
			<< "  --> Gamma = " << Gamma << endl
			<< "  --> br = " << br << endl
			<< "  --> m2 = " << m2 << endl
			<< "  --> m3 = " << m3 << endl << endl;*/

	// clean up
	Delete_resonance_running_sum_vectors();

	return;
}

inline void CorrelationFunction::set_to_zero(double * array, size_t arraylength)
{
	for (size_t arrayidx=0; arrayidx<arraylength; ++arrayidx) array[arrayidx] = 0.0;
}

void CorrelationFunction::Flatten_dN_dypTdpTdphi_moments()
{
	for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
	for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
	{
		// set index for looping
		int qpt_cs_idx = 0;
		
		// flatten arrays, since these are quicker to process
		for (int iqt = 0; iqt < qtnpts; ++iqt)
		for (int iqx = 0; iqx < qxnpts; ++iqx)
		for (int iqy = 0; iqy < qynpts; ++iqy)
		for (int iqz = 0; iqz < qznpts; ++iqz)
		for (int itrig = 0; itrig < 2; ++itrig)
		{
			res_sign_info[ipt][ipphi][qpt_cs_idx] = current_sign_of_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx][iqy][iqz][itrig];
			res_log_info[ipt][ipphi][qpt_cs_idx] = current_ln_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx][iqy][iqz][itrig];
			res_moments_info[ipt][ipphi][qpt_cs_idx] = current_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx][iqy][iqz][itrig];
			++qpt_cs_idx;
		}
	}
}

//End of file
