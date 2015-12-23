#include<iostream>
#include<sstream>
#include<string>
#include<fstream>
#include<cmath>
#include<iomanip>
#include<vector>
#include<stdio.h>

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

void replace_parentheses(std::string & tempstring)
{
	int len = tempstring.length();
	for(unsigned int i = 0; i < len; i++)
	{
		char c = tempstring[i];
		if (c == '(' || c == ')')
			tempstring[i] = '_';
	}
	
	if (tempstring[len - 1] == '_')
		tempstring.erase( len - 1 );
	
	return;
}

void CorrelationFunction::Output_results(int folderindex)
{
	ostringstream filename_stream_HBT;
	filename_stream_HBT << global_path << "/HBTradii_ev" << folderindex << no_df_stem << ".dat";
	ofstream outputHBT;
	outputHBT.open(filename_stream_HBT.str().c_str());
	ostringstream filename_stream_HBTcfs;
	filename_stream_HBTcfs << global_path << "/HBTradii_cfs_ev" << folderindex << no_df_stem << ".dat";
	ofstream outputHBTcoeffs(filename_stream_HBTcfs.str().c_str());

	for(int iKT = 0; iKT < n_localp_T; iKT++)
	{
		for(int Morder=0; Morder<n_order; Morder++)
		{
			outputHBTcoeffs << folderindex << "  " << K_T[iKT] << "  " << Morder
				<< "  " << R2_side_C[iKT][Morder] << "   " << R2_side_S[iKT][Morder] << "  " << R2_out_C[iKT][Morder] << "  " << R2_out_S[iKT][Morder]
				<< "  " << R2_outside_C[iKT][Morder] << "   " << R2_outside_S[iKT][Morder] << "  " << R2_long_C[iKT][Morder] << "  " << R2_long_S[iKT][Morder]
				<< "  " << R2_sidelong_C[iKT][Morder] << "   " << R2_sidelong_S[iKT][Morder] << "  " << R2_outlong_C[iKT][Morder] << "  " << R2_outlong_S[iKT][Morder] << endl;
		}
		for(int iKphi = 0; iKphi < n_localp_phi; iKphi++)
		{
			outputHBT << folderindex << "  " << K_T[iKT] << "  " << K_phi[iKphi]
				<< "  " << R2_side[iKT][iKphi] << "  " << R2_out[iKT][iKphi]
				<< "  " << R2_outside[iKT][iKphi] << "  " << R2_long[iKT][iKphi]
				<< "  " << R2_sidelong[iKT][iKphi] << "  " << R2_outlong[iKT][iKphi] << endl;
		}
	}

	outputHBT.close();
	outputHBTcoeffs.close();

	return;
}

/*
void CorrelationFunction::Output_resonance_spectra(int resonance_pid, int folderindex, double ******* resonance_spectra)
{
	ostringstream filename_stream;
	filename_stream << global_path << "/resonance_" << resonance_pid << "_spectra_ev" << folderindex << no_df_stem << ".dat";
	FILE *out = fopen(filename_stream.str().c_str(), "w");

	for (int iqt = 0; iqt < qnpts; ++iqt)
	{
		double ****** CRS_s1 = resonance_spectra[iqt];
		for (int iqx = 0; iqx < qnpts; ++iqx)
		{
			double ***** CRS_s2 = CRS_s1[iqx];
			for (int iqy = 0; iqy < qnpts; ++iqy)
			{
				double **** CRS_s3 = CRS_s2[iqy];
				for (int iqz = 0; iqz < qnpts; ++iqz)
				{
					double *** CRS_s4 = CRS_s3[iqz];
					double ** CRS_s4a = CRS_s4[0];	//itrig == 0
					double ** CRS_s4b = CRS_s4[1];	//itrig == 1
					for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
					{
						double * CRS_s5a = CRS_s4a[ipt];
						double * CRS_s5b = CRS_s4a[ipt];
						for (int iphi = 0; iphi < n_interp_pphi_pts; ++iphi)
							fprintf(out, "%lf   %lf   ", CRS_s5a[iphi], CRS_s5b[iphi]);
						fprintf(out, "\n");
					}
				}
			}
		}
	}

	return(fclose(out));
}

void CorrelationFunction::Readin_resonance_spectra(int resonance_pid, int folderindex, double ******* resonance_spectra)
{
	ostringstream filename_stream;
	filename_stream << global_path << "/resonance_" << resonance_pid << "_spectra_ev" << folderindex << no_df_stem << ".dat";
	FILE *in = fopen(filename_stream.str().c_str(), "r");

	for (int iqt = 0; iqt < qnpts; ++iqt)
	for (int iqx = 0; iqx < qnpts; ++iqx)
	for (int iqy = 0; iqy < qnpts; ++iqy)
	for (int iqz = 0; iqz < qnpts; ++iqz)
	for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
	for (int iphi = 0; iphi < n_interp_pphi_pts; ++iphi)
		fscanf(in, "%lf", &resonance_spectra[iqt][iqx][iqy][iqz][0][ipt][iphi]
                            &resonance_spectra[iqt][iqx][iqy][iqz][1][ipt][iphi]);

	return(fclose(in));
}

void CorrelationFunction::Readin_daughter_spectra(int resonance_pid, int folderindex, double ******* resonance_spectra)
{
	ostringstream filename_stream;
	filename_stream << global_path << "/resonance_" << resonance_pid << "_spectra_ev" << folderindex << no_df_stem << ".dat";
	FILE *in = fopen(filename_stream.str().c_str(), "r");

	for (int iqt = 0; iqt < qnpts; ++iqt)
	for (int iqx = 0; iqx < qnpts; ++iqx)
	for (int iqy = 0; iqy < qnpts; ++iqy)
	for (int iqz = 0; iqz < qnpts; ++iqz)
	for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
	for (int iphi = 0; iphi < n_interp_pphi_pts; ++iphi)
		fscanf(in, "%lf", &resonance_spectra[iqt][iqx][iqy][iqz][0][ipt][iphi]
                            &resonance_spectra[iqt][iqx][iqy][iqz][1][ipt][iphi]);

	return(fclose(in));
}
*/

void CorrelationFunction::Output_correlationfunction(int folderindex)
{
	ostringstream oCorrFunc_stream;
	string temp_particle_name = particle_name;
	replace_parentheses(temp_particle_name);
	oCorrFunc_stream << global_path << "/correlfunct1D" << "_" << temp_particle_name << ".dat";
	ofstream oCorrFunc;
	oCorrFunc.open(oCorrFunc_stream.str().c_str());

	for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
	for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
	for (int iqo = 0; iqo < qonpts; ++iqo)
	for (int iqs = 0; iqs < qsnpts; ++iqs)
	for (int iql = 0; iql < qlnpts; ++iql)
		oCorrFunc << scientific << setprecision(7) << setw(15)
			<< SPinterp_pT[ipt] << "   " << SPinterp_pphi[ipphi] << "   " << qo_pts[iqo] << "   "
			<< qs_pts[iqs] << "   " << ql_pts[iql] << "   " << CFvals[ipt][ipphi][iqo][iqs][iql] << endl;
				
	return;
}

void CorrelationFunction::Readin_results(int folderindex)
{
double dummy;
	ostringstream filename_stream_HBT;
	filename_stream_HBT << global_path << "/HBTradii_ev" << folderindex << no_df_stem << ".dat";
	ifstream inputHBT(filename_stream_HBT.str().c_str());

for(int iKT = 0; iKT < n_localp_T; iKT++)
{
	for(int iKphi = 0; iKphi < n_localp_phi; iKphi++)
	{
		inputHBT >> dummy;
		inputHBT >> dummy;
        	inputHBT >> dummy;
		inputHBT >> R2_side[iKT][iKphi];
		inputHBT >> R2_out[iKT][iKphi];
		inputHBT >> R2_outside[iKT][iKphi];
		inputHBT >> R2_long[iKT][iKphi];
		inputHBT >> R2_sidelong[iKT][iKphi];
		inputHBT >> R2_outlong[iKT][iKphi];
	}
}

	inputHBT.close();

	return;
}

// **************************************************
// THIS FUNCTION NEEDS TO BE CHECKED AND DEBUGGED!!!
// **************************************************
/*void CorrelationFunction::Output_all_dN_dypTdpTdphi(int folderindex)
{
	ostringstream filename_stream_all_dN_dypTdpTdphi;
	filename_stream_all_dN_dypTdpTdphi << global_path << "/all_res_dN_dypTdpTdphi_ev" << folderindex << no_df_stem << ".dat";
	ofstream output_all_dN_dypTdpTdphi(filename_stream_all_dN_dypTdpTdphi.str().c_str());
	for(int ii = 0; ii < Nparticle; ii++)
	for(int ipphi = 0; ipphi < n_interp_pphi_pts; ipphi++)
	{
		for(int ipt = 0; ipt < n_interp_pT_pts; ipt++)
			output_all_dN_dypTdpTdphi << scientific << setprecision(8) << setw(12) << dN_dypTdpTdphi_moments[ii][ipt][ipphi][0][0][0][0][0] << "   ";
		output_all_dN_dypTdpTdphi << endl;
	}
	output_all_dN_dypTdpTdphi.close();

	return;
}*/

void CorrelationFunction::Output_total_target_dN_dypTdpTdphi(int folderindex)
{
	string local_name = all_particles[target_particle_id].name;
	replace_parentheses(local_name);
	ostringstream filename_stream_target_dN_dypTdpTdphi;
	filename_stream_target_dN_dypTdpTdphi << global_path << "/total_" << local_name << "_dN_dypTdpTdphi_ev" << folderindex << no_df_stem << ".dat";
	ofstream output_target_dN_dypTdpTdphi(filename_stream_target_dN_dypTdpTdphi.str().c_str());

	for(int iphi = 0; iphi < n_interp_pphi_pts; iphi++)
	{
		for(int ipt = 0; ipt < n_interp_pT_pts; ipt++)
			output_target_dN_dypTdpTdphi << scientific << setprecision(8) << setw(12) << spectra[target_particle_id][ipt][iphi] << "   ";
		output_target_dN_dypTdpTdphi << endl;
	}

	output_target_dN_dypTdpTdphi.close();

	return;
}

void CorrelationFunction::Output_total_target_eiqx_dN_dypTdpTdphi(int folderindex)
{
	string local_name = all_particles[target_particle_id].name;
	replace_parentheses(local_name);
	ostringstream filename_stream_target_dN_dypTdpTdphi;
	filename_stream_target_dN_dypTdpTdphi << global_path << "/total_" << local_name << "_eiqx_dN_dypTdpTdphi_ev" << folderindex << no_df_stem << ".dat";
	ofstream output_target_dN_dypTdpTdphi(filename_stream_target_dN_dypTdpTdphi.str().c_str());

	int HDFloadTargetSuccess = Get_resonance_from_HDF_array(target_particle_id, current_dN_dypTdpTdphi_moments);

	for (int iqt = 0; iqt < qnpts; ++iqt)
	for (int iqx = 0; iqx < qnpts; ++iqx)
	for (int iqy = 0; iqy < qnpts; ++iqy)
	for (int iqz = 0; iqz < qnpts; ++iqz)
	for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
	for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
	{
		// addresses NaN issue in sin component when all q^{\mu} == 0
		//if (sqrt(qt_pts[iqt]*qt_pts[iqt]+qx_pts[iqx]*qx_pts[iqx]+qy_pts[iqy]*qy_pts[iqy]+qz_pts[iqz]*qz_pts[iqz]) < 1.e-12)
		//	current_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx][iqy][iqz][1] = 0.0;
		current_dN_dypTdpTdphi_moments[ipt][ipphi][iqt0][iqx0][iqy0][iqz0][1] = 0.0;

		// output all FT'd spectra
		double nonFTd_spectra = spectra[target_particle_id][ipt][ipphi];
		double cos_transf_spectra = current_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx][iqy][iqz][0];
		double sin_transf_spectra = current_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx][iqy][iqz][1];
		output_target_dN_dypTdpTdphi << scientific << setprecision(8) << setw(12)
			<< qt_pts[iqt] << "   " << qx_pts[iqx] << "   " << qy_pts[iqy] << "   " << qz_pts[iqz] << "   "
			<< SPinterp_pT[ipt] << "   " << SPinterp_pphi[ipphi] << "   " << nonFTd_spectra << "   " << cos_transf_spectra << "   " << sin_transf_spectra << "   "
			<< 1. + (cos_transf_spectra*cos_transf_spectra + sin_transf_spectra*sin_transf_spectra)/(nonFTd_spectra*nonFTd_spectra) <<  endl;
	}

	output_target_dN_dypTdpTdphi.close();

	return;
}

void CorrelationFunction::Output_chosen_resonances()
{
	ostringstream filename_stream_crf;
	filename_stream_crf << global_path << "/chosen_resonances.dat";
	ofstream output_crf(filename_stream_crf.str().c_str());

	output_crf << particle_monval << endl;
	for (int icr = 0; icr < (int)chosen_resonances.size(); icr++)
		output_crf << all_particles[chosen_resonances[icr]].monval << endl;

	output_crf.close();

	return;
}

//End of file
