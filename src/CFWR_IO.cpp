#include<iostream>
#include<sstream>
#include<string>
#include<fstream>
#include<cmath>
#include<iomanip>
#include<vector>
#include<stdio.h>

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

void CorrelationFunction::Dump_phases_to_binary(char direction, int ipt, double ** array, const int nd1, const int nd2)
{
	ostringstream filename_stream;
	filename_stream << global_path << "/q" << direction << "_pts_ipt_" << ipt << ".bin";
	ofstream out(filename_stream.str().c_str(), ios::out | ios::binary);

	vector<double> array_copy (nd1*nd2);
	int ac_idx = 0;
	for (int i1 = 0; i1 < nd1; ++i1)
	for (int i2 = 0; i2 < nd2; ++i2)
		array_copy[ac_idx++] = array[i1][i2];

	out.write(reinterpret_cast<const char*>(&array_copy[0]), array_copy.size()*sizeof(double));

	out.close();
	return;
}

void CorrelationFunction::Load_phases_from_binary(char direction, int ipt, double ** array, const int nd1, const int nd2)
{
	ostringstream filename_stream;
	filename_stream << global_path << "/q" << direction << "_pts_ipt_" << ipt << ".bin";
	ifstream in(filename_stream.str().c_str(), ios::in | ios::binary);

	vector<double> array_copy (nd1*nd2);
	in.read(reinterpret_cast<char*>(&array_copy[0]), array_copy.size()*sizeof(double));
	in.close();

	int ac_idx = 0;
	for (int i1 = 0; i1 < nd1; ++i1)
	for (int i2 = 0; i2 < nd2; ++i2)
		array[i1][i2] = array_copy[ac_idx++];

	return;
}

void CorrelationFunction::Output_results(int folderindex)
{
	ostringstream filename_stream_HBT;
	filename_stream_HBT << global_path << "/HBTradii_GF_ev" << folderindex << no_df_stem << ".dat";
	ofstream outputHBT;
	outputHBT.open(filename_stream_HBT.str().c_str());

	for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
	for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
	{
		outputHBT << SPinterp_pT[ipt] << "   " << SPinterp_pphi[ipphi]
			<< "   " << R2_side[ipt][ipphi] << "   " << R2_out[ipt][ipphi]
			<< "   " << R2_outside[ipt][ipphi] << "   " << R2_long[ipt][ipphi]
			<< "   " << R2_sidelong[ipt][ipphi] << "   " << R2_outlong[ipt][ipphi] << endl;
	}

	outputHBT.close();

	return;
}

/*
void CorrelationFunction::Output_resonance_spectra(int resonance_pid, int folderindex, double ******* resonance_spectra)
{
	ostringstream filename_stream;
	filename_stream << global_path << "/resonance_" << resonance_pid << "_spectra_ev" << folderindex << no_df_stem << ".dat";
	FILE *out = fopen(filename_stream.str().c_str(), "w");

	for (int iqt = 0; iqt < qtnpts; ++iqt)
	{
		double ****** CRS_s1 = resonance_spectra[iqt];
		for (int iqx = 0; iqx < qxnpts; ++iqx)
		{
			double ***** CRS_s2 = CRS_s1[iqx];
			for (int iqy = 0; iqy < qynpts; ++iqy)
			{
				double **** CRS_s3 = CRS_s2[iqy];
				for (int iqz = 0; iqz < qznpts; ++iqz)
				{
					double *** CRS_s4 = CRS_s3[iqz];
					double ** CRS_s4a = CRS_s4[0];	//itrig == 0
					double ** CRS_s4b = CRS_s4[1];	//itrig == 1
					for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
					{
						double * CRS_s5a = CRS_s4a[ipt];
						double * CRS_s5b = CRS_s4a[ipt];
						for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
							fprintf(out, "%lf   %lf   ", CRS_s5a[ipphi], CRS_s5b[ipphi]);
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

	for (int iqt = 0; iqt < qtnpts; ++iqt)
	for (int iqx = 0; iqx < qxnpts; ++iqx)
	for (int iqy = 0; iqy < qynpts; ++iqy)
	for (int iqz = 0; iqz < qznpts; ++iqz)
	for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
	for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
		fscanf(in, "%lf", &resonance_spectra[iqt][iqx][iqy][iqz][0][ipt][ipphi]
                            &resonance_spectra[iqt][iqx][iqy][iqz][1][ipt][ipphi]);

	return(fclose(in));
}

void CorrelationFunction::Readin_daughter_spectra(int resonance_pid, int folderindex, double ******* resonance_spectra)
{
	ostringstream filename_stream;
	filename_stream << global_path << "/resonance_" << resonance_pid << "_spectra_ev" << folderindex << no_df_stem << ".dat";
	FILE *in = fopen(filename_stream.str().c_str(), "r");

	for (int iqt = 0; iqt < qtnpts; ++iqt)
	for (int iqx = 0; iqx < qxnpts; ++iqx)
	for (int iqy = 0; iqy < qynpts; ++iqy)
	for (int iqz = 0; iqz < qznpts; ++iqz)
	for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
	for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
		fscanf(in, "%lf", &resonance_spectra[iqt][iqx][iqy][iqz][0][ipt][ipphi]
                            &resonance_spectra[iqt][iqx][iqy][iqz][1][ipt][ipphi]);

	return(fclose(in));
}
*/

void CorrelationFunction::Output_correlationfunction(int folderindex)
{
	ostringstream oCorrFunc_stream;
	string temp_particle_name = particle_name;
	replace_parentheses(temp_particle_name);
	oCorrFunc_stream << global_path << "/correlfunct3D" << "_" << temp_particle_name << ".dat";
	ofstream oCorrFunc;
	oCorrFunc.open(oCorrFunc_stream.str().c_str());

	for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
	for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
	for (int iqx = 0; iqx < qxnpts; ++iqx)
	for (int iqy = 0; iqy < qynpts; ++iqy)
	for (int iqz = 0; iqz < qznpts; ++iqz)
		oCorrFunc << scientific << setprecision(7) << setw(15)
			<< SPinterp_pT[ipt] << "   " << SPinterp_pphi[ipphi] << "   " << qx_PTdep_pts[ipt][iqx] << "   "
			<< qy_PTdep_pts[ipt][iqy] << "   " << qz_PTdep_pts[ipt][iqz] << "   " << CFvals[ipt][ipphi][iqx][iqy][iqz] << endl;

	oCorrFunc.close();
				
	return;
}

void CorrelationFunction::Readin_results(int folderindex)
{
	double dummy;
	ostringstream filename_stream_HBT;
	filename_stream_HBT << global_path << "/HBTradii_GF_ev" << folderindex << no_df_stem << ".dat";
	ifstream inputHBT(filename_stream_HBT.str().c_str());

	for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
	for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
	{
		inputHBT >> dummy;	//pt value
		inputHBT >> dummy;	//pphi value
		inputHBT >> R2_side[ipt][ipphi];
		inputHBT >> R2_out[ipt][ipphi];
		inputHBT >> R2_outside[ipt][ipphi];
		inputHBT >> R2_long[ipt][ipphi];
		inputHBT >> R2_sidelong[ipt][ipphi];
		inputHBT >> R2_outlong[ipt][ipphi];
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

	for(int ipphi = 0; ipphi < n_interp_pphi_pts; ipphi++)
	{
		for(int ipt = 0; ipt < n_interp_pT_pts; ipt++)
			output_target_dN_dypTdpTdphi << scientific << setprecision(8) << setw(12) << spectra[target_particle_id][ipt][ipphi] << "   ";
		output_target_dN_dypTdpTdphi << endl;
	}

	output_target_dN_dypTdpTdphi.close();

	return;
}

/*void CorrelationFunction::Output_current_correlator(int folderindex)
{
	string local_name = all_particles[target_particle_id].name;
	replace_parentheses(local_name);
	ostringstream filename_stream_correlator;
	filename_stream_correlator << global_path << "/correlator_" << local_name << "_vs_resonance_ev" << folderindex << no_df_stem << ".dat";
	ofstream output_correlator(filename_stream_correlator.str().c_str(), ios::app);

	int HDFloadTargetSuccess = Get_resonance_from_HDF_array(target_particle_id, current_dN_dypTdpTdphi_moments);

	for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
	for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
	{
		current_dN_dypTdpTdphi_moments[ipt][ipphi][iqt0][iqx0][iqy0][iqz0][1] = 0.0;
		double nonFTd_spectra = spectra[target_particle_id][ipt][ipphi];

		for (int iqt = 0; iqt < qtnpts; ++iqt)
		for (int iqx = 0; iqx < qxnpts; ++iqx)
		for (int iqy = 0; iqy < qynpts; ++iqy)
		for (int iqz = 0; iqz < qznpts; ++iqz)
		{
			// output all FT'd spectra
			double cos_transf_spectra = current_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx][iqy][iqz][0];
			double sin_transf_spectra = current_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx][iqy][iqz][1];
			output_correlator << scientific << setprecision(8) << setw(12)
				<< current_total_resonance_percentage << "   " << SPinterp_pT[ipt] << "   " << SPinterp_pphi[ipphi]
				<< "   " << qt_pts[iqt] << "   " << qx_pts[iqx] << "   " << qy_pts[iqy] << "   " << qz_pts[iqz]
				<< "   " << 1. + (cos_transf_spectra*cos_transf_spectra + sin_transf_spectra*sin_transf_spectra)/(nonFTd_spectra*nonFTd_spectra) << endl;
		}
	}

	output_correlator.close();

	return;
}

void CorrelationFunction::Readin_current_correlator(int folderindex)
{
	string local_name = all_particles[target_particle_id].name;
	replace_parentheses(local_name);
	ostringstream filename_stream_correlator;
	filename_stream_correlator << global_path << "/correlator_" << local_name << "_vs_resonance_ev" << folderindex << no_df_stem << ".dat";
	ifstream input_correlator(filename_stream_correlator.str().c_str());

	double dummy;


	for (int ir = 0; ir < (int); ++ir)
	for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
	for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
	for (int iqt = 0; iqt < qtnpts; ++iqt)
	for (int iqx = 0; iqx < qxnpts; ++iqx)
	for (int iqy = 0; iqy < qynpts; ++iqy)
	for (int iqz = 0; iqz < qznpts; ++iqz)
	{
		// output all FT'd spectra
		input_correlator >> current_total_resonance_percentage;
		input_correlator >> dummy;
		input_correlator >> dummy;
		input_correlator >> dummy;
		input_correlator >> dummy;
		input_correlator >> dummy;
		input_correlator >> dummy;
		input_correlator >> something;
	}

	input_correlator.close();

	return;
}*/


void CorrelationFunction::Output_total_target_eiqx_dN_dypTdpTdphi(int folderindex)
{
	string local_name = all_particles[target_particle_id].name;
	replace_parentheses(local_name);
	ostringstream filename_stream_target_dN_dypTdpTdphi;
	filename_stream_target_dN_dypTdpTdphi << global_path << "/total_" << local_name << "_eiqx_dN_dypTdpTdphi_ev" << folderindex << no_df_stem << ".dat";
	ofstream output_target_dN_dypTdpTdphi(filename_stream_target_dN_dypTdpTdphi.str().c_str());

	int HDFloadTargetSuccess = Get_resonance_from_HDF_array(target_particle_id, current_dN_dypTdpTdphi_moments);

	double local_epsilon = 1.e-4;

	for (int iqt = 0; iqt < qtnpts; ++iqt)
	for (int iqx = 0; iqx < qxnpts; ++iqx)
	for (int iqy = 0; iqy < qynpts; ++iqy)
	for (int iqz = 0; iqz < qznpts; ++iqz)
	for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
	for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
	{
		// addresses NaN issue in sin component when all q^{\mu} == 0
		//if (sqrt(qt_pts[iqt]*qt_pts[iqt]+qx_pts[iqx]*qx_pts[iqx]+qy_pts[iqy]*qy_pts[iqy]+qz_pts[iqz]*qz_pts[iqz]) < 1.e-12)
		//	current_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx][iqy][iqz][1] = 0.0;
		current_dN_dypTdpTdphi_moments[ipt][ipphi][iqt0][iqx0][iqy0][iqz0][1] = 0.0;

		// output all FT'd spectra
		// with all resonance included
		double nonFTd_spectra = spectra[target_particle_id][ipt][ipphi];
		double cos_transf_spectra = current_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx][iqy][iqz][0];
		double sin_transf_spectra = current_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx][iqy][iqz][1];

		//with no resonances
		double nonFTd_tspectra = thermal_target_dN_dypTdpTdphi_moments[ipt][ipphi][iqt0][iqx0][iqy0][iqz0][0];
		double cos_transf_tspectra = thermal_target_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx][iqy][iqz][0];
		double sin_transf_tspectra = thermal_target_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx][iqy][iqz][1];

		//linearly interpolate between them to get estimate for full resonance calculation
		double num = cos_transf_spectra*cos_transf_spectra + sin_transf_spectra*sin_transf_spectra;
		double den = nonFTd_spectra*nonFTd_spectra;
		double numt = cos_transf_tspectra*cos_transf_tspectra + sin_transf_tspectra*sin_transf_tspectra;
		double dent = nonFTd_tspectra*nonFTd_tspectra;
		//double CF = (num < local_epsilon && den < local_epsilon) ? 1. : 1. + num / den;
		//double thermal_target_CF = (numt < local_epsilon && dent < local_epsilon) ? 1. : 1. + numt / dent;
		double CF = 1. + num / den;
		double thermal_target_CF = 1. + numt / dent;
		double shift = (CF - thermal_target_CF) / fraction_of_resonances;
		if (fraction_of_resonances < 1.e-12)
			shift = 0.0;

		output_target_dN_dypTdpTdphi << scientific << setprecision(8) << setw(12)
			<< qt_PTdep_pts[ipt][iqt] << "   " << qx_PTdep_pts[ipt][iqx] << "   " << qy_PTdep_pts[ipt][iqy] << "   " << qz_PTdep_pts[ipt][iqz] << "   "
			<< SPinterp_pT[ipt] << "   " << SPinterp_pphi[ipphi] << "   " << nonFTd_spectra << "   " << cos_transf_spectra << "   " << sin_transf_spectra << "   "
			<< CF << "   " << thermal_target_CF + shift << endl;
	}

	output_target_dN_dypTdpTdphi.close();

	return;
}

void CorrelationFunction::Readin_total_target_eiqx_dN_dypTdpTdphi(int folderindex)
{
	string local_name = all_particles[target_particle_id].name;
	replace_parentheses(local_name);
	ostringstream filename_stream_target_dN_dypTdpTdphi;
	filename_stream_target_dN_dypTdpTdphi << global_path << "/total_" << local_name << "_eiqx_dN_dypTdpTdphi_ev" << folderindex << no_df_stem << ".dat";
	ifstream input_target_dN_dypTdpTdphi(filename_stream_target_dN_dypTdpTdphi.str().c_str());

	double dummy = 0.0;

	//Need some way to set these guys!  probably write q points out to file and then read them back in as needed...
	//THIS FUNCTION NOT STABLE UNTIL THIS IS FIXED
	//for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
	//	Set_q_pTdep_pts(ipt);

	for (int iqt = 0; iqt < qtnpts; ++iqt)
	for (int iqx = 0; iqx < qxnpts; ++iqx)
	for (int iqy = 0; iqy < qynpts; ++iqy)
	for (int iqz = 0; iqz < qznpts; ++iqz)
	for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
	for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
	{
		//input_target_dN_dypTdpTdphi >> qt_PTdep_pts[ipt][iqt];
		//input_target_dN_dypTdpTdphi >> qx_PTdep_pts[ipt][iqx];
		//input_target_dN_dypTdpTdphi >> qy_PTdep_pts[ipt][iqy];
		//input_target_dN_dypTdpTdphi >> qz_PTdep_pts[ipt][iqz];
		//input_target_dN_dypTdpTdphi >> SPinterp_pT[ipt];
		//input_target_dN_dypTdpTdphi >> SPinterp_pphi[ipphi];
		input_target_dN_dypTdpTdphi >> dummy;
		input_target_dN_dypTdpTdphi >> dummy;
		input_target_dN_dypTdpTdphi >> dummy;
		input_target_dN_dypTdpTdphi >> dummy;
		input_target_dN_dypTdpTdphi >> dummy;
		input_target_dN_dypTdpTdphi >> dummy;
		input_target_dN_dypTdpTdphi >> spectra[target_particle_id][ipt][ipphi];
		input_target_dN_dypTdpTdphi >> current_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx][iqy][iqz][0];
		input_target_dN_dypTdpTdphi >> current_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx][iqy][iqz][1];
		input_target_dN_dypTdpTdphi >> dummy;
		/*cout << scientific << setprecision(8) << setw(12)
			<< qt_PTdep_pts[ipt][iqt] << "   " << qx_PTdep_pts[ipt][iqx] << "   " << qy_PTdep_pts[ipt][iqy] << "   " << qz_PTdep_pts[ipt][iqz] << "   "
			<< SPinterp_pT[ipt] << "   " << SPinterp_pphi[ipphi] << "   "
			<< spectra[target_particle_id][ipt][ipphi] << "   "
			<< current_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx][iqy][iqz][0] << "   "
			<< current_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx][iqy][iqz][1] << "   "
			<< dummy << endl;*/
	}

	input_target_dN_dypTdpTdphi.close();

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
