#include<iostream>
#include<sstream>
#include<string>
#include<fstream>
#include<cmath>
#include<iomanip>
#include<vector>
#include<stdio.h>

#include "CFWR.h"
#include "CFWR_lib.h"
#include "Arsenal.h"
#include "gauss_quadrature.h"
#include "stats.h"

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

//allows possibility of dumping thermal_spectra, spectra, log_spectra, etc...
void CorrelationFunction::Dump_spectra_array(string output_filename, double *** array_to_dump)
{
	ostringstream filename_stream;
	filename_stream << global_path << "/" << output_filename;
	ofstream out(filename_stream.str().c_str());

	for (int ir = 0; ir < Nparticle; ++ir)
	for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
	{
		for (int ipT = 0; ipT < n_interp_pT_pts; ++ipT)
			out << scientific << setprecision(8) << setw(12) << array_to_dump[ir][ipT][ipphi] << "   ";
		out << endl;
	}

	out.close();
}

//allows possibility of reading in thermal_spectra, spectra, log_spectra, etc...
void CorrelationFunction::Load_spectra_array(string input_filename, double *** array_to_read)
{
	ostringstream filename_stream;
	filename_stream << global_path << "/" << input_filename;
	ifstream in(filename_stream.str().c_str());

	for (int ir = 0; ir < Nparticle; ++ir)
	for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
	for (int ipT = 0; ipT < n_interp_pT_pts; ++ipT)
		in >> array_to_read[ir][ipT][ipphi];

	in.close();
}

void CorrelationFunction::Output_results(int folderindex)
{
//debugger(__LINE__, __FILE__);
	ostringstream filename_stream_HBT;
	filename_stream_HBT << global_path << "/HBTradii_GF_ev" << folderindex << no_df_stem << ".dat";
	ofstream outputHBT;
	outputHBT.open(filename_stream_HBT.str().c_str());
	ostringstream filename_stream_HBTcfs;
	filename_stream_HBTcfs << global_path << "/HBTradii_GF_cfs_ev" << folderindex << no_df_stem << ".dat";
	ofstream outputHBTcfs;
	outputHBTcfs.open(filename_stream_HBTcfs.str().c_str());

//debugger(__LINE__, __FILE__);
	for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
	{
//debugger(__LINE__, __FILE__);
		for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
		{
//debugger(__LINE__, __FILE__);
			outputHBT << SPinterp_pT[ipt] << "   " << SPinterp_pphi[ipphi]
				<< "   " << R2_side[ipt][ipphi] << "   " << R2_out[ipt][ipphi]
				<< "   " << R2_outside[ipt][ipphi] << "   " << R2_long[ipt][ipphi]
				<< "   " << R2_sidelong[ipt][ipphi] << "   " << R2_outlong[ipt][ipphi] << endl;
			//ostringstream local_cmd;
			//local_cmd << "wc -l " << global_path << "/HBTradii_GF_ev" << folderindex << no_df_stem << ".dat >> plumberg.out";
			//int cmd_result = system(local_cmd.str().c_str());
//cout << filename_stream_HBT.str() << ": " << SPinterp_pT[ipt] << "   " << SPinterp_pphi[ipphi]
//				<< "   " << R2_side[ipt][ipphi] << "   " << R2_out[ipt][ipphi]
//				<< "   " << R2_outside[ipt][ipphi] << "   " << R2_long[ipt][ipphi]
//				<< "   " << R2_sidelong[ipt][ipphi] << "   " << R2_outlong[ipt][ipphi] << endl;
		}
//debugger(__LINE__, __FILE__);

		//do Fourier transforming here for now...
		double plane_psi = 0.0;
		R2_Fourier_transform(ipt, plane_psi);
		for (int Morder = 0; Morder < n_order; Morder++)
		{
			outputHBTcfs << folderindex << "  " << SPinterp_pT[ipt] << "  " << Morder
				<< "  " << R2_side_C[ipt][Morder] << "   " << R2_side_S[ipt][Morder] << "  " << R2_out_C[ipt][Morder] << "  " << R2_out_S[ipt][Morder]
				<< "  " << R2_outside_C[ipt][Morder] << "   " << R2_outside_S[ipt][Morder] << "  " << R2_long_C[ipt][Morder] << "  " << R2_long_S[ipt][Morder]
				<< "  " << R2_sidelong_C[ipt][Morder] << "   " << R2_sidelong_S[ipt][Morder] << "  " << R2_outlong_C[ipt][Morder] << "  " << R2_outlong_S[ipt][Morder] << endl;
		}
//debugger(__LINE__, __FILE__);
	}

//debugger(__LINE__, __FILE__);
	outputHBT.close();
	outputHBTcfs.close();

	return;
}

void CorrelationFunction::Output_correlationfunction(bool regulated_CF /*==true*/)
{
	ostringstream oCorrFunc_stream;
	string temp_particle_name = particle_name;
	replace_parentheses(temp_particle_name);

	string CF_reg_string = "";
	if (!regulated_CF)
		CF_reg_string = "unregulated_";

	oCorrFunc_stream << global_path << "/correlfunct3D_" << CF_reg_string << temp_particle_name << ".dat";
	ofstream oCorrFunc;
	oCorrFunc.open(oCorrFunc_stream.str().c_str());

	for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
	for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
	for (int iqx = 0; iqx < qxnpts; ++iqx)
	for (int iqy = 0; iqy < qynpts; ++iqy)
	for (int iqz = 0; iqz < qznpts; ++iqz)
	{
		double ckp = cos_SPinterp_pphi[ipphi], skp = sin_SPinterp_pphi[ipphi];
		oCorrFunc << scientific << setprecision(8) << setw(12)
			<< SPinterp_pT[ipt] << "   " << SPinterp_pphi[ipphi] << "   " << qx_pts[iqx] << "   "
			<< qy_pts[iqy] << "   " << qz_pts[iqz] << "   "
			<< qx_pts[iqx] * ckp + qy_pts[iqy] * skp << "   "
			<< -qx_pts[iqx] * skp + qy_pts[iqy] * ckp << "   "
			<< qz_pts[iqz] << "   "
			<< thermalCFvals[ipt][ipphi][iqx][iqy][iqz] << "   " << crosstermCFvals[ipt][ipphi][iqx][iqy][iqz] << "   " << resonancesCFvals[ipt][ipphi][iqx][iqy][iqz] << "   " << CFvals[ipt][ipphi][iqx][iqy][iqz] << endl;
	}

	oCorrFunc.close();
				
	return;
}

void CorrelationFunction::Output_fleshed_out_correlationfunction(int ipt, int ipphi)
{
	ostringstream oCorrFunc_stream;
	string temp_particle_name = particle_name;
	replace_parentheses(temp_particle_name);
	oCorrFunc_stream << global_path << "/correlfunct3D" << "_" << temp_particle_name << "_fleshed_out.dat";
	ofstream oCorrFunc;
	if (ipt==0 && ipphi==0)
		oCorrFunc.open(oCorrFunc_stream.str().c_str());
	else
		oCorrFunc.open(oCorrFunc_stream.str().c_str(), ios::app);

	for (int iqx = 0; iqx < new_nqpts; ++iqx)
	for (int iqy = 0; iqy < new_nqpts; ++iqy)
	for (int iqz = 0; iqz < new_nqpts; ++iqz)
	{
		double ckp = cos_SPinterp_pphi[ipphi], skp = sin_SPinterp_pphi[ipphi];
		oCorrFunc << scientific << setprecision(7) << setw(15)
			<< SPinterp_pT[ipt] << "   " << SPinterp_pphi[ipphi] << "   " << qx_fleshed_out_pts[iqx] << "   "
			<< qy_fleshed_out_pts[iqy] << "   " << qz_fleshed_out_pts[iqz] << "   "
			<< qx_fleshed_out_pts[iqx] * ckp + qy_fleshed_out_pts[iqy] * skp << "   "
			<< -qx_fleshed_out_pts[iqx] * skp + qy_fleshed_out_pts[iqy] * ckp << "   "
			<< qz_fleshed_out_pts[iqz] << "   "
			<< fleshed_out_thermal[iqx][iqy][iqz] << "   " << fleshed_out_crossterm[iqx][iqy][iqz] << "   " << fleshed_out_resonances[iqx][iqy][iqz] << "   " << fleshed_out_CF[iqx][iqy][iqz] << endl;
	}

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

void CorrelationFunction::Output_total_target_eiqx_dN_dypTdpTdphi(int folderindex, double current_fraction /*==-1.0*/)
{
	string local_name = all_particles[target_particle_id].name;
	string current_fraction_string = (current_fraction >= 0.0) ? "_" + patch::to_string(current_fraction) : "";
	replace_parentheses(local_name);
	ostringstream filename_stream_target_dN_dypTdpTdphi;
	filename_stream_target_dN_dypTdpTdphi << global_path << "/total_" << local_name << current_fraction_string << "_eiqx_dN_dypTdpTdphi_ev" << folderindex << no_df_stem << ".dat";
	ofstream output_target_dN_dypTdpTdphi(filename_stream_target_dN_dypTdpTdphi.str().c_str());

	int HDFloadTargetSuccess = Get_resonance_from_HDF_array(target_particle_id, current_dN_dypTdpTdphi_moments);

	// addresses NaN issue in sin component when all q^{\mu} == 0
	if (qtnpts%2==1 && qxnpts%2==1 && qynpts%2==1 && qznpts%2==1)
	{	//if all q-ranges are odd and centered on q=0 ==> q=0 is included!
		int iqt0 = (qtnpts-1)/2;
		int iqx0 = (qxnpts-1)/2;
		int iqy0 = (qynpts-1)/2;
		int iqz0 = (qznpts-1)/2;
		for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
		for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
		{
			current_dN_dypTdpTdphi_moments[ipt][ipphi][iqt0][iqx0][iqy0][iqz0][1] = 0.0;
			thermal_target_dN_dypTdpTdphi_moments[ipt][ipphi][iqt0][iqx0][iqy0][iqz0][1] = 0.0;
		}
	}

	//use this to help look for outliers when CF calculation gets really choppy and noisy
	//Set_target_pphiavgd_CFs();

	for (int iqt = 0; iqt < qtnpts; ++iqt)
	for (int iqx = 0; iqx < qxnpts; ++iqx)
	for (int iqy = 0; iqy < qynpts; ++iqy)
	for (int iqz = 0; iqz < qznpts; ++iqz)
	for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
	for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
	{
		//first, get CF and projected CF
		double CF = get_CF(ipt, ipphi, iqt, iqx, iqy, iqz, false);				//false means don't return projected value
		//double projected_CF = get_CF(ipt, ipphi, iqt, iqx, iqy, iqz, true && !thermal_pions_only);	//true means do return projected value

		//now, regulate results
		//Regulate_CF(ipt, iqt, iqx, iqy, iqz, &CF, &projected_CF);

		//!!!!!!!!!!!!should get projected_CF AFTER regulating CF...!!!!!!!!!!!!
		double projected_CF = get_CF(ipt, ipphi, iqt, iqx, iqy, iqz, true && !thermal_pions_only);	//true means do return projected value

		double nonFTd_spectra = spectra[target_particle_id][ipt][ipphi];
		double cos_transf_spectra = current_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx][iqy][iqz][0];
		double sin_transf_spectra = current_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx][iqy][iqz][1];

		output_target_dN_dypTdpTdphi << scientific << setprecision(8) << setw(12)
			<< qt_pts[iqt] << "   " << qx_pts[iqx] << "   " << qy_pts[iqy] << "   " << qz_pts[iqz] << "   "
			<< SPinterp_pT[ipt] << "   " << SPinterp_pphi[ipphi] << "   "
			<< nonFTd_spectra << "   "																								//non-thermal + thermal
			<< cos_transf_spectra << "   "																							//non-thermal + thermal (cos)
			<< sin_transf_spectra << "   "																							//non-thermal + thermal (sin)
			<< thermal_spectra[target_particle_id][ipt][ipphi] << "   "																//thermal only
			<< thermal_target_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx][iqy][iqz][0] << "   "									//thermal only (cos)
			<< thermal_target_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx][iqy][iqz][1] << "   "									//thermal only (sin)
			<< nonFTd_spectra - thermal_spectra[target_particle_id][ipt][ipphi] << "   "											//non-thermal only
			<< cos_transf_spectra - thermal_target_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx][iqy][iqz][0] << "   "				//non-thermal only (cos)
			<< sin_transf_spectra - thermal_target_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx][iqy][iqz][1] << "   "				//non-thermal only (sin)
			<< CF << "   " << projected_CF << endl;
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
		//input_target_dN_dypTdpTdphi >> qt_pts[iqt];
		//input_target_dN_dypTdpTdphi >> qx_pts[iqx];
		//input_target_dN_dypTdpTdphi >> qy_pts[iqy];
		//input_target_dN_dypTdpTdphi >> qz_pts[iqz];
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
		input_target_dN_dypTdpTdphi >> thermal_spectra[target_particle_id][ipt][ipphi];
		input_target_dN_dypTdpTdphi >> thermal_target_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx][iqy][iqz][0];
		input_target_dN_dypTdpTdphi >> thermal_target_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx][iqy][iqz][1];
		input_target_dN_dypTdpTdphi >> dummy;
		input_target_dN_dypTdpTdphi >> dummy;
		input_target_dN_dypTdpTdphi >> dummy;
		input_target_dN_dypTdpTdphi >> dummy;
		input_target_dN_dypTdpTdphi >> dummy;
		/*cout << scientific << setprecision(8) << setw(12)
			<< qt_pts[iqt] << "   " << qx_pts[iqx] << "   " << qy_pts[iqy] << "   " << qz_pts[iqz] << "   "
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

void CorrelationFunction::Output_resonance_fraction()
{
	ostringstream filename_stream_rf;
	filename_stream_rf << global_path << "/resonance_fraction.dat";
	ofstream output_rf(filename_stream_rf.str().c_str());

	output_rf << fraction_of_resonances << endl;

	output_rf.close();

	return;
}

//End of file
