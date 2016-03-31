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

void CorrelationFunction::Dump_q_pTdep_pts()
{
	ostringstream filename_stream_t;
	filename_stream_t << global_path << "/qt_pTdep_pts.dat";
	ofstream outt(filename_stream_t.str().c_str());
	ostringstream filename_stream_x;
	filename_stream_x << global_path << "/qx_pTdep_pts.dat";
	ofstream outx(filename_stream_x.str().c_str());
	ostringstream filename_stream_y;
	filename_stream_y << global_path << "/qy_pTdep_pts.dat";
	ofstream outy(filename_stream_y.str().c_str());
	ostringstream filename_stream_z;
	filename_stream_z << global_path << "/qz_pTdep_pts.dat";
	ofstream outz(filename_stream_z.str().c_str());

	for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
	{
		for (int iqt = 0; iqt < qtnpts; ++iqt)
			outt << qt_PTdep_pts[ipt][iqt] << "   ";
		outt << endl;
		for (int iqx = 0; iqx < qxnpts; ++iqx)
			outx << qx_PTdep_pts[ipt][iqx] << "   ";
		outx << endl;
		for (int iqy = 0; iqy < qynpts; ++iqy)
			outy << qy_PTdep_pts[ipt][iqy] << "   ";
		outy << endl;
		for (int iqz = 0; iqz < qznpts; ++iqz)
			outz << qz_PTdep_pts[ipt][iqz] << "   ";
		outz << endl;
	}

	outt.close();
	outx.close();
	outy.close();
	outz.close();

	return;
}

void CorrelationFunction::Load_q_pTdep_pts()
{
	ostringstream filename_stream_t;
	filename_stream_t << global_path << "/qt_pTdep_pts.dat";
	ifstream in_t(filename_stream_t.str().c_str());
	ostringstream filename_stream_x;
	filename_stream_x << global_path << "/qx_pTdep_pts.dat";
	ifstream in_x(filename_stream_x.str().c_str());
	ostringstream filename_stream_y;
	filename_stream_y << global_path << "/qy_pTdep_pts.dat";
	ifstream in_y(filename_stream_y.str().c_str());
	ostringstream filename_stream_z;
	filename_stream_z << global_path << "/qz_pTdep_pts.dat";
	ifstream in_z(filename_stream_z.str().c_str());

	for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
	{
		qt_PTdep_pts[ipt] = new double [qtnpts];
		qx_PTdep_pts[ipt] = new double [qxnpts];
		qy_PTdep_pts[ipt] = new double [qynpts];
		qz_PTdep_pts[ipt] = new double [qznpts];

		for (int iqt = 0; iqt < qtnpts; ++iqt)
			in_t >> qt_PTdep_pts[ipt][iqt];
		for (int iqx = 0; iqx < qxnpts; ++iqx)
			in_x >> qx_PTdep_pts[ipt][iqx];
		for (int iqy = 0; iqy < qynpts; ++iqy)
			in_y >> qy_PTdep_pts[ipt][iqy];
		for (int iqz = 0; iqz < qznpts; ++iqz)
			in_z >> qz_PTdep_pts[ipt][iqz];
	}

	in_t.close();
	in_x.close();
	in_y.close();
	in_z.close();

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
	{
		double ckp = cos_SPinterp_pphi[ipphi], skp = sin_SPinterp_pphi[ipphi];
		oCorrFunc << scientific << setprecision(7) << setw(15)
			<< SPinterp_pT[ipt] << "   " << SPinterp_pphi[ipphi] << "   " << qx_PTdep_pts[ipt][iqx] << "   "
			<< qy_PTdep_pts[ipt][iqy] << "   " << qz_PTdep_pts[ipt][iqz] << "   "
			<< qx_PTdep_pts[ipt][iqx] * ckp + qy_PTdep_pts[ipt][iqy] * skp << "   "
			<< -qx_PTdep_pts[ipt][iqx] * skp + qy_PTdep_pts[ipt][iqy] * ckp << "   "
			<< qz_PTdep_pts[ipt][iqz] << "   "
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


/*void CorrelationFunction::Readin_correlationfunction(int folderindex)
{
	ostringstream iCorrFunc_stream;
	string temp_particle_name = particle_name;
	replace_parentheses(temp_particle_name);
	iCorrFunc_stream << global_path << "/correlfunct3D" << "_" << temp_particle_name << ".dat";
	ifstream iCorrFunc;
	iCorrFunc.open(iCorrFunc_stream.str().c_str());

	double dummy = 0.0;

	for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
	for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
	for (int iqx = 0; iqx < qxnpts; ++iqx)
	for (int iqy = 0; iqy < qynpts; ++iqy)
	for (int iqz = 0; iqz < qznpts; ++iqz)
	{
		iCorrFunc 
			>> dummy
			>> dummy
			>> dummy
			>> dummy
			>> dummy
			>> dummy
			>> dummy
			>> dummy
			>> CFvals[ipt][ipphi][iqx][iqy][iqz];
	}

	iCorrFunc.close();
				
	return;
}*/

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

void CorrelationFunction::Regulate_CF(int ipt, int iqt, int iqx, int iqy, int iqz, double * CF, double * projCF)
{
	double mu = target_pphiavgd_CFs[ipt][iqt][iqx][iqy][iqz];
	double sigma = sqrt(target_pphivar_CFs[ipt][iqt][iqx][iqy][iqz]);

	double zthreshhold = 3.0;
	double pTcutoff = 2.0;	//GeV

	double zact = abs(mu - *CF) / sigma;

	if (zact >= zthreshhold || abs(*CF-1.5) > 0.500001)
	{
		//if (SPinterp_pT[ipt] < pTcutoff)
			*global_out_stream_ptr << "WARNING: regulated CF point at pT = " << SPinterp_pT[ipt]
				<< ": (" << *CF << "," << *projCF << ") --> (";
		*CF = mu;	//if it's an outlier, replace with pphi-averaged value
		*projCF = mu;
		//if (SPinterp_pT[ipt] < pTcutoff)
			*global_out_stream_ptr << *CF << "," << *projCF << ")" << endl;
	}

	return;
}

void CorrelationFunction::Regulate_CF_Hampel(int ipt, int iqt, int iqx, int iqy, int iqz,
												double * pphi_CF_slice, double * pphi_CF_slice_term1, double * pphi_CF_slice_term2, double * pphi_CF_slice_term3)
{
	bool * is_outlier = new bool [n_interp_pphi_pts];

	double pphi_median = 0.0;

	find_outliers_Hampel(pphi_CF_slice, n_interp_pphi_pts, is_outlier, &pphi_median);

	for (ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
	if (is_outlier[ipphi] || abs(*pphi_CF_slice[ipphi]-1.5) > 0.500001)
	{
		//if (SPinterp_pT[ipt] < pTcutoff)
			*global_out_stream_ptr << "\t WARNING: regulated CF point at pT = " << SPinterp_pT[ipt]
				<< ": (" << pphi_CF_slice[ipphi] << ", " << pphi_CF_slice_term1[ipphi] << ", " << pphi_CF_slice_term2[ipphi] << ", " << pphi_CF_slice_term3[ipphi] << ") --> (";
		pphi_CF_slice[ipphi] = pphi_median;	//if it's an outlier, replace with pphi-median value
		pphi_CF_slice_term1[ipphi] = get_median(pphi_CF_slice_term1, n_interp_pphi_pts);
		pphi_CF_slice_term2[ipphi] = get_median(pphi_CF_slice_term2, n_interp_pphi_pts);
		pphi_CF_slice_term3[ipphi] = get_median(pphi_CF_slice_term3, n_interp_pphi_pts);
		//if (SPinterp_pT[ipt] < pTcutoff)
			*global_out_stream_ptr << pphi_CF_slice[ipphi] << ", " << pphi_CF_slice_term1[ipphi] << ", " << pphi_CF_slice_term2[ipphi] << ", " << pphi_CF_slice_term3[ipphi] << ")" << endl;
	}

	delete [] is_outlier;

	return;
}


double CorrelationFunction::get_CF(int ipt, int ipphi, int iqt, int iqx, int iqy, int iqz, bool return_projected_value)
{
	double nonFTd_spectra = spectra[target_particle_id][ipt][ipphi];
	double cos_transf_spectra = current_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx][iqy][iqz][0];
	double sin_transf_spectra = current_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx][iqy][iqz][1];

	if (return_projected_value)
	{
		//with no resonances
		double nonFTd_tspectra = thermal_spectra[target_particle_id][ipt][ipphi];
		double cos_transf_tspectra = thermal_target_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx][iqy][iqz][0];
		double sin_transf_tspectra = thermal_target_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx][iqy][iqz][1];

		double projected_nonFTd_spectra = nonFTd_tspectra + (nonFTd_spectra - nonFTd_tspectra) / fraction_of_resonances;
		double projected_cos_transf_spectra = cos_transf_tspectra + (cos_transf_spectra - cos_transf_tspectra) / fraction_of_resonances;
		double projected_sin_transf_spectra = sin_transf_tspectra + (sin_transf_spectra - sin_transf_tspectra) / fraction_of_resonances;

		//projected full value of correlation function by projecting moments *first*
		double projected_num = projected_cos_transf_spectra*projected_cos_transf_spectra + projected_sin_transf_spectra*projected_sin_transf_spectra;
		double projected_den = projected_nonFTd_spectra*projected_nonFTd_spectra;
		return (1. + projected_num / projected_den);
	}
	else
	{
		double num = cos_transf_spectra*cos_transf_spectra + sin_transf_spectra*sin_transf_spectra;
		double den = nonFTd_spectra*nonFTd_spectra;
		return (1. + num / den);
	}
}

void CorrelationFunction::get_CF(double * totalresult, double * thermalresult, double * crosstermresult, double * resonanceresult,
									int ipt, int ipphi, int iqt, int iqx, int iqy, int iqz, bool return_projected_value)
{

	//thermal
	double nonFTd_tspectra = thermal_spectra[target_particle_id][ipt][ipphi];
	double cos_transf_tspectra = thermal_target_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx][iqy][iqz][0];
	double sin_transf_tspectra = thermal_target_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx][iqy][iqz][1];
	//total
	double nonFTd_spectra = spectra[target_particle_id][ipt][ipphi];
	double cos_transf_spectra = current_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx][iqy][iqz][0];
	double sin_transf_spectra = current_dN_dypTdpTdphi_moments[ipt][ipphi][iqt][iqx][iqy][iqz][1];

	if (return_projected_value)
	{
		nonFTd_spectra += (nonFTd_spectra - nonFTd_tspectra) / fraction_of_resonances;
		cos_transf_spectra += (cos_transf_spectra - cos_transf_tspectra) / fraction_of_resonances;
		sin_transf_spectra += (sin_transf_spectra - sin_transf_tspectra) / fraction_of_resonances;
	}

	//non-thermal
	double NT_spectra = nonFTd_spectra - nonFTd_tspectra;
	double cosNT_spectra = cos_transf_spectra - cos_transf_tspectra;
	double sinNT_spectra = sin_transf_spectra - sin_transf_tspectra;

	double num = cos_transf_tspectra*cos_transf_tspectra + sin_transf_tspectra*sin_transf_tspectra;
	double den = nonFTd_spectra*nonFTd_spectra;
	*thermalresult = num / den;

	num = 2.0 * (cosNT_spectra*cos_transf_tspectra+sinNT_spectra*sin_transf_tspectra);			//cross term
	//den = nonFTd_tspectra*nonFTd_tspectra;
	*crosstermresult = num / den;

	num = cosNT_spectra*cosNT_spectra+sinNT_spectra*sinNT_spectra;								//resonances only
	//den = nonFTd_tspectra*nonFTd_tspectra;
	*resonanceresult = num / den;

	num = cos_transf_spectra*cos_transf_spectra + sin_transf_spectra*sin_transf_spectra;
	//den = nonFTd_spectra*nonFTd_spectra;
	*totalresult = num / den;

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

		//output_target_dN_dypTdpTdphi << scientific << setprecision(8) << setw(12)
		//	<< qt_PTdep_pts[ipt][iqt] << "   " << qx_PTdep_pts[ipt][iqx] << "   " << qy_PTdep_pts[ipt][iqy] << "   " << qz_PTdep_pts[ipt][iqz] << "   "
		//	<< SPinterp_pT[ipt] << "   " << SPinterp_pphi[ipphi] << "   " << nonFTd_spectra << "   " << cos_transf_spectra << "   " << sin_transf_spectra << "   "
		//	<< CF << "   " << projected_CF << endl;
		output_target_dN_dypTdpTdphi << scientific << setprecision(8) << setw(12)
			<< qt_PTdep_pts[ipt][iqt] << "   " << qx_PTdep_pts[ipt][iqx] << "   " << qy_PTdep_pts[ipt][iqy] << "   " << qz_PTdep_pts[ipt][iqz] << "   "
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
