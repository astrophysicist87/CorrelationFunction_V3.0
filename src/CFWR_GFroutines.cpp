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

double slopes[4][2];

void CorrelationFunction::Get_GF_HBTradii(FO_surf* FOsurf_ptr, int folderindex)
{
	bool read_in_correlation_function = false;

	*global_out_stream_ptr << "--> Getting HBT radii by Gaussian fit method" << endl;
	for(int iKT = 0; iKT < n_localp_T; iKT++)
	{
      	*global_out_stream_ptr << "   - Calculating K_T = " << K_T[iKT] << " GeV ..." << endl;
		for(int iKphi = 0; iKphi < n_localp_phi; iKphi++)
		{
			*global_out_stream_ptr << "\t\t --> Calculating K_phi = " << K_phi[iKphi] << " ..." << endl;
			if (corrfuncdim == 1)
			{
				if (read_in_correlation_function)
					Read_correlationfunction(iKT, iKphi);
				//else
				//	Cal_correlationfunction();
				//Fit_Correlationfunction1D('o', iKT, iKphi);
				//Fit_Correlationfunction1D('s', iKT, iKphi);
				//Fit_Correlationfunction1D('l', iKT, iKphi);
				//if (!read_in_correlation_function)
				//	Output_correlationfunction(folderindex);
	 		}
			else if (corrfuncdim == 3)
			{
				//Cal_correlationfunction_3D(iKT, iKphi);
				//if (lambdaflag)
				//	Fit_Correlationfunction3D_withlambda(iKT, iKphi);
				//else Fit_Correlationfunction3D(iKT, iKphi);
				//Output_Correlationfunction_3D(iKT, iKphi, folderindex);
				cerr << "corrfuncdim == 3 not currently supported!" << endl;
			}
			else
				cerr << "corrfuncdim == " << corrfuncdim << " not currently supported!" << endl;
		}
		R2_Fourier_transform(iKT, global_plane_psi);
	}
}

double CorrelationFunction::Compute_correlationfunction(int ipt, int ipphi, double * q_interp)
{
	// try using linear-logarithmic interpolation if q-point is within grid
	// otherwise, just use linear interpolation/extrapolation, or throw an exception and do return unity
	double interpolated_result = 0.0;

	double nonFTd_spectra = spectra[target_particle_id][ipt][ipphi];

	int * qidx = new int [4];
	qidx[0] = binarySearch(qt_pts, qnpts, q_interp[0]);
	qidx[1] = binarySearch(qx_pts, qnpts, q_interp[1]);
	qidx[2] = binarySearch(qy_pts, qnpts, q_interp[2]);
	qidx[3] = binarySearch(qz_pts, qnpts, q_interp[3]);
debugger(__LINE__, __FILE__);
	bool q_point_is_outside_grid = ( qidx[0] == -1 ) || ( qidx[1] == -1 ) || ( qidx[2] == -1 ) || ( qidx[3] == -1 );

	if (!q_point_is_outside_grid)
	{
		double * sgnd_q2_interp = new double [4];
		double q_min[4] = {qt_pts[qidx[0]], qx_pts[qidx[1]], qy_pts[qidx[2]], qz_pts[qidx[3]]};
		double q_max[4] = {qt_pts[qidx[0]+1], qx_pts[qidx[1]+1], qy_pts[qidx[2]+1], qz_pts[qidx[3]+1]};
		double * sgnd_q2_min = new double [4];
		double * sgnd_q2_max = new double [4];
debugger(__LINE__, __FILE__);
		double ln_C_at_q[2][2][2][2];
		double ***** current_C_slice = current_dN_dypTdpTdphi_moments[ipt][ipphi];
debugger(__LINE__, __FILE__);
		// set values at all points used for interpolation
		for (int i = 0; i < 2; ++i)
		for (int j = 0; j < 2; ++j)
		for (int k = 0; k < 2; ++k)
		for (int l = 0; l < 2; ++l)
		{
debugger(__LINE__, __FILE__);
			cerr << q_interp[0] << "   " << q_interp[1] << "   " << q_interp[2] << "   " << q_interp[3] << endl;
			cerr << qidx[0]+i << "   " << qidx[1]+j << "   " << qidx[2]+k << "   " << qidx[3]+l << endl;
			double tmp_cos = current_C_slice[qidx[0]+i][qidx[1]+j][qidx[2]+k][qidx[3]+l][0];
debugger(__LINE__, __FILE__);
			double tmp_sin = current_C_slice[qidx[0]+i][qidx[1]+j][qidx[2]+k][qidx[3]+l][1];
debugger(__LINE__, __FILE__);
			ln_C_at_q[i][j][k][l] = log( (tmp_cos*tmp_cos + tmp_sin*tmp_sin)/(nonFTd_spectra*nonFTd_spectra) );
debugger(__LINE__, __FILE__);
		}
debugger(__LINE__, __FILE__);
		// interpolating w.r.t. q^2, not q
		for (int i = 0; i < 4; ++i)
		{
			sgnd_q2_interp[i] = sgn(q_interp[i]) * q_interp[i] * q_interp[i];
			sgnd_q2_min[i] = sgn(q_min[i]) * q_min[i] * q_min[i];
			sgnd_q2_max[i] = sgn(q_max[i]) * q_max[i] * q_max[i];
		}
debugger(__LINE__, __FILE__);

		double interpolated_ln_result = interpolate_4D(sgnd_q2_min, sgnd_q2_max, sgnd_q2_interp, ln_C_at_q);
		interpolated_result = exp(interpolated_ln_result);
debugger(__LINE__, __FILE__);
		// Clean up
		delete [] sgnd_q2_interp;
		delete [] sgnd_q2_min;
		delete [] sgnd_q2_max;
	}
	else
	{
		*global_out_stream_ptr << "Warning: q_interp point was outside of computed grid!" << endl;
		interpolated_result = 1.0;
	}
debugger(__LINE__, __FILE__);
	return (interpolated_result);
}

double CorrelationFunction::interpolate_4D(double * x_min, double * x_max, double * x_interp, double (*vals) [2][2][2])
{
	double total = 0.0;

	for (int i = 0; i < 4; ++i)
	{
		slopes[i][0] = (x_interp[i] - x_min[i]) / (x_max[i] - x_min[i]);
		//slopes[i][1] = (x_max[i] - x_interp[i]) / (x_max[i] - x_min[i]);
		slopes[i][1] = 1.0 - slopes[i][0];
	}

	for (int i = 0; i < 2; ++i)
	for (int j = 0; j < 2; ++j)
	for (int k = 0; k < 2; ++k)
	for (int l = 0; l < 2; ++l)
		total += vals[i][j][k][l] * slopes[0][i] * slopes[1][j] * slopes[2][k] * slopes[3][l];
	
	return (total);
}

void CorrelationFunction::Cal_correlationfunction()
{
	// Can't interpolate if there's only one point in q-space!
	if (qnpts == 1)
		return;

	// chooses the qo, qs, ql points at which to evaluate correlation function
	Set_correlation_function_q_pts();

	int HDFloadTargetSuccess = Get_resonance_from_HDF_array(target_particle_id, current_dN_dypTdpTdphi_moments);

	double * q_interp = new double [4];

	// Then compute full correlation function
	for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
	for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
	for (int iqo = 0; iqo < qonpts; ++iqo)
	for (int iqs = 0; iqs < qsnpts; ++iqs)
	for (int iql = 0; iql < qlnpts; ++iql)
	{
		Get_q_points(qo_pts[iqo], qs_pts[iqs], ql_pts[iql], SPinterp_pT[ipt], SPinterp_pphi[ipphi], q_interp);

		CFvals[ipt][ipphi][iqo][iqs][iql] = Compute_correlationfunction(ipt, ipphi, q_interp);
	}

	return;
}

int CorrelationFunction::Read_correlationfunction(int iKT, int iKphi)
{
   if(fabs(K_y) > 1e-16)
   {
       //cout<<"not support for y not equals 0 yet!" << endl;
       return (0);
   }

   double local_K_T = K_T[iKT];
   double localK_phi = K_phi[iKphi];

ostringstream corrfn_stream;
corrfn_stream << global_path;

if (fabs(local_K_T) < 1e-16) corrfn_stream << "/correlfunct1D" << "_" << particle_name << "_kt_0.0_phi_";
else if (fabs(local_K_T - 1.) < 1e-16) corrfn_stream << "/correlfunct1D" << "_" << particle_name << "_kt_1.0_phi_";
else corrfn_stream << "/correlfunct1D" << "_" << particle_name << "_kt_" << local_K_T << "_phi_";

if (fabs(localK_phi) < 1e-16) corrfn_stream << "0.0.dat";
else corrfn_stream << localK_phi << ".dat";

if ( !fexists( corrfn_stream.str().c_str() ) )
{
	cerr << corrfn_stream.str().c_str() << ": File not found" << endl;
	*global_out_stream_ptr << corrfn_stream.str().c_str() << ": File not found" << endl;
	return (1);
}
  ifstream corrfn(corrfn_stream.str().c_str());

   for(int i = 0; i < qnpts; i++)
   {
            corrfn >> q_out[i];
            corrfn >> Correl_1D_out[i];
            corrfn >> Correl_1D_out_err[i];
            corrfn >> q_side[i];
            corrfn >> Correl_1D_side[i];
            corrfn >> Correl_1D_side_err[i];
            corrfn >> q_long[i];
            corrfn >> Correl_1D_long[i];
            corrfn >> Correl_1D_long_err[i];
   }

   return (0);
}

//End of file
