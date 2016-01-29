#include<iostream>
#include<sstream>
#include<string>
#include<fstream>
#include<cmath>
#include<iomanip>
#include<vector>
#include<stdio.h>

/*#include<gsl/gsl_sf_bessel.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>            // gsl random number generators
#include <gsl/gsl_randist.h>        // gsl random number distributions
#include <gsl/gsl_vector.h>         // gsl vector and matrix definitions
#include <gsl/gsl_blas.h>           // gsl linear algebra stuff
#include <gsl/gsl_multifit_nlin.h>  // gsl multidimensional fitting*/

#include "CFWR.h"
#include "Arsenal.h"
#include "gauss_quadrature.h"

using namespace std;

double slopes[4][2];

void CorrelationFunction::Get_GF_HBTradii(FO_surf* FOsurf_ptr, int folderindex)
{
	*global_out_stream_ptr << "--> Getting HBT radii by Gaussian fit method" << endl;
	for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
	for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
	{
		*global_out_stream_ptr << "   --> Doing pT = " << SPinterp_pT[ipt] << ", pphi = " << SPinterp_pphi[ipphi] << "..." << endl;
		if (USE_LAMBDA)
			Fit_Correlationfunction3D_withlambda( CFvals[ipt][ipphi], ipt, ipphi );
		else
			Fit_Correlationfunction3D( CFvals[ipt][ipphi], ipt, ipphi );
	}

	return;
}

double CorrelationFunction::Compute_correlationfunction(int ipt, int ipphi, double * q_interp)
{
	// try using linear-logarithmic interpolation if q-point is within grid
	// otherwise, just use linear interpolation/extrapolation, or throw an exception and do return 0
	double interpolated_result = 0.0;

	double nonFTd_spectra = spectra[target_particle_id][ipt][ipphi];

	int * qidx = new int [4];
	qidx[0] = binarySearch(qt_pts, qtnpts, q_interp[0]);
	qidx[1] = binarySearch(qx_pts, qxnpts, q_interp[1]);
	qidx[2] = binarySearch(qy_pts, qynpts, q_interp[2]);
	qidx[3] = binarySearch(qz_pts, qznpts, q_interp[3]);

	bool q_point_is_outside_grid = ( qidx[0] == -1 ) || ( qidx[1] == -1 ) || ( qidx[2] == -1 ) || ( qidx[3] == -1 );

	if (!q_point_is_outside_grid)
	{
		double * sgnd_q2_interp = new double [4];
		double q_min[4] = { qt_pts[qidx[0]], qx_pts[qidx[1]], qy_pts[qidx[2]], qz_pts[qidx[3]] };
		double q_max[4] = { qt_pts[qidx[0]+1], qx_pts[qidx[1]+1], qy_pts[qidx[2]+1], qz_pts[qidx[3]+1] };
		double * sgnd_q2_min = new double [4];
		double * sgnd_q2_max = new double [4];

		double ln_C_at_q[2][2][2][2];
		//double C_at_q[2][2][2][2];
		double ***** current_C_slice = current_dN_dypTdpTdphi_moments[ipt][ipphi];

		// set values at all vertices of 4D-hypercube used for interpolation
		for (int i = 0; i < 2; ++i)
		{
			double **** ccs1 = current_C_slice[qidx[0]+i];
			for (int j = 0; j < 2; ++j)
			{
				double *** ccs2 = ccs1[qidx[1]+j];
				for (int k = 0; k < 2; ++k)
				{
					double ** ccs3 = ccs2[qidx[2]+k];
					for (int l = 0; l < 2; ++l)
					{
						double * ccs4 = ccs3[qidx[3]+l];
						double tmp_cos = ccs4[0];
						double tmp_sin = ccs4[1];
						ln_C_at_q[i][j][k][l] = log( (tmp_cos*tmp_cos + tmp_sin*tmp_sin)/(nonFTd_spectra*nonFTd_spectra) );
						//C_at_q[i][j][k][l] = (tmp_cos*tmp_cos + tmp_sin*tmp_sin)/(nonFTd_spectra*nonFTd_spectra);
					}
				}
			}
		}

		// interpolating w.r.t. q^2, not q
		for (int i = 0; i < 4; ++i)
		{
			sgnd_q2_interp[i] = sgn(q_interp[i]) * q_interp[i] * q_interp[i];
			sgnd_q2_min[i] = sgn(q_min[i]) * q_min[i] * q_min[i];
			sgnd_q2_max[i] = sgn(q_max[i]) * q_max[i] * q_max[i];
		}

		double interpolated_ln_result = interpolate_4D(sgnd_q2_min, sgnd_q2_max, sgnd_q2_interp, ln_C_at_q);
		interpolated_result = exp(interpolated_ln_result);
		//double interpolated_result = interpolate_4D(q_min, q_max, q_interp, C_at_q);

		// Clean up
		delete [] sgnd_q2_interp;
		delete [] sgnd_q2_min;
		delete [] sgnd_q2_max;
	}
	else
	{
		//*global_out_stream_ptr << "Warning: q_interp point was outside of computed grid!" << endl
		//						<< "\t q_interp = {" << q_interp[0] << ", " << q_interp[1] << ", " << q_interp[2] << ", " << q_interp[3] << "}" << endl;
		interpolated_result = 0.0;
	}
	return (interpolated_result);
}

double CorrelationFunction::interpolate_4D(double * x_min, double * x_max, double * x_interp, double (*vals) [2][2][2])
{
	double total = 0.0;

	for (int i = 0; i < 4; ++i)
	{
		//slopes[i][0] = (x_interp[i] - x_min[i]) / (x_max[i] - x_min[i]);
		slopes[i][0] = (x_max[i] - x_interp[i]) / (x_max[i] - x_min[i]);
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
	if (qtnpts == 1 || qxnpts == 1 || qynpts == 1 || qznpts == 1)
		return;

	// chooses the qo, qs, ql points at which to evaluate correlation function,
	// and allocates the array to hold correlation function values
	Set_correlation_function_q_pts();
	Allocate_CFvals();

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

//**************************************************************
// Gaussian fit routines below
//**************************************************************

void CorrelationFunction::Fit_Correlationfunction3D(double *** Correl_3D, int ipt, int ipphi)
{
	const size_t data_length = qonpts*qsnpts*qlnpts;  // # of points
	const size_t n_para = 4;  // # of parameters
//debugger(__LINE__, __FILE__);
	// allocate space for a covariance matrix of size p by p
	gsl_matrix *covariance_ptr = gsl_matrix_alloc (n_para, n_para);
//debugger(__LINE__, __FILE__);
	// allocate and setup for generating gaussian distibuted random numbers
	gsl_rng_env_setup ();
	const gsl_rng_type *type = gsl_rng_default;
	gsl_rng *rng_ptr = gsl_rng_alloc (type);
//debugger(__LINE__, __FILE__);
	//set up test data
	struct Correlationfunction3D_data Correlfun3D_data;
	Correlfun3D_data.data_length = data_length;
	Correlfun3D_data.q_o = new double [data_length];
	Correlfun3D_data.q_s = new double [data_length];
	Correlfun3D_data.q_l = new double [data_length];
	Correlfun3D_data.y = new double [data_length];
	Correlfun3D_data.sigma = new double [data_length];
//debugger(__LINE__, __FILE__);
	int idx = 0;
	for(int i = 0; i < qonpts; i++)
	for(int j = 0; j < qsnpts; j++)
	for(int k = 0; k < qlnpts; k++)
	{
		Correlfun3D_data.q_o[idx] = qo_pts[i];
		Correlfun3D_data.q_s[idx] = qs_pts[j];
		Correlfun3D_data.q_l[idx] = ql_pts[k];
		// This sets up the data to be fitted, with gaussian noise added
		// Correlfun3D_data.y[idx] = 1.0*exp( - 0.81*q_out[i]*q_out[i] - 1.21*q_side[j]*q_side[j] - 4.0*q_long[k]*q_long[k] - 0.25*q_out[i]*q_side[j]) + gsl_ran_gaussian(rng_ptr, 	error);
		Correlfun3D_data.y[idx] = Correl_3D[i][j][k];
		Correlfun3D_data.sigma[idx] = Correl_3D_err[i][j][k];
		idx++;
	}
//debugger(__LINE__, __FILE__);
	double para_init[n_para] = { 1.0, 1.0, 1.0, 1.0 };  // initial guesses of parameters
//debugger(__LINE__, __FILE__);
	gsl_vector_view xvec_ptr = gsl_vector_view_array (para_init, n_para);
//debugger(__LINE__, __FILE__);
	// set up the function to be fit 
	gsl_multifit_function_fdf target_func;
	target_func.f = &Fittarget_correlfun3D_f;        // the function of residuals
	target_func.df = &Fittarget_correlfun3D_df;      // the gradient of this function
	target_func.fdf = &Fittarget_correlfun3D_fdf;    // combined function and gradient
	target_func.n = data_length;              // number of points in the data set
	target_func.p = n_para;              // number of parameters in the fit function
	target_func.params = &Correlfun3D_data;  // structure with the data and error bars
//debugger(__LINE__, __FILE__);
	const gsl_multifit_fdfsolver_type *type_ptr = gsl_multifit_fdfsolver_lmsder;
	gsl_multifit_fdfsolver *solver_ptr = gsl_multifit_fdfsolver_alloc (type_ptr, data_length, n_para);
	gsl_multifit_fdfsolver_set (solver_ptr, &target_func, &xvec_ptr.vector);
//debugger(__LINE__, __FILE__);
	size_t iteration = 0;         // initialize iteration counter
	print_fit_state_3D (iteration, solver_ptr);
	int status;  		// return value from gsl function calls (e.g., error)
	do
	{
		iteration++;
      
		// perform a single iteration of the fitting routine
		status = gsl_multifit_fdfsolver_iterate (solver_ptr);

		// print out the status of the fit
		cout << "status = " << gsl_strerror (status) << endl;

		// customized routine to print out current parameters
		print_fit_state_3D (iteration, solver_ptr);

		if (status)    // check for a nonzero status code
		{
			break;  // this should only happen if an error code is returned 
		}

		// test for convergence with an absolute and relative error (see manual)
		status = gsl_multifit_test_delta (solver_ptr->dx, solver_ptr->x, fit_tolerance, fit_tolerance);
	}
	while (status == GSL_CONTINUE && iteration < fit_max_iterations);
//debugger(__LINE__, __FILE__);
	//cerr >> "iteration = " << iteration << endl;
//debugger(__LINE__, __FILE__);
	// calculate the covariance matrix of the best-fit parameters
	gsl_multifit_covar (solver_ptr->J, 0.0, covariance_ptr);

	// print out the covariance matrix using the gsl function (not elegant!)
	cout << endl << "Covariance matrix: " << endl;
	gsl_matrix_fprintf (stdout, covariance_ptr, "%g");
//debugger(__LINE__, __FILE__);
	cout.setf (ios::fixed, ios::floatfield);	// output in fixed format
	cout.precision (5);		                // # of digits in doubles
//debugger(__LINE__, __FILE__);
	int width = 7;		// setw width for output
	cout << endl << "Best fit results:" << endl;
	cout << "R2o = " << setw (width) << get_fit_results (0, solver_ptr)
		<< " +/- " << setw (width) << get_fit_err (0, covariance_ptr) << endl;

	cout << "R2s      = " << setw (width) << get_fit_results (1, solver_ptr)
		<< " +/- " << setw (width) << get_fit_err (1, covariance_ptr) << endl;

	cout << "R2l      = " << setw (width) << get_fit_results (2, solver_ptr)
		<< " +/- " << setw (width) << get_fit_err (2, covariance_ptr) << endl;
  
	cout << "R2os      = " << setw (width) << get_fit_results (3, solver_ptr)
		<< " +/- " << setw (width) << get_fit_err (3, covariance_ptr) << endl;
    
	cout << "status = " << gsl_strerror (status) << endl;
	cout << "--------------------------------------------------------------------" << endl;

	double chi = gsl_blas_dnrm2(solver_ptr->f);
	double dof = data_length - n_para;
	double c = GSL_MAX_DBL(1, chi/sqrt(dof));
//debugger(__LINE__, __FILE__);
	lambda_Correl[ipt][ipphi] = 1.0;
	lambda_Correl_err[ipt][ipphi] = 0.0;
	R2_out[ipt][ipphi] = fabs(get_fit_results(0, solver_ptr))*hbarC*hbarC;
	R2_side[ipt][ipphi] = fabs(get_fit_results(1, solver_ptr))*hbarC*hbarC;
	R2_long[ipt][ipphi] = fabs(get_fit_results(2, solver_ptr))*hbarC*hbarC;
	R2_outside[ipt][ipphi] = fabs(get_fit_results(3, solver_ptr))*hbarC*hbarC;
	R2_out_err[ipt][ipphi] = c*get_fit_err(0, covariance_ptr)*hbarC*hbarC;
	R2_side_err[ipt][ipphi] = c*get_fit_err(1, covariance_ptr)*hbarC*hbarC;
	R2_long_err[ipt][ipphi] = c*get_fit_err(2, covariance_ptr)*hbarC*hbarC;
	R2_outside_err[ipt][ipphi] = c*get_fit_err(3, covariance_ptr)*hbarC*hbarC;
//debugger(__LINE__, __FILE__);
	cout << "final results: " << endl;
	cout << scientific << setw(10) << setprecision(5) 
		<< "chisq/dof = " << chi*chi/dof << endl;
	cout << scientific << setw(10) << setprecision(5) 
		<< " lambda[ipt][ipphi] = " << lambda_Correl[ipt][ipphi] << " +/- " << lambda_Correl_err[ipt][ipphi] << endl;
	cout << " R2_out[ipt][ipphi] = " << R2_out[ipt][ipphi] << " +/- " << R2_out_err[ipt][ipphi] << endl;
	cout << " R2_side[ipt][ipphi] = " << R2_side[ipt][ipphi] << " +/- " << R2_side_err[ipt][ipphi] << endl;
	cout << " R2_long[ipt][ipphi] = " << R2_long[ipt][ipphi] << " +/- " << R2_long_err[ipt][ipphi] << endl;
	cout << " R2_outside[ipt][ipphi] = " << R2_outside[ipt][ipphi] << " +/- " << R2_outside_err[ipt][ipphi] << endl;
//debugger(__LINE__, __FILE__);
	//clean up
	gsl_matrix_free (covariance_ptr);
	gsl_rng_free (rng_ptr);
//debugger(__LINE__, __FILE__);
	delete[] Correlfun3D_data.q_o;
	delete[] Correlfun3D_data.q_s;
	delete[] Correlfun3D_data.q_l;
	delete[] Correlfun3D_data.y;
	delete[] Correlfun3D_data.sigma;
//debugger(__LINE__, __FILE__);
	gsl_multifit_fdfsolver_free (solver_ptr);  // free up the solver

	return;
}

void CorrelationFunction::Fit_Correlationfunction3D_withlambda(double *** Correl_3D, int ipt, int ipphi)
{
	const size_t data_length = qonpts*qsnpts*qlnpts;  // # of points
	const size_t n_para = 5;  // # of parameters

	// allocate space for a covariance matrix of size p by p
	gsl_matrix *covariance_ptr = gsl_matrix_alloc (n_para, n_para);

	// allocate and setup for generating gaussian distibuted random numbers
	gsl_rng_env_setup ();
	const gsl_rng_type *type = gsl_rng_default;
	gsl_rng *rng_ptr = gsl_rng_alloc (type);

	//set up test data
	struct Correlationfunction3D_data Correlfun3D_data;
	Correlfun3D_data.data_length = data_length;
	Correlfun3D_data.q_o = new double [data_length];
	Correlfun3D_data.q_s = new double [data_length];
	Correlfun3D_data.q_l = new double [data_length];
	Correlfun3D_data.y = new double [data_length];
	Correlfun3D_data.sigma = new double [data_length];

	int idx = 0;
	for(int i = 0; i < qonpts; i++)
	for(int j = 0; j < qsnpts; j++)
	for(int k = 0; k < qlnpts; k++)
	{
		Correlfun3D_data.q_o[idx] = qo_pts[i];
		Correlfun3D_data.q_s[idx] = qs_pts[j];
		Correlfun3D_data.q_l[idx] = ql_pts[k];
		// This sets up the data to be fitted, with gaussian noise added
		// Correlfun3D_data.y[idx] = 1.0*exp( - 0.81*q_out[i]*q_out[i] - 1.21*q_side[j]*q_side[j] - 4.0*q_long[k]*q_long[k] - 0.25*q_out[i]*q_side[j]) + gsl_ran_gaussian(rng_ptr, error);
		Correlfun3D_data.y[idx] = Correl_3D[i][j][k];
		Correlfun3D_data.sigma[idx] = Correl_3D_err[i][j][k];
		idx++;
	}
	double para_init[n_para] = { 1.0, 1.0, 1.0, 1.0, 1.0 };  // initial guesses of parameters

	gsl_vector_view xvec_ptr = gsl_vector_view_array (para_init, n_para);
  
	// set up the function to be fit 
	gsl_multifit_function_fdf target_func;
	target_func.f = &Fittarget_correlfun3D_f_withlambda;        // the function of residuals
	target_func.df = &Fittarget_correlfun3D_df_withlambda;      // the gradient of this function
	target_func.fdf = &Fittarget_correlfun3D_fdf_withlambda;    // combined function and gradient
	target_func.n = data_length;              // number of points in the data set
	target_func.p = n_para;              // number of parameters in the fit function
	target_func.params = &Correlfun3D_data;  // structure with the data and error bars

	const gsl_multifit_fdfsolver_type *type_ptr = gsl_multifit_fdfsolver_lmsder;
	gsl_multifit_fdfsolver *solver_ptr = gsl_multifit_fdfsolver_alloc (type_ptr, data_length, n_para);
	gsl_multifit_fdfsolver_set (solver_ptr, &target_func, &xvec_ptr.vector);

	size_t iteration = 0;         // initialize iteration counter
	print_fit_state_3D (iteration, solver_ptr);
	int status;  		// return value from gsl function calls (e.g., error)
	do
	{
		iteration++;
      
		// perform a single iteration of the fitting routine
		status = gsl_multifit_fdfsolver_iterate (solver_ptr);

		// print out the status of the fit
		cout << "status = " << gsl_strerror (status) << endl;

		// customized routine to print out current parameters
		print_fit_state_3D (iteration, solver_ptr);

		if (status)    // check for a nonzero status code
		{
			break;  // this should only happen if an error code is returned 
		}

		// test for convergence with an absolute and relative error (see manual)
		status = gsl_multifit_test_delta (solver_ptr->dx, solver_ptr->x, fit_tolerance, fit_tolerance);
	}
	while (status == GSL_CONTINUE && iteration < fit_max_iterations);

	// calculate the covariance matrix of the best-fit parameters
	gsl_multifit_covar (solver_ptr->J, 0.0, covariance_ptr);

	// print out the covariance matrix using the gsl function (not elegant!)
	cout << endl << "Covariance matrix: " << endl;
	gsl_matrix_fprintf (stdout, covariance_ptr, "%g");

	cout.setf (ios::fixed, ios::floatfield);	// output in fixed format
	cout.precision (5);		                // # of digits in doubles

	int width = 7;		// setw width for output
	cout << endl << "Best fit results:" << endl;
	cout << "lambda      = " << setw (width) << get_fit_results (0, solver_ptr)
		<< " +/- " << setw (width) << get_fit_err (0, covariance_ptr) << endl;

	cout << "R2o = " << setw (width) << get_fit_results (1, solver_ptr)
		<< " +/- " << setw (width) << get_fit_err (1, covariance_ptr) << endl;

	cout << "R2s      = " << setw (width) << get_fit_results (2, solver_ptr)
		<< " +/- " << setw (width) << get_fit_err (2, covariance_ptr) << endl;

	cout << "R2l      = " << setw (width) << get_fit_results (3, solver_ptr)
		<< " +/- " << setw (width) << get_fit_err (3, covariance_ptr) << endl;
  
	cout << "R2os      = " << setw (width) << get_fit_results (4, solver_ptr)
		<< " +/- " << setw (width) << get_fit_err (4, covariance_ptr) << endl;
    
	cout << "status = " << gsl_strerror (status) << endl;
	cout << "--------------------------------------------------------------------" << endl;

	double chi = gsl_blas_dnrm2(solver_ptr->f);
	double dof = data_length - n_para;
	double c = GSL_MAX_DBL(1, chi/sqrt(dof));

	lambda_Correl[ipt][ipphi] = get_fit_results(0, solver_ptr);
	lambda_Correl_err[ipt][ipphi] = c*get_fit_err(0, covariance_ptr);
	R2_out[ipt][ipphi] = fabs(get_fit_results(1, solver_ptr))*hbarC*hbarC;
	R2_side[ipt][ipphi] = fabs(get_fit_results(2, solver_ptr))*hbarC*hbarC;
	R2_long[ipt][ipphi] = fabs(get_fit_results(3, solver_ptr))*hbarC*hbarC;
	R2_outside[ipt][ipphi] = fabs(get_fit_results(4, solver_ptr))*hbarC*hbarC;
	R2_out_err[ipt][ipphi] = c*get_fit_err(1, covariance_ptr)*hbarC*hbarC;
	R2_side_err[ipt][ipphi] = c*get_fit_err(2, covariance_ptr)*hbarC*hbarC;
	R2_long_err[ipt][ipphi] = c*get_fit_err(3, covariance_ptr)*hbarC*hbarC;
	R2_outside_err[ipt][ipphi] = c*get_fit_err(4, covariance_ptr)*hbarC*hbarC;

	cout << "final results: " << endl;
	cout << scientific << setw(10) << setprecision(5) 
		<< "chisq/dof = " << chi*chi/dof << endl;
	cout << scientific << setw(10) << setprecision(5) 
		<< " lambda = " << lambda_Correl << " +/- " << lambda_Correl_err << endl;
	cout << " R2_out[ipt][ipphi] = " << R2_out[ipt][ipphi] << " +/- " << R2_out_err[ipt][ipphi] << endl;
	cout << " R2_side[ipt][ipphi] = " << R2_side[ipt][ipphi] << " +/- " << R2_side_err[ipt][ipphi] << endl;
	cout << " R2_long[ipt][ipphi] = " << R2_long[ipt][ipphi] << " +/- " << R2_long_err[ipt][ipphi] << endl;
	cout << " R2_outside[ipt][ipphi] = " << R2_outside[ipt][ipphi] << " +/- " << R2_outside_err[ipt][ipphi] << endl;

	//clean up
	gsl_matrix_free (covariance_ptr);
	gsl_rng_free (rng_ptr);

	delete[] Correlfun3D_data.q_o;
	delete[] Correlfun3D_data.q_s;
	delete[] Correlfun3D_data.q_l;
	delete[] Correlfun3D_data.y;
	delete[] Correlfun3D_data.sigma;

	gsl_multifit_fdfsolver_free (solver_ptr);  // free up the solver

	return;
}

//*********************************************************************
// 3D case
//*********************************************************************
//  Simple function to print results of each iteration in nice format
int CorrelationFunction::print_fit_state_3D (size_t iteration, gsl_multifit_fdfsolver * solver_ptr)
{
	cout.setf (ios::fixed, ios::floatfield);	// output in fixed format
	cout.precision (5);		// digits in doubles

	int width = 15;		// setw width for output
	cout << scientific
		<< "iteration " << iteration << ": "
		<< "  x = {" << setw (width) << gsl_vector_get (solver_ptr->x, 0)
		<< setw (width) << gsl_vector_get (solver_ptr->x, 1)
		<< setw (width) << gsl_vector_get (solver_ptr->x, 2)
		<< setw (width) << gsl_vector_get (solver_ptr->x, 3)
		<< "}, |f(x)| = " << scientific << gsl_blas_dnrm2 (solver_ptr->f) 
		<< endl << endl;

	return 0;
}
//  Simple function to print results of each iteration in nice format
int CorrelationFunction::print_fit_state_3D_withlambda (size_t iteration, gsl_multifit_fdfsolver * solver_ptr)
{
	cout.setf (ios::fixed, ios::floatfield);	// output in fixed format
	cout.precision (5);		// digits in doubles

	int width = 15;		// setw width for output
	cout << scientific
		<< "iteration " << iteration << ": "
		<< "  x = {" << setw (width) << gsl_vector_get (solver_ptr->x, 0)
		<< setw (width) << gsl_vector_get (solver_ptr->x, 1)
		<< setw (width) << gsl_vector_get (solver_ptr->x, 2)
		<< setw (width) << gsl_vector_get (solver_ptr->x, 3)
		<< setw (width) << gsl_vector_get (solver_ptr->x, 4)
		<< "}, |f(x)| = " << scientific << gsl_blas_dnrm2 (solver_ptr->f) 
		<< endl << endl;

	return 0;
}
//*********************************************************************
//  Function returning the residuals for each point; that is, the 
//  difference of the fit function using the current parameters
//  and the data to be fit.
int Fittarget_correlfun3D_f (const gsl_vector *xvec_ptr, void *params_ptr, gsl_vector *f_ptr)
{
	size_t n = ((struct Correlationfunction3D_data *) params_ptr)->data_length;
	double * q_o = ((struct Correlationfunction3D_data *) params_ptr)->q_o;
	double * q_s = ((struct Correlationfunction3D_data *) params_ptr)->q_s;
	double * q_l = ((struct Correlationfunction3D_data *) params_ptr)->q_l;
	double * y = ((struct Correlationfunction3D_data *) params_ptr)->y;
	double * sigma = ((struct Correlationfunction3D_data *) params_ptr)->sigma;

	//fit parameters
	double R2_o = gsl_vector_get (xvec_ptr, 0);
	double R2_s = gsl_vector_get (xvec_ptr, 1);
	double R2_l = gsl_vector_get (xvec_ptr, 2);
	double R2_os = gsl_vector_get (xvec_ptr, 3);

	size_t i;

	for (i = 0; i < n; i++)
	{
		double Yi = exp(- q_l[i]*q_l[i]*R2_l - q_s[i]*q_s[i]*R2_s - q_o[i]*q_o[i]*R2_o - 2.*q_o[i]*q_s[i]*R2_os);
		gsl_vector_set (f_ptr, i, (Yi - y[i]) / sigma[i]);
	}

	return GSL_SUCCESS;
}

int Fittarget_correlfun3D_f_withlambda (const gsl_vector *xvec_ptr, void *params_ptr, gsl_vector *f_ptr)
{
	size_t n = ((struct Correlationfunction3D_data *) params_ptr)->data_length;
	double * q_o = ((struct Correlationfunction3D_data *) params_ptr)->q_o;
	double * q_s = ((struct Correlationfunction3D_data *) params_ptr)->q_s;
	double * q_l = ((struct Correlationfunction3D_data *) params_ptr)->q_l;
	double * y = ((struct Correlationfunction3D_data *) params_ptr)->y;
	double * sigma = ((struct Correlationfunction3D_data *) params_ptr)->sigma;

	//fit parameters
	double lambda = gsl_vector_get (xvec_ptr, 0);
	double R2_o = gsl_vector_get (xvec_ptr, 1);
	double R2_s = gsl_vector_get (xvec_ptr, 2);
	double R2_l = gsl_vector_get (xvec_ptr, 3);
	double R2_os = gsl_vector_get (xvec_ptr, 4);

	size_t i;

	for (i = 0; i < n; i++)
	{
		//double Yi = lambda*exp(- q_l[i]*q_l[i]*R_l*R_l - q_s[i]*q_s[i]*R_s*R_s
		//             - q_o[i]*q_o[i]*R_o*R_o - q_o[i]*q_s[i]*R_os*R_os);
		double Yi = lambda*exp(- q_l[i]*q_l[i]*R2_l - q_s[i]*q_s[i]*R2_s - q_o[i]*q_o[i]*R2_o - 2.*q_o[i]*q_s[i]*R2_os);
		gsl_vector_set (f_ptr, i, (Yi - y[i]) / sigma[i]);
	}

	return GSL_SUCCESS;
}

//*********************************************************************
//  Function returning the Jacobian of the residual function
int Fittarget_correlfun3D_df (const gsl_vector *xvec_ptr, void *params_ptr,  gsl_matrix *Jacobian_ptr)
{
	size_t n = ((struct Correlationfunction3D_data *) params_ptr)->data_length;
	double * q_o = ((struct Correlationfunction3D_data *) params_ptr)->q_o;
	double * q_s = ((struct Correlationfunction3D_data *) params_ptr)->q_s;
	double * q_l = ((struct Correlationfunction3D_data *) params_ptr)->q_l;
	double * sigma = ((struct Correlationfunction3D_data *) params_ptr)->sigma;

	//fit parameters
	double R2_o = gsl_vector_get (xvec_ptr, 0);
	double R2_s = gsl_vector_get (xvec_ptr, 1);
	double R2_l = gsl_vector_get (xvec_ptr, 2);
	double R2_os = gsl_vector_get (xvec_ptr, 3);

	size_t i;

	for (i = 0; i < n; i++)
	{
		// Jacobian matrix J(i,j) = dfi / dxj, 
		// where fi = (Yi - yi)/sigma[i],      
		//       Yi = A * exp(-lambda * i) + b 
		// and the xj are the parameters (A,lambda,b) 
		double sig = sigma[i];

		// derivatives
		// double common_elemt = exp(- q_l[i]*q_l[i]*R_l*R_l - q_s[i]*q_s[i]*R_s*R_s - q_o[i]*q_o[i]*R_o*R_o - q_o[i]*q_s[i]*R_os*R_os);
		double common_elemt = exp(- q_l[i]*q_l[i]*R2_l - q_s[i]*q_s[i]*R2_s - q_o[i]*q_o[i]*R2_o - 2.*q_o[i]*q_s[i]*R2_os);
      
		gsl_matrix_set (Jacobian_ptr, i, 0, - q_o[i]*q_o[i]*common_elemt/sig);
		gsl_matrix_set (Jacobian_ptr, i, 1, - q_s[i]*q_s[i]*common_elemt/sig);
		gsl_matrix_set (Jacobian_ptr, i, 2, - q_l[i]*q_l[i]*common_elemt/sig);
		gsl_matrix_set (Jacobian_ptr, i, 3, - 2.*q_o[i]*q_s[i]*common_elemt/sig);
	}

	return GSL_SUCCESS;
}

int Fittarget_correlfun3D_df_withlambda (const gsl_vector *xvec_ptr, void *params_ptr,  gsl_matrix *Jacobian_ptr)
{
	size_t n = ((struct Correlationfunction3D_data *) params_ptr)->data_length;
	double * q_o = ((struct Correlationfunction3D_data *) params_ptr)->q_o;
	double * q_s = ((struct Correlationfunction3D_data *) params_ptr)->q_s;
	double * q_l = ((struct Correlationfunction3D_data *) params_ptr)->q_l;
	double * sigma = ((struct Correlationfunction3D_data *) params_ptr)->sigma;

	//fit parameters
	double lambda = gsl_vector_get (xvec_ptr, 0);
	double R2_o = gsl_vector_get (xvec_ptr, 1);
	double R2_s = gsl_vector_get (xvec_ptr, 2);
	double R2_l = gsl_vector_get (xvec_ptr, 3);
	double R2_os = gsl_vector_get (xvec_ptr, 4);

	size_t i;

	for (i = 0; i < n; i++)
	{
		// Jacobian matrix J(i,j) = dfi / dxj, 
		// where fi = (Yi - yi)/sigma[i],      
		//       Yi = A * exp(-lambda * i) + b 
		// and the xj are the parameters (A,lambda,b) 
		double sig = sigma[i];

		//derivatives
		double common_elemt = exp(- q_l[i]*q_l[i]*R2_l - q_s[i]*q_s[i]*R2_s - q_o[i]*q_o[i]*R2_o - 2.*q_o[i]*q_s[i]*R2_os);
      
		gsl_matrix_set (Jacobian_ptr, i, 0, common_elemt/sig);
		gsl_matrix_set (Jacobian_ptr, i, 1, - lambda*q_o[i]*q_o[i]*common_elemt/sig);
		gsl_matrix_set (Jacobian_ptr, i, 2, - lambda*q_s[i]*q_s[i]*common_elemt/sig);
		gsl_matrix_set (Jacobian_ptr, i, 3, - lambda*q_l[i]*q_l[i]*common_elemt/sig);
		gsl_matrix_set (Jacobian_ptr, i, 4, - 2.*lambda*q_o[i]*q_s[i]*common_elemt/sig);
	}

	return GSL_SUCCESS;
}

//*********************************************************************
//  Function combining the residual function and its Jacobian
int Fittarget_correlfun3D_fdf (const gsl_vector* xvec_ptr, void *params_ptr, gsl_vector* f_ptr, gsl_matrix* Jacobian_ptr)
{
	Fittarget_correlfun3D_f(xvec_ptr, params_ptr, f_ptr);
	Fittarget_correlfun3D_df(xvec_ptr, params_ptr, Jacobian_ptr);

	return GSL_SUCCESS;
}

int Fittarget_correlfun3D_fdf_withlambda (const gsl_vector* xvec_ptr, void *params_ptr, gsl_vector* f_ptr, gsl_matrix* Jacobian_ptr)
{
	Fittarget_correlfun3D_f_withlambda(xvec_ptr, params_ptr, f_ptr);
	Fittarget_correlfun3D_df_withlambda(xvec_ptr, params_ptr, Jacobian_ptr);

	return GSL_SUCCESS;
}

//*********************************************************************
//  Function to return the i'th best-fit parameter
inline double CorrelationFunction::get_fit_results(int i, gsl_multifit_fdfsolver * solver_ptr)
{
	return gsl_vector_get (solver_ptr->x, i);
}

//*********************************************************************
//  Function to retrieve the square root of the diagonal elements of
//   the covariance matrix.
inline double CorrelationFunction::get_fit_err (int i, gsl_matrix * covariance_ptr)
{
	return sqrt (gsl_matrix_get (covariance_ptr, i, i));
}

/************************************************************************/

//End of file
