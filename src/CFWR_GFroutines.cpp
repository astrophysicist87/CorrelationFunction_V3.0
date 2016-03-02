#include<iostream>
#include<sstream>
#include<string>
#include<fstream>
#include<cmath>
#include<iomanip>
#include<vector>
#include<stdio.h>

#include "CFWR.h"
#include "CPStopwatch.h"
#include "Arsenal.h"
#include "gauss_quadrature.h"
#include "chebyshev.h"

using namespace std;

void CorrelationFunction::Get_GF_HBTradii(FO_surf* FOsurf_ptr, int folderindex)
{
	if (!VARY_ALPHA)
		*global_out_stream_ptr << "--> Getting HBT radii by Gaussian fit method" << endl;
	else
		*global_out_stream_ptr << "--> Getting HBT radii by Levy-stable fit method" << endl;

	for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
	for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
	{
		*global_out_stream_ptr << "   --> Doing pT = " << SPinterp_pT[ipt] << ", pphi = " << SPinterp_pphi[ipphi] << "..." << endl;
		if (!VARY_ALPHA)
		{
			if (USE_LAMBDA)
				Fit_Correlationfunction3D_withlambda( CFvals[ipt][ipphi], ipt, ipphi );
			else
				Fit_Correlationfunction3D( CFvals[ipt][ipphi], ipt, ipphi );
		}
		else	//varying alpha not yet supported!
		{
			//if (USE_LAMBDA)
			//	Fit_Correlationfunction3D_withlambda( CFvals[ipt][ipphi], ipt, ipphi );
			//else
			//	Fit_Correlationfunction3D( CFvals[ipt][ipphi], ipt, ipphi );
			cerr << "Varying alpha not yet supported!" << endl;
		}

	}

	return;
}

double CorrelationFunction::Compute_correlationfunction(int ipt, int ipphi, int iqx, int iqy, int iqz, double qt_interp, int interp_flag /*==0*/)
{
	// try using linear-logarithmic interpolation if q-point is within grid
	// otherwise, just use linear interpolation/extrapolation, or throw an exception and do return 0
	double interpolated_result = 0.0;

	double nonFTd_spectra = spectra[target_particle_id][ipt][ipphi];

	int qidx = binarySearch(qt_PTdep_pts[ipt], qtnpts, qt_interp);
	double q_min = qt_PTdep_pts[ipt][0] / cos(M_PI / (2.*qtnpts)), q_max = qt_PTdep_pts[ipt][qtnpts-1]/ cos(M_PI / (2.*qtnpts));

	bool q_point_is_outside_grid = ( qidx == -1 && ( qt_interp < q_min || qt_interp > q_max ) );
	bool use_linear_logarithmic_interpolation = false;

	if (!q_point_is_outside_grid)
	{
		if (use_linear_logarithmic_interpolation)
		{
			q_min = qt_PTdep_pts[ipt][qidx];
			q_max = qt_PTdep_pts[ipt][qidx+1];

			double ln_C_at_q[2];
	
			double ***** current_C_slice = current_dN_dypTdpTdphi_moments[ipt][ipphi];
	
			// set values at all vertices of 4D-hypercube used for interpolation
			for (int iqtidx = 0; iqtidx < 1; ++iqtidx)
			{
				double tmp_cos = current_C_slice[qidx+iqtidx][iqx][iqy][iqz][0];
				double tmp_sin = current_C_slice[qidx+iqtidx][iqx][iqy][iqz][1];
				ln_C_at_q[iqtidx] = log( (tmp_cos*tmp_cos + tmp_sin*tmp_sin)/(nonFTd_spectra*nonFTd_spectra) + 1.e-100 );
			}
	
			// interpolating w.r.t. q^2, not q
			double sgnd_qt2_interp = sgn(qt_interp) * qt_interp * qt_interp;
			double * sgnd_qt2_limits = new double [2];
			sgnd_qt2_limits[0] = sgn(q_min) * q_min * q_min;
			sgnd_qt2_limits[1] = sgn(q_max) * q_max * q_max;
	
			double interpolated_ln_result = interpolate1D(sgnd_qt2_limits, ln_C_at_q, sgnd_qt2_interp, 2, 0, true);
			interpolated_result = exp(interpolated_ln_result);
	
			// Clean up
			delete [] sgnd_qt2_limits;
		}
		else
		{
			//q_min = qt_PTdep_pts[ipt][0];
			//q_max = qt_PTdep_pts[ipt][qtnpts-1];
			double C_at_q[qtnpts];
			double ***** current_C_slice = current_dN_dypTdpTdphi_moments[ipt][ipphi];
	
			// set CF values along qt-slice for interpolation
			for (int iqtidx = 0; iqtidx < qtnpts; ++iqtidx)
			{
				double tmp_cos = current_C_slice[iqtidx][iqx][iqy][iqz][0];
				double tmp_sin = current_C_slice[iqtidx][iqx][iqy][iqz][1];
				C_at_q[iqtidx] = (tmp_cos*tmp_cos + tmp_sin*tmp_sin)/(nonFTd_spectra*nonFTd_spectra);
//cerr << "INFODUMP: " << iqtidx << "   " << iqx << "   " << iqy << "   " << iqz << "   " << ipt << "   " << ipphi << "   "
//			<< qt_PTdep_pts[ipt][iqtidx] << "   " << qx_PTdep_pts[ipt][iqx] << "   " << qy_PTdep_pts[ipt][iqy] << "   " << qz_PTdep_pts[ipt][iqz] << "   "
//			<< SPinterp_pT[ipt] << "   " << SPinterp_pphi[ipphi] << "   " << nonFTd_spectra << "   " << tmp_cos << "   " << tmp_sin << endl;
			}

			//assumes qt-grid has already been computed at (adjusted) Chebyshev nodes!!!
			if (QT_POINTS_SPACING == 1 && interp_flag == 0)
			{
				//set up Chebyshev calculation
				int npts_loc[1] = { qtnpts };
				int os[1] = { qtnpts - 1 };
				double lls[1] = { q_min };
				double uls[1] = { q_max };
				double point[1] = { qt_interp };
				int dim_loc = 1;

				Chebyshev cf(C_at_q, npts_loc, os, lls, uls, dim_loc);
	
				interpolated_result = cf.eval(point);
/*if (iqx==2 and iqy==2 and iqz==3)
{
*global_out_stream_ptr << "Setting up Chebyshev approximation (C_at_q = {" << C_at_q[0];
	for (int iqtidx = 1; iqtidx < qtnpts; ++iqtidx)
		*global_out_stream_ptr << ", " << C_at_q[iqtidx];
*global_out_stream_ptr << "}, spectra = " << nonFTd_spectra << "): interpolated_result = " << interpolated_result << endl;
}*/
			}
			else	//if not using Chebyshev nodes in qt-direction, just use straight-up linear(0) or cubic(1) interpolation
			{
*global_out_stream_ptr << "Using cubic interpolation" << endl;
				//interpolated_result = interpolate1D(qt_PTdep_pts[ipt], C_at_q, qt_interp, qtnpts, 0, true);
				interpolated_result = interpolate1D(qt_PTdep_pts[ipt], C_at_q, qt_interp, qtnpts, 1, false);
			}
		}
	}
	else
	{
		*global_out_stream_ptr << "Warning: qt_interp point was outside of computed grid!" << endl
								<< "\t qt_interp = " << qt_interp << " out of {q_min, q_max} = {" << q_min << ", " << q_max << "}" << endl;
		interpolated_result = 0.0;
	}
	return (interpolated_result);
}

//to save time, don't bother making grid large enough to interpolate to OSL
//this function leaves open the option of just interpolating over the qt-direction
void CorrelationFunction::Cal_correlationfunction()
{
	// Can't interpolate if there's only one point in qt-direction!
	if (qtnpts == 1)
		return;

	// store in HDF5 file
	int getHDFresonanceSpectra = Get_resonance_from_HDF_array(target_particle_id, current_dN_dypTdpTdphi_moments);
	if (getHDFresonanceSpectra < 0)
	{
		cerr << "Failed to set this resonance in HDF array!  Exiting..." << endl;
		exit;
	}

	// chooses the qo, qs, ql (or qx, qy, ql) points at which to evaluate correlation function,
	// and allocates the array to hold correlation function values
	Set_correlation_function_q_pts();
	Allocate_CFvals();

	double * q_interp = new double [4];

	// Then compute full correlation function
	*global_out_stream_ptr << "Computing the full correlator in XYZ coordinates..." << endl;
	CPStopwatch sw;
	sw.Start();
	for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
	for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
	for (int iqx = 0; iqx < qxnpts; ++iqx)
	for (int iqy = 0; iqy < qynpts; ++iqy)
	for (int iqz = 0; iqz < qznpts; ++iqz)
	{
		Get_q_points(qx_PTdep_pts[ipt][iqx], qy_PTdep_pts[ipt][iqy], qz_PTdep_pts[ipt][iqz], SPinterp_pT[ipt], SPinterp_pphi[ipphi], q_interp);

		double tmp_spectra = spectra[target_particle_id][ipt][ipphi];

		//CFvals[ipt][ipphi][iqx][iqy][iqz] = Compute_correlationfunction(ipt, ipphi, iqx, iqy, iqz, q_interp[0]);
		double interpcheck = Compute_correlationfunction(ipt, ipphi, iqx, iqy, iqz, q_interp[0], 0);
		//cout << "INTERPOLATION CHECK: " << scientific << setprecision(8) << setw(12) << interpcheck << "   " << SPinterp_pT[ipt] << "   " << SPinterp_pphi[ipphi] << "   "
		//		<< q_interp[0] << "   " << qx_PTdep_pts[ipt][iqx] << "   " << qy_PTdep_pts[ipt][iqy] << "   " << qz_PTdep_pts[ipt][iqz] << "   "
		//		<< Cal_dN_dypTdpTdphi_with_weights_function(current_FOsurf_ptr, target_particle_id, SPinterp_pT[ipt], SPinterp_pphi[ipphi],
		//														q_interp[0], qx_PTdep_pts[ipt][iqx], qy_PTdep_pts[ipt][iqy], qz_PTdep_pts[ipt][iqz]) / (tmp_spectra*tmp_spectra) << endl;
		CFvals[ipt][ipphi][iqx][iqy][iqz] = interpcheck;
	}
	sw.Stop();
	*global_out_stream_ptr << "Finished computing correlator in " << sw.printTime() << " seconds." << endl;

	//put test_interpolator() here...
	//*global_out_stream_ptr << "Testing interpolator..." << endl;
	//sw.Start();
	//test_interpolator();
	//sw.Stop();
	//*global_out_stream_ptr << "Finished testing interpolator in " << sw.printTime() << " seconds." << endl;

	delete [] q_interp;

	//exit(1);

	return;
}

void CorrelationFunction::test_interpolator()
{
	//need grid_points to be Chebyshev-spaced...
	if (QX_POINTS_SPACING != 1 or QY_POINTS_SPACING != 1 or QZ_POINTS_SPACING != 1)
		return;

	int local_qnpts = 25;
	double mtarget = all_particles[target_particle_id].mass;
	for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
	{
		double qxmin = qx_PTdep_pts[ipt][0];
		double qxmax = qx_PTdep_pts[ipt][qxnpts-1];
		double qymin = qy_PTdep_pts[ipt][0];
		double qymax = qy_PTdep_pts[ipt][qynpts-1];
		double qzmin = qz_PTdep_pts[ipt][0];
		double qzmax = qz_PTdep_pts[ipt][qznpts-1];
		double del_qx = (qxmax - qxmin) / double(local_qnpts - 1 + 1.e-100);
		double del_qy = (qymax - qymin) / double(local_qnpts - 1 + 1.e-100);
		double del_qz = (qzmax - qzmin) / double(local_qnpts - 1 + 1.e-100);

		//assumes q(i)-grids has already been computed at (adjusted) Chebyshev nodes!!!
		//set up Chebyshev calculation
		const int dim_loc = 3;
		int npts_loc[dim_loc] = { qxnpts, qynpts, qznpts };
		int os[dim_loc] = { qxnpts-1, qynpts-1, qznpts-1 };
		double lls[dim_loc] = { qxmin, qymin, qzmin };
		double uls[dim_loc] = { qxmax, qymax, qzmax };
		
		for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
		{
			//set up local grid for interpolation
			int iq = 0;
			double C_at_q[qxnpts*qynpts*qznpts];
			double *** current_C_slice = CFvals[ipt][ipphi];
			for (int iqx = 0; iqx < qxnpts; ++iqx)
			for (int iqy = 0; iqy < qynpts; ++iqy)
			for (int iqz = 0; iqz < qznpts; ++iqz)
				C_at_q[iq++] = current_C_slice[iqx][iqy][iqz];
	
			//generate Chebyshev approximation to use below
			Chebyshev cf(C_at_q, npts_loc, os, lls, uls, dim_loc);

			double local_pT = SPinterp_pT[ipt];
			double local_pphi = SPinterp_pphi[ipphi];
			double ckp = cos_SPinterp_pphi[ipphi], skp = sin_SPinterp_pphi[ipphi];
			double local_spectra = spectra[target_particle_id][ipt][ipphi];
			//now use Chebyshev approximator to do interpolation
			for (int iqx = 0; iqx < local_qnpts; ++iqx)
			for (int iqy = 0; iqy < local_qnpts; ++iqy)
			for (int iqz = 0; iqz < local_qnpts; ++iqz)
			{
				//set interpolation point
				double qx_interp = qxmin + (double)iqx * del_qx;
				double qy_interp = qymin + (double)iqy * del_qy;
				double qz_interp = qzmin + (double)iqz * del_qz;
				double xi2 = mtarget*mtarget + local_pT*local_pT + 0.25*(qx_interp*qx_interp+qy_interp*qy_interp+qz_interp*qz_interp);
				double qo = ckp * qx_interp + skp * qy_interp;
				double qt_interp = sqrt(xi2 + qo*local_pT) - sqrt(xi2 - qo*local_pT);	//set qt component
	
				double point[dim_loc] = { qx_interp, qy_interp, qz_interp };
		
				double interpolated_result = cf.eval(point);
				double exact_result = Cal_dN_dypTdpTdphi_with_weights_function(current_FOsurf_ptr, target_particle_id, local_pT, local_pphi,
																				qt_interp, qx_interp, qy_interp, qz_interp) / (local_spectra*local_spectra);

				cout << "test_interpolator(): " << ipt << "   " << ipphi << "   " << iqx << "   " << iqy << "   " << iqz
						<< setprecision(12) << setw(16)
						<< "   " << local_pT << "   " << local_pphi << "   " << qx_interp << "   " << qy_interp << "   " << qz_interp
						<< "   " << interpolated_result << "   " << exact_result << endl;
			}
		}
	}
	return;
}

//**************************************************************
// Gaussian fit routines below
//**************************************************************

void CorrelationFunction::Fit_Correlationfunction3D(double *** Correl_3D, int ipt, int ipphi)
{
	const size_t data_length = q1npts*q2npts*q3npts;  // # of points
	const size_t n_para = 4;  // # of parameters

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
	double ckp = cos_SPinterp_pphi[ipphi], skp = sin_SPinterp_pphi[ipphi];
	for (int i = 0; i < q1npts; i++)
	for (int j = 0; j < q2npts; j++)
	for (int k = 0; k < q3npts; k++)
	{
		Correlfun3D_data.q_o[idx] = qx_PTdep_pts[ipt][i] * ckp + qy_PTdep_pts[ipt][j] * skp;
		Correlfun3D_data.q_s[idx] = -qx_PTdep_pts[ipt][i] * skp + qy_PTdep_pts[ipt][j] * ckp;
		Correlfun3D_data.q_l[idx] = qz_PTdep_pts[ipt][k];
		Correlfun3D_data.y[idx] = Correl_3D[i][j][k];
		Correlfun3D_data.sigma[idx] = Correl_3D_err[i][j][k];
		idx++;
	}

	double para_init[n_para] = { 1.0, 1.0, 1.0, 1.0 };  // initial guesses of parameters

	gsl_vector_view xvec_ptr = gsl_vector_view_array (para_init, n_para);

	// set up the function to be fit 
	gsl_multifit_function_fdf target_func;
	target_func.f = &Fittarget_correlfun3D_f;        // the function of residuals
	target_func.df = &Fittarget_correlfun3D_df;      // the gradient of this function
	target_func.fdf = &Fittarget_correlfun3D_fdf;    // combined function and gradient
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

	//cerr >> "iteration = " << iteration << endl;

	// calculate the covariance matrix of the best-fit parameters
	gsl_multifit_covar (solver_ptr->J, 0.0, covariance_ptr);

	// print out the covariance matrix using the gsl function (not elegant!)
	cout << endl << "Covariance matrix: " << endl;
	gsl_matrix_fprintf (stdout, covariance_ptr, "%g");

	cout.setf (ios::fixed, ios::floatfield);	// output in fixed format
	cout.precision (5);		                // # of digits in doubles

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

	cout << "final results: " << endl;
	cout << scientific << setw(10) << setprecision(5) 
		<< "chisq/dof = " << chi*chi/dof << endl;
	cout << scientific << setw(10) << setprecision(5) 
		<< " lambda[ipt][ipphi] = " << lambda_Correl[ipt][ipphi] << " +/- " << lambda_Correl_err[ipt][ipphi] << endl;
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

void CorrelationFunction::Fit_Correlationfunction3D_withlambda(double *** Correl_3D, int ipt, int ipphi)
{
	const size_t data_length = q1npts*q2npts*q3npts;  // # of points
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
	double ckp = cos_SPinterp_pphi[ipphi], skp = sin_SPinterp_pphi[ipphi];
	for (int i = 0; i < q1npts; i++)
	for (int j = 0; j < q2npts; j++)
	for (int k = 0; k < q3npts; k++)
	{
		Correlfun3D_data.q_o[idx] = qx_PTdep_pts[ipt][i] * ckp + qy_PTdep_pts[ipt][j] * skp;
		Correlfun3D_data.q_s[idx] = -qx_PTdep_pts[ipt][i] * skp + qy_PTdep_pts[ipt][j] * ckp;
		Correlfun3D_data.q_l[idx] = qz_PTdep_pts[ipt][k];
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
