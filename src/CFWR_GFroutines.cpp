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
#include "CPStopwatch.h"
#include "Arsenal.h"
#include "gauss_quadrature.h"
#include "chebyshev.h"

using namespace std;

void CorrelationFunction::Get_GF_HBTradii(int folderindex)
{
	if (!VARY_ALPHA)
		*global_out_stream_ptr << "--> Getting HBT radii by Gaussian fit method" << endl;
	else
		*global_out_stream_ptr << "--> Getting HBT radii by Levy-stable fit method" << endl;

	if (FLESH_OUT_CF)
		Allocate_fleshed_out_CF();

	for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
	for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
	{
		*global_out_stream_ptr << "   --> Doing pT = " << SPinterp_pT[ipt] << ", pphi = " << SPinterp_pphi[ipphi] << "..." << endl;
		
		double *** CF_for_fitting = CFvals[ipt][ipphi];
		if (FLESH_OUT_CF)
		{
			Flesh_out_CF(ipt, ipphi);
			CF_for_fitting = fleshed_out_CF;
		}

		if (!VARY_ALPHA)
		{
			if (USE_LAMBDA)
				Fit_Correlationfunction3D_withlambda( CF_for_fitting, ipt, ipphi, FLESH_OUT_CF );
			else
				Fit_Correlationfunction3D( CF_for_fitting, ipt, ipphi, FLESH_OUT_CF );
		}
		else	//varying alpha not yet supported!
		{
			//if (USE_LAMBDA)
			//	Fit_Correlationfunction3D_withlambda( CF_for_fitting, ipt, ipphi, FLESH_OUT_CF );
			//else
			//	Fit_Correlationfunction3D( CF_for_fitting, ipt, ipphi, FLESH_OUT_CF );
			cerr << "Varying alpha not yet supported!" << endl;
		}

		if (FLESH_OUT_CF && ipt==0 && ipphi==0)
			Output_fleshed_out_correlationfunction(ipt, ipphi);
//if (1) exit(1);
	}

	if (FLESH_OUT_CF)
		Delete_fleshed_out_CF();

	return;
}

double CorrelationFunction::Compute_correlationfunction(int ipt, int ipphi, int iqx, int iqy, int iqz, double qt_interp, int interp_flag /*==0*/)
{
	// try using linear-logarithmic interpolation if q-point is within grid
	// otherwise, just use linear interpolation/extrapolation, or throw an exception and do return 0
	double interpolated_result = 0.0;

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
	
			// set values at all vertices of 4D-hypercube used for interpolation
			for (int iqtidx = 0; iqtidx < 1; ++iqtidx)
				ln_C_at_q[iqtidx] = log( get_CF(ipt, ipphi, qidx+iqtidx, iqx, iqy, iqz, FIT_WITH_PROJECTED_CFVALS && !thermal_pions_only) + 1.e-100 );
	
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
			double C_at_q[qtnpts];	//C - 1
	
			// set CF values along qt-slice for interpolation
			for (int iqtidx = 0; iqtidx < qtnpts; ++iqtidx)
				C_at_q[iqtidx] = get_CF(ipt, ipphi, iqtidx, iqx, iqy, iqz, FIT_WITH_PROJECTED_CFVALS && !thermal_pions_only);

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

void CorrelationFunction::Compute_correlationfunction(double * totalresult, double * thermalresult, double * nonthermalresult,
										int ipt, int ipphi, int iqx, int iqy, int iqz, double qt_interp, int interp_flag /*==0*/)
{
	// try using linear-logarithmic interpolation if q-point is within grid
	// otherwise, just use linear interpolation/extrapolation, or throw an exception and do return 0

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

			double ln_C_at_q[2], ln_Ct_at_q[2], ln_Cnt_at_q[2];
			double tmpC = 0.0, tmpCt = 0.0, tmpCnt = 0.0;
	
			// set values at all vertices of 4D-hypercube used for interpolation
			for (int iqtidx = 0; iqtidx < 1; ++iqtidx)
			{
				//ln_C_at_q[iqtidx] = log( get_CF(ipt, ipphi, qidx+iqtidx, iqx, iqy, iqz, FIT_WITH_PROJECTED_CFVALS && !thermal_pions_only) + 1.e-100 );
				//return C - 1!!!
				get_CF(&tmpC, &tmpCt, &tmpCnt, ipt, ipphi, qidx+iqtidx, iqx, iqy, iqz);
				ln_C_at_q[iqtidx] = log( abs(tmpC) + 1.e-100 );
				ln_Ct_at_q[iqtidx] = log( abs(tmpCt) + 1.e-100 );
				ln_Cnt_at_q[iqtidx] = log( abs(tmpCnt) + 1.e-100 );
			}
	
			// interpolating w.r.t. q^2, not q
			double sgnd_qt2_interp = sgn(qt_interp) * qt_interp * qt_interp;
			double * sgnd_qt2_limits = new double [2];
			sgnd_qt2_limits[0] = sgn(q_min) * q_min * q_min;
			sgnd_qt2_limits[1] = sgn(q_max) * q_max * q_max;
	
			double interpolated_ln_result = interpolate1D(sgnd_qt2_limits, ln_C_at_q, sgnd_qt2_interp, 2, 0, true);
			*totalresult = exp(interpolated_ln_result);
			interpolated_ln_result = interpolate1D(sgnd_qt2_limits, ln_Ct_at_q, sgnd_qt2_interp, 2, 0, true);
			*thermalresult = exp(interpolated_ln_result);
			interpolated_ln_result = interpolate1D(sgnd_qt2_limits, ln_Cnt_at_q, sgnd_qt2_interp, 2, 0, true);
			*nonthermalresult = exp(interpolated_ln_result);
	
			// Clean up
			delete [] sgnd_qt2_limits;
		}
		else
		{
			double C_at_q[qtnpts], Ct_at_q[qtnpts], Cnt_at_q[qtnpts];	//C - 1
			double tmpC = 0.0, tmpCt = 0.0, tmpCnt = 0.0;
	
			// set CF values along qt-slice for interpolation
			for (int iqtidx = 0; iqtidx < qtnpts; ++iqtidx)
			{
				//C_at_q[iqtidx] = get_CF(ipt, ipphi, iqtidx, iqx, iqy, iqz, FIT_WITH_PROJECTED_CFVALS && !thermal_pions_only);
				//return C - 1!!!
				get_CF(&tmpC, &tmpCt, &tmpCnt, ipt, ipphi, iqtidx, iqx, iqy, iqz);
				C_at_q[iqtidx] = tmpC;
				Ct_at_q[iqtidx] = tmpCt;
				//Cnt_at_q[iqtidx] = tmpCnt;
				Cnt_at_q[iqtidx] = tmpC - tmpCt;
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
				Chebyshev cft(Ct_at_q, npts_loc, os, lls, uls, dim_loc);
				Chebyshev cfnt(Cnt_at_q, npts_loc, os, lls, uls, dim_loc);
	
				*totalresult = cf.eval(point);
				*thermalresult = cft.eval(point);
				*nonthermalresult = cfnt.eval(point);
			}
			else	//if not using Chebyshev nodes in qt-direction, just use straight-up linear(0) or cubic(1) interpolation
			{
				*global_out_stream_ptr << "Using cubic interpolation" << endl;
				//interpolated_result = interpolate1D(qt_PTdep_pts[ipt], C_at_q, qt_interp, qtnpts, 0, true);
				*totalresult = interpolate1D(qt_PTdep_pts[ipt], C_at_q, qt_interp, qtnpts, 1, false);
				*thermalresult = interpolate1D(qt_PTdep_pts[ipt], Ct_at_q, qt_interp, qtnpts, 1, false);
				*nonthermalresult = interpolate1D(qt_PTdep_pts[ipt], Cnt_at_q, qt_interp, qtnpts, 1, false);
			}
		}
	}
	else
	{
		*global_out_stream_ptr << "Warning: qt_interp point was outside of computed grid!" << endl
								<< "\t qt_interp = " << qt_interp << " out of {q_min, q_max} = {" << q_min << ", " << q_max << "}" << endl;
		*totalresult = 0.0;
	}
	return;
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

		//double interpcheck = Compute_correlationfunction(ipt, ipphi, iqx, iqy, iqz, q_interp[0], 0);
		//cout << "INTERPOLATION CHECK: " << scientific << setprecision(8) << setw(12) << interpcheck << "   " << SPinterp_pT[ipt] << "   " << SPinterp_pphi[ipphi] << "   "
		//		<< q_interp[0] << "   " << qx_PTdep_pts[ipt][iqx] << "   " << qy_PTdep_pts[ipt][iqy] << "   " << qz_PTdep_pts[ipt][iqz] << "   "
		//		<< Cal_dN_dypTdpTdphi_with_weights_function(current_FOsurf_ptr, target_particle_id, SPinterp_pT[ipt], SPinterp_pphi[ipphi],
		//														q_interp[0], qx_PTdep_pts[ipt][iqx], qy_PTdep_pts[ipt][iqy], qz_PTdep_pts[ipt][iqz]) / (tmp_spectra*tmp_spectra) << endl;
		//CFvals[ipt][ipphi][iqx][iqy][iqz] = interpcheck;
		double tmp1 = 0.0, tmp2 = 0.0, tmp3 = 0.0;
		Compute_correlationfunction(&tmp1, &tmp2, &tmp3, ipt, ipphi, iqx, iqy, iqz, q_interp[0], 0);
		CFvals[ipt][ipphi][iqx][iqy][iqz] = 1.0 + tmp1;			//C
		thermalCFvals[ipt][ipphi][iqx][iqy][iqz] = tmp2;	//Ct-1
		resonancesCFvals[ipt][ipphi][iqx][iqy][iqz] = tmp3;	//Cnt-1
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

/*void CorrelationFunction::test_interpolator_1Dslices()
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
		const int dim_loc = 1;
		int npts_loc_x[dim_loc] = { qxnpts };
		int npts_loc_y[dim_loc] = { qynpts };
		int npts_loc_z[dim_loc] = { qznpts };
		int osx[dim_loc] = { qxnpts-1 };
		int osy[dim_loc] = { qynpts-1 };
		int osz[dim_loc] = { qznpts-1 };
		double llsx[dim_loc] = { qxmin };
		double ulsx[dim_loc] = { qxmax };
		double llsy[dim_loc] = { qymin };
		double ulsy[dim_loc] = { qymax };
		double llsz[dim_loc] = { qzmin };
		double ulsz[dim_loc] = { qzmax };
		
		for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
		{
			//set up local q grids for interpolation
			int iq = 0;
			double C_at_qx[qxnpts];
			double *** current_C_slice = CFvals[ipt][ipphi];
			for (int iqx = 0; iqx < qxnpts; ++iqx)
				C_at_qx[iq++] = current_C_slice[iqx][iqy0][iqz0];
	
			//generate Chebyshev approximation to use below
			Chebyshev cfx(C_at_qx, npts_loc_x, osx, llsx, ulsx, dim_loc);

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
}*/

//******************************************************************
// Routines for refining correlation function grid via interpolation
//******************************************************************

void CorrelationFunction::Allocate_fleshed_out_CF()
{
	fleshed_out_CF = new double ** [new_nqpts];
	fleshed_out_thermal = new double ** [new_nqpts];
	fleshed_out_resonances = new double ** [new_nqpts];
	for (int iqx = 0; iqx < new_nqpts; ++iqx)
	{
		fleshed_out_CF[iqx] = new double * [new_nqpts];
		fleshed_out_thermal[iqx] = new double * [new_nqpts];
		fleshed_out_resonances[iqx] = new double * [new_nqpts];
		for (int iqy = 0; iqy < new_nqpts; ++iqy)
		{
			fleshed_out_CF[iqx][iqy] = new double [new_nqpts];
			fleshed_out_thermal[iqx][iqy] = new double [new_nqpts];
			fleshed_out_resonances[iqx][iqy] = new double [new_nqpts];
			for (int iqz = 0; iqz < new_nqpts; ++iqz)
			{
				fleshed_out_CF[iqx][iqy][iqz] = 0.0;
				fleshed_out_thermal[iqx][iqy][iqz] = 0.0;
				fleshed_out_resonances[iqx][iqy][iqz] = 0.0;
			}
		}
	}

	qx_fleshed_out_pts = new double [new_nqpts];
	qy_fleshed_out_pts = new double [new_nqpts];
	qz_fleshed_out_pts = new double [new_nqpts];

	return;
}

void CorrelationFunction::Delete_fleshed_out_CF()
{
	for (int iqx = 0; iqx < new_nqpts; ++iqx)
	{
		for (int iqy = 0; iqy < new_nqpts; ++iqy)
		{
			delete [] fleshed_out_CF[iqx][iqy];
			delete [] fleshed_out_thermal[iqx][iqy];
			delete [] fleshed_out_resonances[iqx][iqy];
		}
		delete [] fleshed_out_CF[iqx];
		delete [] fleshed_out_thermal[iqx];
		delete [] fleshed_out_resonances[iqx];
	}
	delete [] fleshed_out_CF;
	delete [] fleshed_out_thermal;
	delete [] fleshed_out_resonances;

	delete [] qx_fleshed_out_pts;
	delete [] qy_fleshed_out_pts;
	delete [] qz_fleshed_out_pts;

	return;
}

void CorrelationFunction::Flesh_out_CF(int ipt, int ipphi)
{
	//declare needed quantities here
	double qxmin = 0.9999*qx_PTdep_pts[ipt][0], qxmax = 0.9999*qx_PTdep_pts[ipt][qxnpts-1];
	double qymin = 0.9999*qy_PTdep_pts[ipt][0], qymax = 0.9999*qy_PTdep_pts[ipt][qynpts-1];
	double qzmin = 0.9999*qz_PTdep_pts[ipt][0], qzmax = 0.9999*qz_PTdep_pts[ipt][qznpts-1];

	double new_Del_qx = (qxmax - qxmin)/double(new_nqpts-1);
	double new_Del_qy = (qymax - qymin)/double(new_nqpts-1);
	double new_Del_qz = (qzmax - qzmin)/double(new_nqpts-1);

	//cout << "(qxmin, qxmax, new_Del_qx) = (" << qxmin << ", " << qxmax << ", " << new_Del_qx << ")" << endl;
	//cout << "(qymin, qymax, new_Del_qy) = (" << qymin << ", " << qymax << ", " << new_Del_qy << ")" << endl;
	//cout << "(qzmin, qzmax, new_Del_qz) = (" << qzmin << ", " << qzmax << ", " << new_Del_qz << ")" << endl;

	for (int iqx = 0; iqx < new_nqpts; ++iqx)
	for (int iqy = 0; iqy < new_nqpts; ++iqy)
	for (int iqz = 0; iqz < new_nqpts; ++iqz)
	{
		double qx0 = qxmin + double(iqx) * new_Del_qx;
		double qy0 = qymin + double(iqy) * new_Del_qy;
		double qz0 = qzmin + double(iqz) * new_Del_qz;
		qx_fleshed_out_pts[iqx] = qx0;
		qy_fleshed_out_pts[iqy] = qy0;
		qz_fleshed_out_pts[iqz] = qz0;
		//fleshed_out_CF[iqx][iqy][iqz] = interpolate_CF(CFvals[ipt][ipphi], qx0, qy0, qz0, ipt);
		fleshed_out_thermal[iqx][iqy][iqz] = interpolate_CF(thermalCFvals[ipt][ipphi], qx0, qy0, qz0, ipt, 0);
		fleshed_out_resonances[iqx][iqy][iqz] = interpolate_CF(resonancesCFvals[ipt][ipphi], qx0, qy0, qz0, ipt, 1);
		fleshed_out_CF[iqx][iqy][iqz] = 1.0 + fleshed_out_thermal[iqx][iqy][iqz] + fleshed_out_resonances[iqx][iqy][iqz];
	}

	return;
}

double CorrelationFunction::interpolate_CF(double *** current_C_slice, double qx0, double qy0, double qz0, int ipt, int thermal_or_resonances)
{
	//int thermal_or_resonances - 	0: interpolate thermal CF (assuming basically Gaussian)
	//								1: interpolate resonance contributions (linear near origin, exponential further out)

	int iqx0_loc = binarySearch(qx_PTdep_pts[ipt], qxnpts, qx0, true, true);
	int iqy0_loc = binarySearch(qy_PTdep_pts[ipt], qynpts, qy0, true, true);
	int iqz0_loc = binarySearch(qz_PTdep_pts[ipt], qznpts, qz0, true, true);

	if (iqx0_loc == -1 || iqy0_loc == -1 || iqz0_loc == -1)
	{
		cerr << "Interpolation failed: exiting!" << endl;
		exit(1);
	}

	//interpolate here...
	double qx_0 = qx_PTdep_pts[ipt][iqx0_loc];
	double qx_1 = qx_PTdep_pts[ipt][iqx0_loc+1];
	double qy_0 = qy_PTdep_pts[ipt][iqy0_loc];
	double qy_1 = qy_PTdep_pts[ipt][iqy0_loc+1];
	double qz_0 = qz_PTdep_pts[ipt][iqz0_loc];
	double qz_1 = qz_PTdep_pts[ipt][iqz0_loc+1];

	double fxiyizi = 0.0;
	if (thermal_or_resonances == 0)
	{
		double sqx02 = sgn(qx0)*qx0*qx0;
		double sqy02 = sgn(qy0)*qy0*qy0;
		double sqz02 = sgn(qz0)*qz0*qz0;
		double sqx_02 = sgn(qx_0)*qx_0*qx_0;
		double sqy_02 = sgn(qy_0)*qy_0*qy_0;
		double sqz_02 = sgn(qz_0)*qz_0*qz_0;
		double sqx_12 = sgn(qx_1)*qx_1*qx_1;
		double sqy_12 = sgn(qy_1)*qy_1*qy_1;
		double sqz_12 = sgn(qz_1)*qz_1*qz_1;

		//interpolate over qz-points first
		///////////////////////
		//interpolate over first pair of qz-points
		double fx0y0z0 = current_C_slice[iqx0_loc][iqy0_loc][iqz0_loc];
		double fx0y0z1 = current_C_slice[iqx0_loc][iqy0_loc][iqz0_loc+1];
		double fx0y0zi = interpolate_qi(sqz02, sqz_02, sqz_12, fx0y0z0, fx0y0z1, false);
	
		//interpolate over second pair of qz-points
		double fx0y1z0 = current_C_slice[iqx0_loc][iqy0_loc+1][iqz0_loc];
		double fx0y1z1 = current_C_slice[iqx0_loc][iqy0_loc+1][iqz0_loc+1];
		double fx0y1zi = interpolate_qi(sqz02, sqz_02, sqz_12, fx0y1z0, fx0y1z1, false);    
	
		//interpolate over third pair of qz-points
		double fx1y0z0 = current_C_slice[iqx0_loc+1][iqy0_loc][iqz0_loc];
		double fx1y0z1 = current_C_slice[iqx0_loc+1][iqy0_loc][iqz0_loc+1];
		double fx1y0zi = interpolate_qi(sqz02, sqz_02, sqz_12, fx1y0z0, fx1y0z1, false);
	
		//interpolate over fourth pair of qz-points
		double fx1y1z0 = current_C_slice[iqx0_loc+1][iqy0_loc+1][iqz0_loc];
		double fx1y1z1 = current_C_slice[iqx0_loc+1][iqy0_loc+1][iqz0_loc+1];
		double fx1y1zi = interpolate_qi(sqz02, sqz_02, sqz_12, fx1y1z0, fx1y1z1, false);
		///////////////////////
	
		//interpolate over qy-points next
		double fx0yizi = interpolate_qi(sqy02, sqy_02, sqy_12, fx0y0zi, fx0y1zi, false);
		double fx1yizi = interpolate_qi(sqy02, sqy_02, sqy_12, fx1y0zi, fx1y1zi, false);
	
		//finally, interpolate over qx-points
		fxiyizi = interpolate_qi(sqx02, sqx_02, sqx_12, fx0yizi, fx1yizi, false);
	}
	else if (thermal_or_resonances == 1)
	{
		int cqx = (qxnpts - 1)/2;           //central qx index
		int cqy = (qynpts - 1)/2;
		int cqz = (qznpts - 1)/2;
		int iqxh0 = (qxnpts - cqx - 1)/2;   //halfway between cqx and nqpts-1
		int iqyh0 = (qynpts - cqy - 1)/2;
		int iqzh0 = (qznpts - cqz - 1)/2;

		//if any of these are false, use logarithmic interpolation in that direction
		//(assuming approximately exponential fall-off)
		bool use_linear_qx = (iqx0_loc <= cqx + iqxh0) && (iqx0_loc >= cqx - iqxh0);
		bool use_linear_qy = (iqy0_loc <= cqy + iqyh0) && (iqy0_loc >= cqy - iqyh0);
		bool use_linear_qz = (iqz0_loc <= cqz + iqzh0) && (iqz0_loc >= cqz - iqzh0);

		//interpolate over qz-points first
		///////////////////////
		//interpolate over first pair of qz-points
		double fx0y0z0 = current_C_slice[iqx0_loc][iqy0_loc][iqz0_loc];
		double fx0y0z1 = current_C_slice[iqx0_loc][iqy0_loc][iqz0_loc+1];
		double fx0y0zi = interpolate_qi(qz0, qz_0, qz_1, fx0y0z0, fx0y0z1, use_linear_qz);
	
		//interpolate over second pair of qz-points
		double fx0y1z0 = current_C_slice[iqx0_loc][iqy0_loc+1][iqz0_loc];
		double fx0y1z1 = current_C_slice[iqx0_loc][iqy0_loc+1][iqz0_loc+1];
		double fx0y1zi = interpolate_qi(qz0, qz_0, qz_1, fx0y1z0, fx0y1z1, use_linear_qz);    
	
		//interpolate over third pair of qz-points
		double fx1y0z0 = current_C_slice[iqx0_loc+1][iqy0_loc][iqz0_loc];
		double fx1y0z1 = current_C_slice[iqx0_loc+1][iqy0_loc][iqz0_loc+1];
		double fx1y0zi = interpolate_qi(qz0, qz_0, qz_1, fx1y0z0, fx1y0z1, use_linear_qz);
	
		//interpolate over fourth pair of qz-points
		double fx1y1z0 = current_C_slice[iqx0_loc+1][iqy0_loc+1][iqz0_loc];
		double fx1y1z1 = current_C_slice[iqx0_loc+1][iqy0_loc+1][iqz0_loc+1];
		double fx1y1zi = interpolate_qi(qz0, qz_0, qz_1, fx1y1z0, fx1y1z1, use_linear_qz);
		///////////////////////
	
		//interpolate over qy-points next
		double fx0yizi = interpolate_qi(qy0, qy_0, qy_1, fx0y0zi, fx0y1zi, use_linear_qy);
		double fx1yizi = interpolate_qi(qy0, qy_0, qy_1, fx1y0zi, fx1y1zi, use_linear_qy);
	
		//finally, interpolate over qx-points
		fxiyizi = interpolate_qi(qx0, qx_0, qx_1, fx0yizi, fx1yizi, use_linear_qx);
	}

	return (fxiyizi);
}

double CorrelationFunction::interpolate_qi(double q0, double qi0, double qi1, double f1, double f2, bool use_linear)
{
	double if1 = f1;
	double if2 = f2;
	bool use_log = (!use_linear) && (f1 > 0.0) && (f2 > 0.0);
	if (use_log)
	{
		if1 = log(f1+1.e-100);
		if2 = log(f2+1.e-100);
	}
	double tmp_result = lin_int(q0 - qi0, 1./(qi1-qi0), if1, if2);
	if (use_log)
		tmp_result = exp(tmp_result);

	return (tmp_result);
}

//**************************************************************
// Gaussian fit routines below
//**************************************************************

void CorrelationFunction::Fit_Correlationfunction3D(double *** Correl_3D, int ipt, int ipphi, bool fleshing_out_CF /*== true*/)
{
	double * q1pts = qx_PTdep_pts[ipt];
	double * q2pts = qy_PTdep_pts[ipt];
	double * q3pts = qz_PTdep_pts[ipt];
	if (fleshing_out_CF)
	{
		q1npts = new_nqpts;
		q2npts = new_nqpts;
		q3npts = new_nqpts;
		q1pts = qx_fleshed_out_pts;
		q2pts = qy_fleshed_out_pts;
		q3pts = qz_fleshed_out_pts;
	}
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
		Correlfun3D_data.q_o[idx] = q1pts[i] * ckp + q2pts[j] * skp;
		Correlfun3D_data.q_s[idx] = -q1pts[i] * skp + q2pts[j] * ckp;
		Correlfun3D_data.q_l[idx] = q3pts[k];
		Correlfun3D_data.y[idx] = Correl_3D[i][j][k];
		//Correlfun3D_data.sigma[idx] = Correl_3D_err[i][j][k];
		Correlfun3D_data.sigma[idx] = 1.e-2;
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
	/*cout << endl << "Best fit results:" << endl;
	cout << "R2o = " << setw (width) << get_fit_results (0, solver_ptr)
		<< " +/- " << setw (width) << get_fit_err (0, covariance_ptr) << endl;

	cout << "R2s      = " << setw (width) << get_fit_results (1, solver_ptr)
		<< " +/- " << setw (width) << get_fit_err (1, covariance_ptr) << endl;

	cout << "R2l      = " << setw (width) << get_fit_results (2, solver_ptr)
		<< " +/- " << setw (width) << get_fit_err (2, covariance_ptr) << endl;
  
	cout << "R2os      = " << setw (width) << get_fit_results (3, solver_ptr)
		<< " +/- " << setw (width) << get_fit_err (3, covariance_ptr) << endl;
    
	cout << "status = " << gsl_strerror (status) << endl;
	cout << "--------------------------------------------------------------------" << endl;*/

	double chi = gsl_blas_dnrm2(solver_ptr->f);
	double dof = data_length - n_para;
	double c = GSL_MAX_DBL(1, chi/sqrt(dof));

	lambda_Correl[ipt][ipphi] = 1.0;
	lambda_Correl_err[ipt][ipphi] = 0.0;
	R2_out[ipt][ipphi] = fabs(get_fit_results(0, solver_ptr))*hbarC*hbarC;
	R2_side[ipt][ipphi] = fabs(get_fit_results(1, solver_ptr))*hbarC*hbarC;
	R2_long[ipt][ipphi] = fabs(get_fit_results(2, solver_ptr))*hbarC*hbarC;
	R2_outside[ipt][ipphi] = get_fit_results(3, solver_ptr)*hbarC*hbarC;
	R2_out_err[ipt][ipphi] = c*get_fit_err(0, covariance_ptr)*hbarC*hbarC;
	R2_side_err[ipt][ipphi] = c*get_fit_err(1, covariance_ptr)*hbarC*hbarC;
	R2_long_err[ipt][ipphi] = c*get_fit_err(2, covariance_ptr)*hbarC*hbarC;
	R2_outside_err[ipt][ipphi] = c*get_fit_err(3, covariance_ptr)*hbarC*hbarC;

	cout << "final results: " << endl;
	cout << scientific << setw(10) << setprecision(5) 
		<< "chisq/dof = " << chi*chi/dof << endl;
	cout << scientific << setw(10) << setprecision(5) 
		<< " lambda[ipt=" << ipt << "][ipphi=" << ipphi << "] = " << lambda_Correl[ipt][ipphi] << " +/- " << lambda_Correl_err[ipt][ipphi] << endl;
	cout << " R2_out[ipt=" << ipt << "][ipphi=" << ipphi << "] = " << R2_out[ipt][ipphi] << " +/- " << R2_out_err[ipt][ipphi] << endl;
	cout << " R2_side[ipt=" << ipt << "][ipphi=" << ipphi << "] = " << R2_side[ipt][ipphi] << " +/- " << R2_side_err[ipt][ipphi] << endl;
	cout << " R2_long[ipt=" << ipt << "][ipphi=" << ipphi << "] = " << R2_long[ipt][ipphi] << " +/- " << R2_long_err[ipt][ipphi] << endl;
	cout << " R2_outside[ipt=" << ipt << "][ipphi=" << ipphi << "] = " << R2_outside[ipt][ipphi] << " +/- " << R2_outside_err[ipt][ipphi] << endl;

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

void CorrelationFunction::Fit_Correlationfunction3D_withlambda(double *** Correl_3D, int ipt, int ipphi, bool fleshing_out_CF /*== true*/)
{
	double * q1pts = qx_PTdep_pts[ipt];
	double * q2pts = qy_PTdep_pts[ipt];
	double * q3pts = qz_PTdep_pts[ipt];
	if (fleshing_out_CF)
	{
		q1npts = new_nqpts;
		q2npts = new_nqpts;
		q3npts = new_nqpts;
		q1pts = qx_fleshed_out_pts;
		q2pts = qy_fleshed_out_pts;
		q3pts = qz_fleshed_out_pts;
	}
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
		Correlfun3D_data.q_o[idx] = q1pts[i] * ckp + q2pts[j] * skp;
		Correlfun3D_data.q_s[idx] = -q1pts[i] * skp + q2pts[j] * ckp;
		Correlfun3D_data.q_l[idx] = q3pts[k];
		Correlfun3D_data.y[idx] = Correl_3D[i][j][k];
		//Correlfun3D_data.sigma[idx] = Correl_3D_err[i][j][k];
		Correlfun3D_data.sigma[idx] = 1.e-3;
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
	/*cout << endl << "Best fit results:" << endl;
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
	cout << "--------------------------------------------------------------------" << endl;*/

	double chi = gsl_blas_dnrm2(solver_ptr->f);
	double dof = data_length - n_para;
	double c = GSL_MAX_DBL(1, chi/sqrt(dof));

	lambda_Correl[ipt][ipphi] = get_fit_results(0, solver_ptr);
	lambda_Correl_err[ipt][ipphi] = c*get_fit_err(0, covariance_ptr);
	R2_out[ipt][ipphi] = fabs(get_fit_results(1, solver_ptr))*hbarC*hbarC;
	R2_side[ipt][ipphi] = fabs(get_fit_results(2, solver_ptr))*hbarC*hbarC;
	R2_long[ipt][ipphi] = fabs(get_fit_results(3, solver_ptr))*hbarC*hbarC;
	R2_outside[ipt][ipphi] = get_fit_results(4, solver_ptr)*hbarC*hbarC;
	R2_out_err[ipt][ipphi] = c*get_fit_err(1, covariance_ptr)*hbarC*hbarC;
	R2_side_err[ipt][ipphi] = c*get_fit_err(2, covariance_ptr)*hbarC*hbarC;
	R2_long_err[ipt][ipphi] = c*get_fit_err(3, covariance_ptr)*hbarC*hbarC;
	R2_outside_err[ipt][ipphi] = c*get_fit_err(4, covariance_ptr)*hbarC*hbarC;

	cout << "final results: " << endl;
	cout << scientific << setw(10) << setprecision(5) 
		<< "chisq/dof = " << chi*chi/dof << endl;
	cout << scientific << setw(10) << setprecision(5) 
		<< " lambda = " << lambda_Correl[ipt][ipphi] << " +/- " << lambda_Correl_err[ipt][ipphi] << endl;
	cout << " R2_out[ipt=" << ipt << "][ipphi=" << ipphi << "] = " << R2_out[ipt][ipphi] << " +/- " << R2_out_err[ipt][ipphi] << endl;
	cout << " R2_side[ipt=" << ipt << "][ipphi=" << ipphi << "] = " << R2_side[ipt][ipphi] << " +/- " << R2_side_err[ipt][ipphi] << endl;
	cout << " R2_long[ipt=" << ipt << "][ipphi=" << ipphi << "] = " << R2_long[ipt][ipphi] << " +/- " << R2_long_err[ipt][ipphi] << endl;
	cout << " R2_outside[ipt=" << ipt << "][ipphi=" << ipphi << "] = " << R2_outside[ipt][ipphi] << " +/- " << R2_outside_err[ipt][ipphi] << endl;

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
		double Yi = 1.0 + exp(- q_l[i]*q_l[i]*R2_l - q_s[i]*q_s[i]*R2_s - q_o[i]*q_o[i]*R2_o - 2.*q_o[i]*q_s[i]*R2_os);
		gsl_vector_set (f_ptr, i, (Yi - y[i]) / sigma[i]);
cout << "i = " << i << ": " << y[i] << "   " << Yi << "   " << (Yi - y[i]) / sigma[i]
		<< "   " << q_o[i] << "   " << q_s[i] << "   " << q_l[i] << "   " << R2_o << "   " << R2_s << "   " << R2_l << "   " << R2_os << endl;
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
		double Yi = 1.0 + lambda*exp(- q_l[i]*q_l[i]*R2_l - q_s[i]*q_s[i]*R2_s - q_o[i]*q_o[i]*R2_o - 2.*q_o[i]*q_s[i]*R2_os);
		gsl_vector_set (f_ptr, i, (Yi - y[i]) / sigma[i]);
cout << "i = " << i << ": " << y[i] << "   " << Yi << "   " << (Yi - y[i]) / sigma[i] << "   " << lambda
		<< "   " << q_o[i] << "   " << q_s[i] << "   " << q_l[i] << "   " << R2_o << "   " << R2_s << "   " << R2_l << "   " << R2_os << endl;
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

/*
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
*/

/************************************************************************/

//End of file
