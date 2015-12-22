#ifndef CFWR_H
#define CFWR_H

#include<iostream>
#include<sstream>
#include<fstream>
#include<cmath>
#include<iomanip>
#include<string>
#include<fstream>
#include<vector>
#include<set>
#include<queue>

#include<gsl/gsl_sf_bessel.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>// gsl random number generators
#include <gsl/gsl_randist.h>  // gsl random number distributions
#include <gsl/gsl_vector.h>   // gsl vector and matrix definitions
#include <gsl/gsl_blas.h>     // gsl linear algebra stuff
#include <gsl/gsl_multifit_nlin.h>  // gsl multidimensional fitting

#include "H5Cpp.h"

#include "readindata.h"
#include "parameters.h"
#include "Arsenal.h"
#include "gauss_quadrature.h"

using namespace std;

typedef struct
{
   int resonance_particle_id;		// keeps track of current resonance's index in all_particles array
   int resonance_idx;			// keeps track of current resonance's index in chosen_resonances vector
   int nbody;
   int resonance_sign;
   double resonance_mass;
   double resonance_mu;
   double resonance_gspin;
   double resonance_Gamma;
   double resonance_total_br;
   double resonance_direct_br;
   double * resonance_decay_masses;
   int * resonance_decay_monvals;
   double * resonance_decay_Gammas;
   string resonance_name;
   bool include_channel;
}decay_info;

struct Correlationfunction1D_data
{
  size_t data_length;
  double * q;
  double * y;
  double * sigma;
};

int Fittarget_correlfun1D_f (const gsl_vector *xvec_ptr, void *params_ptr, gsl_vector *f_ptr);
int Fittarget_correlfun1D_df (const gsl_vector *xvec_ptr, void *params_ptr,  gsl_matrix *Jacobian_ptr);
int Fittarget_correlfun1D_fdf (const gsl_vector* xvec_ptr, void *params_ptr, gsl_vector* f_ptr, gsl_matrix* Jacobian_ptr);


class CorrelationFunction
{
	private:
		//particle information 
		string particle_name;
		double particle_mass;
		int particle_monval;
		int particle_id;     //particle id
		double particle_sign;   //+/- 1 for Fermi/Bose statistics for baryon/meson
		double particle_gspin;  //particle degeneracy 
		double particle_mu;
		particle_info * all_particles;
		vector<int> chosen_resonances;
		bool thermal_pions_only;
		int Nparticle;
		int target_particle_id;		//the particle whose spectra (with resonance contributions) you want to compute
		int current_level_of_output;
		int qspace_cs_slice_length;
		int full_FO_length;
		int nFO_cutoff;

		int n_zeta_pts, n_v_pts, n_s_pts;
		double v_min, v_max, zeta_min, zeta_max, s_min, s_max;
		
		//integrated space-time moments vars.
		double S_r_to_part, S_xs_r_to_part, S_xo_r_to_part, S_xl_r_to_part, S_t_r_to_part;
		double S_xs2_r_to_part, S_xo2_r_to_part, S_xl2_r_to_part, S_t2_r_to_part;
		double S_xs_xo_r_to_part, S_xs_xl_r_to_part, S_xs_t_r_to_part, S_xo_xl_r_to_part, S_xo_t_r_to_part, S_xl_t_r_to_part;
		
		//array to hold previous and current resonance info
		decay_info * decay_channels;
		int current_decay_channel_idx, current_resonance_particle_id, previous_resonance_idx, current_resonance_idx, current_reso_nbody;
		double current_resonance_mu, current_resonance_mass, current_resonance_Gamma, current_m2_Gamma, current_m3_Gamma;
		double current_resonance_total_br, current_resonance_direct_br, current_daughter_mass, current_daughter_Gamma;
		double * current_resonance_decay_masses, * P_eval, * alpha_mu;
		int previous_decay_channel_idx, previous_resonance_particle_id, previous_reso_nbody;
		double previous_resonance_mu, previous_resonance_mass, previous_resonance_Gamma, previous_m2_Gamma, previous_m3_Gamma;
		double previous_resonance_total_br, previous_resonance_direct_br, previous_daughter_mass, previous_daughter_Gamma;
		double * previous_resonance_decay_masses;
		
		//arrays to hold results of resonance phase-space integrations
		double ******* current_dN_dypTdpTdphi_moments;
		double ******* current_ln_dN_dypTdpTdphi_moments;
		double ******* current_sign_of_dN_dypTdpTdphi_moments;
		double ******** current_daughters_dN_dypTdpTdphi_moments;
		double ******** current_daughters_ln_dN_dypTdpTdphi_moments;
		double ******** current_daughters_sign_of_dN_dypTdpTdphi_moments;

		// needed these to avoid too many trigonometric evaluations
		double **** osc0, *** osc1, *** osc2, **** osc3;
	
		//needed for resonance calculations
		//kinematic info
		double pstar, Estar, Yp, Ym, DeltaY;
		double MTbar, DeltaMT, MTp, MTm;
		double Qfunc;
		//pair momentum info, currently assumes pT != 0
		double p_y, pT, pphi, mT, mass;
		double ch_p_y, sh_p_y;
		//resonance momentum info
		double P_Y, PT, PPhi, MT, Mres, PPhip, PPhim;
		double m2, m3, Gamma, br, m2Gamma, m3Gamma;
		double * Pp;
		double * Pm;

		//SP momentum arrays for interpolation grid
		double * SPinterp_pT;
		double * SPinterp_pphi;
		double * sin_SPinterp_pphi;
		double * cos_SPinterp_pphi;
		double ** SPinterp_p0, ** SPinterp_p0_Transpose;
		double ** SPinterp_pz, ** SPinterp_pz_Transpose;

		FO_surf* current_FOsurf_ptr;
		//FO surface info that is constant - to save time
		double Tdec, Edec, Pdec, muRES, signRES, gRES;
		double S_prefactor;
	
		//single particle spectra for plane angle determination
		//int n_order;
		double SP_p_y, global_plane_psi, mean_pT;
		double * SP_pT, * SP_pphi, * SP_pT_weight, * SP_pphi_weight, * dN_dypTdpT;
		double * dN_dydphi, * pTdN_dydphi, * plane_angle;
		double ** dN_dypTdpTdphi, ** cosine_iorder, ** sine_iorder;
		size_t *** most_important_FOcells;
		double * giant_array, ** giant_array_slice;
		int ** number_of_FOcells_above_cutoff_array;

		//pair momentum
		double K_y, ch_K_y, sh_K_y;
		double current_K_phi, cos_cKphi, sin_cKphi;
		double beta_perp, beta_l;
		double * K_T, * K_phi, * K_phi_weight;
		    
		//spatial rapidity grid
		double * eta_s, * ch_eta_s, * sh_eta_s, * eta_s_weight;

		//points and weights for resonance integrals
		double * zeta_pts, * v_pts, ** s_pts, * NEW_s_pts;
		double * zeta_wts, * v_wts, ** s_wts, * NEW_s_wts;

		//some arrays to save unnecessary multiple calculations for resonances
		//use these for n_body = 2
		double VEC_n2_spt, VEC_n2_pstar, VEC_n2_Estar, VEC_n2_psBmT, VEC_n2_DeltaY, VEC_n2_Yp, VEC_n2_Ym;
		double * VEC_n2_P_Y, * VEC_n2_MTbar, * VEC_n2_DeltaMT, * VEC_n2_MTp, * VEC_n2_MTm;
		double ** VEC_n2_MT, ** VEC_n2_PPhi_tilde, ** VEC_n2_PPhi_tildeFLIP, ** VEC_n2_PT;
		double VEC_n2_s_factor;
		double * VEC_n2_v_factor;
		double ** VEC_n2_zeta_factor;
		double VEC_n2_g_s;
		//use these for n_body = 3
		double * VEC_pstar, * VEC_Estar, * VEC_DeltaY, * VEC_Yp, * VEC_Ym, * VEC_s_factor, * VEC_g_s;
		double ** VEC_P_Y, ** VEC_MTbar, ** VEC_DeltaMT, ** VEC_MTp, ** VEC_MTm, ** VEC_v_factor;
		double *** VEC_MT, *** VEC_PPhi_tilde, *** VEC_PPhi_tildeFLIP, *** VEC_PT, *** VEC_zeta_factor;
        	double * ssum_vec, * vsum_vec, * zetasum_vec, * Csum_vec;
		
		//Emission function
		int FO_length;
		int Emissionfunction_length;
		
		double *** spectra, *** abs_spectra;
		double **** CFvals;
		
		double * qo_pts, * qs_pts, * ql_pts, * q_pts, * q_axes, * qt_pts, * qx_pts, * qy_pts, * qz_pts;
		double * q_out, * q_side, * q_long;
		
		//store correlation functions
		double * Correl_1D_out;
		double * Correl_1D_out_err;
		double * Correl_1D_side;
		double * Correl_1D_side_err;
		double * Correl_1D_long;
		double * Correl_1D_long_err;
		double ** lambda_Correl, ** lambda_Correl_err;

		//HBT radii coefficients
		double **R2_side, **R2_out, **R2_long, **R2_outside, **R2_sidelong, **R2_outlong;
		double **R2_side_err, **R2_out_err, **R2_long_err, **R2_outside_err, **R2_sidelong_err, **R2_outlong_err;
		double **R2_side_C, **R2_side_S;
		double **R2_out_C, **R2_out_S;
		double **R2_long_C, **R2_long_S;
		double **R2_outside_C, **R2_outside_S;
		double **R2_outlong_C, **R2_outlong_S;
		double **R2_sidelong_C, **R2_sidelong_S;

		double *** res_sign_info, *** res_log_info, *** res_moments_info;
		double ***** local_temp_moments, ******* temp_moments_array;

		double **** S_p_withweight_array;
		//double ** S_pt_pphi_withweight_array;
		bool **** zero_FOcell_flag;
		
		//miscellaneous
		ofstream * global_out_stream_ptr;
		int global_folderindex;
		string global_path;
		string global_runfolder;
		string global_resultsfolder_stem;
		string no_df_stem;
		int n_resonance, n_decay_channels;
		int n_body;
		set<int> daughter_resonance_indices;

		//some private methods		
		bool particles_are_the_same(int idx1, int idx2);
		bool Search_for_similar_particle(int dc_idx, int * result);

	public:
		void Determine_plane_angle(FO_surf* FOsurf_ptr, int dc_idx);
		void Compute_correlation_function(FO_surf* FOsurf_ptr);
		void Update_sourcefunction(particle_info* particle, int FOarray_length, int particle_idx);
		bool fexists(const char *filename);

		// HDF routines
		int Set_giant_HDF_array();
		int Get_small_array_from_giant_HDF_array(int isurf, int ieta, double * small_array);
		int Set_giant_chunked_HDF_array();
		int Get_chunk(int isurf, double * small_array);
		int Clean_up_HDF_miscellany();
		int Get_resonance_from_HDF_array(int local_pid, double ******* resonance_array_to_fill);
		int Set_resonance_in_HDF_array(int local_pid, double ******* resonance_array_to_use);
		int Initialize_resonance_HDF_array();
		int Copy_chunk(int current_resonance_index, int reso_idx_to_be_copied);
		int Get_2D_resonance_from_HDF_array(int local_pid, double ******* resonance_array_to_fill);
		int Set_2D_resonance_in_HDF_array(int local_pid, double ******* resonance_array_to_use);
		int Initialize_2D_resonance_HDF_array();
		int Copy_2D_chunk(int current_resonance_index, int reso_idx_to_be_copied);

		void Set_giant_array_slice(int iqt, int iqx, int iqy, int iqz);
		void addElementToQueue(priority_queue<pair<double, size_t> >& p, pair<double, size_t> elem, size_t max_size);

		void Set_dN_dypTdpTdphi_moments(FO_surf* FOsurf_ptr, int dc_idx);
		void Cal_dN_dypTdpTdphi(double** SP_p0, double** SP_px, double** SP_py, double** SP_pz, FO_surf* FOsurf_ptr);
		void Cal_dN_dypTdpTdphi_heap(FO_surf* FOsurf_ptr, int local_pid);
		void Cal_dN_dypTdpTdphi_with_weights(FO_surf* FOsurf_ptr, int local_pid, double cutoff);
		double Cal_dN_dypTdpTdphi_function(FO_surf* FOsurf_ptr, int local_pid, double pT, double pphi);
		void Do_resonance_integrals(int iKT, int iKphi, int dc_idx);
		void Flatten_dN_dypTdpTdphi_moments();
		void get_rapidity_dependence(double * rap_indep_vector, double * rap_dep_vector, double rap_val);
		void combine_sourcevariances(double * output, double * input, double * alpha_vec);
		void Set_current_daughter_info(int dc_idx, int daughter_idx);
		void Set_current_particle_info(int dc_idx);
		bool Do_this_decay_channel(int dc_idx);
		bool Do_this_daughter_particle(int dc_idx, int daughter_idx, int * daughter_resonance_pid);
		void Get_spacetime_moments(FO_surf* FOsurf_ptr, int dc_idx);
		void Recycle_spacetime_moments();
		void Load_resonance_and_daughter_spectra(int local_pid);
		void Update_daughter_spectra();
		void Set_current_resonance_logs_and_signs();
		void Set_current_daughters_resonance_logs_and_signs(int n_daughters);
		void Allocate_decay_channel_info();
		void Load_decay_channel_info(int dc_idx, double K_T_local, double K_phi_local);
		void Delete_decay_channel_info();
		void Compute_source_variances(int iKT, int iKphi);
		int Set_daughter_list(int parent_resonance_index);

		void form_trig_sign_z(int isurf, int ieta, int iqt, int iqx, int iqy, int iqz, int ii, double * results);

		void Get_source_variances(int, int);
		void Calculate_R2_side(int, int);
		void Calculate_R2_out(int, int);
		void Calculate_R2_long(int, int);
		void Calculate_R2_outside(int, int);
		void Calculate_R2_sidelong(int, int);
		void Calculate_R2_outlong(int, int);
		void R2_Fourier_transform(int iKT, double plane_psi);

		//miscellaneous
		void Set_ofstream(ofstream& myout);
		void Set_path(string path);
		void Set_resultsfolder_stem(string usrdef_stem);
		void Set_runfolder(string runfolder);
		void Set_use_delta_f(bool usrdef_usedeltaf);
		void Set_particle_mass(double usrdef_particle_mass);
		void Set_current_FOsurf_ptr(FO_surf* FOsurf_ptr);
		double get_Q();
		double g(double s);
		inline void set_to_zero(double * array, size_t arraylength);
		void adaptive_simpson_integration(void (*f) (double, double *), double a, double b, double * results);
		double place_in_range(double phi, double min, double max);
		void Get_current_decay_string(int dc_idx, string * decay_string);
		int lookup_resonance_idx_from_particle_id(int particle_id);
		int list_daughters(int parent_resonance_index, set<int> * daughter_resonance_indices_ptr, particle_info * particle, int Nparticle);
		static inline double lin_int(double x_m_x1, double one_by_x2_m_x1, double f1, double f2);
		void Edndp3(double ptr, double phir, double * results);
		void Set_q_points();
		void Get_q_points(double qo, double qs, double ql, double KT, double Kphi, double * qgridpts);
		void Allocate_resonance_running_sum_vectors();
		void Delete_resonance_running_sum_vectors();
		void Zero_resonance_running_sum_vector(double * vec);
		void Setup_temp_arrays(double ***** local_temp_moments, double ******* temp_moments_array);
		void Teardown_temp_arrays(double ***** local_temp_moments, double ******* temp_moments_array);
		void Setup_current_daughters_dN_dypTdpTdphi_moments(int n_daughter);
		void Cleanup_current_daughters_dN_dypTdpTdphi_moments(int n_daughter);
		void Delete_S_p_withweight_array();
		void Reset_zero_FOcell_flag_array();

		// Gaussian fit / correlation function routines
		void Get_GF_HBTradii(FO_surf* FOsurf_ptr, int folderindex);
		void Cal_correlationfunction();
		void Fit_Correlationfunction1D(char osl_switch, int iKT, int iKphi);
		int print_fit_state_1D (size_t iteration, gsl_multifit_fdfsolver * solver_ptr);
		int Read_correlationfunction(int iKT, int iKphi);
		inline double get_fit_results(int i, gsl_multifit_fdfsolver * solver_ptr);
		inline double get_fit_err (int i, gsl_matrix * covariance_ptr);

		// input and output function prototypes
		void Output_dN_dypTdpTdphi(int folderindex);
		void Output_dN_dypTdpT(int folderindex);
		void Output_all_dN_dypTdpTdphi(int folderindex);
		void Output_total_target_dN_dypTdpTdphi(int folderindex);
		void Output_total_target_eiqx_dN_dypTdpTdphi(int folderindex);
		void Output_results(int folderindex);
		void Readin_results(int folderindex);
		void Read_in_all_dN_dypTdpTdphi(int folderindex);
		void Output_chosen_resonances();
		void Output_correlationfunction(int folderindex);

		//parameters that the user is free to define
		double plumberg_test_variable;
		bool use_delta_f;
		bool append_output;
		int n_events;
		int initial_event, currentfolderindex;
		bool read_in_all_dN_dypTdpTdphi, output_all_dN_dypTdpTdphi;
		double fraction_of_resonances;
		double * SPinterp_pT_public, * pTdep_fractions_of_resonances;

		// need to hold giant array stuff
		H5::DataSpace * dataspace, * memspace;
		H5::H5File * file;
		H5::DataSet * dataset;
		H5::DataSpace * SPAdataspace, * SPAmemspace;
		H5::H5File * SPAfile;
		H5::DataSet * SPAdataset;
		H5::DataSpace * resonance_dataspace, * resonance_memspace;
		H5::H5File * resonance_file;
		H5::DataSet * resonance_dataset;

		CorrelationFunction(particle_info* particle, particle_info* all_particles_in, int Nparticle,
				FO_surf* FOsurf_ptr, vector<int> chosen_resonances, int particle_idx, ofstream& myout);
		~CorrelationFunction();

};

#endif
