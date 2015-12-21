#include<iostream>
#include<sstream>
#include<string>
#include<fstream>
#include<cmath>
#include<iomanip>
#include<vector>
#include<stdio.h>
#include<time.h>

#include "H5Cpp.h"

#include "CFWR.h"

using namespace std;

const int SPA_RANK = 2;
const int giant_FOslice_array_size = eta_s_npts * qnpts * qnpts * qnpts * qnpts * ntrig;
const int chunk_size = n_interp_pT_pts * n_interp_pphi_pts * qnpts * qnpts * qnpts * qnpts * ntrig;
const int small_array_size = qnpts * qnpts * qnpts * qnpts * ntrig;

int CorrelationFunction::Set_giant_HDF_array()
{
	ostringstream filename_stream_ga;
	filename_stream_ga << global_path << "/giant_array.h5";
	const H5std_string FILE_NAME(filename_stream_ga.str().c_str());
	const H5std_string DATASET_NAME("ga");

	try
	{
		Exception::dontPrint();
	
		file = new H5::H5File(FILE_NAME, H5F_ACC_TRUNC);
	
		hsize_t dims[RANK] = {FO_length * giant_FOslice_array_size};               // dataset dimensions
		dataspace = new H5::DataSpace (RANK, dims);
	
		dataset = new H5::DataSet( file->createDataSet(DATASET_NAME, PredType::NATIVE_DOUBLE, *dataspace) );

		hsize_t dimsm[RANK] = {giant_FOslice_array_size};
		hsize_t offset[RANK] = {0};
		hsize_t count[RANK] = {giant_FOslice_array_size};
	
		memspace = new H5::DataSpace (RANK, dimsm, NULL);

		double * tmp_results_ii0 = new double [ntrig];
		double * tmp_results_ii1 = new double [ntrig];
		for (int isurf = 0; isurf < FO_length; ++isurf)
		{
			int igFOsa = 0;
			offset[0] = isurf;
			dataspace->selectHyperslab(H5S_SELECT_SET, count, offset);
			double * giant_FOslice_array = new double [giant_FOslice_array_size];
			for (int ieta = 0; ieta < eta_s_npts; ++ieta)
			for (int iqt = 0; iqt < qnpts; ++iqt)
			for (int iqx = 0; iqx < qnpts; ++iqx)
			for (int iqy = 0; iqy < qnpts; ++iqy)
			for (int iqz = 0; iqz < qnpts; ++iqz)
			{
				form_trig_sign_z(isurf, ieta, iqt, iqx, iqy, iqz, 0, tmp_results_ii0);
				form_trig_sign_z(isurf, ieta, iqt, iqx, iqy, iqz, 1, tmp_results_ii1);
				giant_FOslice_array[igFOsa] = tmp_results_ii0[0] + tmp_results_ii1[0];
				giant_FOslice_array[igFOsa+1] = tmp_results_ii0[1] + tmp_results_ii1[1];
				igFOsa += 2;
			}
			dataset->write(giant_FOslice_array, PredType::NATIVE_DOUBLE, *memspace, *dataspace);
			delete giant_FOslice_array;
		}
	
		delete [] tmp_results_ii0;
		delete [] tmp_results_ii1;

		memspace->close();
		dataset->close();
		file->close();
		delete memspace;
		delete file;
		delete dataset;
		// Now set-up for remainder of calculations
		file = new H5::H5File(FILE_NAME, H5F_ACC_RDONLY);
		dataset = new H5::DataSet( file->openDataSet( DATASET_NAME ) );
		dimsm[0] = small_array_size;
		memspace = new H5::DataSpace (RANK, dimsm, NULL);
    }

    catch(FileIException error)
    {
		error.printError();
		return -1;
    }

    catch(H5::DataSetIException error)
    {
		error.printError();
		return -1;
    }

    catch(H5::DataSpaceIException error)
    {
		error.printError();
		return -1;
    }

	return (0);
}

int CorrelationFunction::Get_small_array_from_giant_HDF_array(int isurf, int ieta, double * small_array)
{
	try
	{
		Exception::dontPrint();

		hsize_t offset[RANK] = {isurf * eta_s_npts + ieta};
		hsize_t count[RANK] = {small_array_size};

		dataspace->selectHyperslab(H5S_SELECT_SET, count, offset);

		dataset->read(small_array, PredType::NATIVE_DOUBLE, *memspace, *dataspace);
    }

    catch(FileIException error)
    {
		error.printError();
		return -1;
    }

    catch(H5::DataSetIException error)
    {
		error.printError();
		return -1;
    }

    catch(H5::DataSpaceIException error)
    {
		error.printError();
		return -1;
    }

	return (0);
}

int CorrelationFunction::Clean_up_HDF_miscellany()
{
	try
	{
		Exception::dontPrint();

		delete memspace;
		delete file;
		delete dataspace;
		delete dataset;
    }

    catch(FileIException error)
    {
		error.printError();
		return -1;
    }

    catch(H5::DataSetIException error)
    {
		error.printError();
		return -1;
    }

    catch(H5::DataSpaceIException error)
    {
		error.printError();
		return -1;
    }

	return (0);
}

//*************************************************************
// Chunked versions here
//*************************************************************

int CorrelationFunction::Set_giant_chunked_HDF_array()
{
	double local_small_array[giant_FOslice_array_size];

	ostringstream filename_stream_ga;
	filename_stream_ga << global_path << "/giant_array.h5";
	const H5std_string FILE_NAME(filename_stream_ga.str().c_str());
	const H5std_string DATASET_NAME("ga");

	try
    {
		Exception::dontPrint();
	
		file = new H5::H5File(FILE_NAME, H5F_ACC_TRUNC);
	
		DSetCreatPropList cparms;
		hsize_t chunk_dims[RANK] = {giant_FOslice_array_size};
		cparms.setChunk( RANK, chunk_dims );
	
		hsize_t dims[RANK] = {FO_length * giant_FOslice_array_size};	// == FO_length * giant_FOslice_array_size
		dataspace = new H5::DataSpace (RANK, dims);
	
		dataset = new H5::DataSet( file->createDataSet(DATASET_NAME, PredType::NATIVE_DOUBLE, *dataspace, cparms) );

		hsize_t offset[RANK];
		hsize_t count[RANK] = {giant_FOslice_array_size};
		hsize_t dimsm[RANK] = {giant_FOslice_array_size};

		memspace = new H5::DataSpace (RANK, dimsm, NULL);

		double * tmp_results_ii0 = new double [ntrig];
		double * tmp_results_ii1 = new double [ntrig];
		for (int isurf = 0; isurf < FO_length; ++isurf)
		{
			int igFOsa = 0;
			offset[0] = isurf * giant_FOslice_array_size;
			dataspace->selectHyperslab(H5S_SELECT_SET, count, offset);
			for (int ieta = 0; ieta < eta_s_npts; ++ieta)
			for (int iqt = 0; iqt < qnpts; ++iqt)
			for (int iqx = 0; iqx < qnpts; ++iqx)
			for (int iqy = 0; iqy < qnpts; ++iqy)
			for (int iqz = 0; iqz < qnpts; ++iqz)
			{
				form_trig_sign_z(isurf, ieta, iqt, iqx, iqy, iqz, 0, tmp_results_ii0);
				form_trig_sign_z(isurf, ieta, iqt, iqx, iqy, iqz, 1, tmp_results_ii1);
				local_small_array[igFOsa] = tmp_results_ii0[0] + tmp_results_ii1[0];
				local_small_array[igFOsa+1] = tmp_results_ii0[1] + tmp_results_ii1[1];
				igFOsa += 2;
			}
			dataset->write(local_small_array, PredType::NATIVE_DOUBLE, *memspace, *dataspace);
		}
	
		delete [] tmp_results_ii0;
		delete [] tmp_results_ii1;

		memspace->close();
		dataset->close();
		file->close();
		delete memspace;
		delete file;
		delete dataset;
		file = new H5::H5File(FILE_NAME, H5F_ACC_RDONLY);
		dataset = new H5::DataSet( file->openDataSet( DATASET_NAME ) );
		memspace = new H5::DataSpace (RANK, dimsm, NULL);
    }

    catch(FileIException error)
    {
		error.printError();
		cerr << "FileIException error!" << endl;
		return -1;
    }

    catch(H5::DataSetIException error)
    {
		error.printError();
		cerr << "H5::DataSetIException error!" << endl;
		return -2;
    }

    catch(H5::DataSpaceIException error)
    {
		error.printError();
		cerr << "H5::DataSpaceIException error!" << endl;
		return -3;
    }

	return (0);
}

int CorrelationFunction::Get_chunk(int isurf, double * local_small_array)
{
	try
	{
		Exception::dontPrint();

		hsize_t offset[RANK] = {isurf * giant_FOslice_array_size};
		hsize_t count[RANK] = {giant_FOslice_array_size};				// == chunk_dims

		dataspace->selectHyperslab(H5S_SELECT_SET, count, offset);

		dataset->read(local_small_array, PredType::NATIVE_DOUBLE, *memspace, *dataspace);
    }

    catch(FileIException error)
    {
		error.printError();
		cerr << "FileIException error!" << endl;
		return -1;
    }

    catch(H5::DataSetIException error)
    {
		error.printError();
		cerr << "H5::DataSetIException error!" << endl;
		return -2;
    }

    catch(H5::DataSpaceIException error)
    {
		error.printError();
		cerr << "H5::DataSpaceIException error!" << endl;
		return -3;
    }

	return (0);
}

int CorrelationFunction::Initialize_resonance_HDF_array()
{
	double * resonance_chunk = new double [chunk_size];

	ostringstream filename_stream_ra;
	filename_stream_ra << global_path << "/resonance_spectra.h5";
	H5std_string RESONANCE_FILE_NAME(filename_stream_ra.str().c_str());
	H5std_string RESONANCE_DATASET_NAME("ra");

	try
    {
		Exception::dontPrint();
	
		resonance_file = new H5::H5File(RESONANCE_FILE_NAME, H5F_ACC_TRUNC);

		DSetCreatPropList cparms;
		hsize_t chunk_dims[RANK] = {chunk_size};
		cparms.setChunk( RANK, chunk_dims );

		hsize_t dims[RANK] = {Nparticle * chunk_size};
		resonance_dataspace = new H5::DataSpace (RANK, dims);

		resonance_dataset = new H5::DataSet( resonance_file->createDataSet(RESONANCE_DATASET_NAME, PredType::NATIVE_DOUBLE, *resonance_dataspace, cparms) );

		hsize_t count[RANK] = {chunk_size};
		hsize_t dimsm[RANK] = {chunk_size};
		hsize_t offset[RANK] = {0};

		resonance_memspace = new H5::DataSpace (RANK, dimsm, NULL);
		for (int ir = 0; ir < Nparticle; ++ir)
		{
			offset[0] = ir * chunk_size;
			resonance_dataspace->selectHyperslab(H5S_SELECT_SET, count, offset);

			for (int iidx = 0; iidx < chunk_size; ++iidx)
				resonance_chunk[iidx] = 0.0;

			//initialize everything with zeros
			resonance_dataset->write(resonance_chunk, PredType::NATIVE_DOUBLE, *resonance_memspace, *resonance_dataspace);
		}
		resonance_memspace->close();
		resonance_dataset->close();
		resonance_file->close();
		delete resonance_memspace;
		delete resonance_file;
		delete resonance_dataset;
		resonance_file = new H5::H5File(RESONANCE_FILE_NAME, H5F_ACC_RDWR);
		resonance_dataset = new H5::DataSet( resonance_file->openDataSet( RESONANCE_DATASET_NAME ) );
		resonance_memspace = new H5::DataSpace (RANK, dimsm, NULL);
    }

    catch(FileIException error)
    {
		error.printError();
		cerr << "FileIException error!" << endl;
		return -1;
    }

    catch(H5::DataSetIException error)
    {
		error.printError();
		cerr << "H5::DataSetIException error!" << endl;
		return -2;
    }

    catch(H5::DataSpaceIException error)
    {
		error.printError();
		cerr << "H5::DataSpaceIException error!" << endl;
		return -3;
    }

	delete [] resonance_chunk;

	return (0);
}

int CorrelationFunction::Get_resonance_from_HDF_array(int local_pid, double ******* resonance_array_to_fill)
{
	double * resonance_chunk = new double [chunk_size];

	ostringstream filename_stream_ra;
	filename_stream_ra << global_path << "/resonance_spectra.h5";
	H5std_string RESONANCE_FILE_NAME(filename_stream_ra.str().c_str());
	H5std_string RESONANCE_DATASET_NAME("ra");

	try
    {
		Exception::dontPrint();
	
		hsize_t offset[RANK] = {local_pid * chunk_size};
		hsize_t count[RANK] = {chunk_size};				// == chunk_dims
		resonance_dataspace->selectHyperslab(H5S_SELECT_SET, count, offset);

		resonance_dataset->read(resonance_chunk, PredType::NATIVE_DOUBLE, *resonance_memspace, *resonance_dataspace);

		// use loaded chunk to fill resonance_array_to_fill
		int iidx = 0;
		for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
		for (int iphi = 0; iphi < n_interp_pphi_pts; ++iphi)
		for (int iqt = 0; iqt < qnpts; ++iqt)
		for (int iqx = 0; iqx < qnpts; ++iqx)
		for (int iqy = 0; iqy < qnpts; ++iqy)
		for (int iqz = 0; iqz < qnpts; ++iqz)
		for (int itrig = 0; itrig < ntrig; ++itrig)
		{
			double temp = resonance_chunk[iidx];
			resonance_array_to_fill[ipt][iphi][iqt][iqx][iqy][iqz][itrig] = temp;
			++iidx;
		}

   }

    catch(FileIException error)
    {
		error.printError();
		cerr << "FileIException error!" << endl;
		return -1;
    }

    catch(H5::DataSetIException error)
    {
		error.printError();
		cerr << "H5::DataSetIException error!" << endl;
		return -2;
    }

    catch(H5::DataSpaceIException error)
    {
		error.printError();
		cerr << "H5::DataSpaceIException error!" << endl;
		return -3;
    }

	delete [] resonance_chunk;

	return (0);
}

int CorrelationFunction::Set_resonance_in_HDF_array(int local_pid, double ******* resonance_array_to_use)
{
	double * resonance_chunk = new double [chunk_size];

	ostringstream filename_stream_ra;
	filename_stream_ra << global_path << "/resonance_spectra.h5";
	H5std_string RESONANCE_FILE_NAME(filename_stream_ra.str().c_str());
	H5std_string RESONANCE_DATASET_NAME("ra");

	try
    {
		Exception::dontPrint();
	
		hsize_t offset[RANK] = {local_pid * chunk_size};
		hsize_t count[RANK] = {chunk_size};				// == chunk_dims
		resonance_dataspace->selectHyperslab(H5S_SELECT_SET, count, offset);

		// use loaded chunk to fill resonance_array_to_fill
		int iidx = 0;
		for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
		for (int iphi = 0; iphi < n_interp_pphi_pts; ++iphi)
		for (int iqt = 0; iqt < qnpts; ++iqt)
		for (int iqx = 0; iqx < qnpts; ++iqx)
		for (int iqy = 0; iqy < qnpts; ++iqy)
		for (int iqz = 0; iqz < qnpts; ++iqz)
		for (int itrig = 0; itrig < ntrig; ++itrig)
		{
			double temp = resonance_array_to_use[ipt][iphi][iqt][iqx][iqy][iqz][itrig];
			resonance_chunk[iidx] = temp;
			++iidx;
		}

		resonance_dataset->write(resonance_chunk, PredType::NATIVE_DOUBLE, *resonance_memspace, *resonance_dataspace);
   }

    catch(FileIException error)
    {
		error.printError();
		cerr << "FileIException error!" << endl;
		return -1;
    }

    catch(H5::DataSetIException error)
    {
		error.printError();
		cerr << "H5::DataSetIException error!" << endl;
		return -2;
    }

    catch(H5::DataSpaceIException error)
    {
		error.printError();
		cerr << "H5::DataSpaceIException error!" << endl;
		return -3;
    }

	delete [] resonance_chunk;

	return (0);
}

int CorrelationFunction::Copy_chunk(int current_resonance_index, int reso_idx_to_be_copied)
{
	double * resonance_chunk = new double [chunk_size];

	ostringstream filename_stream_ra;
	filename_stream_ra << global_path << "/resonance_spectra.h5";
	H5std_string RESONANCE_FILE_NAME(filename_stream_ra.str().c_str());
	H5std_string RESONANCE_DATASET_NAME("ra");

	try
    {
		Exception::dontPrint();
	
		hsize_t offset[RANK] = {reso_idx_to_be_copied * chunk_size};
		hsize_t count[RANK] = {chunk_size};
		resonance_dataspace->selectHyperslab(H5S_SELECT_SET, count, offset);

		// load resonance to be copied first
		resonance_dataset->read(resonance_chunk, PredType::NATIVE_DOUBLE, *resonance_memspace, *resonance_dataspace);

		// now set this to current resonance
		offset[0] = current_resonance_index * chunk_size;
		resonance_dataspace->selectHyperslab(H5S_SELECT_SET, count, offset);
		resonance_dataset->write(resonance_chunk, PredType::NATIVE_DOUBLE, *resonance_memspace, *resonance_dataspace);
   }

    catch(FileIException error)
    {
		error.printError();
		cerr << "FileIException error!" << endl;
		return -1;
    }

    catch(H5::DataSetIException error)
    {
		error.printError();
		cerr << "H5::DataSetIException error!" << endl;
		return -2;
    }

    catch(H5::DataSpaceIException error)
    {
		error.printError();
		cerr << "H5::DataSpaceIException error!" << endl;
		return -3;
    }

	delete [] resonance_chunk;

	return (0);
}

int CorrelationFunction::Create_S_p_withweight_HDFarray()
{
	const int full_FO_length = FO_length * eta_s_npts;
	const int momentum_grid_size = n_interp_pT_pts * n_interp_pphi_pts;

	ostringstream filename_stream_ga;
	filename_stream_ga << global_path << "/S_p_withweight_array.h5";
	const H5std_string SPA_FILE_NAME(filename_stream_ga.str().c_str());
	const H5std_string SPA_DATASET_NAME("spwwa");

	double local_array [1][full_FO_length];
	for (int ifFO = 0; ifFO < full_FO_length; ++ifFO)
	{
		local_array[0][ifFO] = 0.0;
	}

	try
    {
		Exception::dontPrint();
	
		SPAfile = new H5::H5File(SPA_FILE_NAME, H5F_ACC_TRUNC);
	
		DSetCreatPropList cparms;
		hsize_t chunk_dims[SPA_RANK] = {1, full_FO_length};
		cparms.setChunk( SPA_RANK, chunk_dims );
	
		hsize_t dims[SPA_RANK] = {momentum_grid_size, full_FO_length};
		SPAdataspace = new H5::DataSpace (SPA_RANK, dims);
	
		SPAdataset = new H5::DataSet( SPAfile->createDataSet(SPA_DATASET_NAME, PredType::NATIVE_DOUBLE, *SPAdataspace, cparms) );

		hsize_t offset[SPA_RANK] = {0, 0};
		hsize_t count[SPA_RANK] = {1, full_FO_length};
		hsize_t dimsm[SPA_RANK] = {1, full_FO_length};

		SPAmemspace = new H5::DataSpace (SPA_RANK, dimsm, NULL);

		for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
		for (int iphi = 0; iphi < n_interp_pphi_pts; ++iphi)
		{
			offset[0] = ipt * n_interp_pphi_pts + iphi;

			// select appropriate slice
			SPAdataspace->selectHyperslab(H5S_SELECT_SET, count, offset);

			// just initialize with zeros
			SPAdataset->write(local_array, PredType::NATIVE_DOUBLE, *SPAmemspace, *SPAdataspace);
		}
		
		SPAmemspace->close();
		SPAdataset->close();
		SPAfile->close();
		delete SPAmemspace;
		delete SPAfile;
		delete SPAdataset;
		SPAfile = new H5::H5File(SPA_FILE_NAME, H5F_ACC_RDONLY);
		SPAdataset = new H5::DataSet( SPAfile->openDataSet( SPA_DATASET_NAME ) );
		SPAmemspace = new H5::DataSpace (SPA_RANK, dimsm, NULL);
    }

    catch(FileIException error)
    {
		error.printError();
		cerr << "FileIException error!" << endl;
		return -1;
    }

    catch(H5::DataSetIException error)
    {
		error.printError();
		cerr << "DataSetIException error!" << endl;
		return -2;
    }

    catch(H5::DataSpaceIException error)
    {
		error.printError();
		cerr << "DataSpaceIException error!" << endl;
		return -3;
    }

	return (0);
}

int CorrelationFunction::Get_S_p_withweight_chunk(int ipt, int iphi, double * local_array)
{
	const int full_FO_length = FO_length * eta_s_npts;

	try
	{
		Exception::dontPrint();

		hsize_t offset[SPA_RANK] = {ipt * n_interp_pphi_pts + iphi, 0};
		hsize_t count[SPA_RANK] = {1, full_FO_length};				// == chunk_dims

		SPAdataspace->selectHyperslab(H5S_SELECT_SET, count, offset);

		SPAdataset->read(local_array, PredType::NATIVE_DOUBLE, *SPAmemspace, *SPAdataspace);
    }

    catch(FileIException error)
    {
		error.printError();
		cerr << "FileIException error!" << endl;
		return -1;
    }

    catch(H5::DataSetIException error)
    {
		error.printError();
		cerr << "DataSetIException error!" << endl;
		return -2;
    }

    catch(H5::DataSpaceIException error)
    {
		error.printError();
		cerr << "DataSpaceIException error!" << endl;
		return -3;
    }

	return (0);
}

int CorrelationFunction::Set_S_p_withweight_chunk(int ipt, int iphi, double * local_array)
{
	try
	{
		Exception::dontPrint();

		hsize_t offset[SPA_RANK] = {ipt * n_interp_pphi_pts + iphi, 0};
		hsize_t count[SPA_RANK] = {1, full_FO_length};				// == chunk_dims

		SPAdataspace->selectHyperslab(H5S_SELECT_SET, count, offset);

		SPAdataset->write(local_array, PredType::NATIVE_DOUBLE, *SPAmemspace, *SPAdataspace);
    }

    catch(FileIException error)
    {
		error.printError();
		cerr << "FileIException error!" << endl;
		return -1;
    }

    catch(H5::DataSetIException error)
    {
		error.printError();
		cerr << "DataSetIException error!" << endl;
		return -2;
    }

    catch(H5::DataSpaceIException error)
    {
		error.printError();
		cerr << "DataSpaceIException error!" << endl;
		return -3;
    }

	return (0);
}

//End of file
