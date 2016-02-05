#include<iostream>
#include<sstream>
#include<string>
#include<fstream>
#include<cmath>
#include<iomanip>
#include<vector>
#include<stdio.h>
#include<time.h>
#include<algorithm>

#include "H5Cpp.h"

#include "CFWR.h"

using namespace std;

const int SPA_RANK = 2;
const int RANK2D = 2;
const int giant_FOslice_array_size = eta_s_npts * qtnpts * qxnpts * qynpts * qznpts * ntrig;
const int chunk_size = n_interp_pT_pts * n_interp_pphi_pts * qtnpts * qxnpts * qynpts * qznpts * ntrig;
const int small_array_size = qtnpts * qxnpts * qynpts * qznpts * ntrig;
const int q_space_size = qtnpts * qxnpts * qynpts * qznpts;

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
			for (int iqt = 0; iqt < qtnpts; ++iqt)
			for (int iqx = 0; iqx < qxnpts; ++iqx)
			for (int iqy = 0; iqy < qynpts; ++iqy)
			for (int iqz = 0; iqz < qznpts; ++iqz)
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
		cerr << "FileIException error!" << endl;
debugger(__LINE__, __FILE__);
		return -1;
    }

    catch(H5::DataSetIException error)
    {
		error.printError();
		cerr << "DataSetIException error!" << endl;
debugger(__LINE__, __FILE__);
		return -2;
    }

    catch(H5::DataSpaceIException error)
    {
		error.printError();
		cerr << "DataSpaceIException error!" << endl;
debugger(__LINE__, __FILE__);
		return -3;
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
		cerr << "FileIException error!" << endl;
debugger(__LINE__, __FILE__);
		return -1;
    }

    catch(H5::DataSetIException error)
    {
		error.printError();
		cerr << "DataSetIException error!" << endl;
debugger(__LINE__, __FILE__);
		return -2;
    }

    catch(H5::DataSpaceIException error)
    {
		error.printError();
		cerr << "DataSpaceIException error!" << endl;
debugger(__LINE__, __FILE__);
		return -3;
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
		cerr << "FileIException error!" << endl;
debugger(__LINE__, __FILE__);
		return -1;
    }

    catch(H5::DataSetIException error)
    {
		error.printError();
		cerr << "DataSetIException error!" << endl;
debugger(__LINE__, __FILE__);
		return -2;
    }

    catch(H5::DataSpaceIException error)
    {
		error.printError();
		cerr << "DataSpaceIException error!" << endl;
debugger(__LINE__, __FILE__);
		return -3;
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
			for (int iqt = 0; iqt < qtnpts; ++iqt)
			for (int iqx = 0; iqx < qxnpts; ++iqx)
			for (int iqy = 0; iqy < qynpts; ++iqy)
			for (int iqz = 0; iqz < qznpts; ++iqz)
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
debugger(__LINE__, __FILE__);
		return -1;
    }

    catch(H5::DataSetIException error)
    {
		error.printError();
		cerr << "DataSetIException error!" << endl;
debugger(__LINE__, __FILE__);
		return -2;
    }

    catch(H5::DataSpaceIException error)
    {
		error.printError();
		cerr << "DataSpaceIException error!" << endl;
debugger(__LINE__, __FILE__);
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
debugger(__LINE__, __FILE__);
		return -1;
    }

    catch(H5::DataSetIException error)
    {
		error.printError();
		cerr << "DataSetIException error!" << endl;
debugger(__LINE__, __FILE__);
		return -2;
    }

    catch(H5::DataSpaceIException error)
    {
		error.printError();
		cerr << "DataSpaceIException error!" << endl;
debugger(__LINE__, __FILE__);
		return -3;
    }

	return (0);
}

//*******************************************
// HDF array for resonances
//*******************************************

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
		hsize_t chunk_dims[RANK2D] = {1, chunk_size};
		cparms.setChunk( RANK2D, chunk_dims );

		hsize_t dims[RANK2D] = {Nparticle, chunk_size};
		resonance_dataspace = new H5::DataSpace (RANK2D, dims);

		resonance_dataset = new H5::DataSet( resonance_file->createDataSet(RESONANCE_DATASET_NAME, PredType::NATIVE_DOUBLE, *resonance_dataspace, cparms) );

		hsize_t count[RANK2D] = {1, chunk_size};
		hsize_t dimsm[RANK2D] = {1, chunk_size};
		hsize_t offset[RANK2D] = {0, 0};

		resonance_memspace = new H5::DataSpace (RANK2D, dimsm, NULL);
		for (int ir = 0; ir < Nparticle; ++ir)
		{
			offset[0] = ir;
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
		resonance_memspace = new H5::DataSpace (RANK2D, dimsm, NULL);
    }

    catch(FileIException error)
    {
		error.printError();
		cerr << "FileIException error!" << endl;
debugger(__LINE__, __FILE__);
		return -1;
    }

    catch(H5::DataSetIException error)
    {
		error.printError();
		cerr << "DataSetIException error!" << endl;
debugger(__LINE__, __FILE__);
		return -2;
    }

    catch(H5::DataSpaceIException error)
    {
		error.printError();
		cerr << "DataSpaceIException error!" << endl;
debugger(__LINE__, __FILE__);
		return -3;
    }

	delete [] resonance_chunk;

	return (0);
}

int CorrelationFunction::Open_resonance_HDF_array()
{
	ostringstream filename_stream_ra;
	filename_stream_ra << global_path << "/resonance_spectra.h5";
	H5std_string RESONANCE_FILE_NAME(filename_stream_ra.str().c_str());
	H5std_string RESONANCE_DATASET_NAME("ra");

	try
    {
		Exception::dontPrint();
	
		hsize_t dimsm[RANK2D] = {1, chunk_size};
		hsize_t dims[RANK2D] = {Nparticle, chunk_size};

		resonance_dataspace = new H5::DataSpace (RANK2D, dims);
		resonance_file = new H5::H5File(RESONANCE_FILE_NAME, H5F_ACC_RDWR);
		resonance_dataset = new H5::DataSet( resonance_file->openDataSet( RESONANCE_DATASET_NAME ) );
		resonance_memspace = new H5::DataSpace (RANK2D, dimsm, NULL);
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
	
		hsize_t offset[RANK2D] = {local_pid, 0};
		hsize_t count[RANK2D] = {1, chunk_size};				// == chunk_dims
		resonance_dataspace->selectHyperslab(H5S_SELECT_SET, count, offset);

		// use loaded chunk to fill resonance_array_to_fill
		int iidx = 0;
		for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
		for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
		for (int iqt = 0; iqt < qtnpts; ++iqt)
		for (int iqx = 0; iqx < qxnpts; ++iqx)
		for (int iqy = 0; iqy < qynpts; ++iqy)
		for (int iqz = 0; iqz < qznpts; ++iqz)
		for (int itrig = 0; itrig < ntrig; ++itrig)
		{
			double temp = resonance_array_to_use[ipt][ipphi][iqt][iqx][iqy][iqz][itrig];
			resonance_chunk[iidx] = temp;
			++iidx;
		}

		resonance_dataset->write(resonance_chunk, PredType::NATIVE_DOUBLE, *resonance_memspace, *resonance_dataspace);
   }

    catch(FileIException error)
    {
		error.printError();
		cerr << "FileIException error!" << endl;
debugger(__LINE__, __FILE__);
		return -1;
    }

    catch(H5::DataSetIException error)
    {
		error.printError();
		cerr << "DataSetIException error!" << endl;
debugger(__LINE__, __FILE__);
		return -2;
    }

    catch(H5::DataSpaceIException error)
    {
		error.printError();
		cerr << "DataSpaceIException error!" << endl;
debugger(__LINE__, __FILE__);
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
	
		hsize_t offset[RANK2D] = {local_pid, 0};
		hsize_t count[RANK2D] = {1, chunk_size};				// == chunk_dims
		resonance_dataspace->selectHyperslab(H5S_SELECT_SET, count, offset);

		resonance_dataset->read(resonance_chunk, PredType::NATIVE_DOUBLE, *resonance_memspace, *resonance_dataspace);

		// use loaded chunk to fill resonance_array_to_fill
		int iidx = 0;
		for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
		for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
		for (int iqt = 0; iqt < qtnpts; ++iqt)
		for (int iqx = 0; iqx < qxnpts; ++iqx)
		for (int iqy = 0; iqy < qynpts; ++iqy)
		for (int iqz = 0; iqz < qznpts; ++iqz)
		for (int itrig = 0; itrig < ntrig; ++itrig)
		{
			double temp = resonance_chunk[iidx];
			//cerr << "INFODUMP: " << local_pid << "   " << ipt << "   " << ipphi << "   " << iqt << "   "
			//	<< iqx << "   " << iqy << "   " << iqz << "   " << itrig << "   " << temp << endl;
			resonance_array_to_fill[ipt][ipphi][iqt][iqx][iqy][iqz][itrig] = temp;
			++iidx;
		}

   }

    catch(FileIException error)
    {
		error.printError();
		cerr << "FileIException error!" << endl;
debugger(__LINE__, __FILE__);
		return -1;
    }

    catch(H5::DataSetIException error)
    {
		error.printError();
		cerr << "DataSetIException error!" << endl;
debugger(__LINE__, __FILE__);
		return -2;
    }

    catch(H5::DataSpaceIException error)
    {
		error.printError();
		cerr << "DataSpaceIException error!" << endl;
debugger(__LINE__, __FILE__);
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
	
		hsize_t offset[RANK2D] = {reso_idx_to_be_copied, 0};
		hsize_t count[RANK2D] = {1, chunk_size};
		resonance_dataspace->selectHyperslab(H5S_SELECT_SET, count, offset);

		// load resonance to be copied first
		resonance_dataset->read(resonance_chunk, PredType::NATIVE_DOUBLE, *resonance_memspace, *resonance_dataspace);

		// now set this to current resonance
		offset[0] = current_resonance_index;
		resonance_dataspace->selectHyperslab(H5S_SELECT_SET, count, offset);
		resonance_dataset->write(resonance_chunk, PredType::NATIVE_DOUBLE, *resonance_memspace, *resonance_dataspace);
   }

    catch(FileIException error)
    {
		error.printError();
		cerr << "FileIException error!" << endl;
debugger(__LINE__, __FILE__);
		return -1;
    }

    catch(H5::DataSetIException error)
    {
		error.printError();
		cerr << "DataSetIException error!" << endl;
debugger(__LINE__, __FILE__);
		return -2;
    }

    catch(H5::DataSpaceIException error)
    {
		error.printError();
		cerr << "DataSpaceIException error!" << endl;
debugger(__LINE__, __FILE__);
		return -3;
    }

	delete [] resonance_chunk;

	return (0);
}

int CorrelationFunction::Dump_resonance_HDF_array_spectra(string output_filename, double ******* resonance_array_to_use)
{
	double * resonance_chunk = new double [chunk_size];

	ostringstream filename_stream_ra;
	filename_stream_ra << global_path << "/resonance_spectra.h5";
	H5std_string RESONANCE_FILE_NAME(filename_stream_ra.str().c_str());
	H5std_string RESONANCE_DATASET_NAME("ra");

	ostringstream filename_stream;
	filename_stream << global_path << "/" << output_filename;
	ofstream out;
	out.open(filename_stream.str().c_str());

	try
    {
		Exception::dontPrint();
	
		// use loaded chunk to fill resonance_array_to_fill
		for(int ir = 0; ir < Nparticle; ir++)
		{
			hsize_t offset[RANK2D] = {ir, 0};
			hsize_t count[RANK2D] = {1, chunk_size};				// == chunk_dims
			resonance_dataspace->selectHyperslab(H5S_SELECT_SET, count, offset);
	
			resonance_dataset->read(resonance_chunk, PredType::NATIVE_DOUBLE, *resonance_memspace, *resonance_dataspace);

			int iidx = 0;
			for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
			for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
			for (int iqt = 0; iqt < qtnpts; ++iqt)
			for (int iqx = 0; iqx < qxnpts; ++iqx)
			for (int iqy = 0; iqy < qynpts; ++iqy)
			for (int iqz = 0; iqz < qznpts; ++iqz)
			for (int itrig = 0; itrig < ntrig; ++itrig)
				resonance_array_to_use[ipt][ipphi][iqt][iqx][iqy][iqz][itrig] = resonance_chunk[iidx++];

			for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
			{
				for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
					out << scientific << setprecision(8) << setw(12) << resonance_array_to_use[ipt][ipphi][iqt0][iqx0][iqy0][iqz0][0] << "   ";
				out << endl;
			}
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
		cerr << "DataSetIException error!" << endl;
		return -2;
    }

    catch(H5::DataSpaceIException error)
    {
		error.printError();
		cerr << "DataSpaceIException error!" << endl;
		return -3;
    }

	delete [] resonance_chunk;
	out.close();

	return (0);

}

//*******************************************
// HDF array for correlator snapshots
//*******************************************

int CorrelationFunction::Initialize_snapshot_HDF_array()
{
	double * snapshot_chunk = new double [q_space_size];

	ostringstream filename_stream_ra;
	filename_stream_ra << global_path << "/correlator_snapshots.h5";
	H5std_string SNAPSHOT_FILE_NAME(filename_stream_ra.str().c_str());
	H5std_string SNAPSHOT_DATASET_NAME("cs");

	try
    {
		Exception::dontPrint();
//debugger(__LINE__, __FILE__);
		snapshot_file = new H5::H5File(SNAPSHOT_FILE_NAME, H5F_ACC_TRUNC);
//debugger(__LINE__, __FILE__);
		DSetCreatPropList cparms;
		hsize_t chunk_dims[RANK2D] = {1, q_space_size};
		cparms.setChunk( RANK2D, chunk_dims );
//debugger(__LINE__, __FILE__);
		cout << "Initializing with dimensions: " << chosen_resonances.size() * n_interp_pT_pts * n_interp_pphi_pts << " x " << q_space_size << endl;
		hsize_t dims[RANK2D] = {(chosen_resonances.size()+1) * n_interp_pT_pts * n_interp_pphi_pts, q_space_size};
		snapshot_dataspace = new H5::DataSpace (RANK2D, dims);
//debugger(__LINE__, __FILE__);
		snapshot_dataset = new H5::DataSet( snapshot_file->createDataSet(SNAPSHOT_DATASET_NAME, PredType::NATIVE_DOUBLE, *snapshot_dataspace, cparms) );
//debugger(__LINE__, __FILE__);
		hsize_t count[RANK2D] = {1, q_space_size};
		hsize_t dimsm[RANK2D] = {1, q_space_size};
		hsize_t offset[RANK2D] = {0, 0};
//debugger(__LINE__, __FILE__);
		snapshot_memspace = new H5::DataSpace (RANK2D, dimsm, NULL);
		int snapshot_idx = 0;
//debugger(__LINE__, __FILE__);
		for (int ir = 0; ir < chosen_resonances.size()+1; ++ir)
		for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
		for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
		{
			offset[0] = snapshot_idx;
			snapshot_dataspace->selectHyperslab(H5S_SELECT_SET, count, offset);

			for (int iidx = 0; iidx < q_space_size; ++iidx)
				snapshot_chunk[iidx] = 0.0;

			//initialize everything with zeros
			snapshot_dataset->write(snapshot_chunk, PredType::NATIVE_DOUBLE, *snapshot_memspace, *snapshot_dataspace);

			++snapshot_idx;
		}
//debugger(__LINE__, __FILE__);
		snapshot_memspace->close();
		snapshot_dataset->close();
		snapshot_file->close();
//debugger(__LINE__, __FILE__);
		delete snapshot_memspace;
		delete snapshot_file;
		delete snapshot_dataset;
//debugger(__LINE__, __FILE__);
		snapshot_file = new H5::H5File(SNAPSHOT_FILE_NAME, H5F_ACC_RDWR);
		snapshot_dataset = new H5::DataSet( snapshot_file->openDataSet( SNAPSHOT_DATASET_NAME ) );
		snapshot_memspace = new H5::DataSpace (RANK2D, dimsm, NULL);
//debugger(__LINE__, __FILE__);
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

	delete [] snapshot_chunk;

	return (0);
}

int CorrelationFunction::Set_correlator_snapshot(int icr, double ******* snapshot_array_to_use)
{
	double * snapshot_chunk = new double [q_space_size];
string tmpstring = (icr == 0) ?
					all_particles[target_particle_id].name :
					all_particles[chosen_resonances[icr-1]].name;
//if (VERBOSE > 0) *global_out_stream_ptr << "Particle(chosen_resonances = " << chosen_resonances[icr-1] << ") = " << tmpstring << endl;

	ostringstream filename_stream_ra;
	filename_stream_ra << global_path << "/correlator_snapshots.h5";
	H5std_string SNAPSHOT_FILE_NAME(filename_stream_ra.str().c_str());
	H5std_string SNAPSHOT_DATASET_NAME("cs");

	try
    {
		Exception::dontPrint();
	
		hsize_t count[RANK2D] = {1, q_space_size};				// == chunk_dims

		// use loaded chunk to fill snapshot_array_to_fill
		int snapshot_idx = icr * n_interp_pT_pts * n_interp_pphi_pts;
		for (int ipt = 0; ipt < n_interp_pT_pts; ++ipt)
		for (int ipphi = 0; ipphi < n_interp_pphi_pts; ++ipphi)
		{
			snapshot_array_to_use[ipt][ipphi][iqt0][iqx0][iqy0][iqz0][1] = 0.0;
			double nonFTd_spectra = spectra[target_particle_id][ipt][ipphi];
			//double nonFTd_spectra = snapshot_array_to_use[ipt][ipphi][iqt0][iqx0][iqy0][iqz0][0];

			hsize_t offset[RANK2D] = {snapshot_idx, 0};
			snapshot_dataspace->selectHyperslab(H5S_SELECT_SET, count, offset);
			
			int iidx = 0;
			for (int iqt = 0; iqt < qtnpts; ++iqt)
			for (int iqx = 0; iqx < qxnpts; ++iqx)
			for (int iqy = 0; iqy < qynpts; ++iqy)
			for (int iqz = 0; iqz < qznpts; ++iqz)
			{
				// output all FT'd spectra
				double cos_transf_spectra = snapshot_array_to_use[ipt][ipphi][iqt][iqx][iqy][iqz][0];
				double sin_transf_spectra = snapshot_array_to_use[ipt][ipphi][iqt][iqx][iqy][iqz][1];
				snapshot_chunk[iidx] = 1. + (cos_transf_spectra*cos_transf_spectra + sin_transf_spectra*sin_transf_spectra)/(nonFTd_spectra*nonFTd_spectra);
				//if (icr==0)
				//	cerr << icr << "   " << iqt << "   " << iqx << "   " << iqy << "   " << iqz << "   "
				//		<< nonFTd_spectra << "   " << cos_transf_spectra << "   " << sin_transf_spectra << "   " << snapshot_chunk[iidx] << endl;
	
				//initialize everything with zeros
				snapshot_dataset->write(snapshot_chunk, PredType::NATIVE_DOUBLE, *snapshot_memspace, *snapshot_dataspace);
				++iidx;
			}
			++snapshot_idx;
		}

		snapshot_dataset->write(snapshot_chunk, PredType::NATIVE_DOUBLE, *snapshot_memspace, *snapshot_dataspace);
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

	delete [] snapshot_chunk;

	return (0);
}

int CorrelationFunction::Extrapolate_over_snapshot_HDF_array()
{
	double * snapshot_chunk = new double [q_space_size];
	double ** pT_pphi_snapshots  = new double * [chosen_resonances.size()+1];
	for (int icr = 0; icr < chosen_resonances.size()+1; ++icr)
		pT_pphi_snapshots[icr] = new double [q_space_size];


	ostringstream filename_stream_ra;
	filename_stream_ra << global_path << "/correlator_snapshots.h5";
	H5std_string SNAPSHOT_FILE_NAME(filename_stream_ra.str().c_str());
	H5std_string SNAPSHOT_DATASET_NAME("cs");

	try
    {
		Exception::dontPrint();

		//put pT-pphi loops in momentarily...
		int ipt = 0;
		int ipphi = 0;
	
		hsize_t count[RANK2D] = {1, q_space_size};				// == chunk_dims

		//read in all resonance snapshots for fixed pT-pphi combo.
		for (int icr = 0; icr < chosen_resonances.size()+1; ++icr)
		{
debugger(__LINE__, __FILE__);
			int snapshot_idx = icr * n_interp_pT_pts * n_interp_pphi_pts;
	
			hsize_t offset[RANK2D] = {snapshot_idx, 0};

			snapshot_dataspace->selectHyperslab(H5S_SELECT_SET, count, offset);
	
			snapshot_dataset->read(snapshot_chunk, PredType::NATIVE_DOUBLE, *snapshot_memspace, *snapshot_dataspace);
			for (int iq = 0; iq < q_space_size; ++iq)
				pT_pphi_snapshots[icr][iq] = snapshot_chunk[iq];
		}

		// with pT_pphi_snapshots fully loaded, do extrapolation over each point in q-space
		for (int iq = 0; iq < q_space_size; ++iq)
		{
			vector<double> resonance_pc_markers;
			vector<double> pT_pphi_qpt_snapshots;
			//vector<double> projected_results;
			double chisq = 0.0;
			int polynomial_fit_order = min(4, (int)chosen_resonances.size());

			resonance_pc_markers.push_back(snapshot_fractions[0]);
			pT_pphi_qpt_snapshots.push_back( pT_pphi_snapshots[0][iq] );
			
			for (int icr = chosen_resonances.size(); icr > 0; --icr)
			{
				resonance_pc_markers.push_back( snapshot_fractions[icr] - snapshot_fractions[icr-1] + resonance_pc_markers.back() );
				cout << "RE debug: " << snapshot_fractions[icr] << "   " << snapshot_fractions[icr-1] << endl;
				pT_pphi_qpt_snapshots.push_back( pT_pphi_snapshots[icr][iq] - pT_pphi_snapshots[icr-1][iq] + pT_pphi_qpt_snapshots.back() );
				cout << "RE2 debug: " << pT_pphi_snapshots[icr][iq] << "   " << pT_pphi_snapshots[icr-1][iq] << endl;
			}

			for (int icr = 0; icr < chosen_resonances.size()+1; ++icr)
				cout << "Resonance extrapolation: " << iq << "   " << icr << "   " << resonance_pc_markers[icr] << "   " << pT_pphi_qpt_snapshots[icr] << endl;

			cout << "iq = " << iq << ": Projected result = " << gsl_polynomial_fit(resonance_pc_markers, pT_pphi_qpt_snapshots, polynomial_fit_order, chisq) << endl;
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
		cerr << "DataSetIException error!" << endl;
		return -2;
    }

    catch(H5::DataSpaceIException error)
    {
		error.printError();
		cerr << "DataSpaceIException error!" << endl;
		return -3;
    }

	delete [] snapshot_chunk;
	for (int icr = 0; icr < chosen_resonances.size()+1; ++icr)
		delete [] pT_pphi_snapshots[icr];
	delete [] pT_pphi_snapshots;

	return (0);
}

//End of file
