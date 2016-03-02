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

	//bool file_does_not_already_exist = !fexists(filename_stream_ra.str().c_str());
	bool file_does_not_already_exist = true;	//force full initialization for timebeing...

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
		if (file_does_not_already_exist)
		{
			*global_out_stream_ptr << "HDF resonance file doesn't exist!  Initializing to zero..." << endl;

			for (int ir = 0; ir < Nparticle; ++ir)
			{
				//debugger(__LINE__, __FILE__);
				//print_now();
	
				offset[0] = ir;
				resonance_dataspace->selectHyperslab(H5S_SELECT_SET, count, offset);
	
				for (int iidx = 0; iidx < chunk_size; ++iidx)
					resonance_chunk[iidx] = 0.0;
	
				//initialize everything with zeros
				resonance_dataset->write(resonance_chunk, PredType::NATIVE_DOUBLE, *resonance_memspace, *resonance_dataspace);
			}
		}
		/*resonance_memspace->close();
		resonance_dataset->close();
		resonance_file->close();
		delete resonance_memspace;
		delete resonance_file;
		delete resonance_dataset;
		resonance_file = new H5::H5File(RESONANCE_FILE_NAME, H5F_ACC_RDWR);
		resonance_dataset = new H5::DataSet( resonance_file->openDataSet( RESONANCE_DATASET_NAME ) );
		resonance_memspace = new H5::DataSpace (RANK2D, dimsm, NULL);*/
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

//End of file
