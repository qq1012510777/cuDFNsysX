#pragma once

#include <cuda_runtime.h>
#include <thrust/host_vector.h>
#include <string>
#include "H5Cpp.h"

namespace cuDFNsys
{

    void WriteDouble3ToHDF5(const std::string &filename,
                            const std::string &datasetName,
                            const double3 *data, const int &n, const bool &AddDataSet)
    {
        try
        {
            // Define dimensions
            // hsize_t n = data.size();  // Number of elements in the vector
            hsize_t dims[2] = {(hsize_t)n, 3}; // n x 3 dataset

            // Convert vector<double3> to a simple array for HDF5
            std::vector<double> flatData;
            flatData.resize(n * 3);
            for (int i = 0; i < n * 3; i += 3)
            {
                flatData[i] = data[(i / 3)].x;
                flatData[i + 1] = data[(i / 3)].y;
                flatData[i + 2] = data[(i / 3)].z;
            }

            // Create HDF5 file
            H5::H5File file(filename, (AddDataSet == false ? H5F_ACC_TRUNC : H5F_ACC_RDWR));

            // Create dataspace with specified dimensions
            H5::DataSpace dataspace(2, dims);

            // Create dataset for double data type
            H5::DataSet dataset = file.createDataSet(datasetName, H5::PredType::NATIVE_DOUBLE, dataspace);

            // Write data to dataset
            dataset.write(flatData.data(), H5::PredType::NATIVE_DOUBLE);

            file.close();
        }
        catch (H5::FileIException &error)
        {
            error.printErrorStack();
        }
        catch (H5::DataSetIException &error)
        {
            error.printErrorStack();
        }
        catch (H5::DataSpaceIException &error)
        {
            error.printErrorStack();
        }
    };

    void WriteInt4ToHDF5(const std::string &filename,
                         const std::string &datasetName,
                         const int4 *data, const int &n, const bool &AddDataSet)
    {
        try
        {
            // Define dimensions
            // hsize_t n = data.size();  // Number of elements in the vector
            hsize_t dims[2] = {(hsize_t)n, 4}; // n x 3 dataset

            // Convert vector<double3> to a simple array for HDF5
            std::vector<int> flatData;
            flatData.resize(n * 4);
            for (int i = 0; i < n * 4; i += 4)
            {
                flatData[i] = data[(i / 4)].x;
                flatData[i + 1] = data[(i / 4)].y;
                flatData[i + 2] = data[(i / 4)].z;
                flatData[i + 3] = data[(i / 4)].w;
            }

            // Create HDF5 file
            H5::H5File file(filename, (AddDataSet == false ? H5F_ACC_TRUNC : H5F_ACC_RDWR));

            // Create dataspace with specified dimensions
            H5::DataSpace dataspace(2, dims);

            // Create dataset for double data type
            H5::DataSet dataset = file.createDataSet(datasetName, H5::PredType::NATIVE_INT, dataspace);

            // Write data to dataset
            dataset.write(flatData.data(), H5::PredType::NATIVE_INT);

            file.close();
        }
        catch (H5::FileIException &error)
        {
            error.printErrorStack();
        }
        catch (H5::DataSetIException &error)
        {
            error.printErrorStack();
        }
        catch (H5::DataSpaceIException &error)
        {
            error.printErrorStack();
        }
    };

    void WriteDoubleToHDF5(const std::string &filename,
                           const std::string &datasetName,
                           const double *data, const int &n, const bool &AddDataSet)
    {
        try
        {
            // Define dimensions
            //hsize_t n = data.size();  // Number of elements in the vector
            hsize_t dims[2] = {(hsize_t)n, 1}; // n x 3 dataset

            // Convert vector<double3> to a simple array for HDF5
            // std::vector<double> flatData(data.begin(), data.end());

            // Create HDF5 file
            H5::H5File file(filename, (AddDataSet == false ? H5F_ACC_TRUNC : H5F_ACC_RDWR));

            // Create dataspace with specified dimensions
            H5::DataSpace dataspace(2, dims);

            // Create dataset for double data type
            H5::DataSet dataset = file.createDataSet(datasetName, H5::PredType::NATIVE_DOUBLE, dataspace);

            // Write data to dataset
            dataset.write(data, H5::PredType::NATIVE_DOUBLE);

            file.close();
        }
        catch (H5::FileIException &error)
        {
            error.printErrorStack();
        }
        catch (H5::DataSetIException &error)
        {
            error.printErrorStack();
        }
        catch (H5::DataSpaceIException &error)
        {
            error.printErrorStack();
        }
    };

    void WriteIntToHDF5(const std::string &filename,
                        const std::string &datasetName,
                        const int *data, const int &n, const bool &AddDataSet)
    {
        try
        {
            // Define dimensions
            //hsize_t n = data.size();  // Number of elements in the vector
            hsize_t dims[2] = {(hsize_t)n, 1}; // n x 3 dataset

            // Convert vector<double3> to a simple array for HDF5
            // std::vector<int> flatData(data.begin(), data.end());

            // Create HDF5 file
            H5::H5File file(filename, (AddDataSet == false ? H5F_ACC_TRUNC : H5F_ACC_RDWR));

            // Create dataspace with specified dimensions
            H5::DataSpace dataspace(2, dims);

            // Create dataset for double data type
            H5::DataSet dataset = file.createDataSet(datasetName, H5::PredType::NATIVE_INT, dataspace);

            // Write data to dataset
            dataset.write(data, H5::PredType::NATIVE_INT);

            file.close();
        }
        catch (H5::FileIException &error)
        {
            error.printErrorStack();
        }
        catch (H5::DataSetIException &error)
        {
            error.printErrorStack();
        }
        catch (H5::DataSpaceIException &error)
        {
            error.printErrorStack();
        }
    };

};