#pragma once

#include <cuda_runtime.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <stdexcept>
#include <thrust/random.h>
#include <string>
#include <fstream>

#include "PDF.cuh"
#include "Operators.cuh"
#include "Quaternion.cuh"
#include "FunctionGenFractureAttributes.cuh"
#include "HDF5API.cuh"

__global__ void GenerateVertices() {

};

namespace cuDFNsys
{
    class DFNGen
    {
    public:
        // attributes
        thrust::host_vector<double> Conductivity;
        thrust::host_vector<double> Aperture;
        // thrust::host_vector<double> Radius;
        thrust::host_vector<double3> Location;
        thrust::host_vector<double3> NormalVec;
        thrust::host_vector<double3> FractureVertices;
        thrust::host_vector<double3> TruncatedFractureVertices;
        thrust::host_vector<int4> VerticesIndex;

        unsigned long RandomSeed = 0;

        // should be defined before gen
        int NumFractures = 0;
        // thrust::host_vector<int> NumVerticesEachFracture;

    public:
        DFNGen() {};

        double3 (*FractureLocationDistributionFunctionPointer)(const double *rand_0_1, const double *Para);
        void GenerateFractureLocations(const double *Para);

        double3 (*FractureNormalDistributionFunctionPointer)(const double *rand_0_1, const double *Para);
        void GenerateFractureNormals(const double *Para);

        int (*FractureVerticesFunctionPointer)(const double *rand_0_1, const double *Para, double3 *vertices);
        void GenerateFractureVertices(const double *Para);

        void OutputDFNGen(const std::string HDF5filename);
        void OutputXMF(const std::string XMFfilename, const std::string HDF5filename);
    };

    void cuDFNsys::DFNGen::GenerateFractureLocations(const double *Para)
    {
        if (this->NumFractures == 0)
            throw std::runtime_error("From `cuDFNsys::DFNGen::GenerateFractureLocation`, the number of fractures is zero");
        this->Location.resize(this->NumFractures);

        thrust::default_random_engine rng(RandomSeed);

        thrust::uniform_real_distribution<double> dist(0.0, 1.0);

        if (!FractureLocationDistributionFunctionPointer)
            throw std::runtime_error("From `cuDFNsys::DFNGen::GenerateFractureLocation`, the `FractureLocationDistributionFunctionPointer` is not specified");

        for (size_t i = 0; i < this->NumFractures; ++i)
        {
            double rand_0_1[3] = {dist(rng), dist(rng), dist(rng)};
            this->Location[i] = this->FractureLocationDistributionFunctionPointer(rand_0_1, Para);
        }
    };

    void cuDFNsys::DFNGen::GenerateFractureNormals(const double *Para)
    {
        if (this->NumFractures == 0)
            throw std::runtime_error("From `cuDFNsys::DFNGen::GenerateFractureNormals`, the number of fractures is zero");
        this->NormalVec.resize(this->NumFractures);

        thrust::default_random_engine rng(RandomSeed);

        thrust::uniform_real_distribution<double> dist(0.0, 1.0);

        if (!FractureNormalDistributionFunctionPointer)
            throw std::runtime_error("From `cuDFNsys::DFNGen::GenerateFractureNormals`, the `FractureNormalDistributionFunctionPointer` is not specified");

        for (size_t i = 0; i < this->NumFractures; ++i)
        {
            double rand_0_1[3] = {dist(rng), dist(rng), dist(rng)};
            this->NormalVec[i] = this->FractureNormalDistributionFunctionPointer(rand_0_1, Para);

            this->NormalVec[i] = this->NormalVec[i] / Double3Norm(this->NormalVec[i]);
            if (this->NormalVec[i].z < 0)
                this->NormalVec[i] = this->NormalVec[i] * -1.;
        }
    };

    void cuDFNsys::DFNGen::GenerateFractureVertices(const double *Para)
    {
        if (this->NumFractures == 0)
            throw std::runtime_error("From `cuDFNsys::DFNGen::GenerateFractureVertices`, the number of fractures is zero");

        if (!FractureVerticesFunctionPointer)
            throw std::runtime_error("From `cuDFNsys::DFNGen::GenerateFractureVertices`, the `FractureVerticesFunctionPointer` is not specified");

        // generate 2D polygons first
        thrust::default_random_engine rng(RandomSeed);

        thrust::uniform_real_distribution<double> dist(0.0, 1.0);

        int NumVertices;
        {
            double Rand_0_1[1] = {dist(rng)};
            double3 Vertices_1[1000];
            NumVertices = FractureVerticesFunctionPointer(Rand_0_1, Para, Vertices_1);
        }

        this->FractureVertices.resize(NumVertices * this->NumFractures);

        // now just 2D polygons
        cuDFNsys::Quaternion ff;
        cuDFNsys::Quaternion ff2;

        VerticesIndex.resize(this->NumFractures);

        for (int i = 0; i < this->NumFractures; ++i)
        {
            double Rand_0_1[1] = {dist(rng)};
            double3 *Vertices_1 = new double3[NumVertices];
            NumVertices = FractureVerticesFunctionPointer(Rand_0_1, Para, Vertices_1);

            double Range_[2] = {0, 2 * M_PI};
            double RotatAngle = cuDFNsys::RandomUniform(dist(rng), Range_);
            // printf("%f,\n", RotatAngle);
            ff = ff.DescribeRotation(make_double3(0, 0, 1), RotatAngle);
            // rotate to normal
            double3 RotationAxis = Double3CrossProduct(make_double3(0, 0, 1), this->NormalVec[i]);
            RotationAxis = RotationAxis / Double3Norm(RotationAxis);
            ff2 = ff2.DescribeRotation(RotationAxis, this->NormalVec[i].z);

            for (int j = 0; j < NumVertices; ++j)
            {
                // now just 2D polygons for each Vertices_1[j]
                // rotate a random angle
                this->FractureVertices[i * NumVertices + j] = Vertices_1[j];
                this->FractureVertices[i * NumVertices + j] = ff.Rotate(this->FractureVertices[i * NumVertices + j]);
                this->FractureVertices[i * NumVertices + j] = ff2.Rotate(this->FractureVertices[i * NumVertices + j]);
                this->FractureVertices[i * NumVertices + j] = this->FractureVertices[i * NumVertices + j] + this->Location[i];
            }
            delete[] Vertices_1;

            VerticesIndex[i].x = i * NumVertices;
            VerticesIndex[i].y = VerticesIndex[i].x + NumVertices - 1;
        }
        this->Conductivity.resize(this->NumFractures);
        this->Aperture.resize(this->NumFractures);
    };

    void cuDFNsys::DFNGen::OutputDFNGen(const std::string HDF5filename)
    {
        cuDFNsys::WriteDoubleToHDF5(HDF5filename, "Conductivity", Conductivity.data(), Conductivity.size(), false);
        cuDFNsys::WriteDoubleToHDF5(HDF5filename, "Aperture", Aperture.data(), Aperture.size(), true);
        // cuDFNsys::WriteDoubleToHDF5(HDF5filename, "Radius", Radius.data(), Radius.size(), true);

        cuDFNsys::WriteDouble3ToHDF5(HDF5filename, "Location", Location.data(), Location.size(), true);
        cuDFNsys::WriteDouble3ToHDF5(HDF5filename, "NormalVec", NormalVec.data(), NormalVec.size(), true);
        cuDFNsys::WriteDouble3ToHDF5(HDF5filename, "FractureVertices", FractureVertices.data(), FractureVertices.size(), true);
        cuDFNsys::WriteDouble3ToHDF5(HDF5filename, "TruncatedFractureVertices", TruncatedFractureVertices.data(), TruncatedFractureVertices.size(), true);

        cuDFNsys::WriteInt4ToHDF5(HDF5filename, "VerticesIndex", VerticesIndex.data(), VerticesIndex.size(), true);

        int Tem[1] = {(int)RandomSeed};
        cuDFNsys::WriteIntToHDF5(HDF5filename, "RandomSeed", Tem, 1, true);
        Tem[0] = NumFractures;
        cuDFNsys::WriteIntToHDF5(HDF5filename, "NumFractures", Tem, 1, true);
    };

    void cuDFNsys::DFNGen::OutputXMF(const std::string XMFfilename, const std::string HDF5filename)
    {

        std::ofstream xmfFile(XMFfilename);
        if (!xmfFile.is_open())
        {
            xmfFile.close();
            throw std::runtime_error("From `cuDFNsys::DFNGen::OutputXMF`, failed to open file: " + XMFfilename + "\n");
        }
        int Dimension_allF = 0;
        for (int i = 0; i < this->VerticesIndex.size(); ++i)
            Dimension_allF = Dimension_allF + (this->VerticesIndex[i].y - this->VerticesIndex[i].x) + 3;

        xmfFile << "<?xml version=\"1.0\" ?>\n";
        xmfFile << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
        xmfFile << "<Xdmf Version=\"2.0\">\n";
        xmfFile << "  <Domain>\n";
        xmfFile << "   <Grid Name=\"DFN_Fractures\">\n";
        xmfFile << "      <Geometry GeometryType=\"XYZ\">\n";
        xmfFile << "        <DataItem Format=\"HDF\" DataType=\"Float\" Precision=\"8\" Dimensions=\"" << this->FractureVertices.size() << " 3\">\n";
        xmfFile << "          " << HDF5filename << ":/FractureVertices\n";
        xmfFile << "        </DataItem>\n";
        xmfFile << "      </Geometry>\n";
        xmfFile << "\n";

        xmfFile << "      <!-- Topology: Mixed polygons -->\n";
        xmfFile << "      <Topology TopologyType=\"Mixed\" NumberOfElements=\"" << this->NumFractures << "\">\n";
        xmfFile << "        <DataItem Dimensions=\"" << Dimension_allF << "\" Format=\"XML\">\n";
        for (int i = 0; i < this->VerticesIndex.size(); ++i)
        {
            xmfFile << "          3 3 ";
            for (int j = 0;; j++)
            {
                xmfFile << this->VerticesIndex[i].x + j;
                if (this->VerticesIndex[i].x + j == this->VerticesIndex[i].y)
                    break;
                xmfFile << " ";
                // std::cout << this->VerticesIndex[i].x + j << ", " << this->VerticesIndex[i].y << std::endl;
            }
            xmfFile << "\n";
        }
        xmfFile << "        </DataItem>\n";
        xmfFile << "      </Topology>\n";

        xmfFile << "      <Attribute Name=\"Conductivity\" AttributeType=\"Scalar\" Center=\"Cell\">\n";
        xmfFile << "        <DataItem Dimensions=\"" << this->NumFractures << "\" Format=\"HDF\" NumberType=\"Float\" Precision=\"8\">\n";
        xmfFile << "          " << HDF5filename << ":/Conductivity\n";
        xmfFile << "        </DataItem>\n";
        xmfFile << "      </Attribute>\n";

        xmfFile << "      <Attribute Name=\"Aperture\" AttributeType=\"Scalar\" Center=\"Cell\">\n";
        xmfFile << "        <DataItem Dimensions=\"" << this->NumFractures << "\" Format=\"HDF\" NumberType=\"Float\" Precision=\"8\">\n";
        xmfFile << "          " << HDF5filename << ":/Aperture\n";
        xmfFile << "        </DataItem>\n";
        xmfFile << "      </Attribute>\n";

        xmfFile << "   </Grid>\n";
        xmfFile << "  </Domain>\n";
        xmfFile << "</Xdmf>\n";

        xmfFile.close();
    }
};