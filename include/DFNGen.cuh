#pragma once

#include <cuda_runtime.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <stdexcept>
#include <thrust/random.h>
#include <string>
#include <fstream>
#include <omp.h>
#include <thrust/remove.h>
#include <thrust/functional.h>

#include "CPUSecond.cuh"
#include "PDF.cuh"
#include "Operators.cuh"
#include "Quaternion.cuh"
#include "FunctionGenFractureAttributes.cuh"
#include "HDF5API.cuh"
#include "Graph.cuh"

#include "OCCTAPI.cuh"
#include "DeviceIfFracturesIntersect.cuh"
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
        thrust::host_vector<int4> VerticesIndex;             // vertice label from 0
        thrust::host_vector<int2> IntersectionFracturePairs; // fracture label from 0
        thrust::host_vector<thrust::host_vector<int>> FractureClusters;

        unsigned long RandomSeed = 0;

        // should be defined before gen
        int NumFractures = 0;
        // thrust::host_vector<int> NumVerticesEachFracture;
        int NumFracturePairsCheckAtOnce = 256;
        int Nproc = 10;

        thrust::default_random_engine rng;

        bool IfUsingGPUCheckFracturesIntersect = true;

    public:
        DFNGen() {};

        void InitRandomEngine();

        double3 (*FractureLocationDistributionFunctionPointer)(const double *rand_0_1, const double *Para);
        void GenerateFractureLocations(const double *Para);

        double3 (*FractureNormalDistributionFunctionPointer)(const double *rand_0_1, const double *Para);
        void GenerateFractureNormals(const double *Para);

        int (*FractureVerticesFunctionPointer)(const double *rand_0_1, const double *Para, double3 *vertices);
        void GenerateFractureVertices(const double *Para);

        void IdentifyIntersectedFractures();
        void IdentifyFractureClusters();

        void TruncateFracturesByStlDomain(const std::string &stlname)
        {
            TopoDS_Solid Domain_stl = cuDFNsys::LoadAStlDomain(stlname);

            for (int i = 0; i < this->NumFractures; ++i)
            {
                int NumVertices = this->VerticesIndex[i].y - this->VerticesIndex[i].x + 1;
                double3 *V[NumVertices];
                for (int j = 0; j < NumVertices; ++j)
                    V[j] = &(this->FractureVertices[this->VerticesIndex[i].x + j]);

                std::cout << "Frac i = " << i << "\n";
                for (int j = 0; j < NumVertices; ++j)
                    std::cout << *V[j] << " ";
                std::cout << "\n";

                TopoDS_Face Polygon = cuDFNsys::FractureToTopoDS_Face(*V, NumVertices);
                thrust::host_vector<double3> NewFrac = cuDFNsys::TruncatePolygon(Polygon, Domain_stl);

                std::cout << "Frac i after truncation = " << i << "\n";
                for (auto e : NewFrac)
                    std::cout << e << " ";
                std::cout << "\n\n";

                this->VerticesIndex[i].z = this->TruncatedFractureVertices.size();
                this->VerticesIndex[i].w = this->TruncatedFractureVertices.size() + NewFrac.size() - 1;
                this->TruncatedFractureVertices.insert(this->TruncatedFractureVertices.end(), NewFrac.begin(), NewFrac.end());
                
            }
        };

        void OutputDFNGen(const std::string HDF5filename);
        void OutputXMF(const std::string XMFfilename, const std::string HDF5filename);
    };

    cuDFNsys::DFNGen operator+(const cuDFNsys::DFNGen &a, const cuDFNsys::DFNGen &b)
    {
        return a;
    };

    void cuDFNsys::DFNGen::InitRandomEngine()
    {
        this->rng = thrust::default_random_engine(RandomSeed);
    };

    void cuDFNsys::DFNGen::GenerateFractureLocations(const double *Para)
    {
        if (this->NumFractures == 0)
            throw std::runtime_error("From `cuDFNsys::DFNGen::GenerateFractureLocation`, the number of fractures is zero");
        this->Location.resize(this->NumFractures);

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

        thrust::uniform_real_distribution<double> dist(0.0, 1.0);

        if (!FractureNormalDistributionFunctionPointer)
            throw std::runtime_error("From `cuDFNsys::DFNGen::GenerateFractureNormals`, the `FractureNormalDistributionFunctionPointer` is not specified");

        for (size_t i = 0; i < this->NumFractures; ++i)
        {
            double rand_0_1[3] = {dist(rng), dist(rng), dist(rng)};
            // std::cout << rand_0_1[0] << ", " << rand_0_1[1] << ", " << rand_0_1[2] << std::endl;
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
            ff2 = ff2.DescribeRotation(RotationAxis, acos(this->NormalVec[i].z));

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

    void cuDFNsys::DFNGen::IdentifyIntersectedFractures()
    {

        if (!this->IfUsingGPUCheckFracturesIntersect)
        { // Generate Pairs
            thrust::host_vector<int2> TotalPair(
                this->NumFractures * floor((this->NumFractures - 1) / 2) + (this->NumFractures - 1) % 2 * this->NumFractures * 0.5);

            for (int i = 0; i < TotalPair.size(); ++i)
            {
                TotalPair[i].x = floor((pow(2 * (i + 1), 0.5) + 1 / 2.0));
                TotalPair[i].y = i - 0.5 * TotalPair[i].x * (TotalPair[i].x - 1);
            }

            int NumFracturePairTotal = this->NumFractures * floor((this->NumFractures - 1) / 2) + (this->NumFractures - 1) % 2 * this->NumFractures * 0.5;

            for (int i = 0; i < NumFracturePairTotal; i += NumFracturePairsCheckAtOnce)
            {
                int HowManyPairsCheck = (i + NumFracturePairsCheckAtOnce >= NumFracturePairTotal ? NumFracturePairTotal : i + NumFracturePairsCheckAtOnce) - i;

                thrust::host_vector<int2> SubPair(HowManyPairsCheck, make_int2(-1, -1));
#pragma omp parallel for schedule(static) num_threads(this->Nproc)
                for (int j = i; j < i + HowManyPairsCheck; ++j)
                {
                    int x_ = floor((pow(2 * (j + 1), 0.5) + 1 / 2.0));
                    int y_ = j - 0.5 * x_ * (x_ - 1);

                    thrust::host_vector<double3> FractureV1(this->FractureVertices.begin() + this->VerticesIndex[x_].x,
                                                            this->FractureVertices.begin() + this->VerticesIndex[x_].y + 1);
                    thrust::host_vector<double3> FractureV2(this->FractureVertices.begin() + this->VerticesIndex[y_].x,
                                                            this->FractureVertices.begin() + this->VerticesIndex[y_].y + 1);

                    thrust::host_vector<double3> IntersectionEdge;
                    bool If_intersect = cuDFNsys::FractureIntersectionCheckOCCT(FractureV1.data(), FractureV2.data(), this->VerticesIndex[x_].y - this->VerticesIndex[x_].x + 1,
                                                                                this->VerticesIndex[y_].y - this->VerticesIndex[y_].x + 1, false, IntersectionEdge);
                    if (If_intersect)
                        SubPair[j - i] = make_int2(x_, y_);

                    // std::cout << x_ << ", " << y_ << ", " << If_intersect << "\n";
                    // for (auto e : FractureV1)
                    //     std::cout << e << " ";
                    // std::cout << "\n";
                    // for (auto e : FractureV2)
                    //     std::cout << e << " ";
                    // std::cout << "\n";
                }

                auto new_end = thrust::remove_if(
                    SubPair.begin(),
                    SubPair.end(),
                    [](const int2 &pair)
                    {
                        return pair.x == -1; // Predicate to remove elements where .x == -1
                    });
                SubPair.erase(new_end, SubPair.end());
                IntersectionFracturePairs.insert(IntersectionFracturePairs.end(), SubPair.begin(), SubPair.end());
            }

            // for (auto e : IntersectionFracturePairs)
            //     std::cout << e << " ";
            // std::cout << "\n";
        }
        else
        {
            int NumFracturePairTotal = this->NumFractures * floor((this->NumFractures - 1) / 2) + (this->NumFractures - 1) % 2 * this->NumFractures * 0.5;

            thrust::device_vector<double3> FractureVertices_DEV = this->FractureVertices;
            thrust::device_vector<int4> VerticesIndex_DEV = this->VerticesIndex;
            thrust::device_vector<int2> IntersectionFracturePairs_DEV(NumFracturePairTotal, make_int2(-1, -1));

            double3 *FractureVertices_DEV_ptr = thrust::raw_pointer_cast(FractureVertices_DEV.data());
            int4 *VerticesIndex_DEV_ptr = thrust::raw_pointer_cast(VerticesIndex_DEV.data());
            int2 *IntersectionFracturePairs_DEV_ptr = thrust::raw_pointer_cast(IntersectionFracturePairs_DEV.data());

            cuDFNsys::DeviceIfFracturesIntersect<<<NumFracturePairTotal / 256 + 1, 256>>>(NumFracturePairTotal, FractureVertices_DEV_ptr, VerticesIndex_DEV_ptr, IntersectionFracturePairs_DEV_ptr);
            cudaDeviceSynchronize();

            auto new_end = thrust::remove_if(
                IntersectionFracturePairs_DEV.begin(),
                IntersectionFracturePairs_DEV.end(),
                [] __host__ __device__(const int2 &pair)
                {
                    return pair.x == -1; // Inline predicate
                });
            // // Resize the vector to remove the "removed" elements
            IntersectionFracturePairs_DEV.erase(new_end, IntersectionFracturePairs_DEV.end());
            this->IntersectionFracturePairs = IntersectionFracturePairs_DEV;

            // for (auto e : IntersectionFracturePairs)
            //     std::cout << e << " ";
            // std::cout << "\n";
        }
    };

    void cuDFNsys::DFNGen::IdentifyFractureClusters()
    {
        cuDFNsys::Graph G(this->NumFractures, this->IntersectionFracturePairs);
        G.UseDFS(this->FractureClusters);
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
        cuDFNsys::WriteInt2ToHDF5(HDF5filename, "IntersectionFracturePairs", IntersectionFracturePairs.data(),
                                  IntersectionFracturePairs.size(), true);

        for (int i = 0; i < this->FractureClusters.size(); ++i)
            cuDFNsys::WriteIntToHDF5(HDF5filename, "FractureCluster" + std::to_string(i), this->FractureClusters[i].data(),
                                     this->FractureClusters[i].size(), true);

        int Tem[1] = {(int)RandomSeed};
        cuDFNsys::WriteIntToHDF5(HDF5filename, "RandomSeed", Tem, 1, true);
        Tem[0] = NumFractures;
        cuDFNsys::WriteIntToHDF5(HDF5filename, "NumFractures", Tem, 1, true);

        Tem[0] = NumFracturePairsCheckAtOnce;
        cuDFNsys::WriteIntToHDF5(HDF5filename, "NumFracturePairsCheckAtOnce", Tem, 1, true);
        Tem[0] = Nproc;
        cuDFNsys::WriteIntToHDF5(HDF5filename, "Nproc", Tem, 1, true);
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

        thrust::host_vector<int> ScalarValue_cluster(this->NumFractures);
        for (int i = 0; i < this->FractureClusters.size(); ++i)
            for (int j = 0; j < this->FractureClusters[i].size(); ++j)
                ScalarValue_cluster[this->FractureClusters[i][j]] = i;

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
        xmfFile << "      <Attribute Name=\"Fracture normal vector\" AttributeType=\"Vector\" Center=\"Cell\">\n";
        xmfFile << "        <DataItem Format=\"HDF\" DataType=\"Float\" Precision=\"8\" Dimensions=\"" << this->NumFractures << " 3\">\n";
        xmfFile << "          " << HDF5filename << ":/NormalVec\n";
        xmfFile << "        </DataItem>\n";
        xmfFile << "      </Attribute>\n";

        xmfFile << "      <Attribute Name=\"Cluster\" AttributeType=\"Scalar\" Center=\"Cell\">\n";
        xmfFile << "        <DataItem Dimensions=\"" << this->NumFractures << "\" Format=\"XML\" NumberType=\"Float\" Precision=\"8\">\n";
        xmfFile << "          ";
        for (int i = 0; i < ScalarValue_cluster.size(); ++i)
            xmfFile << ScalarValue_cluster[i] << (i == ScalarValue_cluster.size() - 1 ? "\n" : " ");
        xmfFile << "        </DataItem>\n";
        xmfFile << "      </Attribute>\n";

        xmfFile << "   </Grid>\n";
        xmfFile << "  </Domain>\n";
        xmfFile << "</Xdmf>\n";

        xmfFile.close();
    }

};