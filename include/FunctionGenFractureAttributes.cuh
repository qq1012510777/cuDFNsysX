#pragma once
#include <cuda_runtime.h>
#include "PDF.cuh"

namespace cuDFNsys
{

    __device__ __host__ double3 FractureLocationGenerationFunctionUniform(const double *rand_0_1, const double *Para)
    {
        double a1 = Para[0];
        double b1 = Para[1];

        double a2 = Para[2];
        double b2 = Para[3];

        double a3 = Para[4];
        double b3 = Para[5];

        double3 random = make_double3(a1 + (b1 - a1) * rand_0_1[0],
                                      a2 + (b2 - a2) * rand_0_1[1],
                                      a3 + (b3 - a3) * rand_0_1[2]);
        return random;
    }

    __device__ __host__ double3 FractureNormalVectorGenerationFunctionUniform(const double *rand_0_1, const double *Para)
    {
        double a1 = Para[0];
        double b1 = Para[1];

        double a2 = Para[2];
        double b2 = Para[3];

        double a3 = Para[4];
        double b3 = Para[5];

        double3 random = make_double3(a1 + (b1 - a1) * rand_0_1[0],
                                      a2 + (b2 - a2) * rand_0_1[1],
                                      a3 + (b3 - a3) * rand_0_1[2]);
        return random;
    }

    __device__ __host__ int FractureVerticesGenerationFunctionRegularTriangleUniformRadius(const double *rand_0_1, const double *Para, double3 *vertices)
    {
        // the radius of the
        double R = RandomUniform(rand_0_1[0], Para);

        double angles[3] = {0, 2 * M_PI / 3, 4 * M_PI / 3};

        for (int i = 0; i < 3; ++i)
        {
            vertices[i].x = R * cos(angles[i]); // x-coordinate
            vertices[i].y = R * sin(angles[i]); // y-coordinate
            vertices[i].z = 0;
        }

        return 3; // return number of vertices
    }

}