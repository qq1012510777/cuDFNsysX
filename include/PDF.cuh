#pragma once
#include <cuda_runtime.h>

namespace cuDFNsys
{
    __device__ __host__ double RandomPowerlaw(const double &rand_0_1, const double *Para)
    {
        double x0 = Para[0];
        double x1 = Para[1];
        double alpha_g = Para[2];

        double x_g = (pow(x1, 1.0f - alpha_g) - pow(x0, 1.0f - alpha_g)) * rand_0_1 + pow(x0, 1.0f - alpha_g);
        x_g = pow(x_g, 1.0f / (1.0f - alpha_g));
        return x_g;
    }; // RandomPowerlaw

    __device__ __host__ double RandomUniform(const double &rand_0_1, const double *Para)
    {
        double a = Para[0];
        double b = Para[1];
        double random = a + (b - a) * rand_0_1;
        return random;
    }; // RandomUniform

};