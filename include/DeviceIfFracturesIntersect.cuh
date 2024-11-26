#pragma once

#include <cuda_runtime.h>
#include <math.h>
#include "Operators.cuh"
#include "GeometryAlgorithm.cuh"

namespace cuDFNsys
{
    __global__ void DeviceIfFracturesIntersect(const int NumFracPairs,
                                               const double3 *FractureVertices_DEV_ptr,
                                               const int4 *VerticesIndex_DEV_ptr,
                                               int2 *IntersectionFracturePairs_DEV_ptr)
    {
        int i = threadIdx.x + blockIdx.x * blockDim.x;
        if (i > NumFracPairs - 1)
            return;

        int x_ = floor((pow(2 * (i + 1), 0.5) + 1 / 2.0));
        int y_ = i - 0.5 * x_ * (x_ - 1);

        int NumV1 = VerticesIndex_DEV_ptr[x_].y - VerticesIndex_DEV_ptr[x_].x + 1;
        int NumV2 = VerticesIndex_DEV_ptr[y_].y - VerticesIndex_DEV_ptr[y_].x + 1;

        const double3 *V1[100];
        for (int j = 0; j < NumV1; ++j)
            V1[j] = &FractureVertices_DEV_ptr[VerticesIndex_DEV_ptr[x_].x] + j;

        const double3 *V2[100];
        for (int j = 0; j < NumV2; ++j)
            V2[j] = &FractureVertices_DEV_ptr[VerticesIndex_DEV_ptr[y_].x] + j;   

        // for (int j = 0; j < NumV1; ++j)
        //     printf("(%f, %f, %f) ", (*V1[j]).x, (*V1[j]).y, (*V1[j]).z);
        // printf("\n");
        // for (int j = 0; j < NumV2; ++j)
        //     printf("(%f, %f, %f) ", (*V2[j]).x, (*V2[j]).y, (*V2[j]).z);
        // printf("\n");

        if(PolygonsIntersect(*V1, *V2, NumV1, NumV2))
            IntersectionFracturePairs_DEV_ptr[i] = make_int2(x_, y_);
            
    };
};
