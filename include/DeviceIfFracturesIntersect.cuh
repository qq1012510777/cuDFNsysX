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
                                               int2 *IntersectionFracturePairs_DEV_ptr,
                                               const bool IfConsiderTruncatedFractures = false)
    {
        int i = threadIdx.x + blockIdx.x * blockDim.x;

        if (i > NumFracPairs - 1)
            return;

        int x_ = floor((pow(2 * (i + 1), 0.5) + 1 / 2.0));
        int y_ = i - 0.5 * x_ * (x_ - 1);

        // printf("(%d, %d)\n", x_, y_);

        int NumV1 = 0;
        int NumV2 = 0;

        if (!IfConsiderTruncatedFractures)
        {
            // printf("OOO\n");
            NumV1 = VerticesIndex_DEV_ptr[x_].y - VerticesIndex_DEV_ptr[x_].x + 1,
            NumV2 = VerticesIndex_DEV_ptr[y_].y - VerticesIndex_DEV_ptr[y_].x + 1;
        }
        else
        {
            // printf("TTT\n");
            if (VerticesIndex_DEV_ptr[x_].w != -1 && VerticesIndex_DEV_ptr[y_].w != -1)
                NumV1 = VerticesIndex_DEV_ptr[x_].w - VerticesIndex_DEV_ptr[x_].z + 1,
                NumV2 = VerticesIndex_DEV_ptr[y_].w - VerticesIndex_DEV_ptr[y_].z + 1;
            else
                return;
        }
        printf("(%d, %d)\n", NumV1, NumV2);

        int B1 = (!IfConsiderTruncatedFractures ? VerticesIndex_DEV_ptr[x_].x : VerticesIndex_DEV_ptr[x_].z);
        double3 V1[50];
        for (int j = 0; j < NumV1; ++j)
            V1[j] = FractureVertices_DEV_ptr[B1 + j];

        int B2 = (!IfConsiderTruncatedFractures ? VerticesIndex_DEV_ptr[y_].x : VerticesIndex_DEV_ptr[y_].z);
        double3 V2[50];
        for (int j = 0; j < NumV2; ++j)
            V2[j] = FractureVertices_DEV_ptr[B2 + j];

        if (NumV1 > 50 || NumV2 > 50)
            printf("In DeviceIfFracturesIntersect, the number of vertices is too high  NumV1 = %d, NumV2 = %d\n", NumV1, NumV2);

        // for (int j = 0; j < NumV1; ++j)
        //     printf("(%f, %f, %f) ", (V1[j]).x, (V1[j]).y, (V1[j]).z);
        // printf("\n");
        // for (int j = 0; j < NumV2; ++j)
        //     printf("(%f, %f, %f) ", (V2[j]).x, (V2[j]).y, (V2[j]).z);
        // printf("\n");

        if (PolygonsIntersect(V1, V2, NumV1, NumV2))
        {
            // printf("Intersecting pair: %d, %d\n", x_, y_);
            IntersectionFracturePairs_DEV_ptr[i] = make_int2(x_, y_);
        }
    };
};
