#pragma once
#include <iostream>
#include <cuda_runtime.h>
#include "Operators.cuh"

namespace cuDFNsys
{

    // __host__ __device__ bool IsPointInTriangle(const double3 &p, const double3 &a, const double3 &b, const double3 &c)
    // {
    //     double3 ab = b - a;
    //     double3 ac = c - a;
    //     double3 ap = p - a;
    //     double3 cross1 = Double3CrossProduct(ab, ap);
    //     double3 cross2 = Double3CrossProduct(ap, ac);
    //     double3 cross3 = Double3CrossProduct(ac, ab);
    //     double dot1 = Double3Dot(cross1, cross3);
    //     double dot2 = Double3Dot(cross2, cross3);
    //     printf("%f, %f, %f, %f\n", dot1, dot2, dot1+dot2, Double3Dot(cross3, cross3));
    //     return (dot1 >= 0 && dot2 >= 0 && dot1 + dot2 <= Double3Dot(cross3, cross3));
    // }
    // __host__ __device__ bool IsPointInTriangle(const double3 &p, const double3 *Triangle)
    // {
    //     return IsPointInTriangle(p, Triangle[0], Triangle[1], Triangle[2]);
    // }

    __host__ __device__ bool RayIntersectsTriangle(const double3 &rayOrigin, const double3 &rayDir,
                                                   const double3 &v0, const double3 &v1, const double3 &v2, double &t)
    {
        const double EPSILON = 1e-8;
        double3 edge1 = v1 - v0;
        double3 edge2 = v2 - v0;
        double3 h = Double3CrossProduct(rayDir, edge2);
        double a = Double3Dot(edge1, h);

        if (fabs(a) < EPSILON) // Parallel ray
            return false;

        double f = 1.0 / a;
        double3 s = rayOrigin - v0;
        double u = f * Double3Dot(s, h);
        if (u < 0.0 || u > 1.0)
            return false;

        double3 q = Double3CrossProduct(s, edge1);
        double v = f * Double3Dot(rayDir, q);
        if (v < 0.0 || u + v > 1.0)
            return false;

        t = f * Double3Dot(edge2, q);
        return t > EPSILON; // Intersection detected
    };

    __host__ __device__ bool TrianglesIntersect(const double3 *triangle1, const double3 *triangle2)
    {
        // Unpack triangle vertices
        const double3 &p1 = triangle1[0];
        const double3 &q1 = triangle1[1];
        const double3 &r1 = triangle1[2];

        const double3 &p2 = triangle2[0];
        const double3 &q2 = triangle2[1];
        const double3 &r2 = triangle2[2];

        // Compute triangle normals
        double3 normal1 = Double3CrossProduct(q1 - p1, r1 - p1);
        double3 normal2 = Double3CrossProduct(q2 - p2, r2 - p2);

        // Compute signed distances of triangle2 vertices relative to triangle1's plane
        double dp2 = Double3Dot(normal1, p2 - p1);
        double dq2 = Double3Dot(normal1, q2 - p1);
        double dr2 = Double3Dot(normal1, r2 - p1);

        // If all points of triangle2 are on the same side of triangle1's plane, no intersection
        if ((dp2 > 0 && dq2 > 0 && dr2 > 0) || (dp2 < 0 && dq2 < 0 && dr2 < 0))
            return false;

        // Compute signed distances of triangle1 vertices relative to triangle2's plane
        double dp1 = Double3Dot(normal2, p1 - p2);
        double dq1 = Double3Dot(normal2, q1 - p2);
        double dr1 = Double3Dot(normal2, r1 - p2);

        // If all points of triangle1 are on the same side of triangle2's plane, no intersection
        if ((dp1 > 0 && dq1 > 0 && dr1 > 0) || (dp1 < 0 && dq1 < 0 && dr1 < 0))
            return false;

        // Edge-edge intersection tests
        const double3 edges1[] = {q1 - p1, r1 - q1, p1 - r1};
        const double3 edges2[] = {q2 - p2, r2 - q2, p2 - r2};

        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                double3 axis = Double3CrossProduct(edges1[i], edges2[j]);

                // Skip degenerate axes
                if (Double3Norm(axis) < 1e-6)
                    continue;

                double min1 = INFINITY, max1 = -INFINITY;
                double min2 = INFINITY, max2 = -INFINITY;

                // Project triangle1 vertices onto the axis
                for (const double3 &v : {p1, q1, r1})
                {
                    double proj = Double3Dot(v, axis);
                    min1 = fmin(min1, proj);
                    max1 = fmax(max1, proj);
                }

                // Project triangle2 vertices onto the axis
                for (const double3 &v : {p2, q2, r2})
                {
                    double proj = Double3Dot(v, axis);
                    min2 = fmin(min2, proj);
                    max2 = fmax(max2, proj);
                }

                // Check if projections overlap
                if (max1 < min2 || max2 < min1)
                    return false;
            }
        }

        // All tests passed; triangles intersect
        return true;
    }

    // Main Polygon-Polygon Intersection Test
    // cannot guarantee correctness for concave polygons
    __host__ __device__ bool PolygonsIntersect(
        const double3 *Vertices_Polygon1, const double3 *Vertices_Polygon2,
        const int NumVertices_Polygon1, const int NumVertices_Polygon2)
    {
        // Iterate over all triangle pairs from the two polygons
        for (int i = 0; i < NumVertices_Polygon1 - 2; ++i)
        {
            for (int j = 0; j < NumVertices_Polygon2 - 2; ++j)
            {
                const double3 *V1[3] = {&Vertices_Polygon1[0], &Vertices_Polygon1[i + 1], &Vertices_Polygon1[i + 2]};
                const double3 *V2[3] = {&Vertices_Polygon2[0], &Vertices_Polygon2[j + 1], &Vertices_Polygon2[j + 2]};

                // double3 V1[3] = {make_double3(-11.194467308199275, 30.016399834222312, 88.73795816690222),
                //                  make_double3(2.424293150903748, 8.642098266568603, 80.16985476223483),
                //                  make_double3(10.238139329403761, 21.4950350078845, 102.29426718207047)};
                // double3 V2[3] = {
                //     make_double3(31.080304930514355, 96.77713966099772, 29.648916007562498),
                //     make_double3(19.31127382807052, 91.81921927419917, 6.401495649195706),
                //     make_double3(45.21811656200058, 87.01105530631769, 9.442988668768649)};
                if (cuDFNsys::TrianglesIntersect(
                        *V1, *V2))
                    return true; // Intersection found
            }
        }
        return false; // No intersection
    }
};