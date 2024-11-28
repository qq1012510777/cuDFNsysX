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

    __host__ __device__ void CalculatePlaneCoefficients(
        const double3 &P1, const double3 &P2, const double3 &P3,
        double &A, double &B, double &C, double &D)
    {
        // Compute vectors
        double3 V1 = P2 - P1;
        double3 V2 = P3 - P1;
        // V1 = V1 / Double3Norm(V1);
        // V2 = V2 / Double3Norm(V2);

        // Compute the normal vector (A, B, C)
        double3 N = Double3CrossProduct(V1, V2);
        // N = N / Double3Norm(N);

        A = N.x;
        B = N.y;
        C = N.z;

        // Compute D using P1
        D = -(A * P1.x + B * P1.y + C * P1.z);
    }

    __host__ __device__ bool TrianglesIntersect(const double3 *triangle1, const double3 *triangle2)
    {
        for (int i = 0; i < 3; ++i)
            printf("%.40f, %.40f, %.40f\n", triangle1[i].x, triangle1[i].y, triangle1[i].z);
        for (int i = 0; i < 3; ++i)
            printf("%.40f, %.40f, %.40f\n", triangle2[i].x, triangle2[i].y, triangle2[i].z);
        printf("-------------------------\n\n");
        double3 e1 = triangle1[1] - triangle1[0]; // Edge 1 of triangle1
        double3 e2 = triangle1[2] - triangle1[0]; // Edge 2 of triangle1
        double3 n1 = Double3CrossProduct(e1, e2); // Normal vector of triangle1

        double3 e3 = triangle2[1] - triangle2[0]; // Edge 1 of triangle2
        double3 e4 = triangle2[2] - triangle2[0]; // Edge 2 of triangle2
        double3 n2 = Double3CrossProduct(e3, e4); // Normal vector of triangle2

        // Plane equation tests
        double d1 = Double3Dot(n1, triangle1[0]);
        double d2 = Double3Dot(n2, triangle2[0]);

        // Test triangle1 vertices against triangle2's plane
        double dist0 = Double3Dot(n2, triangle1[0]) - d2;
        double dist1 = Double3Dot(n2, triangle1[1]) - d2;
        double dist2 = Double3Dot(n2, triangle1[2]) - d2;

        if (dist0 * dist1 > 0.0 && dist0 * dist2 > 0.0)
            return false;

        // Test triangle2 vertices against triangle1's plane
        double dist3 = Double3Dot(n1, triangle2[0]) - d1;
        double dist4 = Double3Dot(n1, triangle2[1]) - d1;
        double dist5 = Double3Dot(n1, triangle2[2]) - d1;

        // printf("%.20f, %.20f, %20f ; %.20f, %.20f, %20f\n", dist0,
        //        dist1,
        //        dist2,
        //        dist3,
        //        dist4,
        //        dist5);

        if (dist3 * dist4 > 0.0 && dist3 * dist5 > 0.0)
            return false;
        // Compute intersection line direction
        double3 dir = Double3CrossProduct(n1, n2);
        dir = dir / Double3Norm(dir);

        // Project triangle1 vertices onto the intersection line
        double proj1[3];
        for (int i = 0; i < 3; ++i)
        {
            proj1[i] = Double3Dot(dir, triangle1[i]);
        }

        // Project triangle2 vertices onto the intersection line
        double proj2[3];
        for (int i = 0; i < 3; ++i)
        {
            proj2[i] = Double3Dot(dir, triangle2[i]);
        }

        // Find projection intervals
        double min1 = fmin(fmin(proj1[0], proj1[1]), proj1[2]);
        double max1 = fmax(fmax(proj1[0], proj1[1]), proj1[2]);

        double min2 = fmin(fmin(proj2[0], proj2[1]), proj2[2]);
        double max2 = fmax(fmax(proj2[0], proj2[1]), proj2[2]);

        // Check if projection intervals overlap
        return (fmax(min1, max1) >= fmin(min2, max2) && fmax(min2, max2) >= fmin(min1, max1));
    }

    __host__ __device__ bool ArePointsCollinear(const double3 &A, const double3 &B, const double3 &C, double epsilon = 1e-8)
    {
        double3 AB = B - A;
        double3 AC = C - A;
        double3 cross = Double3CrossProduct(AB, AC);
        return Double3Norm(cross) < epsilon; // Check if the cross product's magnitude is near zero
    }

    // Main Polygon-Polygon Intersection Test
    // cannot guarantee correctness for concave polygons
    __host__ __device__ bool PolygonsIntersect(
        const double3 *Vertices_Polygon1, const double3 *Vertices_Polygon2,
        const int NumVertices_Polygon1, const int NumVertices_Polygon2)
    {
        double CollinearEPSILON = 1e-3;
        // Iterate over all triangle pairs from the two polygons
        for (int i = 0; i < NumVertices_Polygon1 - 2; ++i)
        {
            if (ArePointsCollinear(Vertices_Polygon1[0], Vertices_Polygon1[i + 1], Vertices_Polygon1[i + 2], CollinearEPSILON))
                continue;

            const double3 V1[3] = {Vertices_Polygon1[0], Vertices_Polygon1[i + 1], Vertices_Polygon1[i + 2]};

            for (int j = 0; j < NumVertices_Polygon2 - 2; ++j)
            {

                if (ArePointsCollinear(Vertices_Polygon2[0], Vertices_Polygon2[j + 1], Vertices_Polygon2[j + 2], CollinearEPSILON))
                    continue;

                const double3 V2[3] = {Vertices_Polygon2[0], Vertices_Polygon2[j + 1], Vertices_Polygon2[j + 2]};
                // printf("%d -  %d\n", i, j);

                // double3 V1[3] = {make_double3(1.000000000000000000000000000000, 0.443320452273296028433691162718, 0.000000000000000000000000000000),
                //                  make_double3(0.669065646357200338734116940032, 1.000000000000000000000000000000, 0.330934353642799661265883059968),
                //                  make_double3(0.316041686360868157024128777266, 1.000000000000000000000000000000, 0.000000000000000000000000000000)};
                // double3 V2[3] = {
                //     make_double3(0.690740882962605162731506425189, 0.690740882962605162731506425189, 0.000000000000000000000000000000),
                //     make_double3(0.755442582729919487327663318865, 1.000000000000000000000000000000, 0.000000000000000000000000000000),
                //     make_double3(0.000000000000000000000000000000, 1.000000000000000000000000000000, 0.684407326723559861214596367063)};
                // printf("~%d\n", cuDFNsys::TrianglesIntersect(
                //                     V1, V2));
                if (cuDFNsys::TrianglesIntersect(
                        V1, V2))
                {
                    printf("Intsected~~~~~\n\n");
                    return true;
                } // Intersection found
                else
                {
                    printf("Not !! Intsected~~~~~\n\n");
                }
            }
            // printf("\n\n\n\n");
        }
        return false; // No intersection
    }
};