#pragma once
#include <thrust/host_vector.h>
#include <cmath>
#include <algorithm>
#include <vector>
#include "Operators.cuh"

__host__ __device__ bool IsEqual(const double3 &a, const double3 &b, double tolerance = 0)
{
    return fabs(a.x - b.x) < tolerance &&
           fabs(a.y - b.y) < tolerance &&
           fabs(a.z - b.z) < tolerance;
}

__host__ double3 ComputeCentroid(const thrust::host_vector<double3> &points)
{
    double3 centroid = make_double3(0.0, 0.0, 0.0);
    for (const auto &p : points)
    {
        centroid = centroid + p;
    }
    return centroid / static_cast<double>(points.size());
}

__host__ thrust::host_vector<double3> CorrectPolygon(const thrust::host_vector<double3> &points, double tolerance = 0)
{
    // Step 1: Remove duplicates
    thrust::host_vector<double3> uniquePoints;
    for (const auto &p : points)
    {
        if (std::none_of(uniquePoints.begin(), uniquePoints.end(),
                         [&](const double3 &q)
                         { return IsEqual(p, q, tolerance); }))
        {
            uniquePoints.push_back(p);
        }
    }

    // Step 2: Compute centroid
    double3 centroid = ComputeCentroid(uniquePoints);

    // Step 3: Sort vertices based on their angular order
    std::sort(uniquePoints.begin(), uniquePoints.end(),
              [&](const double3 &a, const double3 &b)
              {
                  double3 vecA = a - centroid;
                  double3 vecB = b - centroid;
                  double crossZ = Double3CrossProduct(vecA, vecB).z;
                  if (fabs(crossZ) < tolerance)
                  {
                      // Collinear points, sort by distance to centroid
                      return Double3Norm(vecA) < Double3Norm(vecB);
                  }
                  return crossZ > 0; // Counter-clockwise order
              });

    return uniquePoints;
}
