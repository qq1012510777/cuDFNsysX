#pragma once
#include <iostream>
#include <cuda_runtime.h>

// double3
std::ostream &operator<<(std::ostream &os, const double3 &v)
{
    os << v.x << ", " << v.y << ", " << v.z << "; ";
    return os;
}
__host__ __device__ double3 Double3CrossProduct(const double3 &u, const double3 &v)
{
    return make_double3(
        u.y * v.z - u.z * v.y,
        u.z * v.x - u.x * v.z,
        u.x * v.y - u.y * v.x);
}
__host__ __device__ double Double3Norm(const double3 &v)
{
    return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}
__host__ __device__ double Double3Dot(const double3 &v, const double3 &w)
{
    return (v.x * w.x + v.y * w.y + v.z * w.z);
}
__host__ __device__ double3 operator/(const double3 &vec, double scalar)
{
    return make_double3(vec.x / scalar, vec.y / scalar, vec.z / scalar);
}
__host__ __device__ double3 operator*(const double3 &vec, double scalar)
{
    return make_double3(vec.x * scalar, vec.y * scalar, vec.z * scalar);
}
__host__ __device__ double3 operator-(const double3 &vec, const double3 &vec2)
{
    return make_double3(vec.x - vec2.x, vec.y - vec2.y, vec.z - vec2.z);
}
__host__ __device__ double3 operator+(const double3 &vec, const double3 &vec2)
{
    return make_double3(vec.x + vec2.x, vec.y + vec2.y, vec.z + vec2.z);
}

// int2
std::ostream &operator<<(std::ostream &os, const int2 &v)
{
    os << "(" << v.x << ", " << v.y << ")";
    return os;
}