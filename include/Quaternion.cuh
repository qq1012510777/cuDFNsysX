#pragma once
#include <cuda_runtime.h>
#include "Operators.cuh"

namespace cuDFNsys
{
    struct Quaternion
    {
    protected:
        // 1,i,j,k
        double4 QuaternionNum;

    public:
        // describe quaternion
        __device__ __host__ Quaternion DescribeRotation(const double3 v, const double angle)
        {
            // v = v / Double3Norm(v);
            double sina_2 = sin(angle * 0.5000);
            double cosa_2 = cos(angle * 0.5000);
            Quaternion result;
            result.QuaternionNum.x = cosa_2,
            result.QuaternionNum.y = sina_2 * v.x,
            result.QuaternionNum.z = sina_2 * v.y,
            result.QuaternionNum.w = sina_2 * v.z;
            double norm_r = pow(result.QuaternionNum.x * result.QuaternionNum.x + result.QuaternionNum.y * result.QuaternionNum.y +
                                    result.QuaternionNum.z * result.QuaternionNum.z + result.QuaternionNum.w * result.QuaternionNum.w,
                                0.5);
            result.QuaternionNum.x /= norm_r;
            result.QuaternionNum.y /= norm_r;
            result.QuaternionNum.z /= norm_r;
            result.QuaternionNum.w /= norm_r;
            // printf("tt %f, %f, %f, %f\n", result.QuaternionNum.w, result.QuaternionNum.x,
            //        result.QuaternionNum.y, result.QuaternionNum.z);

            return result;
        }; // Quaternion::DescribeRotation;

        // rotate
        __device__ __host__ double3 Rotate(const double3 v)
        {
            // printf("* %f, %f, %f, %f\n", QuaternionNum.w, QuaternionNum.x, QuaternionNum.y, QuaternionNum.z);

            double t2 = QuaternionNum.x * QuaternionNum.y,
                   t3 = QuaternionNum.x * QuaternionNum.z,
                   t4 = QuaternionNum.x * QuaternionNum.w,
                   t5 = -QuaternionNum.y * QuaternionNum.y,
                   t6 = QuaternionNum.y * QuaternionNum.z,
                   t7 = QuaternionNum.y * QuaternionNum.w,
                   t8 = -QuaternionNum.z * QuaternionNum.z,
                   t9 = QuaternionNum.z * QuaternionNum.w,
                   t10 = -QuaternionNum.w * QuaternionNum.w;

            double3 DF;

            DF.x = 2.0 * ((t8 + t10) * v.x + (t6 - t4) * v.y + (t3 + t7) * v.z) + v.x,
            DF.y = 2.0 * ((t4 + t6) * v.x + (t5 + t10) * v.y + (t9 - t2) * v.z) + v.y,
            DF.z = 2.0 * ((t7 - t3) * v.x + (t2 + t9) * v.y + (t5 + t8) * v.z) + v.z;

            return DF;
        };

        // get QuaternionNum
        __device__ __host__ double4 GetQuaternionNum()
        {
            double4 f = QuaternionNum;
            return f;
        };

        __device__ __host__ void SetQuaternionNum(const double4 &f)
        {
            this->QuaternionNum = f;
        };

        // rotation back
        __device__ __host__ Quaternion InvertRotationAxis()
        {
            Quaternion result;
            result.QuaternionNum.x = QuaternionNum.x;
            result.QuaternionNum.y = QuaternionNum.y * -1.;
            result.QuaternionNum.z = QuaternionNum.z * -1.;
            result.QuaternionNum.w = QuaternionNum.w * -1.;
            return result;
        }

        // rotate a 3D triangle to 2D
        //__device__ __host__ Quaternion Rotate3DTriangleto2DTriangle(double3 &A, double3 &B, double3 &C)
        //{
        //    double3 U = C - A, V = B - A;
        //    U = U / NormDouble3(U);
        //    V = V / NormDouble3(V);
        //    double3 N_f = CrossDouble3(U, V);
        //    N_f = N_f / NormDouble3(N_f);
        //    if (N_f.z < 0)
        //        N_f.x *= -1, N_f.y *= -1, N_f.z *= -1;
        //};
    };
};