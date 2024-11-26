#pragma once
#include <sys/time.h> // For gettimeofday
#include <iostream>

namespace cuDFNsys
{
    double CPUSecond()
    {
        struct timeval tp;
        gettimeofday(&tp, NULL);
        return ((double)tp.tv_sec + (double)tp.tv_usec * 1e-6);
    };
};