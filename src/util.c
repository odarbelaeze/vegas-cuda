#include <gsl/gsl_rng.h>

#include "util.h"


float rand_gaussian ()
{
    float u = (float) rand() / (float) RAND_MAX;
    float v = (float) rand() / (float) RAND_MAX;
    return sqrt(- 2.0f * log(1 - u)) * cos(2.0f * M_PI * v);
}


float norm_float3 (const float x, const float y, const float z)
{
    return sqrt(x * x + y * y + z * z);
}


void index_to_xy (uint idx, uint width, uint * x, uint * y)
{
    *x = idx % width;
    *y = idx / width;
}


void xy_to_index (uint * idx, uint width, uint x, uint y)
{
    *idx = y * width + x;
}


int compare_uint(const void * a, const void * b)
{
    return * ((uint *) a) - * ((uint *) b);
}


void rot(uint width, uint * x, uint * y, uint rx, uint ry)
{
    if (ry == 0U)
    {
        if (rx == 1U)
        {
            *x = width - 1U - *x;
            *y = width - 1U - *y;
        }
        uint t = *x;
        *x = *y;
        *y = t;
    }
}


uint xy_to_d (uint width, uint x, uint y)
{
    uint rx;
    uint ry;
    uint d = 0U;
    for (uint s = width / 2U; s > 0U; s /= 2U)
    {
        rx = (x & s) > 0U;
        ry = (y & s) > 0U;
        d += s * s * ((3U * rx) ^ ry);
        rot(s, &x, &y, rx, ry);
    }
    return d;
}


void d_to_xy (uint width, uint d, uint * x, uint * y)
{
    *x = *y = 0U;
    uint t = d;
    uint rx;
    uint ry;
    for (uint s = 1U; s < width; s *= 2U)
    {
        rx = 1U & (t / 2U);
        ry = 1U & (t ^ rx);
        rot(s, x, y, rx, ry);
        *x += s * rx;
        *y += s * ry;
        t /= 4U;
    }
}

