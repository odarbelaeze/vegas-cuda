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

void index_to_xy (unsigned int idx, unsigned int width, unsigned int * x, unsigned int * y)
{
    *x = idx % width;
    *y = idx / width;
}

void xy_to_index (unsigned int * idx, unsigned int width, unsigned int x, unsigned int y)
{
    *idx = y * width + x;
}

int compare_uint(const void * a, const void * b)
{
    return * ((unsigned int *) a) - * ((unsigned int *) b);
}

