#ifndef VEGAS_UTIL_H_
#define VEGAS_UTIL_H_

#include <math.h>
#include <stdlib.h>

#include <gsl/gsl_rng.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef unsigned int uint;

float rand_gaussian();
float norm_float3(const float x, const float y, const float z);
void index_to_xy (uint, uint, uint *, uint *);
void xy_to_index (uint *, uint, uint, uint);
int compare_uint(const void *, const void *);
void rot(uint width, uint * x, uint * y, uint rx, uint ry);
uint xy_to_d (uint width, uint x, uint y);
void d_to_xy (uint width, uint d, uint * x, uint * y);

#endif
