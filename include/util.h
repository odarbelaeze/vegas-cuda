#ifndef VEGAS_UTIL_H_
#define VEGAS_UTIL_H_

#include <math.h>
#include <stdlib.h>

#include <gsl/gsl_rng.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

float rand_gaussian();
float norm_float3(const float x, const float y, const float z);
void index_to_xy (unsigned int, unsigned int, unsigned int *, unsigned int *);
void xy_to_index (unsigned int *, unsigned int, unsigned int, unsigned int);
int compare_uint(const void *, const void *);

#endif
