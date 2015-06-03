#ifndef MC_H_
#define MC_H_

#include "system.h"

struct _perturbation {
    uint site;
    spin_t pert;
    float prob;
};

typedef struct _perturbation perturbation_t;

void mc_integrate (
        const system_t * const sys,
        const perturbation_t * const perts,
        uint time, float temp,
        float * const __restrict enes,
        float * const __restrict mags);

void random_preturbations_gsl (
        perturbation_t * const __restrict perturbations,
        const gsl_rng * __restrict rng,
        uint n_sites,
        uint count);

#endif
