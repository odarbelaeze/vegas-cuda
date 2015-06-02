#ifndef MC_H_
#define MC_H_

#include "system.h"

struct _perturbation {
    uint site;
    spin_t pert;
    float prob;
};

typedef struct _perturbation perturbation_t;

void mc_integrate (system_t * sys, perturbation_t * perts,
                   uint time, float temp,
                   float * enes, float * mags);

void random_preturbations_gsl (
        perturbation_t * perturbations,
        gsl_rng * rng,
        uint n_sites,
        uint count);

#endif
