#ifndef SYSTEM_H_
#define SYSTEM_H_

#include <math.h>
#include <stdio.h>
#include <string.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include "util.h"

struct _spin {
    float x;
    float y;
    float z;
};

typedef struct _spin spin_t;

struct _system {
    spin_t * spins;
    unsigned int n_sites;
    unsigned int n_links;
    unsigned int * limits;
    unsigned int * neighbors;
};

typedef struct _system system_t;

float exchange_interaction(const spin_t *, const spin_t *, float);
float spin_norm(const spin_t * const);
void random_spins(spin_t * const, unsigned int);
void random_spin_gsl (spin_t * const, const gsl_rng *, unsigned int);
system_t * create_lattice(unsigned int, unsigned short);
void print_spins_csv(const system_t * const);
void dipolar_field(
        spin_t * const __restrict,
        const system_t * const __restrict,
        uint);
void spin_delta(
        const spin_t * const __restrict,
        const spin_t * const __restrict,
        spin_t * const __restrict);
float compute_magnetization(const system_t * const);
float compute_energy(const system_t * const);

#endif
