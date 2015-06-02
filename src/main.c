#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "system.h"

struct _perturbation {
    uint site;
    spin_t pert;
    float prob;
};

typedef struct _perturbation perturbation_t;

void random_preturbations_gsl (
        perturbation_t * perturbations,
        gsl_rng * rng,
        uint n_sites,
        uint count)
{
    for (uint i = 0U; i < count; i++)
    {
        perturbations[i].site = gsl_rng_uniform_int(rng, n_sites);
        random_spin_gsl(&(perturbations[i].pert), rng, 1);
        perturbations[i].prob = gsl_rng_uniform(rng);
    }
}

spin_t * dipolar_field(system_t * sys, uint site)
{
    spin_t * field = (spin_t *) malloc (sizeof(spin_t));
    field->x = 0.0f; field->y = 0.0f; field->z = 0.0f;
    uint behind = (site == 0U) ? 0U : sys->limits[site - 1U];
    for (uint i = behind; i < sys->limits[site]; i++)
    {
        field->x += sys->spins[sys->neighbors[i]].x;
        field->y += sys->spins[sys->neighbors[i]].y;
        field->z += sys->spins[sys->neighbors[i]].z;
    }
    return field;
}

spin_t * spin_delta(spin_t * a, spin_t * b)
{
    spin_t * delta = (spin_t *) malloc (sizeof(spin_t));
    delta->x = b->x - a->x;
    delta->y = b->y - a->y;
    delta->z = b->z - a->z;
    return delta;
}

float compute_magnetization(system_t * sys)
{
    double acum_x = 0.0;
    double acum_y = 0.0;
    double acum_z = 0.0;

    for (uint i = 0U; i < sys->n_sites; i++)
    {
        acum_x += sys->spins[i].x;
        acum_y += sys->spins[i].y;
        acum_z += sys->spins[i].z;
    }

    return norm_float3(acum_x, acum_y, acum_z);
}

float compute_energy(system_t * sys)
{
    double ene_accum = 0.0;

    for (uint i = 0U; i < sys->n_sites; i++)
    {
        spin_t * field = dipolar_field (sys, i);
        ene_accum = exchange_interaction(sys->spins + 1, field, 1.0f);
        free (field);
    }

    return (float) (0.5 * ene_accum);
}

void mc_integrate(
        system_t * sys, perturbation_t * perts, float temp,
        float * enes, float * mags, uint time)
{
    for (uint i = 0U; i < time; i++)
    {
        for (uint j = 0U; j < sys->n_sites; j++)
        {
            uint prog = i * sys->n_sites + j;
            uint site = perts[prog].site;
            spin_t * delta;
            spin_t * field;

            field = dipolar_field (sys, site);
            delta = spin_delta(sys->spins + site, &(perts[prog].pert));

            float delta_e = exchange_interaction(delta, field, 1.0f);

            free(field);
            free(delta);

            if (delta_e < 0.0)
            {
                memcpy ((void *) &(perts[prog].pert),
                        (void *) sys->spins + site,
                        sizeof(spin_t));
            }
            else
            {
                if (perts[prog].prob < exp(- delta_e / temp))
                {
                    memcpy ((void *) &(perts[prog].pert),
                            (void *) sys->spins + site,
                            sizeof(spin_t));
                }
            }
        }

        mags[i] = compute_magnetization(sys);
        enes[i] = compute_energy(sys);
    }
}

int main(int argc, char **argv)
{
    gsl_rng * rng = gsl_rng_alloc (gsl_rng_mt19937);

    uint size = 512U;
    uint steps = 100U;

    system_t * sys = create_lattice(256U, 1U);
    random_spin_gsl (sys->spins, rng, sys->n_sites);

    float * energies;
    float * magnetizations;

    energies = (float *) malloc (steps * sizeof(float));
    magnetizations = (float *) malloc (steps * sizeof(float));

    uint n_attempts = sys->n_sites * steps;

    perturbation_t * perturbations;
    perturbations = (perturbation_t *) malloc (n_attempts * sizeof(perturbation_t));
    random_preturbations_gsl(perturbations, rng, sys->n_sites, n_attempts);

    mc_integrate(sys, perturbations, 0.5, energies, magnetizations, steps);

    free(energies);
    free(magnetizations);
    free(perturbations);

    for (uint i = 1U; i < sys->n_sites; i++)
    {
        uint behind = sys->limits[i - 1U];
        uint front = sys->limits[i];
        if (front < behind) {
            printf("**ERROR** %u should be greater than %u", front, behind);
        }
    }

    for (uint i = 0U; i < sys->n_links; i++)
    {
        unsigned short nbhidx = sys->neighbors[i];
        if (nbhidx >= sys->n_sites)
        {
            printf("**ERROR** %u should be no greater than %u", nbhidx, sys->n_sites);
        }
    }

    printf ("%u should be %u\n",
            sys->n_links, sys->limits[sys->n_sites - 1U]);

    gsl_rng_free (rng);

    free(sys->spins);
    free(sys->limits);
    free(sys->neighbors);
    free(sys);

    return 0;
}
