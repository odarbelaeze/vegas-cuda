#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "system.h"

spin_t * dipolar_field(system_t * sys, unsigned int site)
{
    spin_t * field = (spin_t *) malloc (sizeof(spin_t));
    field->x = 0.0f; field->y = 0.0f; field->z = 0.0f;
    unsigned int behind = (site == 0U) ? 0U : sys->limits[site - 1U];
    for (unsigned int i = behind; i < sys->limits[site]; i++)
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

void mc_integrate(
        system_t * sys, spin_t * pert, unsigned int * sites,
        float * enes, float * mags, unsigned int time)
{
    for (unsigned int i = 0U; i < time; i++)
    {
        for (unsigned int j = 0U; j < sys->n_sites; j++)
        {
            unsigned int site = sites[i];
            spin_t * delta;
            spin_t * field;

            field = dipolar_field (sys, site);
            delta = spin_delta(sys->spins + site, pert + i * sys->n_sites + j);

            float delta_e = exchange_interaction(delta, field, 1.0);

            free(field);
            free(delta);

            if (delta_e < 0.0)
            {
            }
        }
    }
}

int main(int argc, char **argv)
{
    system_t * sys = create_lattice(1024U, 1U);
    gsl_rng * rng = gsl_rng_alloc (gsl_rng_mt19937);

    random_spin_gsl (sys->spins, rng, sys->n_sites);

    printf("%u should be %u\n",
           sys->n_links, sys->limits[sys->n_sites - 1U]);

    free(sys->spins);
    free(sys->limits);
    free(sys->neighbors);
    free(sys);

    return 0;
}
