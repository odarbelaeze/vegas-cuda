#include "system.h"
#include "mc.h"


int main(int argc, char **argv)
{
    gsl_rng * rng = gsl_rng_alloc (gsl_rng_mt19937);

    uint size = 1024U;
    uint steps = 100U;

    system_t * sys = create_lattice(size, 1U);
    random_spin_gsl (sys->spins, rng, sys->n_sites);

    float * energies;
    float * magnetizations;

    energies = (float *) malloc (steps * sizeof(float));
    magnetizations = (float *) malloc (steps * sizeof(float));

    uint n_attempts = sys->n_sites * steps;

    perturbation_t * perturbations;
    perturbations = (perturbation_t *) malloc (n_attempts * sizeof(perturbation_t));
    random_preturbations_gsl(perturbations, rng, sys->n_sites, n_attempts);

    mc_integrate(sys, perturbations, steps, 0.5, energies, magnetizations);

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

    destroy_system(sys);

    return 0;
}
