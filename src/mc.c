#include "mc.h"


void mc_integrate (
        const system_t * const __restrict sys,
        perturbation_t const * const __restrict perts,
        uint time, float temp,
        float * const __restrict enes,
        float * const __restrict mags)
{
    spin_t field;
    spin_t delta;

    for (uint i = 0U; i < time; i++)
    {
        for (uint j = 0U; j < sys->n_sites; j++)
        {
            uint prog = i * sys->n_sites + j;
            uint site = perts[prog].site;

            dipolar_field (&field, sys, site);
            spin_delta(sys->spins + site, &(perts[prog].pert), &delta);

            float delta_e = exchange_interaction(&delta, &field, 1.0f);


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


void random_preturbations_gsl (
        perturbation_t * const __restrict perturbations,
        const gsl_rng * __restrict rng,
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
