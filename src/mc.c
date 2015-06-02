#include "mc.h"


void mc_integrate (system_t * sys, perturbation_t * perts,
                   uint time, float temp,
                   float * enes, float * mags)
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
