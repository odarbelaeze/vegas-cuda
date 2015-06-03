#include "system.h"

float exchange_interaction(const spin_t * a, const spin_t * b, float exchange)
{
    float energy = sqrt(a->x * b->x + a->y * b->y + a->z * b->z);
    return - exchange * energy;
}

float spin_norm(const spin_t * const spin)
{
    return sqrt(spin->x * spin->x + spin->y * spin->y + spin->z * spin->z);
}

void random_spins(spin_t * const spins, uint count)
{
    for (uint i = 0U; i < count; i++)
    {
        spin_t * spin = spins + i;
        spin->x = rand_gaussian();
        spin->y = rand_gaussian();
        spin->z = rand_gaussian();
        float norm = spin_norm(spin);
        spin->x /= norm;
        spin->y /= norm;
        spin->z /= norm;
    }
}

void random_spin_gsl (spin_t * const spins, const gsl_rng * rng, uint count)
{
    for (uint i = 0U; i < count; i++)
    {
        spin_t * spin = spins + i;
        spin->x = gsl_ran_gaussian(rng, 1.0);
        spin->y = gsl_ran_gaussian(rng, 1.0);
        spin->z = gsl_ran_gaussian(rng, 1.0);
        float norm = spin_norm(spin);
        spin->x /= norm;
        spin->y /= norm;
        spin->z /= norm;
    }
}

system_t * create_lattice(uint width, unsigned short pbc)
{
    system_t * sys = (system_t *) malloc (sizeof(system_t));
    sys->n_sites = width * width;
    sys->n_links = sys->n_sites * 4U;

    if (pbc == 0U)
    {
        sys->n_links -= 4U * width;
    }

    sys->spins = (spin_t *) malloc (sys->n_sites * sizeof(spin_t));
    sys->limits = (uint *) malloc (sys->n_sites * sizeof(uint));
    sys->neighbors = (uint *) malloc (sys->n_links * sizeof(uint));

    if (pbc == 0U)
    {
        for (uint i = 0U; i < sys->n_sites; i++)
        {
            uint n_behind = (i != 0U)? sys->limits[i - 1] : 0U;

            uint x, y;
            uint nbh_count = 0U;

            index_to_xy(i, width, &x, &y);

            if (x > 0) sys->neighbors[n_behind + nbh_count++] = i - 1;
            if (x < (width - 1)) sys->neighbors[n_behind + nbh_count++] = i + 1;
            if (y > 0) sys->neighbors[n_behind + nbh_count++] = i - width;
            if (y < (width - 1)) sys->neighbors[n_behind + nbh_count++] = i + width;

            sys->limits[i] = n_behind + nbh_count;
        }
    }
    else
    {
        for (uint i = 0U; i < sys->n_sites; i++)
        {
            uint n_behind = (i != 0U)? sys->limits[i - 1] : 0U;

            uint x, y;

            index_to_xy(i, width, &x, &y);

            sys->neighbors[n_behind + 0U] =
                (x > 0)? i - 1 : i + width - 1;
            sys->neighbors[n_behind + 1U] =
                (x < (width - 1))? i + 1 : i - width + 1;
            sys->neighbors[n_behind + 2U] =
                (y > 0)? i - width : i + (width * (width - 1U));
            sys->neighbors[n_behind + 3U] =
                (y < (width - 1))? i + width : i - (width * (width - 1U));

            sys->limits[i] = n_behind + 4U;
        }
    }


    uint * keys = (uint *) malloc (sys->n_sites * sizeof(uint));
    uint * nw_lims = (uint *) malloc (sys->n_sites * sizeof(uint));
    uint * nw_neig = (uint *) malloc (sys->n_links * sizeof(uint));

    for (uint i = 0U; i < sys->n_sites; i++)
    {
        uint * x = (uint *) malloc (sizeof(uint));
        uint * y = (uint *) malloc (sizeof(uint));
        index_to_xy(i, width, x, y);
        keys[i] = xy_to_d(width, *x, *y);
        uint n_behind = (i != 0U)? sys->limits[i - 1] : 0U;
        uint n_neighs = sys->limits[i] - n_behind;
        nw_lims[keys[i]] = n_neighs;
        free(x);
        free(y);
    }


    for (uint i = 1U; i < sys->n_sites; i++)
    {
        nw_lims[i] += nw_lims[i - 1];
    }

    for (uint i = 0U; i < sys->n_sites; i++)
    {
        uint n_behind = (i != 0U)? sys->limits[i - 1] : 0U;
        uint n_neighs = sys->limits[i] - n_behind;
        uint site = keys[i];
        uint start_point = (site != 0U)? nw_lims[site - 1] : 0U;
        memcpy((void *) nw_neig + start_point,
                (void *) sys->neighbors + n_behind,
                n_neighs * sizeof(uint));
    }

    for (uint i = 0U; i < sys->n_links; i++)
    {
        nw_neig[i] = keys[nw_neig[i]];
    }

    free(sys->limits);
    free(sys->neighbors);

    sys->limits = nw_lims;
    sys->neighbors = nw_neig;

    free(keys);

    for (uint i = 0U; i < sys->n_sites; i++)
    {
        uint n_behind = (i != 0U)? sys->limits[i - 1] : 0U;
        uint n_items = sys->limits[i] - n_behind;
        qsort(sys->neighbors + n_behind, n_items, sizeof(uint), compare_uint);
    }

    return sys;
}

void print_spins_csv(const system_t * const sys)
{
    for (uint i = 0U; i < sys->n_sites; i++)
    {
        printf("%2.8f %2.8f %2.8f\n",
                sys->spins[i].x,
                sys->spins[i].y,
                sys->spins[i].z);
    }
}

void dipolar_field(
        spin_t * const __restrict field,
        const system_t * const __restrict sys,
        uint site)
{
    field->x = 0.0f; field->y = 0.0f; field->z = 0.0f;
    uint behind = (site != 0U) ? sys->limits[site - 1U] : 0U;
    for (uint i = behind; i < sys->limits[site]; i++)
    {
        field->x += sys->spins[sys->neighbors[i]].x;
        field->y += sys->spins[sys->neighbors[i]].y;
        field->z += sys->spins[sys->neighbors[i]].z;
    }
}

void spin_delta(
        const spin_t * const __restrict a,
        const spin_t * const __restrict b,
        spin_t * const __restrict delta)
{
    delta->x = b->x - a->x;
    delta->y = b->y - a->y;
    delta->z = b->z - a->z;
}

float compute_magnetization (const system_t * const sys)
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

float compute_energy (const system_t * const sys)
{
    double ene_accum = 0.0;
    spin_t field;

    for (uint i = 0U; i < sys->n_sites; i++)
    {
        dipolar_field (&field, sys, i);
        ene_accum = exchange_interaction(sys->spins + 1, &field, 1.0f);
    }

    return (float) (0.5 * ene_accum);
}
