#include "system.h"

float exchange_interaction(const spin_t * a, const spin_t * b, float exchange)
{
    float energy = sqrt(a->x * b->x + a->y * b->y + a->z * b->z);
    return - exchange * energy;
}

float spin_norm(const spin_t * spin)
{
    return sqrt(spin->x * spin->x + spin->y * spin->y + spin->z * spin->z);
}

void random_spins(spin_t * spins, uint count)
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

void random_spin_gsl (spin_t * spins, const gsl_rng * rng, uint count)
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

    /* Now order the sites according to keys */

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


    /* Integrate the new limits */
    for (uint i = 0U; i < sys->n_sites; i++)
    {
        if (i != 0U)
        {
            nw_lims[i] += nw_lims[i - 1];
        }
    }

    /* Copy the new neighbors */
    for (uint i = 0U; i < sys->n_sites; i++)
    {
        uint n_behind = (i != 0U)? sys->limits[i - 1] : 0U;
        uint n_neighs = sys->limits[i] - n_behind;
        uint site = keys[i];
        uint start_point = (site != 0U)? nw_lims[site - 1] : 0U;
        memcpy((void *) nw_neig+ start_point,
                (void *) sys->neighbors + n_behind,
                n_neighs * sizeof(uint));
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

void print_spins_csv(const system_t * sys)
{
    for (uint i = 0U; i < sys->n_sites; i++)
    {
        printf("%2.8f %2.8f %2.8f\n",
                sys->spins[i].x,
                sys->spins[i].y,
                sys->spins[i].z);
    }
}

