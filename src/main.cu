#include "system_cu.h"

#include <curand_kernel.h>

typedef unsigned int uint;

struct _result {
    float3 spin;
    uint site;
};

typedef struct _result result_t;

struct _site {
    float3 spin;
    uint limit;
};

typedef struct _site site_t;

/**
  * Initializes random number generator states for each thread
  * assuming that there is allocated memory to do so.
  */
__global__ void intit_rng (curandState * const state, const uint seed)
{
    uint tid = blockIdx.x * blockDim.x + threadIdx.x;
    curand_init (seed, tid, 0, &(state[tid]));
}

/**
  * Create a lot of moves and store what would the tread
  * do if it successes.
  */
__global__ void compute_move (
        result_t * const moves, curandState * const state,
        site_t * const sites, uint * const neighbors,
        uint n_sites, const float temp)
{
    uint tid = blockIdx.x * blockDim.x + threadIdx.x;
    uint site_idx = curand(state + tid) % n_sites;
    site_t site = sites[site_idx];

    float3 nw_spin;
    nw_spin.x = curand_uniform(state + tid);
    nw_spin.y = curand_uniform(state + tid);
    nw_spin.z = curand_uniform(state + tid);

    float3 delta;
    delta.x = site.spin.x - nw_spin.x;
    delta.y = site.spin.y - nw_spin.y;
    delta.z = site.spin.z - nw_spin.z;

    float3 field;
    field.x = 0.0f;
    field.y = 0.0f;
    field.z = 0.0f;

    uint start = (site_idx > 0)? sites[site_idx - 1].limit : 0;
    for (uint i = start; i < site.limit; i++)
    {
        site_t neighbor = sites[i];
        field.x += neighbor.spin.x;
        field.y += neighbor.spin.y;
        field.z += neighbor.spin.z;
    }

    float delta_e = - (
            delta.x * field.x +
            delta.y + field.y +
            delta.z + field.z
            );

    if (delta_e < 0.0)
    {
        moves[tid].site = site_idx;
        moves[tid].spin = nw_spin;
    }
    else
    {
        float r = curand_uniform(state + tid);
        if (r < expf(- delta_e / temp))
        {
            moves[tid].site = site_idx;
            moves[tid].spin = nw_spin;
        }
        else
        {
            moves[tid].site = n_sites;
        }
    }
}

int main(int argc, char ** argv)
{
    uint blockSize = 256;
    uint gridSize = 100;
    uint n_threads = gridSize * blockSize;
    uint seed = 131;
    uint size = 1024U;
    uint steps = 100U;

    curandState * d_state;
    cudaMalloc(&d_state, n_threads * sizeof(curandState));
    intit_rng <<< gridSize, blockSize >>> (d_state, seed);

    system_t * sys = create_lattice(size, 1U);
    uint n_sites = sys->n_sites;
    random_spins(sys->spins, n_sites);

    /* Move sites to device */

    site_t * h_sites;
    site_t * d_sites;

    h_sites = (site_t *) malloc (n_sites * sizeof(site_t));

    for (unsigned int i = 0; i < n_sites; ++i)
    {
        h_sites[i].spin.x = sys->spins[i].x;
        h_sites[i].spin.y = sys->spins[i].y;
        h_sites[i].spin.z = sys->spins[i].z;
        h_sites[i].limit = sys->limits[i];
    }

    cudaMalloc(&d_sites, n_sites * sizeof(site_t));

    cudaMemcpy(d_sites, h_sites, n_sites * sizeof(site_t), cudaMemcpyHostToDevice);

    free(h_sites);

    /* Move neighbor references to device */

    uint * d_neighbors;
    cudaMalloc(&d_neighbors, sys->n_links * sizeof(uint));
    cudaMemcpy(d_neighbors, sys->neighbors, sys->n_links * sizeof(uint), cudaMemcpyHostToDevice);

    /* The system in the host is not needed anymore */

    destroy_system(sys);

    /* Allocate memory for the results */

    result_t * d_results;
    cudaMalloc(&d_results, n_threads * sizeof(result_t));

    for (unsigned int i = 0U; i < steps; ++i)
    {
        for (unsigned int j = 0U; j < n_sites; j += n_threads)
        {
            compute_move <<< gridSize, blockSize >>> (
                    d_results, d_state, d_sites, d_neighbors,
                    n_sites, 0.5f);
        }
    }

    /* Free memory on device */
    cudaFree(d_sites);
    cudaFree(d_state);
    cudaFree(d_results);
}
