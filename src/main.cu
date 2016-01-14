#include "../include/system_cu.h"

#include <curand_kernel.h>
#include <thrust/device_ptr.h>
#include <thrust/sort.h>

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

__device__ float3 rand_spin (curandState * const state)
{
    float3 nw_spin;
    nw_spin.x = curand_normal(state);
    nw_spin.y = curand_normal(state);
    nw_spin.z = curand_normal(state);
    float norm = sqrtf(
            nw_spin.x * nw_spin.x +
            nw_spin.y * nw_spin.y +
            nw_spin.z * nw_spin.z);
    nw_spin.x /= norm;
    nw_spin.y /= norm;
    nw_spin.z /= norm;
    return nw_spin;
}

/**
  * Create a lot of moves and store what would the tread
  * do if it successes.
  */
__global__ void compute_move (
        result_t * const moves, uint * const active,
        curandState * const state,
        site_t * const sites, uint * const neighbors,
        uint n_sites, const float temp)
{
    uint tid = blockIdx.x * blockDim.x + threadIdx.x;
    uint site_idx = curand(state + tid) % n_sites;
    site_t site = sites[site_idx];

    float3 nw_spin = rand_spin (state + tid);

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
        site_t neighbor = sites[neighbors[i]];
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
        active[tid] = site_idx;
    }
    else
    {
        float r = curand_uniform(state + tid);
        if (r < expf(- delta_e / temp))
        {
            moves[tid].site = site_idx;
            moves[tid].spin = nw_spin;
            active[tid] = site_idx;
        }
        else
        {
            moves[tid].site = n_sites;
            active[tid] = n_sites;
        }
    }
}

__global__ void perform_moves (
        result_t * const moves, uint * const active,
        site_t * const sites, uint * const neighbors,
        uint n_sites, uint n_threads)
{
    uint tid = blockIdx.x * blockDim.x + threadIdx.x;
    result_t my_move = moves[tid];
    uint site_idx = my_move.site;
    if (site_idx < n_sites)
    {
        /* detect conflicts */
        bool conflicts_found = false;
        uint start = (site_idx > 0)? sites[site_idx - 1].limit : 0;

        for (uint i = start; i < sites[site_idx].limit; i++)
        {
            /*
               uint imin = 0U;
               uint imax = n_threads;
               uint imid;
               uint key = neighbors[i];
            */

            bool found = false;
            for (uint j = 1; j < n_threads; j *= 2)
                found = (j == tid);

            /*
                while (imax >= imin)
                {
                    imid = (imax + imin) / 2U;
                    if (active[imid] == key)
                    {
                        found = true;
                        break;
                    }
                    else if (active[imid] < key) imin = imid + 1;
                    else imax = imid - 1;
                }
            */

            if (found) {
                conflicts_found = true;
                break;
            }
        }

        if (!conflicts_found)
        {
            sites[site_idx].spin = my_move.spin;
        }
    }
}

inline void fail_on_cuda_error (int line)
{
#ifdef CUDA_DEBUG
    cudaError_t error = cudaGetLastError();
    if (error != cudaSuccess)
    {
        printf("%s on line: %i\n", cudaGetErrorString(error), line);
        exit(-1);
    }
#endif
}

inline void fail_on_bad_device_alloc(int line)
{
#ifdef CUDA_DEBUG
    printf("bad device alloc! in line: %i", line);
    exit(-1);
#endif
}

int main(int argc, char ** argv)
{
    uint blockSize = 256;
    uint gridSize = 64;
    uint n_threads = gridSize * blockSize;
    uint seed = 131;
    uint size = 1024U;
    uint steps = 100U;

    curandState * d_state;
    cudaMalloc((void **) &d_state, n_threads * sizeof(curandState));
    if (d_state == NULL)
    {
        fail_on_bad_device_alloc(__LINE__);
    }
    intit_rng <<< gridSize, blockSize >>> (d_state, seed);
    fail_on_cuda_error(__LINE__);

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

    cudaMalloc((void **) &d_sites, n_sites * sizeof(site_t));
    if (d_sites == NULL)
    {
        fail_on_bad_device_alloc(__LINE__);
    }

    cudaMemcpy(d_sites, h_sites, n_sites * sizeof(site_t), cudaMemcpyHostToDevice);

    free(h_sites);

    /* Move neighbor references to device */

    uint * d_neighbors;
    cudaMalloc((void **) &d_neighbors, sys->n_links * sizeof(uint));
    if (d_sites == NULL)
    {
        fail_on_bad_device_alloc(__LINE__);
    }
    cudaMemcpy(d_neighbors, sys->neighbors, sys->n_links * sizeof(uint), cudaMemcpyHostToDevice);

    /* The system in the host is not needed anymore */

    destroy_system(sys);

    /* Allocate memory for the results */

    result_t * d_results;
    uint * d_active;
    cudaMalloc((void **) &d_results, n_threads * sizeof(result_t));
    if (d_results == NULL)
    {
        fail_on_bad_device_alloc(__LINE__);
    }
    cudaMalloc((void **) &d_active, n_threads * sizeof(uint));
    if (d_results == NULL)
    {
        fail_on_bad_device_alloc(__LINE__);
    }

    /* Launch kernels */


    for (unsigned int i = 0U; i < steps; ++i)
    {
        for (unsigned int j = 0U; j < n_sites; j += n_threads)
        {
            compute_move <<< gridSize, blockSize >>> (
                    d_results, d_active, d_state, d_sites, d_neighbors,
                    n_sites, 0.5f);
            fail_on_cuda_error(__LINE__);

            thrust::device_ptr<uint> dt_active = thrust::device_pointer_cast(d_active);
            thrust::sort(dt_active, dt_active + n_threads);
            fail_on_cuda_error(__LINE__);

            perform_moves <<< gridSize, blockSize >>> (
                    d_results, d_active, d_sites, d_neighbors, n_sites, n_threads);
            fail_on_cuda_error(__LINE__);
        }
    }

    /* Free memory on device */
    cudaFree(d_sites);
    cudaFree(d_state);
    cudaFree(d_results);
    cudaFree(d_active);
}

