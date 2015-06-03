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
    dim3 blockSize(1, 1, 1);
    dim3 gridSize(1, 1, 1);
}
