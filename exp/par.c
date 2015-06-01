#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

void generate_rand_omp(float * harvest, unsigned int size)
{
    unsigned int state = 13U;
#pragma omp parallel for schedule(static) firstprivate(state)
    for (unsigned int i = 0; i < size; i++)
    {
        harvest[i] = (float) rand_r(&state) / (float) RAND_MAX;
    }
}

int main(int argc, char *argv[])
{
    float * harvest;
    unsigned int size = 10;

    harvest = (float *) malloc (sizeof(float) * size);
    generate_rand_omp(harvest, size);

    for (unsigned int i = 0; i < size; i++)
    {
        printf("%f\n", harvest[i]);
    }

    unsigned int state = 10;
#pragma omp parallel firstprivate(state)
    printf ("%u out of %u: %u at %p\n", omp_get_thread_num(), omp_get_num_threads(),
            state, &state);

    free(harvest);
    return 0;
}
