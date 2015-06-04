#ifndef SYSTEM_CU_H_
#define SYSTEM_CU_H_

struct _spin {
    float x;
    float y;
    float z;
};

typedef struct _spin spin_t;

struct _system {
    spin_t * spins;
    unsigned int n_sites;
    unsigned int n_links;
    unsigned int * limits;
    unsigned int * neighbors;
};

typedef struct _system system_t;

extern "C"
void random_spins(spin_t * const, unsigned int);

extern "C"
system_t * create_lattice(unsigned int, unsigned short);

extern "C"
void destroy_system(system_t *);

#endif
