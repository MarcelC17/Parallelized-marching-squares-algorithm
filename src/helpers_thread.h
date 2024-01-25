#include "helpers.h"
#include <stdlib.h>
#include <stdint.h>

typedef struct {
    int step_x, step_y, P, Id;
    ppm_image *image, *unscaled_image;
    unsigned char ** grid;
    ppm_image** contour_map;
    pthread_barrier_t* barrier;
} thread_mem;

void Create_thread_mem(thread_mem* mem,int step_x, int step_y, int P, int Id, ppm_image* new_image, ppm_image* image,
unsigned char** grid, pthread_barrier_t* barrier, ppm_image** contour_map);

unsigned char** Init_grid(ppm_image* scaled_image, int step_x, int step_y);
