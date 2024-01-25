#define _XOPEN_SOURCE 600

#include "helpers_thread.h"
#include <stdio.h>
#include <stdlib.h>

//Initialization of threaad function argument.
//mem - structure passed to pthread_create as argument
//new_image - scaled image
//image - input image
//P - thread number
//Id - thread id
void Create_thread_mem(thread_mem* mem,int step_x, int step_y, int P, int Id, ppm_image* new_image, 
ppm_image* image, unsigned char** grid, pthread_barrier_t* barrier, ppm_image** contour_map){
    mem->step_x = step_x;
    mem->step_y = step_y;
    mem->P = P;
    mem->Id = Id;
    //Checks if the image needs to scale. If scaling is necessary 
    //sets output image to the scaled image, otherwise output is 
    // set to initial image.
    if (image->x > RESCALE_X && image->y > RESCALE_Y) {
        mem->image = new_image;
        mem->unscaled_image = image;
    }else{
        mem->image = image;
        mem->unscaled_image = new_image;
    }
    
    mem->grid = grid;
    mem->barrier = barrier;
    mem->contour_map = contour_map;
}

//Initialization of grid
unsigned char ** Init_grid(ppm_image* scaled_image, int step_x, int step_y){
    int p = scaled_image->x / step_x;
    int q = scaled_image->y / step_y;


    unsigned char **grid = (unsigned char **)malloc((p + 1) * sizeof(unsigned char*));
    if (!grid) {
        fprintf(stderr, "Unable to allocate memory\n");
        exit(1);
    }

    for (int i = 0; i <= p; i++) {
        grid[i] = (unsigned char *)malloc((q + 1) * sizeof(unsigned char));
        if (!grid[i]) {
            fprintf(stderr, "Unable to allocate memory\n");
            exit(1);
        }
    }

    return grid;
}
