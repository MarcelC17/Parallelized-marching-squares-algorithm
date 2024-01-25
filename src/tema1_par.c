#define _XOPEN_SOURCE 600

#include "helpers.h"
#include "helpers_thread.h"

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <pthread.h>
#include <math.h>

#define CONTOUR_CONFIG_COUNT    16
#define FILENAME_MAX_SIZE       50
#define STEP                    8
#define SIGMA                   200
#define RESCALE_X               2048
#define RESCALE_Y               2048

#define CLAMP(v, min, max) if(v < min) { v = min; } else if(v > max) { v = max; }


// Creates a map between the binary configuration (e.g. 0110_2) and the corresponding pixels
// that need to be set on the output image. An array is used for this map since the keys are
// binary numbers in 0-15. Contour images are located in the './contours' directory.
ppm_image **init_contour_map() {
    ppm_image **map = (ppm_image **)malloc(CONTOUR_CONFIG_COUNT * sizeof(ppm_image *));
    if (!map) {
        fprintf(stderr, "Unable to allocate memory\n");
        exit(1);
    }

    for (int i = 0; i < CONTOUR_CONFIG_COUNT; i++) {
        char filename[FILENAME_MAX_SIZE];
        sprintf(filename, "./contours/%d.ppm", i);
        map[i] = read_ppm(filename);
    }

    return map;
}

// Updates a particular section of an image with the corresponding contour pixels.
// Used to create the complete contour image.
void update_image(ppm_image *image, ppm_image *contour, int x, int y) {
    for (int i = 0; i < contour->x; i++) {
        for (int j = 0; j < contour->y; j++) {
            int contour_pixel_index = contour->x * i + j;
            int image_pixel_index = (x + i) * image->y + y + j;

            image->data[image_pixel_index].red = contour->data[contour_pixel_index].red;
            image->data[image_pixel_index].green = contour->data[contour_pixel_index].green;
            image->data[image_pixel_index].blue = contour->data[contour_pixel_index].blue;
        }
    }
}

// Calls `free` method on the utilized resources.
void free_resources(ppm_image *image, ppm_image* new_image, ppm_image **contour_map,
unsigned char **grid, int step_x, pthread_barrier_t* barrier) {
    for (int i = 0; i < CONTOUR_CONFIG_COUNT; i++) {
        free(contour_map[i]->data);
        free(contour_map[i]);
    }
    free(contour_map);

    for (int i = 0; i <= image->x / step_x; i++) {
        free(grid[i]);
    }
    free(grid);

    free(image->data);
    free(image);

    free(new_image->data);
    free(new_image);

    pthread_barrier_destroy(barrier);
}

//Rescales image using bicubic interpolation
void rescale_image(ppm_image *image, ppm_image *scaled_image,
 int Id, int P) {
    uint8_t sample[3];
    //Parallel itteration boundaries x dimension
    int start_x = Id * (double)scaled_image->x / P;
 	int end_x = fmin(scaled_image->x, (Id + 1) * (double)scaled_image->x / P);    

    //Bicubic interpolation for scaling
    for (int i = start_x; i < end_x; i++) {
        for (int j = 0; j < scaled_image->y; j++) {
            float u = (float)i / (float)(scaled_image->x - 1);
            float v = (float)j / (float)(scaled_image->y - 1);
            sample_bicubic(image, u, v, sample);

            scaled_image->data[i * scaled_image->y + j].red = sample[0];
            scaled_image->data[i * scaled_image->y + j].green = sample[1];
            scaled_image->data[i * scaled_image->y + j].blue = sample[2];
        }
    }

}

// Corresponds to step 1 of the marching squares algorithm, which focuses on sampling the image.
// Builds a p x q grid of points with values which can be either 0 or 1, depending on how the
// pixel values compare to the `sigma` reference value. The points are taken at equal distances
// in the original image, based on the `step_x` and `step_y` arguments.
void sample_grid(int p, int q, int start_p, int end_p, int start_q, int end_q, ppm_image* image,
 int step_x, int step_y, unsigned char ** grid, int Id){
    
    //The iteration of matrix elements is executed in parallel, each thread itterates start_q - end_q 
    // elements of a pariticular row("the row is divided").
    for (int i = 0; i < p; i++) {
        for (int j = start_q; j < end_q; j++) {
            ppm_pixel curr_pixel = image->data[i * step_x * image->y + j * step_y];

            unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

            if (curr_color > SIGMA) {
                grid[i][j] = 0;
            } else {
                grid[i][j] = 1;
            }
        }
    }
    
    //Right-bottom corner set to 0. Conditioned to avoid race condition
    if (Id == 0){
        grid[p][q] = 0;
    }

    //The iteration through rows is executed in parallel.
    for (int i = start_p; i < end_p; i++) {

        ppm_pixel curr_pixel = image->data[i * step_x * image->y + image->x - 1];

        unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

        if (curr_color > SIGMA) {
            grid[i][q] = 0;
        } else {
            grid[i][q] = 1;
        }
        
    }

    //The iteration through columns is executed in parallel.
    for (int j = start_q; j < end_q; j++) {
        ppm_pixel curr_pixel = image->data[(image->x - 1) * image->y + j * step_y];

        unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

        if (curr_color > SIGMA) {
            grid[p][j] = 0;
        } else {
            grid[p][j] = 1;
        }
    }
    
}

// Corresponds to step 2 of the marching squares algorithm, which focuses on identifying the
// type of contour which corresponds to each subgrid. It determines the binary value of each
// sample fragment of the original image and replaces the pixels in the original image with
// the pixels of the corresponding contour image accordingly.
void march(int p, unsigned char ** grid, ppm_image* image, ppm_image ** contour_map, int step_x,
int step_y, int start_q, int end_q){
    
    for (int i = 0; i < p; i++) {
        for (int j = start_q; j < end_q; j++) {
            unsigned char k = 8 * grid[i][j] + 4 * grid[i][j + 1] + 2 * grid[i + 1][j + 1] + 1 * grid[i + 1][j];
            update_image(image, contour_map[k], i * step_x, j * step_y);
        }
    }
}

//Thread function
void *f(void *arg)
{   
    //Reads structure passed as function argument
    thread_mem data = *(thread_mem *)arg;

    //Checks if scaling is necessary. Scales the image.
    if (data.unscaled_image->x > RESCALE_X && data.unscaled_image->y > RESCALE_Y) {
        rescale_image(data.unscaled_image, data.image, data.Id, data.P);
    }

    //Grid dimensions
    int p = data.image->x / data.step_x;
    int q = data.image->y / data.step_y;

    //Itteration boundaries :
    //1. for p dimmension
    int start_p = data.Id * (double)p / data.P;
 	int end_p = fmin(p, (data.Id + 1) * (double)p / data.P);

    //2. for q dimmension
    int start_q = data.Id * (double)q / data.P;
 	int end_q = fmin(q, (data.Id + 1) * (double)q / data.P);
    
    //Waits for scaling step to finish
    pthread_barrier_wait(data.barrier);

    //1st step
    sample_grid(p, q, start_p, end_p, start_q, end_q, data.image, data.step_x, data.step_y, data.grid, data.Id);

    //Waits for sample_grid to finish
    pthread_barrier_wait(data.barrier);
    //2nd step
    march( p, data.grid, data.image, data.contour_map, data.step_x, data.step_y, start_q, end_q);

	pthread_exit(NULL);
}

//Alocation of rescaled image
ppm_image* alloc_rescale(){
    ppm_image *new_image = (ppm_image *)malloc(sizeof(ppm_image));
    if (!new_image) {
        fprintf(stderr, "Unable to allocate memory\n");
        exit(1);
    }
    new_image->x = RESCALE_X;
    new_image->y = RESCALE_Y;

    new_image->data = (ppm_pixel*)malloc(new_image->x * new_image->y * sizeof(ppm_pixel));
    if (!new_image) {
        fprintf(stderr, "Unable to allocate memory\n");
        exit(1);
    }

    return new_image;
}


int main(int argc, char *argv[]) {
    if (argc < 4) {
        fprintf(stderr, "Usage: ./tema1 <in_file> <out_file> <P>\n");
        return 1;
    }

    
    ppm_image *image = read_ppm(argv[1]);
    int step_x = STEP;
    int step_y = STEP;
    //Reads number of threads
    int THREAD_NUMBER = *(char*)argv[3] - '0';

    //Itteration and error
    int i,r;

	void *status;
	pthread_t threads[THREAD_NUMBER];

    //Data used as argument by thread function
    thread_mem data[THREAD_NUMBER];

    //Barrier declaration and initialization
    pthread_barrier_t barrier;
    pthread_barrier_init(&barrier, NULL, THREAD_NUMBER);


   
    ppm_image **contour_map = init_contour_map();

    //Rescaled image allocation
    ppm_image* new_image = alloc_rescale();

    //Initialize grid of points
    unsigned char **grid = Init_grid(new_image, step_x, step_y);

    


    //Start threads
    for (i = 0; i < THREAD_NUMBER; i++) {
        
        //Initialize of data 
        Create_thread_mem(&data[i], step_x, step_y, THREAD_NUMBER, i, new_image, image, grid, &barrier, contour_map);

		r = pthread_create(&threads[i], NULL, f, &data[i]);

		if (r) {
			printf("Eroare la crearea thread-ului %d\n", i);
			exit(-1);
		}
	}

    //Join threads
	for (i = 0; i < THREAD_NUMBER; i++) {
		r = pthread_join(threads[i], &status);

		if (r) {
			printf("Eroare la asteptarea thread-ului %d\n", i);
			exit(-1);
		}
	}

    //Write output
    write_ppm(data[0].image, argv[2]);


    free_resources(data[0].image, data[0].unscaled_image, contour_map, grid, step_x, &barrier);

    return 0;
}