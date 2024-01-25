# Parallelization of Marching Squares algorithm

------------------------------------------------------------------------

## Table of Contents:
- Thread Functions and Structure
- Sample Grid Parallelization
- March Parallelization
- Rescale Image Parallelization
------------------------------------------------------------------------

## Thread Functions and Structure


### Create_thread_mem

Thread function argument is structured using the `thread_mem` data structure. This structure is used to pass essential data to the thread function. Its key components are:

- `step_x`, `step_y`:     Horizontal and vertical step sizes used for grid sampling.
- `P`:     The total number of threads or thread count.
- `Id`:     The unique identifier for each thread.
- `image`, `unscaled_image`:     Represent scaled and unscaled images . The choice between these two images depends on whether scaling is required.
- `grid`:     Two-dimensional array used for grid sampling.
- `barrier`:     Pointer to a pthread barrier, which is used to synchronize thread execution.
- `contour_map`:     Array of contour images.

### Init_grid

Creates and initializes the grid used for image sampling.

- `scaled_image`:     Image for which the grid is being created.
- `step_x` and `step_y`:     Step sizes used for grid sampling.

The function calculates the dimensions of the grid (`p` and `q`) based on the scaled image's dimensions. It then allocates memory for the grid and returns a pointer to the initialized grid.

------------------------------------------------------------------------


## Rescale Image Parallelization

The `rescale_image` function is responsible for rescaling the image using bicubic interpolation. This function may be computationally expensive, especially for large images. Each thread independently performs bicubic interpolation on its segment, resulting in faster rescaling. The parallelizationn is performed on `x` dimension of the scaled matrix.

`int start_x = Id * (double)scaled_image->x / P
`int end_x = fmin(scaled_image->x, (Id + 1) * (double)scaled_image->x / P)



------------------------------------------------------------------------

## Sample Grid Parallelization


The `sample_grid` function is responsible for sampling the image and creating a grid of binary values based on a threshold (`SIGMA`). To speed up this process, we've parallelized the `sample_grid` function by dividing the image into smaller segments:

`int start_p = data.Id * (double)p / data.P`
`int end_p = fmin(p, (data.Id + 1) * (double)p / data.P)

`int start_q = data.Id * (double)q / data.P;
`int end_q = fmin(q, (data.Id + 1) * (double)q / data.P);

Each segment to a separate thread.

------------------------------------------------------------------------


## March Parallelization

The `march` function is used to identify the type of contour that corresponds to each subgrid.

Similar to the `sample_grid` function, we've parallelized the `march` function by dividing the grid into smaller segments based on `q` dimension of the matrix. 
