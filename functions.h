#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "helpers.h"
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
#define uchar                   unsigned char

#define CLAMP(v, min, max) if(v < min) { v = min; } else if(v > max) { v = max; }

// The structure that is used to store the data the functions need
typedef struct {
    int thread_id;
    int size;
    ppm_image *image_org, *image_rscl, **contour_map;
    pthread_barrier_t *barrier;
    uchar **grid;
} ppm_thread;

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
void free_resources(ppm_image *image, ppm_image **contour_map, unsigned char **grid, int step_x) {
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
}

// The thread function that is used to parallelize
// This function incapsulates the original functions of this algorithm such as
// rescale_image, sample_grid, march and init_contour_map
void *thread_func(void *data)
{
    // Get the information from the thread structure
    ppm_thread *info = (ppm_thread *) data;
    ppm_image *image = info->image_org;
    ppm_image *image_rscl = info->image_rscl;
    int size = info->size;
    int thread_id = info->thread_id;
    pthread_barrier_t *barrier = info->barrier;
    uchar **grid = info->grid;
    ppm_image **contour_map = info->contour_map;

    // Verifies if the original image and the rescaled image are not the same
    if (image != image_rscl) {
        
        // Creates the start point and the end point for the rescale
        int start_rscl = thread_id * image_rscl->y / size;
        int end_rscl = fmin((thread_id + 1) * image_rscl->y / size, image_rscl->y); 
        uint8_t sample[3];
        
        // Uses bicubic interpolation for scaling
        for (int i = 0; i < image_rscl->x; i++) {
            for (int j = start_rscl; j < end_rscl; j++) {
                float u = (float)i / (float)(image_rscl->x - 1);
                float v = (float)j / (float)(image_rscl->y - 1);
                sample_bicubic(image, u, v, sample);

                image_rscl->data[i * image_rscl->y + j].red = sample[0];
                image_rscl->data[i * image_rscl->y + j].green = sample[1];
                image_rscl->data[i * image_rscl->y + j].blue = sample[2];
            }
        }

        // Calls the barrier
        pthread_barrier_wait(barrier);
    }

    int p = image_rscl->x / STEP;
    int q = image_rscl->y / STEP;

    // Creates the start points and the end points for sample
    int start_grid1 = thread_id * q / size;
    int end_grid1 = fmin((thread_id + 1) * q / size, q);

    int start_grid2 = thread_id * p / size;
    int end_grid2 = fmin((thread_id + 1) * p / size, p);

    // Corresponds to step 1 of the marching squares algorithm, which focuses on sampling the image.
    // Builds a p x q grid of points with values which can be either 0 or 1, depending on how the
    // pixel values compare to the `sigma` reference value. The points are taken at equal distances
    // in the original image, based on the `step_x` and `step_y` arguments
    for (int i = 0; i < p; i++) {
        for (int j = start_grid1; j < end_grid1; j++) {
            ppm_pixel curr_pixel = image_rscl->data[i * STEP * image_rscl->y + j * STEP];

            uchar curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

            if (curr_color > SIGMA) {
                grid[i][j] = 0;
            } else {
                grid[i][j] = 1;
            }
        }
    }
    
    grid[p][q] = 0;

    // Last sample points have no neighbors below / to the right, so we use pixels on the
    // last row / column of the input image for them
    for (int i = start_grid2; i < end_grid2; i++) {
        ppm_pixel curr_pixel = image_rscl->data[i * STEP * image_rscl->y + image_rscl->x - 1];

        uchar curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

        if (curr_color > SIGMA) {
            grid[i][q] = 0;
        } else {
            grid[i][q] = 1;
        }
    }
    for (int j = start_grid1; j < end_grid1; j++) {
        ppm_pixel curr_pixel = image_rscl->data[(image_rscl->x - 1) * image_rscl->y + j * STEP];

        uchar curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

        if (curr_color > SIGMA) {
            grid[p][j] = 0;
        } else {
            grid[p][j] = 1;
        }
    }

    // Calls the barrier
    pthread_barrier_wait(barrier);

    // Creates the start point and point for map_init
    int start_map = thread_id * CONTOUR_CONFIG_COUNT / size;
    int end_map = fmin((thread_id + 1) * CONTOUR_CONFIG_COUNT / size, CONTOUR_CONFIG_COUNT);

    // Creates a map between the binary configuration (e.g. 0110_2) and the corresponding pixels
    // that need to be set on the output image. An array is used for this map since the keys are
    // binary numbers in 0-15. Contour images are located in the './contours' directory
    for (int i = start_map; i < end_map; i++) {
        char filename[FILENAME_MAX_SIZE];
        sprintf(filename, "./contours/%d.ppm", i);
        contour_map[i] = read_ppm(filename);
    }

    // Calls the barrier
    pthread_barrier_wait(barrier);

    // Corresponds to step 2 of the marching squares algorithm, which focuses on identifying the
    // type of contour which corresponds to each subgrid. It determines the binary value of each
    // sample fragment of the original image and replaces the pixels in the original image with
    // the pixels of the corresponding contour image accordingly
    for (int i = 0; i < p; i++) {
        for (int j = start_grid1; j < end_grid1; j++) {
            unsigned char k = 8 * grid[i][j] + 4 * grid[i][j + 1] + 2 * grid[i + 1][j + 1] + 1 * grid[i + 1][j];
            update_image(image_rscl, contour_map[k], i * STEP, j * STEP);
        }
    }

    pthread_exit(NULL);
}

#endif