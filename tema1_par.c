// Author: APD team, except where source was noted

#include "helpers.h"
#include "functions.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <pthread.h>
#include <math.h>

// The main function
int main(int argc, char *argv[]) {
    if (argc < 4) {
        fprintf(stderr, "Usage: ./tema1 <in_file> <out_file> <P>\n");
        return 1;
    }

    ppm_image *image = read_ppm(argv[1]);
    ppm_image *img_aux;
    int P = atoi(argv[3]);

    pthread_t threads[P];
    ppm_thread thread[P];
    pthread_barrier_t barrier;

    if (image->x <= RESCALE_X && image->y <= RESCALE_Y) {
        img_aux = image;
    }
    else {
        // Allocates memory for image
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
        img_aux = new_image;
    }

    int p = img_aux->x / STEP;
    int q = img_aux->y / STEP;

    // Allocates memory for grid
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

    // Allocates memory for map
    ppm_image **map = (ppm_image **)malloc(CONTOUR_CONFIG_COUNT * sizeof(ppm_image *));
    if (!map) {
        fprintf(stderr, "Unable to allocate memory\n");
        exit(1);
    }

    // Initializes the barrier
    pthread_barrier_init(&barrier, NULL, P);

    // Creates the threads
    for (int i = 0; i < P; i++)
    {
        thread[i].thread_id = i;
        thread[i].size = P;
        thread[i].image_org = image;
        thread[i].image_rscl = img_aux; 
        thread[i].barrier = &barrier;
        thread[i].grid = grid;
        thread[i].contour_map = map;
        pthread_create(&threads[i], NULL, thread_func, &thread[i]);
    }

    // Waits for the threads
    for (int i = 0; i < P; i++)
    {
        pthread_join(threads[i], NULL);
    }
    
    // Destroys the barrier
    pthread_barrier_destroy(&barrier);

    write_ppm(img_aux, argv[2]);
    
    // Verifies if the original image and the rescaled one are not same and if it is so,
    // we also free the memory that was allocated for the original one
    if (img_aux != image)
    {
        free(image->data);
        free(image);
    }
    
    free_resources(img_aux, map, grid, STEP);

    return 0;
}
