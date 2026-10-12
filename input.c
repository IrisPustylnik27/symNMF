#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "input.h"
#include "matrix.h"

static double **read_vectors_from_file(FILE *file, int *n, int *d);

static int resize_if_needed(double ***vectors, int *capacity_vec_list, int n);

static int is_blank_line(const char *line);

static double *parse_vector(const char *line, int *d);

static int resize_matrix(double ***matrix, size_t new_rows);

static int resize_vector(double **vector, size_t new_size);


/* Loads vectors from a file, sets their shape, and returns NULL on failure. */
double **load_vectors(char *file_name, int *n, int *d)
{
    FILE *file = NULL;
    double **vectors = NULL;
    if (!(file = fopen(file_name, "r")))
    {
        fprintf(stderr, "An Error Has Occurred\n");
        return NULL;
    }
    if (!(vectors = read_vectors_from_file(file, n, d)))
    {
        fprintf(stderr, "An Error Has Occurred\n");
        fclose(file);
        return NULL;
    }
    fclose(file);
    return vectors;
}

/* Reads all nonblank vector rows and returns their matrix, or NULL on failure. */
static double **read_vectors_from_file(FILE *file, int *n, int *d)
{
    size_t len = 0;
    char *line = NULL;
    int capacity_vec_list = 1;
    double *vector = NULL;
    double **vectors = NULL;

    if (!(vectors = malloc(capacity_vec_list * sizeof(double *))))
        return NULL;

    while ((getline(&line, &len, file)) != -1)
    {
        if (is_blank_line(line))
            continue;
        if (!(vector = parse_vector(line, d)))
            goto error;
        if (!resize_if_needed(&vectors, &capacity_vec_list, *n))
            goto error;
        vectors[*n] = vector;
        vector = NULL;
        (*n)++;
    }

    if (ferror(file) || *n == 0)
        goto error;

    free(line);
    return vectors;

error:
    free_matrix(vectors, *n);
    free(vector);
    free(line);
    return NULL;
}

/* Expands the vector-pointer array when full and reports whether it succeeded. */
static int resize_if_needed(double ***vectors, int *capacity_vec_list, int n)
{
    int new_capacity;

    if (n < *capacity_vec_list)
        return 1;

    new_capacity = 2 * *capacity_vec_list;
    if (!(resize_matrix(vectors, new_capacity)))
        return 0;
    *capacity_vec_list = new_capacity;

    return 1;
}

/* Returns whether a line contains only supported whitespace characters. */
static int is_blank_line(const char *line)
{
    int i = 0;

    if (line == NULL)
        return 1;

    while (line[i] != '\0')
    {
        if (line[i] != ' ' &&
            line[i] != '\t' &&
            line[i] != '\n' &&
            line[i] != '\r')
        {
            return 0;
        }
        i++;
    }
    return 1;
}

/* Parses one comma-separated vector and returns NULL on allocation failure. */
static double *parse_vector(const char *line, int *d)
{
    double *vector = NULL;
    char *token, *line_copy;
    int capacity = 1, current_d = 0;

    if (!(line_copy = strdup(line)))
        return NULL;

    if (!(vector = malloc(capacity * sizeof(double))))
        goto error;

    token = strtok(line_copy, ",");
    while (token != NULL)
    {
        if (current_d == capacity)
        {
            capacity *= 2;
            if (!(resize_vector(&vector, capacity)))
                goto error;
        }
        vector[current_d] = atof(token);
        current_d++;
        token = strtok(NULL, ",");
    }

    if (*d == 0)
        *d = current_d;
    goto no_error;

error:
    free(vector);
    vector = NULL;

no_error:
    free(line_copy);
    return vector;
}

/* Resizes a matrix pointer array and reports whether reallocation succeeded. */
static int resize_matrix(double ***matrix, size_t new_rows)
{
    double **new_matrix = realloc(*matrix, new_rows * sizeof(double *));
    if (new_matrix == NULL)
        return 0;
    *matrix = new_matrix;
    return 1;
}

/* Resizes a vector and reports whether reallocation succeeded. */
static int resize_vector(double **vector, size_t new_size)
{
    double *new_vector = realloc(*vector, new_size * sizeof(double));
    if (new_vector == NULL)
        return 0;
    *vector = new_vector;
    return 1;
}
