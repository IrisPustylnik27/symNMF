#define _GNU_SOURCE
#define ITER 300
#define EPSILON 1e-4
#define BETA 0.5
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "symnmf.h"
#include "matrix.h"

static void print_matrix(double **matrix, int rows, int cols);

static double **calculate_new_H_matrix(double **pH_old, double **norm_matrix, int n, int K);

static int is_valid_goal(const char *goal);

static int run_goal(const char *goal, double **vectors, int n, int d);

static double squared_euclidean_distance(double *vec1, double *vec2, int d);

static double **read_vectors_from_file(FILE *file, int *n, int *d);

static int resize_vector(double **vector, size_t new_size);

static int resize_matrix(double ***matrix, size_t new_rows);

static double *parse_vector(const char *line, int *d);

static double **load_vectors(char *file_name, int *n, int *d);

static void fill_new_H(double **H_new, double **H_old, double **numerator, double **denominator, int n, int k);

static int is_blank_line(const char *line);

static int resize_if_needed(double ***vectors, int *capacity_vec_list, int n);

int main(int argc, char *argv[])
{
    char *goal = NULL, *file_name = NULL;
    double **vectors = NULL;
    int n = 0, d = 0;

    if (argc != 3)
    {
        fprintf(stderr, "An Error Has Occurred\n");
        return 1;
    }
    goal = argv[1];
    file_name = argv[2];

    if (!is_valid_goal(goal))
    {
        fprintf(stderr, "An Error Has Occurred\n");
        return 1;
    }

    if (!(vectors = load_vectors(file_name, &n, &d)))
    {
        return 1;
    }

    if (!run_goal(goal, vectors, n, d))
    {
        free_matrix(vectors, n);
        fprintf(stderr, "An Error Has Occurred\n");
        return 1;
    }

    free_matrix(vectors, n);
    return 0;
}

static double **load_vectors(char *file_name, int *n, int *d)
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

static int is_valid_goal(const char *goal)
{
    return strcmp(goal, "sym") == 0 || strcmp(goal, "ddg") == 0 || strcmp(goal, "norm") == 0;
}

static int run_goal(const char *goal, double **vectors, int n, int d)
{
    double **A = NULL, **D = NULL, **W = NULL;
    int success = 1;

    if (!(A = similarity_matrix(vectors, n, d)))
        goto error;
    if (strcmp(goal, "sym") == 0)
    {
        print_matrix(A, n, n);
        goto cleanup;
    }
    if (!(D = diagonal_matrix(A, n)))
        goto error;
    if (strcmp(goal, "ddg") == 0)
    {
        print_matrix(D, n, n);
        goto cleanup;
    }
    if (!(W = normalized_similarity_matrix(D, A, n)))
        goto error;
    print_matrix(W, n, n);
    goto cleanup;

error:
    success = 0;

cleanup:
    free_matrix(A, n);
    free_matrix(D, n);
    free_matrix(W, n);

    return success;
}

static void print_matrix(double **matrix, int rows, int cols)
{
    int i, j;
    for (i = 0; i < rows; i++)
    {
        for (j = 0; j < cols; j++)
        {
            printf("%.4f", matrix[i][j]);
            if (j < cols - 1)
            {
                printf(",");
            }
        }
        printf("\n");
    }
}

double **similarity_matrix(double **vectors, int n, int d)
{
    int i, j;
    double distance;
    double **A;

    if (!(A = allocate_matrix(n, n, 0)))
        return NULL;

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            if (i == j)
            {
                A[i][j] = 0.0;
            }
            else
            {
                distance = squared_euclidean_distance(vectors[i], vectors[j], d);
                A[i][j] = exp(-distance / 2.0);
            }
        }
    }
    return A;
}

double **diagonal_matrix(double **A, int n)
{
    int i, j;
    double **D;

    if (!(D = allocate_matrix(n, n, 1)))
        return NULL;

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            D[i][i] += A[i][j];
        }
    }

    return D;
}

double **normalized_similarity_matrix(double **D, double **A, int n)
{
    int i, j;
    double **W;
    double denominator;

    if (!(W = allocate_matrix(n, n, 1)))
        return NULL;

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            denominator = sqrt(D[i][i] * D[j][j]);
            if (denominator == 0.0)
            {
                free_matrix(W, n);
                return NULL;
            }
            W[i][j] = A[i][j] / denominator;
        }
    }
    return W;
}

double **optimize_H_matrix(double **H, double **W, int n, int k)
{
    int t;
    double diff = 0.0;
    double **H_new, **H_old;

    if (!(H_old = allocate_matrix(n, k, 1)))
        return NULL;

    copy_matrix(H, H_old, n, k);

    for (t = 0; t < ITER; t++)
    {
        H_new = calculate_new_H_matrix(H_old, W, n, k);
        if (H_new == NULL)
        {
            free_matrix(H_old, n);
            return NULL;
        }

        diff = squared_frobenius_norm(H_new, H_old, n, k);
        free_matrix(H_old, n);
        H_old = H_new;

        if (diff < EPSILON)
            break;
    }

    return H_old;
}



static double **calculate_new_H_matrix(double **H_old, double **W, int n, int k)
{
    double **H_new = NULL, **H_old_transpose = NULL, **numerator = NULL, **denominator = NULL, **temp = NULL;

    if (!(H_new = allocate_matrix(n, k, 1)))
        return NULL;

    if (!(H_old_transpose = transpose_matrix(H_old, n, k)))
        goto error;

    if (!(numerator = multiply_matrices(W, H_old, n, n, k)))
        goto error;

    if (!(temp = multiply_matrices(H_old, H_old_transpose, n, k, n)))
        goto error;

    if (!(denominator = multiply_matrices(temp, H_old, n, n, k)))
        goto error;

    fill_new_H(H_new, H_old, numerator, denominator, n, k);

    goto cleanup;

error:
    free_matrix(H_new, n);
    H_new = NULL;

cleanup:
    free_matrix(numerator, n);
    free_matrix(denominator, n);
    free_matrix(temp, n);
    free_matrix(H_old_transpose, k);

    return H_new;
}

static void fill_new_H(double **H_new, double **H_old, double **numerator, double **denominator, int n, int k)
{
    int i, j;
    double ratio;

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < k; j++)
        {
            ratio = 0.0;
            if (denominator[i][j] != 0.0)
            {
                ratio = numerator[i][j] / denominator[i][j];
                H_new[i][j] = H_old[i][j] * (1 - BETA + BETA * ratio);
            }
        }
    }
}

static double squared_euclidean_distance(double *vec1, double *vec2, int d)
{
    int i;
    double sum = 0.0;

    for (i = 0; i < d; i++)
    {
        sum += (vec1[i] - vec2[i]) * (vec1[i] - vec2[i]);
    }
    return sum;
}

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

static int resize_matrix(double ***matrix, size_t new_rows)
{
    double **new_matrix = realloc(*matrix, new_rows * sizeof(double *));
    if (new_matrix == NULL)
        return 0;
    *matrix = new_matrix;
    return 1;
}

static int resize_vector(double **vector, size_t new_size)
{
    double *new_vector = realloc(*vector, new_size * sizeof(double));
    if (new_vector == NULL)
        return 0;
    *vector = new_vector;
    return 1;
}
