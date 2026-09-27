#define _GNU_SOURCE
#define ITER 300
#define EPSILON 1e-4
#define BETA 0.5
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>
#include "symnmf.h"


int main( int argc, char *argv[])
{
    char *goal = NULL;
    char *file_name = NULL;
    double **pListOfVec = NULL;
    double **pSimilarityMatrix, **pDiagonalMatrix, **pLaplacianMatrix;
    File *file = NULL;
    int N = 0, d = 0;

    if (!parse_arguments(argc, argv, &goal, &file_name))
    {
        return 1;
    }

    file = fopen(file_name, "r");
    if (file == NULL)
    {
        fprintf(stderr, "An Error Has Occured\n");
        return 1;
    }       

    pListOfVec = read_vectors_from_file(file, &N, &d);

    if(goal == NULL || 
        (strcmp(goal, "sym") != 0 &&
        strcmp(goal, "ddg") != 0 && 
        strcmp(goal, "norm") != 0))
    {
        fprintf(stderr, "An Error Has Occured\n");
        free_matrix(pListOfVec, N);
        return 1;
    }

    pSimilarityMatrix = similarity_matrix(pListOfVec, N, d);
    pDiagonalMatrix = diagonal_matrix(pSimilarityMatrix, N);
    pLaplacianMatrix = laplacian_matrix(pDiagonalMatrix, pSimilarityMatrix, N);

    if(strcmp(goal, "sym") == 0)
    {
        print_matrix(pSimilarityMatrix, N, N);
    }
    if(strcmp(goal, "ddg") == 0)

    {
        print_matrix(pDiagonalMatrix, N, N);
    }
    if(strcmp(goal, "norm") == 0)
    {
        print_matrix(pLaplacianMatrix, N, N);
    }

    free_matrix(pListOfVec, N);
    free_matrix(pSimilarityMatrix, N);
    free_matrix(pDiagonalMatrix, N);
    free_matrix(pLaplacianMatrix, N);
    return 0;
}

void print_matrix(double **matrix, int rows, int cols)
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

double **similarity_matrix(double **pListOfVec, int N, int d)
{
    int i, j;
    double **pSimilarityMatrix = malloc(N * sizeof(double *));
    for (i = 0; i < N; i++)
    {
        pSimilarityMatrix[i] = malloc(N * sizeof(double));
        for (j = 0; j < N; j++)
        {
            if (i == j)
            {
                pSimilarityMatrix[i][j] = 0.0;
            }
            else
            {
                double distance = calc_square_euclidian_distance(pListOfVec[i], pListOfVec[j], d);
                pSimilarityMatrix[i][j] = exp(-distance/2.0);
            }
        }
    }
    return pSimilarityMatrix;
}

double **diagonal_matrix(double **pSimilarityMatrix, int N)
{
    int i, k, j;
    double **pDiagonalMatrix = malloc(N * sizeof(double *));

    for (i = 0; i < N; i++)
    {
        pDiagonalMatrix[i] = malloc(N * sizeof(double));
        for (j = 0; j < N; j++)
        {
            if (i == j)
            {
                double sum = 0.0;
                for (k = 0; k < N; k++)
                {
                    sum += pSimilarityMatrix[i][k];
                }
                pDiagonalMatrix[i][j] = sum;
            }
            else
            {
                pDiagonalMatrix[i][j] = 0.0;
            }
        }
    }
    return pDiagonalMatrix;
}

double **norm_matrix(double **pDiagonalMatrix, double **pSimilarityMatrix, int N)
{
    int i, j;
    double **pLaplacianMatrix = malloc(N * sizeof(double *));

    for (i = 0; i < N; i++)
    {
        pLaplacianMatrix[i] = malloc(N * sizeof(double));
        for (j = 0; j < N; j++)
        {
            pLaplacianMatrix[i][j] = pSimilarityMatrix[i][j]/sqrt(pDiagonalMatrix[i][i] * pDiagonalMatrix[j][j]);
        }
    }
    return pLaplacianMatrix;
}

double **optimize_H_matrix(double **pH, double **pLaplacianMatrix, int N, int K)
{
    int i, j, t;
    double ratio;
    double diff = 0.0;
    double **pH_new, **pH_old_transpose, **pNumerator, **pDenominator, **temp;
    double **pH_old = malloc(N * sizeof(double *));

    for (i = 0; i < N; i++)
    {
        pH_old[i] = malloc(K * sizeof(double));
        for (j = 0; j < K; j++)
        {
            pH_old[i][j] = pH[i][j];
        }
    }

    for (t = 0; t < ITER; t++)
    {
        pH_new = malloc(N * sizeof(double *));
        for (i = 0; i < N; i++)
        {
            pH_new[i] = malloc(K * sizeof(double));
        }

        pH_old_transpose = transpose_matrix(pH_old, N, K);
        pNumerator = multiply_matrices(pLaplacianMatrix, pH_old, N, N, K);
        temp = multiply_matrices(pH_old, pH_old_transpose, N, K, N);
        pDenominator = multiply_matrices(pDenominator, pH_old, N, N, K);

        

        for (i = 0; i < N; i++)
        {
            for (j = 0; j < K; j++)
            {
                ratio = 0.0;
                if (pDenominator[i][j] != 0.0){
                    ratio = pNumerator[i][j] / pDenominator[i][j];
                    pH_new[i][j] = pH_old[i][j] * (1 - BETA + BETA * ratio);
                }
            }
        }

        diff = calc_square_frobenius_norm(pH_new, pH_old, N, K);

        free_matrix(pNumerator, N);
        free_matrix(pDenominator, N);
        free_matrix(temp, N);

        free_matrix(pH_old_transpose, K);
        free_matrix(pH_old, N);
        pH_old = pH_new;

        if (diff < EPSILON || t == ITER - 1)
        {
            break;
        }
    }

    return pH_old;
}

double **multiply_matrices(double **A, double **B, int rowsA, int colsA, int colsB)
{
    int i, j, k;
    double **C = malloc(rowsA * sizeof(double *));

    for (i = 0; i < rowsA; i++)
    {
        C[i] = malloc(colsB * sizeof(double));
        for (j = 0; j < colsB; j++)
        {
            C[i][j] = 0.0;
            for (k = 0; k < colsA; k++)
            {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return C;
}

double **transpose_matrix(double **A, int rows, int cols)
{
    int i, j;
    double **B = malloc(cols * sizeof(double *));

    for (i = 0; i < cols; i++)
    {
        B[i] = malloc(rows * sizeof(double));
        for (j = 0; j < rows; j++)
        {
            B[i][j] = A[j][i];
        }
    }
    return B;
}

double calc_square_frobenius_norm(double **A, double **B, int rows, int cols)
{
    int i, j;
    double sum = 0.0;

    for (i = 0; i < rows; i++)
    {
        for (j = 0; j < cols; j++)
        {
            sum += (A[i][j] - B[i][j]) * (A[i][j] - B[i][j]);
        }
    }
    return sum;
}

double calc__square_euclidian_distance(double *vec1, double *vec2, int d)
{
    int i;
    double sum = 0.0;

    for (i = 0; i < d; i++)
    {
        sum += (vec1[i] - vec2[i]) * (vec1[i] - vec2[i]);
    }
    return sum;
}

double **read_vectors_from_file(FILE *file, int *N, int *d)
{
    size_t len = 0;
    char *pLine = NULL, *token;
    int capacity = 1, capacityVecList = 1, current_d = 0;
    double **pVector = NULL, **temp;
    double **pListOfVec = malloc(capacityVecList * sizeof(double *));

    while ((getline(&pLine, &len, file)) != -1)
    {
        current_d = 0;
        pVector = malloc(capacity * sizeof(double));
        token = strtok(pLine, ",");
        temp = NULL;

        /* initialize the vector */
        while (token != NULL)
        {
            if (current_d == capacity)
            {
                capacity *= 2;
                pVector = safe_realloc(pVector, capacity * sizeof(double));
            }
            pVector[current_d] = atof(token);
            current_d++;
            token = strtok(NULL, ",");
        }

        if (*d == 0)
        {
            *d = current_d;
        }
        if (*N == capacityVecList)
        {
            capacityVecList *= 2;
            temp = safe_realloc(pListOfVec, capacityVecList * sizeof(double *));
            pListOfVec = temp;
        }
        pListOfVec[*N] = safe_realloc(pVector, *d * sizeof(double));
        (*N)++;
    }

    free(pLine);
    free(token);
    return pListOfVec;
}

void *safe_realloc(void *pointer, size_t new_size)
{
    void *new_pointer = realloc(pointer, new_size);

    if (new_pointer == NULL)
    {
        fprintf(stderr, "An Error Has Occured\n");
        free(pointer);
        exit(1);
        return NULL;
    }
    return new_pointer;
}

void free_matrix(double **matrix, int n){
    int i;

    if (matrix == NULL){
        return;
    }
    for (i = 0; i < n; i++)
        {
            free(matrix[i]);
        }
    free(matrix);
}