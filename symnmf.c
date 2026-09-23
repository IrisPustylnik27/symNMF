#define _GNU_SOURCE
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
    int N = 0, d = 0;
    const int ITER = 300;
    const double EPSILON = 1e-4;
    const double BETA = 0.5;

    if (!parse_arguments(argc, argv, &goal, &file_name))
    {
        return 1;
    }

    FILE *file = fopen(file_name, "r");
    if (file == NULL)
    {
        fprintf(stderr, "An Error Has Occured\n");
        return 1;
    }       

    pListOfVec = read_vectors_from_file(file, &N, &d);

    if(goal == NULL || strcmp(goal, "spk") != 0)
    {
        fprintf(stderr, "An Error Has Occured\n");
        return 1;
    }

    double **pSimilarityMatrix = similarity_matrix(pListOfVec, N, d);
    double **pDiagonalMatrix = diagonal_matrix(pSimilarityMatrix, N);
    double **pLaplacianMatrix = laplacian_matrix(pDiagonalMatrix, pSimilarityMatrix, N);

    if(goal=="sym")
    {
        print_matrix(pSimilarityMatrix, N, N);
    }
    if(goal=="ddg")

    {
        print_matrix(pDiagonalMatrix, N, N);
    }
    if(goal=="norm")
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
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
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
                double distance = calc_euclidian_distance(pListOfVec[i], pListOfVec[j], d);
                pSimilarityMatrix[i][j] = exp(-distance/2.0);
            }
        }
    }
    return pSimilarityMatrix;
}

double **diagonal_matrix(double **pSimilarityMatrix, int N)
{
    int i,k,j;
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

double **laplacian_matrix(double **pDiagonalMatrix, double **pSimilarityMatrix, int N)
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
    int i, j, k;
    double **pH_old;
    double **pNumerator;
    double **pDenominator;
    double **pH_new = malloc(N * sizeof(double *));
    for (i = 0; i < N; i++)
    {
        pH_new[i] = malloc(K * sizeof(double));
        for (j = 0; j < K; j++)
        {
            pH_new[i][j] = pH[i][j];
        }
    }

    for (int t = 0; t < ITER; t++)
    {
        pH_old = pH_new;
        pNumerator = multiply_matrices(pLaplacianMatrix, pH_old, N, N, K);
        pDenominator = multiply_matrices(pH_old, transpose_matrix(pH_old, N, K), N, K, N);
        pDenominator = multiply_matrices(pDenominator, pH_old, N, N, K);

        for (i = 0; i < N; i++)
        {
            for (j = 0; j < K; j++)
            {
                if (pDenominator[i][j] != 0)
                {
                    pH_new[i][j] *= pH_old[i][j] * (1 - BETA + BETA*(pNumerator[i][j] / pDenominator[i][j]));
                }
            }
        }

        double max_diff = 0.0;
        for (i = 0; i < N; i++)
        {
            for (j = 0; j < K; j++)
            {
                double diff = calc_frobenius_norm(&pH_new[i][j], &pH_old[i][j], 1, 1);
                if (diff > max_diff)
                {
                    max_diff = diff;
                }
            }
        }

        if (max_diff < EPSILON || t == iter - 1)
        {
            break;
        }
    }

    return pH_new;
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

double calc_frobenius_norm(double **A, double **B, int rows, int cols)
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
    return sqrt(sum);
}

double calc_euclidian_distance(double *vec1, double *vec2, int d)
{
    int i;
    double sum = 0.0;
    for (i = 0; i < d; i++)
    {
        sum += (vec1[i] - vec2[i]) * (vec1[i] - vec2[i]);
    }
    return sqrt(sum);
}

double **read_vectors_from_file(FILE *file, int *N, int *d)
{
    size_t len = 0;
    char *pLine = NULL;
    int capacity = 1, capacityVecList = 1;
    double **pListOfVec = malloc(capacityVecList * sizeof(double *));

    while ((getline(&pLine, &len, file)) != -1)
    {
        int current_d = 0;
        double *pVector = malloc(capacity * sizeof(double));
        char *token = strtok(pLine, ",");
        double **temp;
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
    return pListOfVec;
}

void cleanup(double **pListOfVec, double **pCentroids, double *pCentroid, int *assignmentsToClusters, int N, int K, char *pLine, double *pVector)
{
    int i;
    if (pLine != NULL){
        free(pLine);
    }
    if (pVector != NULL)
    {
        free(pVector);
    }
    if (assignmentsToClusters != NULL)
    {
        free(assignmentsToClusters);
    }
    if (pListOfVec != NULL)
    {
        for (i = 0; i < N; i++)
        {
            free(pListOfVec[i]);
        }
        free(pListOfVec);
    }
    if (pCentroids != NULL)
    {
        for (i = 0; i < K; i++)
        {
            free(pCentroids[i]);
        }
        free(pCentroids);
    }
    if (pCentroid != NULL)
    {
        free(pCentroid);
    }
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