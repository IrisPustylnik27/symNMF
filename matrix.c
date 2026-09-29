#include <stdlib.h>
#include "matrix.h"

double **allocate_matrix(int rows, int cols, int zero_initialize)
{
    double **matrix;
    int i;

    if (rows <= 0 || cols <= 0)
        return NULL;
    if (!(matrix = malloc(rows * sizeof(double *))))
        return NULL;

    for (i = 0; i < rows; i++)
    {
        if (zero_initialize)
        {
            matrix[i] = calloc(cols, sizeof(double));
        }
        else
        {
            matrix[i] = malloc(cols * sizeof(double));
        }

        if (matrix[i] == NULL)
        {
            free_matrix(matrix, i);
            return NULL;
        }
    }
    return matrix;
}

void free_matrix(double **matrix, int n)
{
    int i;

    if (matrix == NULL)
    {
        return;
    }
    for (i = 0; i < n; i++)
    {
        free(matrix[i]);
    }
    free(matrix);
}

double **transpose_matrix(double **A, int rows, int cols)
{
    int i, j;
    double **B;

    if (!(B = allocate_matrix(cols, rows, 0)))
        return NULL;

    for (i = 0; i < cols; i++)
    {
        for (j = 0; j < rows; j++)
        {
            B[i][j] = A[j][i];
        }
    }
    return B;
}

double **multiply_matrices(double **A, double **B, int rowsA, int colsA, int colsB)
{
    int i, j, k;
    double **C;

    if (!(C = allocate_matrix(rowsA, colsB, 1)))
        return NULL;

    for (i = 0; i < rowsA; i++)
    {
        for (j = 0; j < colsB; j++)
        {
            for (k = 0; k < colsA; k++)
            {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return C;
}

void copy_matrix(double **old_matrix, double **new_matrix, int n, int k)
{
    int i, j;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < k; j++)
        {
            new_matrix[i][j] = old_matrix[i][j];
        }
    }
}

double squared_frobenius_norm(double **A, double **B, int rows, int cols)
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
