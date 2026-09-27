#ifndef SYMNMF_H
#define SYMNMF_H

#include <stdio.h>
#include <stdlib.h>

void print_matrix(double **matrix, int rows, int cols);

void free_matrix(double **matrix, int n);

double **similarity_matrix(double **pListOfVec, int N, int d);

double **diagonal_matrix(double **pSimilarityMatrix, int N);

double **laplacian_matrix(double **pDiagonalMatrix, double **pSimilarityMatrix, int N);

double **optimize_H_matrix(double **pH, double **pLaplacianMatrix, int N, int K);

double **multiply_matrices(double **A, double **B, int rowsA, int colsA, int colsB);

double **transpose_matrix(double **A, int rows, int cols);

double calc_frobenius_norm(double **A, double **B, int rows, int cols);

double calc__square_euclidian_distance(double *vec1, double *vec2, int d);

double **read_vectors_from_file(FILE *file, int *N, int *d);

void *safe_realloc(void *pointer, size_t new_size);

#endif