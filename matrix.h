#ifndef MATRIX_H
#define MATRIX_H

double **multiply_matrices(double **A, double **B, int rowsA, int colsA, int colsB);

double **transpose_matrix(double **A, int rows, int cols);

double squared_frobenius_norm(double **A, double **B, int rows, int cols);

void copy_matrix(double **old_matrix, double **new_matrix, int n, int k);

void free_matrix(double **matrix, int n);

double **allocate_matrix(int rows, int cols, int zero_initialize);

#endif
