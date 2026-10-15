#ifndef MATRIX_H
#define MATRIX_H

/* Multiplies a rowsA-by-colsA matrix A by a colsA-by-colsB matrix B. */
double **multiply_matrices(double **A, double **B, int rowsA, int colsA, int colsB);

/* Returns a newly allocated transpose of the given rows-by-cols matrix. */
double **transpose_matrix(double **A, int rows, int cols);

/* Returns the squared Frobenius norm of the difference between A and B. */
double squared_frobenius_norm(double **A, double **B, int rows, int cols);

/* Copies all entries of an n-by-k matrix into an allocated destination matrix. */
void copy_matrix(double **old_matrix, double **new_matrix, int n, int k);

/* Frees a matrix containing n allocated rows; accepts NULL. */
void free_matrix(double **matrix, int n);

/* Allocates a rows-by-cols matrix, optionally initializing all entries to zero. */
double **allocate_matrix(int rows, int cols, int zero_initialize);

#endif
