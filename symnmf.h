#ifndef SYMNMF_H
#define SYMNMF_H

void free_matrix(double **matrix, int n);

double **similarity_matrix(double **vectors, int n, int d);

double **diagonal_matrix(double **A, int n);

double **normalized_similarity_matrix(double **D, double **A, int n);

double **optimize_H_matrix(double **H, double **W, int n, int k);

double **allocate_matrix(int rows, int cols, int zero_initialize);

#endif
