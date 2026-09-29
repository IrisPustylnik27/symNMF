#ifndef SYMNMF_H
#define SYMNMF_H

double **similarity_matrix(double **vectors, int n, int d);

double **diagonal_matrix(double **A, int n);

double **normalized_similarity_matrix(double **D, double **A, int n);

double **optimize_H_matrix(double **H, double **W, int n, int k);

#endif
