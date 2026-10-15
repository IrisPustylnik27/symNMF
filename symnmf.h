#ifndef SYMNMF_H
#define SYMNMF_H

/* Builds the n-by-n similarity matrix for the given d-dimensional vectors. */
double **similarity_matrix(double **vectors, int n, int d);

/* Builds the n-by-n diagonal degree matrix of the similarity matrix A. */
double **diagonal_matrix(double **A, int n);

/* Builds the n-by-n normalized similarity matrix from D and A. */
double **normalized_similarity_matrix(double **D, double **A, int n);

/* Optimizes the n-by-k matrix H against W and returns a newly allocated result. */
double **optimize_H_matrix(double **H, double **W, int n, int k);

#endif
