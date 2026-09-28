#ifndef SYMNMF_H
#define SYMNMF_H

#include <stdio.h>
#include <stdlib.h>

void free_matrix(double **matrix, int n);

double **similarity_matrix(double **pListOfVec, int n, int d);

double **diagonal_matrix(double **pSimilarityMatrix, int n);

double **normalized_similarity_matrix(double **pDiagonalMatrix, double **pSimilarityMatrix, int N);

double **optimize_H_matrix(double **pH, double **pLaplacianMatrix, int n, int k);

double **allocate_matrix(int rows, int cols, int zero_initialize);

#endif
