# Symmetric Non-negative Matrix Factorization with Python and C

This project implements Symmetric Non-negative Matrix Factorization (SymNMF) for clustering numerical data.

The project combines Python and C: Python handles input processing, matrix initialization, and output formatting, while the main matrix computations and optimization algorithm are implemented in C and exposed through a custom Python extension.

## Overview

SymNMF is a clustering method that approximates a normalized similarity matrix `W` as the product of a non-negative matrix `H` and its transpose:

```text
W ≈ H · Hᵀ
```

Each row of `H` represents a data point's association with the clusters. A point is assigned to the cluster corresponding to the largest value in its row.

The project also includes an analysis program that compares SymNMF with a custom K-Means implementation using silhouette scores and the Adjusted Rand Index.

## Features

* SymNMF clustering with a C computational backend
* Similarity matrix calculation
* Diagonal degree matrix calculation
* Normalized similarity matrix calculation
* Custom Python C extension module
* Random initialization with a fixed seed for reproducible results
* Standalone C executable for matrix calculations
* Comparison with K-Means using silhouette scores
* Adjusted Rand Index for comparing the two cluster assignments
* Matrix output with 4 digits after the decimal point
* Dynamic memory allocation and cleanup in C

## Technologies

* Python
* C
* Python C API
* NumPy
* scikit-learn
* setuptools
* GCC
* Make

## Main Concepts Practiced

* Symmetric Non-negative Matrix Factorization
* Similarity-based clustering
* Matrix multiplication and normalization
* Iterative optimization and convergence checks
* K-Means clustering
* Clustering evaluation metrics
* Python and C integration
* Converting Python lists to C matrices and back
* Manual memory management
* Input validation and error handling

## Installation

Install the Python dependencies:

```bash
python3 -m pip install numpy scikit-learn setuptools
```

Build the C extension module:

```bash
python3 setup.py build_ext --inplace
```

Building the extension requires a C compiler and Python development headers.

## Usage

### Python Program

Run the program with:

```bash
python3 symnmf.py k goal file_name
```

Where:

* `k` is the number of clusters and must satisfy `1 < k < N`
* `N` is the number of input vectors
* `goal` selects the calculation to perform
* `file_name` is the path to the input file
* The `k` argument is required for every goal

Supported goals:

| Goal | Description |
| --- | --- |
| `sym` | Calculate the similarity matrix |
| `ddg` | Calculate the diagonal degree matrix |
| `norm` | Calculate the normalized similarity matrix |
| `symnmf` | Calculate the final association matrix `H` |

Examples:

```bash
python3 symnmf.py 2 sym input.txt
python3 symnmf.py 2 ddg input.txt
python3 symnmf.py 2 norm input.txt
python3 symnmf.py 2 symnmf input.txt
```

### Standalone C Program

Compile the program:

```bash
make
```

Run the program with:

```bash
./symnmf goal file_name
```

The standalone C program supports `sym`, `ddg`, and `norm`.

Example:

```bash
./symnmf norm input.txt
```

### Clustering Analysis

Compare SymNMF and K-Means:

```bash
python3 analysis.py k file_name
```

Example:

```bash
python3 analysis.py 2 input.txt
```

The program prints:

* `nmf` — silhouette score for SymNMF clustering
* `kmeans` — silhouette score for K-Means clustering
* `ari` — Adjusted Rand Index between the two cluster assignments

The ARI compares the results of the two algorithms with each other.

## Input Format

The input file should contain numerical vectors, one vector per line, with coordinates separated by commas.

All vectors must have the same number of coordinates. The file should not contain a header or an index column.

Example:

```text
1.0,2.0
1.5,2.2
8.0,9.0
8.5,9.2
```

## Output Format

The matrix programs print comma-separated rows with 4 digits after the decimal point.

* `sym`, `ddg`, and `norm` produce an `N × N` matrix
* `symnmf` produces an `N × k` association matrix

The analysis program prints three labeled scores, each with 4 digits after the decimal point.

## Algorithm Flow

1. Parse command-line arguments and read the input vectors
2. Calculate the similarity matrix from pairwise squared Euclidean distances
3. Calculate the diagonal degree matrix
4. Calculate the normalized similarity matrix `W`
5. Initialize a non-negative `N × k` matrix `H` using a fixed random seed
6. Update `H` iteratively to approximate `W ≈ H · Hᵀ`
7. Stop when the convergence threshold is reached or after 300 iterations
8. Print the resulting association matrix

For clustering analysis, each point is assigned to the cluster with the largest value in its row of `H`. The resulting assignments are then compared with K-Means.

## Project Structure

| File | Purpose |
| --- | --- |
| `symnmf.py` | Python command-line interface and initialization of `H` |
| `analysis.py` | K-Means implementation and clustering comparison |
| `symnmfmodule.c` | Python C extension interface |
| `symnmf.c` | SymNMF computations and standalone C entry point |
| `symnmf.h` | SymNMF function declarations |
| `matrix.c` / `matrix.h` | Matrix utilities and memory management |
| `input.c` / `input.h` | C input handling |
| `setup.py` | Build configuration for the Python C extension |
| `Makefile` | Build configuration for the standalone C executable |
| `tests/run_c_tests.py` | C test runner |

## Error Handling

The programs handle invalid arguments, input-reading failures, and computation errors.

Possible error messages include:

```text
Incorrect number of clusters!
An Error Has Occurred
```

## Team

This project was developed by:

* Ksenia Iaremenko
* Iris Pustylnik

## Notes

This academic project was developed by Iris Pustylnik and Ksenia Iaremenko as part of the Software Project course at Tel Aviv University.
