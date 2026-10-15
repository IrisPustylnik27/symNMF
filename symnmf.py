import sys
from math import sqrt

import numpy as np

import symnmfmodule


def main():
    """Parse command-line input, run the selected goal, and print its result."""
    k, N, goal, vectors = parse_arguments()
    ans_matrix = run_goal(vectors, k, N, goal)
    print_matrix(ans_matrix)


def run_goal(vectors, k, N, goal):
    """Run the requested SymNMF goal and return its resulting matrix."""
    if goal == "symnmf":
        try:
            laplacian = symnmfmodule.norm(vectors)
        except Exception:
            exception_handler()
        decomposition_matrix = initialize_H(laplacian, k, N)

        try:
            ans_matrix = symnmfmodule.symnmf(decomposition_matrix, laplacian)
        except Exception:
            exception_handler()
    elif goal == "sym":
        try:
            ans_matrix = symnmfmodule.sym(vectors)
        except Exception:
            exception_handler()
    elif goal == "ddg":
        try:
            ans_matrix = symnmfmodule.ddg(vectors)
        except Exception:
            exception_handler()
    elif goal == "norm":
        try:
            ans_matrix = symnmfmodule.norm(vectors)
        except Exception:
            exception_handler()
    else:
        exception_handler()

    return ans_matrix


def initialize_H(laplacian, k, N):
    """Create the randomized N-by-k initial H matrix from normalized similarity."""
    np.random.seed(1234)
    avg = np.array(laplacian).mean()
    low, high = 0, 2 * sqrt(avg / k)
    return np.random.uniform(low, high, size=(N, k)).tolist()


def exception_handler():
    """Print the standard error message and terminate with a failure status."""
    print("An Error Has Occurred")
    sys.exit(1)


def parse_arguments():
    """Validate CLI arguments and return k, N, goal, and the input vectors."""
    k, N, goal, file_name = 0, 0, None, None

    if len(sys.argv) != 4:
        exception_handler()

    try:
        k = int(sys.argv[1])
    except ValueError:
        print("Incorrect number of clusters!")
        sys.exit(1)
    goal = sys.argv[2]
    file_name = sys.argv[3]

    try:
        vectors = np.loadtxt(file_name, delimiter=",", ndmin=2).tolist()
        N = len(vectors)
        if N == 0:
            exception_handler()
    except Exception:
        exception_handler()

    if k <= 1 or k >= N:
        print("Incorrect number of clusters!")
        sys.exit(1)

    return k, N, goal, vectors


def fit(X, k):
    """Compute and return the SymNMF association matrix for X and k."""
    vectors = np.array(X, dtype=float, ndmin=2).tolist()
    N = len(vectors)
    ans_matrix = run_goal(vectors, k, N, "symnmf")
    return ans_matrix


def print_matrix(matrix):
    """Print a matrix as comma-separated rows with four decimal places."""
    for row in matrix:
        print(",".join(f"{val:.4f}" for val in row))


if __name__ == "__main__":
    main()
