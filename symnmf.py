from math import sqrt
import pandas as pd
import symnmfmodule
import numpy as np
import sys

np.random.seed(1234)

def main():
    k, N, goal, file_name = 0, 0, None, None
    
    # taking args
    if len(sys.argv) != 4:
        print("An Error Has Occurred")
        sys.exit(1)  

    try: 
        k = int(sys.argv[1])
        goal = str(sys.argv[2])
        file_name = sys.argv[3]

    except ValueError:
        print("An Error Has Occurred")
        sys.exit(1)

    try: 
        # read data from file
        df = pd.read_csv(file_name, header=None)
        vectors = df.values.tolist()
        N = len(vectors)
        if N == 0:
            print("An Error Has Occurred")
            sys.exit(1)
    except Exception:
            print("An Error Has Occurred")
            sys.exit(1)

    if k < 1 or k >= N:
        print("An Error Has Occurred")
        sys.exit(1)

    if goal == "symnmf":
        laplacian = symnmfmodule.norm(vectors)
        avg = np.array(laplacian).mean()
        low, high = 0, 2*sqrt(avg/k)
        decomposition_matrix = np.random.uniform(low, high, size=(N, k)).tolist()
        ans_matrix = symnmfmodule.symnmf(decomposition_matrix, laplacian)
        
    elif goal == "sym":
        ans_matrix = symnmfmodule.sym(vectors)

    elif goal == "ddg":
        ans_matrix = symnmfmodule.ddg(vectors)

    elif goal == "norm":
        ans_matrix = symnmfmodule.norm(vectors)

    else:
        print("An Error Has Occurred")
        sys.exit(1)

    print_matrix(ans_matrix)
    

def fit(X, k):
    N = len(X)
    laplacian_matrix = symnmfmodule.norm(X)
    avg = np.array(laplacian_matrix).mean()
    low, high = 0, 2*sqrt(avg/k)
    decomposition_matrix = np.random.uniform(low, high, size=(N, k)).tolist()
    ans_matrix = symnmfmodule.symnmf(decomposition_matrix, laplacian_matrix)
    return ans_matrix


def print_matrix(matrix):
    for row in matrix:
        print(",".join(f"{val:.4f}" for val in row))

if __name__ == "__main__":
    main()