# generate 10000 random matrices and generate their inverses

import numpy as np
import time

def generate_random_matrix(n):
    res = np.random.rand(n, n)
    while np.linalg.matrix_rank(res) != n:
        res = np.random.rand(n, n)
    return res

def generate_random_matrices(n, num_matrices):
    matrices = []
    for i in range(num_matrices):
        matrices.append(generate_random_matrix(n))
    return matrices

def generate_inverses(matrices):
    inverses = []
    for matrix in matrices:
        inverses.append(np.linalg.inv(matrix))
    return inverses

def main():
    n = 50
    num_matrices = n
    matrices = generate_random_matrices(n, num_matrices)
    inverses = generate_inverses(matrices)
    # conctenate matrices and inverses
    concataneted = None
    for i in range(num_matrices):
        if i == 0:
            concatanated = np.concatenate((matrices[i], inverses[i]), axis=1)
        else:
            concatanated = np.concatenate((concatanated, np.concatenate((matrices[i], inverses[i]), axis=1)), axis=0)
    np.savetxt("./test/test_inverse.csv", concatanated, delimiter=",", header=str(n))
    print("Done")

if __name__ == "__main__":
    main()