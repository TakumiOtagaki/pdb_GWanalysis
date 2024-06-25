import numpy as np
from scipy.special import logsumexp
# import ot
def my_sinkhorn(a, b, C, epsilon, max_iter=100000):
    n, m = C.shape
    f = np.ones(n)
    g = np.ones(m)
    for t in range(max_iter):
        if t % 1000 == 0:
            print(f"Iteration {t} / {max_iter}")
        # f は g から
        f = epsilon * np.log(a) - epsilon * logsumexp((-C + g.reshape(1, m)) / epsilon, axis=1) 
        # g は f から
        #         g = -epsilon * torch.logsumexp((-C + f.reshape(n, 1)) / epsilon, dim=0) + epsilon * torch.log(q)
        g = epsilon * np.log(b) - epsilon * logsumexp((-C + f.reshape(n, 1)) / epsilon, axis=0)
    P = np.exp((-C + f.reshape(n, 1) + g.reshape(1, m)) / epsilon)
    return P

def check_convergence(P, C, small_value = 1e-3):
    return np.linalg.norm(P.sum(axis=1) - 1) < small_value and np.linalg.norm(P.sum(axis=0) - 1) < small_value

def my_gromov_wasserstein_distance2(matrix1, matrix2, max_iter=100000, epsilon=1e-4, lambda_=1e1):
    n, m = matrix1.shape[0], matrix2.shape[0]
    a = np.ones(n) / n
    b = np.ones(m) / m

    P = np.ones((n, m)) / (n * m)
    for t in range(max_iter):
        C = -4 * matrix1 @ P @ matrix2.T + (epsilon - lambda_) * np.log(P)
        P = my_sinkhorn(a, b, C, lambda_)
        if check_convergence(P, C, epsilon):
            break
    return np.sum(P * C)