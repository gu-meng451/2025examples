import numpy as np


class Quadrule:
    def __init__(self, label, n, dim, ξ, w):
        self.label = label
        self.n = n
        self.dim = dim
        self.ξ = ξ
        self.w = w
        self.iterator = list(zip(ξ, w))  # Creating an iterable of nodes and weights


def gauss_legendre_1d(n):
    """Computes the 1D Gauss-Legendre quadrature rule."""
    beta = 0.5 / np.sqrt(1 - (2 * np.arange(1.0, n)) ** -2)  # Compute beta coefficients
    T = np.diag(beta, -1) + np.diag(beta, 1)  # Construct symmetric tridiagonal matrix
    λ, V = np.linalg.eigh(T)  # Compute eigenvalues and eigenvectors
    p = np.argsort(λ)  # Sort eigenvalues

    ξ = λ[p]  # Sorted nodes
    w = 2 * (V[0, p] ** 2)  # Compute weights

    return Quadrule("1D GL", n, 1, ξ, w)

