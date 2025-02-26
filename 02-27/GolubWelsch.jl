## Gauss-Legendre Rules from Golub-Welsch Algorithm
# 
# This is a very simple implementation.

using LinearAlgebra

function gauss_legendre(n)
    # Construct the Jacobi matrix
    beta = 0.5 ./ sqrt.(1 .- (2 .* (1:n-1)) .^ (-2))  # Subdiagonal elements
    J = SymTridiagonal(zeros(n), beta)  # Jacobi matrix

    # Compute eigenvalues and eigenvectors
    vals, vecs = eigen(J)

    # Nodes are the eigenvalues
    nodes = vals

    # Weights are given by the square of the first row of the eigenvectors
    weights = 2 * (vecs[1, :]) .^ 2

    return nodes, weights
end

# Example usage
n = 5;  # Number of quadrature points
nodes, weights = gauss_legendre(n);

println("Nodes:", nodes)
println("Weights:", weights)
