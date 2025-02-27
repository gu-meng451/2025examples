## Gauss-Legendre Rules from Golub-Welsch Algorithm
# 
# This is a very simple implementation.

using LinearAlgebra

function gauss_legendre(n)
    # Construct the Jacobi matrix
    β = @. 0.5 / sqrt(1 - (2 * (1:n-1)) ^ (-2))  # Subdiagonal elements
    T = SymTridiagonal(zeros(n), β)  # Jacobi matrix

    # Compute eigenvalues and eigenvectors
    λ, V = eigen(T)
    p = sortperm(λ)

    # Nodes are the eigenvalues
    nodes = λ[p]

    # Weights are given by the square of the first row of the eigenvectors
    weights = 2V[1, p] .^ 2

    return nodes, weights
end

# Example usage
n = 30;  # Number of quadrature points
nodes, weights = gauss_legendre(n);

println("Nodes:", nodes)
println("Weights:", weights)
