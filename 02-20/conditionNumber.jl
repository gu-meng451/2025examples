using LinearAlgebra
using Printf

function Hilbert(n::Int)
    return [1//(i+j-1) for i in 1:n, j in 1:n]
end;

n = 4;
A = Hilbert(n)

b = rand(n)
κ₂ = cond(A);
@printf("κ₂ = %.3e\n", κ₂)

δb = 1e-4
δA = 1e-4

κ₂ * (δA / norm(A) + δb / norm(b))

x = A\b

A*x-b |> norm
