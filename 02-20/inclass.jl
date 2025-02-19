using LinearAlgebra
using Printf

include("gd.jl")
include("cg.jl")

A = [20 -1 0.0;
     -1 20 -10;
     0 -10 20];
eigvals(A)

b = [0; 1; 2];

x0 = [0; 0; 0]

x, iter = gd(A, b, x0, maxIter=1000)
A*x - b |> norm

x, iter = cg(A, b, x0)
A*x - b |> norm
