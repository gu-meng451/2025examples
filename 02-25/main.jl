using Pkg;
Pkg.activate(".");
using Plots

include("Spline.jl")
using .Spline

X = [-1, 2, 3, 4, 7]
f(x) = sin(x)
Y = f.(X)

S = spinterp(X, Y)
S(0.1)

scatter(X, Y, marker=3, label="control points",
    xlabel="x", ylabel="y(x)")
plot!(S, -1, 7, label="Spline")
plot!(f, -1, 7, label="Original function")
