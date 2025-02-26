using Pkg;
Pkg.activate(".");
using Plots

include("Spline.jl")
include("LagrangeInterp.jl")
using .Spline

X = [-1, 2, 3, 4, 7]
f(x) = sin(x)
Y = f.(X)

S = spinterp(X, Y)
S(0.1)

L = LagrangeInterp.interp(X,Y)

scatter(X, Y, marker=3, label="control points",
    xlabel="x", ylabel="y(x)")
plot!(S, -1, 7, label="Spline", lw=5)
plot!(L, -1, 7, label="Lagrange", lw=5)
plot!(f, -1, 7, label="Original function", lw=5)
