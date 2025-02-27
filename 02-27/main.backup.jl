# Integration
using Pkg
Pkg.activate(".")
using LinearAlgebra
using Plots

include("QuadRules.backup.jl")
# using .QuadRules

### Integrand:
f(x) = sin(x)
a = 0.0
b = pi

I_true = -cos(b) + cos(a)

## Trap
n_list = 2 .^ (1:8);
I_list = [QuadRules.trap(f, a, b, n) for n in n_list]
err1 = @. abs(I_list - I_true) / I_true

plt = plot(n_list, err1, lw=2,
    label="Trap",
    xlabel="n points", ylabel="Normalized Error",
    xscale=:log10, yscale=:log10)
plot!(plt, n -> n^-2, lw=2, 
    linestyle=:dash, label="O(n^-2)")

## Simpson 1/3
n_list = 2 .^ (1:8) .+ 1;
I_list = [QuadRules.simpson13(f, a, b, n) for n in n_list]
err2 = @. abs(I_list - I_true) / I_true
plot!(plt, n_list, err2, lw=2,
    label="Simpson 1/3")
plot!(plt, n -> n^-4, 
    linestyle=:dash, lw=2, label="O(n^-4)")

## Gauss-Legendre
n_list = 2 .^ (1:8);
I_list = [QuadRules.gl(f, a, b, n) for n in n_list]
err3 = @. abs(I_list - I_true) / I_true + eps()

plot!(plt, n_list, err3, lw=2,
    label="GL")