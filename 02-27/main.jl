# Integration
using Pkg
Pkg.activate(".")
using LinearAlgebra
using Plots

include("QuadRules.jl")
using .QuadRules

### Integrand:
f(x) = sin(x)
a = 0.0
b = pi

I_true = -cos(b) + cos(a)

## Trap
n_list = 2 .^ (1:8);
I_list = [trap(f, a, b, n) for n in n_list]
err1 = @. abs(I_list - I_true) / I_true

plt = plot(n_list, err1, lw=2,
    label="Trap",
    xlabel="n points", ylabel="Normalized Error",
    xscale=:log10, yscale=:log10)
plot!(plt, n -> n^-2, lw=2, 
    linestyle=:dash, label="O(n^-2)")
