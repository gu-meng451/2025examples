## simple bar
#= 
Consider a uniform bar that is fixed to the wall at x = 0, and free 
at x= L.  It has a constant loading f(x) = f0.
We want to find the displacement and stress along the bar.
=#

##
using Pkg
Pkg.activate(".")

include("Preprocess.jl")
include("Quadrature.jl")
include("Bar1D.jl")
include("SimpleVisualization.jl")
using GLMakie

quad_rules = Dict("line" => Quadrature.gauss_legendre_1d(2))

## Properties of the uniform body
prop = (E=1.0, # Young's Modulus [F/m^2]
    A=1.0, # cross-sectional area [L^2]
    L=1.0, # length [L]
    f0=1.0, # loading [F/L]
) 

## mesh connectivity
# IEN(e,a)
nnp = 5 # number of nodes
nel = nnp - 1 # number of elements
nee = 2 # number of equations per element
IEN = Dict("line" => zeros(Int, nel, nee))
# let's use the same ordering from class
IEN["line"][1,:] = [1, 3]
IEN["line"][2,:] = [3, 4]
IEN["line"][3,:] = [4, 5]
IEN["line"][4,:] = [5, 2]

# make the list of node locations
x = LinRange(0, prop.L, nnp) |> collect
# fix them to match where we put the nodes
x[:] = x[[1,5,2,3,4]]

## Alternate: using mode nodes
nnp = 300 # number of nodes
nel = nnp - 1 # number of elements
nee = 2 # number of equations per element
IEN = Dict("line" => zeros(Int, nel, nee))
for e in 1:nel
    IEN["line"][e, :] = [e, e + 1]
end
x = LinRange(0, prop.L, nnp) |> collect

## Essential Boundary Conditions: [i,A]
BC_fix_list = zeros(Bool, 1, nnp)
BC_fix_list[1, 1] = true

# BC values
BC_g_list = zeros(1, nnp)
BC_g_list[1, 1] = 0.0

## Loading force
f(x) = prop.f0

## Build mesh
m = Preprocess.build_mesh(x, [], [], IEN, 1, BC_fix_list, BC_g_list)

## Assemble the global matrices
K = Bar1D.assemble_stiffness(m, prop)
F = Bar1D.assemble_rhs(m, f, quad_rules)

## Solve the system
q = zeros(m.nnp * m.ned)

# apply essential BCs in q[r2]
idx = findall(BC_fix_list)
q[ m.ID[idx] ] = BC_g_list[idx]  

# solve
r1 = m.free_range
r2 = m.freefix_range
q[r1] = K[r1, r1] \ (F[r1] - K[r1, r2] * q[r2])

## Plot the result
fig, ax = SimpleVisualization.init_plot()
SimpleVisualization.draw_element(ax, m, q, linewidth=3, linestyle=:dash, color=:red)
u_true(x) = prop.f0 / (2 * prop.E*prop.A) * (2*prop.L - x)*x
lines!(ax, LinRange(0, prop.L, 20), u_true)
fig

## Plot the stress
σ_true(x) = prop.f0/prop.A*(prop.L-x)
fig, ax = SimpleVisualization.init_plot(ylabel="Stress")
SimpleVisualization.draw_element_stress(ax, m, q, prop)
lines!(ax, LinRange(0, prop.L, 20), σ_true)
fig

# show sparse matrix of K, flip the y axis
fig = Figure()
ax = Axis(fig[1,1], yreversed=true)
spy!(ax,K)
fig
