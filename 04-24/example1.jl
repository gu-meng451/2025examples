using Pkg
Pkg.activate(".")
using LinearAlgebra
using GLMakie

## Simple 3 dof example
M = diagm(0 => [1.0, 1.0, 1.0])
K = diagm(0 => [2.0, 2.0, 2.0], -1 => [-1.0, -1.0], 1 => [-1.0, -1.0])
κ1 = 0.01
κ2 = 0.01
C = κ1 * M + κ2 * K
Fext(t) = [0.0, 0.0, 0.0]

## Initial conditions
d0 = [0.0, 1.0, 2.0]
v0 = [0.0, 0.0, 0.0]

## Time step
Δt = 0.01

## Integrate the system
function integrate_step(dn, vn, an, tn, Δt, Fext, M, K, C, P1, P2, P3, P4)

    R = αf * Fext(tn) - αm * M * an - αf * C * vn - αf * K * dn
    rhs = (1 - αf) * Fext(tn + Δt) + R + P2 * an + P3 * vn + P4 * dn
    dn1 = P1 \ rhs
    an1 = (dn1 - dn - Δt * vn - Δt^2 * (1 / 2 - β) * an) / (Δt^2 * β)
    vn1 = vn + Δt * (1 - γ) * an + Δt * γ * an1
    return dn1, vn1, an1

end

function integrate_system(d0, v0, tfinal, Δt, ρ∞, Fext, M, K, C)
    dn = d0
    vn = v0

    t = 0:Δt:tfinal

    a0 = M \ (Fext(0) - C * v0 - K * d0)
    an = a0

    # let's save everything
    Dn = zeros(length(d0), length(t))
    Vn = zeros(length(d0), length(t))
    An = zeros(length(d0), length(t))

    ## Build the system matrices
    ## Build the generalized-alpha method
    ρ∞ = 0.9
    αf = ρ∞ / (1 + ρ∞)
    αm = (2ρ∞ - 1) / (1 + ρ∞)
    γ = 1 / 2 + (1 - ρ∞) / (1 + ρ∞)
    β = 1 / (1 + ρ∞)^2
    P1 = (1 - αm) / β / Δt^2 * M + (1 - αf) * γ / β / Δt * C + (1 - αf) * K
    P2 = (1 - αm) * (1 - 2β) / 2 / β * M + (1 - αf) * (γ - 2β) / 2 / β * Δt * C
    P3 = (1 - αm) / β / Δt * M + (1 - αf) * (γ - β) / β * C
    P4 = (1 - αm) / β / Δt^2 * M + (1 - αf) * γ / β / Δt * C

    for (i, tn) in enumerate(t)
        dn, vn, an = integrate_step(dn, vn, an, tn, Δt, Fext, M, K, C, P1, P2, P3, P4)
        Dn[:, i] = dn
        Vn[:, i] = vn
        An[:, i] = an
    end
    return Dn, Vn, An
end

## Run the integration
tfinal = 10.0
Dn, Vn, An = integrate_system(d0, v0, tfinal, Δt, ρ∞, Fext, M, K, C)

## Plot the results
fig = Figure()
ax = Axis(fig[1, 1], title = "Displacement", xlabel = "Time", ylabel = "Displacement")
lines!(ax, 0:Δt:tfinal, Dn[1, :], label = "DOF 1")
lines!(ax, 0:Δt:tfinal, Dn[2, :], label = "DOF 2")
lines!(ax, 0:Δt:tfinal, Dn[3, :], label = "DOF 3")
# add the legend
Legend(fig, ax)