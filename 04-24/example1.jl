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
Δt = 0.1
tfinal = 100.0
ρ∞ = 0.5

## Trap Rule
function trap_step(dn, vn, an, tn, Δt, Fext, M, K, C)
    rhs = Fext(tn+Δt) + C*(2/Δt*dn + vn) + M*(4/Δt^2*dn + 4/Δt*vn + an)
    P1 = 4/Δt^2 * M + 2/Δt * C + K
    dn1 = P1 \ rhs
    vn1 = 2/Δt * (dn1 - dn) - vn
    an1 = 4/Δt^2 * (dn1 - dn) - 4/Δt * vn - an
    return dn1, vn1, an1
end

function trap(d0, v0, tfinal, Δt, Fext, M, K, C)
    dn = d0
    vn = v0

    t = 0:Δt:tfinal

    a0 = M \ (Fext(0) - C * v0 - K * d0)
    an = a0

    # let's save everything
    Dn = zeros(length(d0), length(t))
    Vn = zeros(length(d0), length(t))
    An = zeros(length(d0), length(t))

    for (i, tn) in enumerate(t)
        dn, vn, an = trap_step(dn, vn, an, tn, Δt, Fext, M, K, C)

        Dn[:, i] = dn
        Vn[:, i] = vn
        An[:, i] = an
    end
    return Dn, Vn, An
end



## Generalized alpha method
function genα_step(dn, vn, an, tn, Δt, Fext, M, K, C, 
    P1, P2, P3, P4, αf, αm, γ , β )

    R = αf * Fext(tn) - αm * M * an - αf * C * vn - αf * K * dn
    rhs = (1 - αf) * Fext(tn + Δt) + R + P2 * an + P3 * vn + P4 * dn
    dn1 = P1 \ rhs
    an1 = (dn1 - dn - Δt * vn - Δt^2 * (1 / 2 - β) * an) / (Δt^2 * β)
    vn1 = vn + Δt * (1 - γ) * an + Δt * γ * an1
    return dn1, vn1, an1

end

function genα(d0, v0, tfinal, Δt, ρ∞, Fext, M, K, C)
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
    αf = ρ∞ / (1 + ρ∞)
    αm = (2ρ∞ - 1) / (1 + ρ∞)
    γ = 1 / 2 + (1 - ρ∞) / (1 + ρ∞)
    β = 1 / (1 + ρ∞)^2
    P1 = (1 - αm) / β / Δt^2 * M + (1 - αf) * γ / β / Δt * C + (1 - αf) * K
    P2 = (1 - αm) * (1 - 2β) / 2 / β * M + (1 - αf) * (γ - 2β) / 2 / β * Δt * C
    P3 = (1 - αm) / β / Δt * M + (1 - αf) * (γ - β) / β * C
    P4 = (1 - αm) / β / Δt^2 * M + (1 - αf) * γ / β / Δt * C

    for (i, tn) in enumerate(t)
        dn, vn, an = genα_step(dn, vn, an, tn, Δt, Fext, M, K, C, P1, P2, P3, P4, αf, αm, γ , β )

        Dn[:, i] = dn
        Vn[:, i] = vn
        An[:, i] = an
    end
    return Dn, Vn, An
end

## Run the integration
Dn_trap, Vn_trap, An_trap = trap(d0, v0, tfinal, Δt, Fext, M, K, C)
Dn_ga, Vn_ga, An_ga = genα(d0, v0, tfinal, Δt, ρ∞, Fext, M, K, C)

## Plot the results
fig = Figure()
ax = Axis(fig[1, 1], title = "Displacement", xlabel = "Time", ylabel = "Displacement")

lines!(ax, 0:Δt:tfinal, Dn_trap[1, :], label="trap DOF 1")
lines!(ax, 0:Δt:tfinal, Dn_trap[2, :], label="trap DOF 2")
lines!(ax, 0:Δt:tfinal, Dn_trap[3, :], label="trap DOF 3")


lines!(ax, 0:Δt:tfinal, Dn_ga[1, :], label="gen α DOF 1")
lines!(ax, 0:Δt:tfinal, Dn_ga[2, :], label="gen α DOF 2")
lines!(ax, 0:Δt:tfinal, Dn_ga[3, :], label="gen α DOF 3")
# add the legend
Legend(fig, ax)