module QuadRules

using LinearAlgebra

export trap, simpson13, gl

function trap(f, a, b, n_points::Int)

    n_seg = n_points - 1
    x = LinRange(a, b, n_points)
    F = f.(x)
    h = (b - a) / n_seg

    I = 0.0

    I += F[1] * h / 2
    I += sum(F[2:n_points-1] * h)
    I += F[n_points] * h / 2

    return I
end

function simpson13(f::Function, a::Real, b::Real, n_points::Int)
    # ensure we have an odd number of n_points
    @assert isodd(n_points) && n_points > 2

    #
    n_seg = n_points - 1
    h = (b-a)/n_seg;
    x(i) = a + i*h;

    I = 0.0
    I += h/3*f(x(0))

    # odd steps: i = 1,3,5...
    idx = 1:round(Int, n_seg/2)
    for i in idx
        I += 4/3*h*f(x(2i-1))
    end

    # even steps
    idx = 1:round(Int, n_seg / 2 - 1)
    for i in idx
        I += 2 / 3 * h * f(x(2i))
    end

    # last step
    I += h/3*f(b)

    return I

end

function gauss_legendre_rule(n)
    # Construct the Jacobi matrix
    β = @. 0.5 / sqrt(1 - (2 * (1:n-1))^(-2))  # Subdiagonal elements
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

function gl(f::Function, a::Real, b::Real, n_points::Int)
    
    Ξ, W = gauss_legendre_rule(n_points);
    I = 0.0;
    J = (b-a)/2;

    for (ξ,w) in zip(Ξ,W)
        t = (b-a)/2*ξ + (a+b)/2;
        I += f(t)*w*J;
    end
    return I
end


end