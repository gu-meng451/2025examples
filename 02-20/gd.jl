"""
Gradient Descent Ax=b
A must be symmetric
https://en.wikipedia.org/wiki/Gradient_descent
"""
function gd(A, b, x0; tol=1e-6, maxIter=2 * length(b)^2, verbose=true)

    # should include check to ensure A is symmetric

    iter = 0
    r = b - A * x0

    if norm(r) <= tol
        return (x0, 0)
    end

    x = copy(x0)

    flag = 0
    while flag == 0
        iter += 1

        Ar = A * r

        γ = (r' * r) / (r' * Ar)
        x += γ * r
        r1 = r - γ * Ar

        if verbose
            @printf("iter=%3d, |res|=%.3e\n", iter, norm(r1))
        end

        if norm(r1) <= tol
            return x, iter
        elseif iter >= maxIter
            error("Failed to converge")
        end

        r = copy(r1)

    end

end