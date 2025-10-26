using LinearAlgebra
using Polynomials

function minimax(f, df, ddf, a, b, n; maxiter = 5, atol = 1e-14)
    x0 = chebnodes(a, b, n + 2)
    
    coeffs = minimaxiterate(f, x0)
    for _ in 1:maxiter
        p = Polynomial(coeffs[1:end-1])
        dp = derivative(p)
        ddp = derivative(dp)

        g(x) = f(x) - p(x)
        dg(x) = df(x) - dp(x)
        ddg(x) = ddf(x) - ddp(x)

        zinner = [newtonroot(g, dg, x0[k], x0[k+1]) for k ∈ 1:n+1]
        z = vcat(a, zinner, b)

        xold = copy(x0)
        x0 = [newtonopt(g, dg, ddg, z[k], z[k+1]; atol = atol, maxiter = maxiter) for k ∈ 1:n+2]

        coeffs = minimaxiterate(f, x0)

        if norm(xold - x0) / norm(x0) < atol
            println("Converged early!")
            break
        end
    end

    return Polynomial(coeffs[1:end-1])
end

function chebnodes(a, b, n)
    N = (n-1):-1:0
    x = cospi.((N .+ 1/2) / n)

    return ((a + b) .+ (b - a) * x) / 2
end

function newtonroot(f, df, a, b; atol = 1e-14, maxiter = 20)
    x0 = (a + b) / 2

    for _ ∈ 1:maxiter
        if abs(f(x0)) < atol
            return x0
        end

        x0 -= f(x0) / df(x0)

        if x0 < a || b < x0
            x0 = abs(f(a)) < abs(f(b)) ? a : b
            return x0
        end
    end

    return x0
end

function newtonopt(f, df, ddf, a, b; atol = 1e-14, maxiter = 20)
    x0 = (a + b) / 2

    for _ in 1:maxiter
        xprev = x0
        x0 -= df(x0) / ddf(x0)

        if x0 < a || b < x0
            return abs(f(a)) > abs(f(b)) ? a : b
        end

        if abs(xprev - x0) < atol
            return x0
        end
    end

    return x0
end

function minimaxiterate(f, x0)
    n = length(x0) - 2

    A0 = x0.^(transpose(0:n))
    A1 = (-1).^(0:n+1)

    A = hcat(A0, A1)
    b = f.(x0)

    return A \ b
end
