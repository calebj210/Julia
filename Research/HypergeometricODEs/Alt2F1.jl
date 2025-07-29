#=
# 2F1 computation via the other linearly independent solution
#
# Author: Caleb Jacobs
# DLM: July 16, 2025
=#

include("pFq.jl")
using OrdinaryDiffEq
using FastGaussQuadrature

function ivp_2f1(a, b, c, z)
    z0, u0 = initialize(a, b, c, z, .5)
    z0 = [z0]
    u0 = [first(u0)]
    if z0 == z
        return u0[1]
    end
    tspan = (0, abs(z - z0[1]))

    Z(t) = z0[1] + sign(z - z0[1]) * t

    W(z) = (1 - c) * z^(-c) * (1 - z)^(c - a - b - 1)

    f2(z) = z^(1 - c) * taylor_2f1(a - c + 1, b - c + 1, 2 - c, z)
    df2(z) = (1 - c) / z * f2(z) + z^(1 - c) * (a - c + 1) * (b - c + 1) / (2 - c) * taylor_2f1(a - c + 2, b - c + 2, 3 - c, z)

    P(z) = df2(z) / f2(z)
    Q(z) = -W(z) / f2(z)

    function f(du, u, p, t)
        du[1] = (P(Z(t)) * u[1] + Q(Z(t))) / conj(sign(z - z0[1]))
    end

    prob = ODEProblem(f, u0, tspan)
    return solve(prob, Tsit5();
                 save_on = false, 
                 save_start = false, 
                 save_end = true, 
                 reltol = eps(),
     ).u[1][1]
end

glweights = gausslegendre(100)

function integral_2f1(a, b, c, z; N = 10)
    z0, f1 = initialize(a, b, c, z, .5)
    f1 = f1[1]
    if z0 == z
        return f1
    end

    Z(t) = (z0 * (1 - t) + z * (t + 1)) / 2
    dz = (z - z0) / 2
    W(z) = (1 - c) * z^(-c) * (1 - z)^(c - a - b - 1)
    f2(z) = z^(1 - c) * taylor_2f1(a - c + 1, b - c + 1, 2 - c, z)

    F(z) = -W(z) / (f2(z))^2

    t, w = glweights

    int = sum(w .* F.(Z.(t)))

    return f2(z) * (f1 / f2(z0) + int * dz)
end
