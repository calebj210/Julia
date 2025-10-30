#=
#   Transformations for 2F1
#
# Author: Caleb Jacobs
# DLM: October 30, 2025
=#

using SpecialFunctions
include("ConformalBase.jl")

function z_2f1(a, b, c, z, ord)
    p = a + b - c

    if real(p) < 0 
        return conformal_2f1(a, b, c, z, ord)
    else
        return (1 - z)^(-p) * conformal_2f1(c - a, c - b, c, z, ord)
    end
end

function zoverzminusone_2f1(a, b, c, z, ord)
    if real(a) < real(b)
        g = (1 - z)^(-a)
        b = c - b
    else
        g = (1 - z)^(-b)
        a = c - a
    end

    f = z_2f1(a, b, c, z / (z - 1), ord)

    return g * f
end

function oneminusoneoverz_2f1(a, b, c, z, ord)
    g1 = gamma(c - a - b) * gamma(c) / gamma(c - a) / gamma(c - b) * z^(-a)
    g2 = gamma(a + b - c) * gamma(c) / gamma(a) / gamma(b) * (1 - z)^(c - a - b) * z^(a - c)

    f1 = z_2f1(a, a - c + 1, a + b - c + 1, 1 - 1 / z, ord)
    f2 = z_2f1(c - a, 1 - a, c - a - b + 1, 1 - 1 / z, ord)

    return g1 * f1 + g2 * f2
end

function oneminusz_2f1(a, b, c, z, ord)
    g1 = gamma(c - a - b) * gamma(c) / gamma(c - a) / gamma(c - b)
    g2 = gamma(a + b - c) * gamma(c) / gamma(a) / gamma(b) * (1 - z)^(c - a - b)

    f1 = z_2f1(a, b,         a + b - c + 1, 1 - z, ord)
    f2 = z_2f1(c - a, c - b, c - a - b + 1, 1 - z, ord)
    
    return g1 * f1 + g2 * f2
end

function oneoverz_2f1(a, b, c, z, ord)
    g1 = gamma(b - a) / gamma(b) * gamma(c) / gamma(c - a) * (-z)^(-a)
    g2 = gamma(a - b) / gamma(a) * gamma(c) / gamma(c - b) * (-z)^(-b)

    f1 = z_2f1(a, a - c + 1, a - b + 1, 1 / z, ord)
    f2 = z_2f1(b, b - c + 1, b - a + 1, 1 / z, ord)

    return g1 * f1 + g2 * f2
end

function oneoveroneminusz_2f1(a, b, c, z, ord)
    g1 = gamma(b - a) / gamma(b) * gamma(c) / gamma(c - a) * (1 - z)^(-a)
    g2 = gamma(a - b) / gamma(a) * gamma(c) / gamma(c - b) * (1 - z)^(-b)

    f1 = z_2f1(a, c - b, a - b + 1, 1 / (1 - z), ord)
    f2 = z_2f1(b, c - a, b - a + 1, 1 / (1 - z), ord)

    return g1 * f1 + g2 * f2
end

transformations = (;
    z = z_2f1,

    zoverzminusone = zoverzminusone_2f1,

    oneminusz = oneminusz_2f1,
    oneminusoneoverz = oneminusoneoverz_2f1,

    oneoverz = oneoverz_2f1,
    oneoveroneminusz = oneoveroneminusz_2f1,
)
