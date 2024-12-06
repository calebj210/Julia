using SpecialFunctions

"""
    U(a, b, z,; n, m)
Compute the Krummer U hypergeometric function using the exponentially-improved expansion given in https://dlmf.nist.gov/13.7.E10.
"""
function U(a, b, z; n = 10, m = 10) 
    as  = a .+ (0:n - 2)
    bs  = a - b + 1 .+ (0:n - 2)
    den = 1:n - 1
    zs  = repeat([-1 / z], n - 1)

    partial_sum = 1 / z^a * (1 + sum(cumprod(as ./ den .* bs .* zs)))

    return partial_sum + Rn(n, m, a, b, z)
end

function Rn(n, m, a, b, z)
    as  = 1 - a .+ (0:m - 2)
    bs  = a - b .+ (0:m - 2)
    den = 1:m - 1
    zs  = repeat([-1 / z], m - 1)
    Gs  = Gp.(n + 2a - b .- (1:m - 1), z)

    partial_sum = Gp(n + 2a - b, z) + sum(cumprod(as ./ den .* bs .* zs .* Gs)) # + remainder

    return (-1)^n * 2π * z^(a - b) / gamma(a) / gamma(a - b + 1) * partial_sum
end

Gp(p::Number, z::Number) = exp(z) / (2π) * gamma(p) * gamma(1 - p, z)
