#=
# Surace tensor elements
#
# Author: Caleb Jacobs
# DLM: 2022-06-10
=#

using ForwardDiff
using LinearAlgebra

# Basic differential operators
∇   = ForwardDiff.gradient
jac = ForwardDiff.jacobian
hes = ForwardDiff.hessian

# Basic surface elements
RS(s::Vector, z) = vcat(s, z(s))                    # Surface position vector
function MS(s, z)                                   # Covariant metric tensor
    d = length(s)
    J = jac(x -> RS(x, z), s)

    return [J[:, i] ⋅ J[:, j] for i ∈ 1:d, j ∈ 1:d]
end
MSS(s::Vector, z) = MS(s, z)^(-1)                   # Contravariant metric tensor
A(s::Vector, z) = det(MS(s, z))                     # Area element

# Surface operators
function Δ(s::Vector, z, f)                         # Laplace-Beltrami operator
    d = length(s)
    return A(s, z)^(-1/2) * sum(α -> ∇(x -> sqrt(A(x, z)) * 
                            sum(β -> MSS(x, z)[α, β] * ∇(f, x)[β], 
                            1:d), s)[α], 1:d)
end

Δ²(s, z, f) = Δ(s, z, x -> Δ(x, z, f))              # Biharmonic operator

function ∇Γ(s::Vector, z, f)                        # Surface gradient
    d = length(s)
    g = MSS(s, z)
    J = jac(x -> RS(x, z), s)
    ∇F = ∇(f, s)
    
    return sum(i -> sum(j -> g[i, j] * J[:, i] * ∇F[j], 1:d), 1:d)
end

function ∇Γdot(s::Vector, z, F)                     # Surface divergence
    d = length(s)
    g = MSS(s, z)
    JS = jac(x -> RS(x, z), s)
    JF = jac(F, s)
    
    return sum(i -> sum(j -> g[i, j] * JS[:, i] ⋅ JF[:, j], 1:d), 1:d)
end


# 3D ambient space specific formulas
function normal(s::Vector, z)                       # Normal to surface
    J = jac(x -> RS(x, z), s)
    
    return normalize(J[:, 1] × J[:, 2])
end

function H(s::Vector, z)                             # Mean curvature
    ∂z  = ∇(z, s)   # First partials of z
    ∂²z = hes(z, s) # Second partials of z

    return  ((1 + ∂z[1]^2) * ∂²z[2, 2] - 2 * ∂z[1] * ∂z[2] * ∂²z[1, 2] + (1 + ∂z[2]^2) * ∂²z[1, 1]) /
            (2 * (1 + ∂z[1]^2 + ∂z[2]^2)^(3/2))
end
