#=
# Elements of the ambient space
#
# Author Caleb Jacobs
# DLM: 2022-06-09
=#

using ForwardDiff
using LinearAlgebra

# Simplify ForwardDiff notation
∂   = ForwardDiff.derivative
∇   = ForwardDiff.gradient
jac = ForwardDiff.jacobian

# Begin defining ambient space elements
d = 3                                # Dimension of ambient space
R(z::Vector)   = [z[1] * cos(z[2]), 
                  z[1] * sin(z[2]), 
                  z[3]]  # Position vector
J(z::Vector)   = jac(R, z)           # Coordinate jacobian J_{i'}^i 
JJ(z::Vector)  = J(z)^(-1)           # inverse of ^
Z(z::Vector)   = [J(z)[:, i] for i ∈ 1:d]
MZ(z::Vector)  = [Z(z)[i] ⋅ Z(z)[j] for i ∈ 1:d, j ∈ 1:d]
MZZ(z::Vector) = MZ(z)^(-1)
ZZ(z::Vector)  = [sum(j -> MZZ(z)[i, j] * Z(z)[j], 1:d) for i ∈ 1:d]
V(z::Vector)   = det(MZ(z))
