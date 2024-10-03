#=
# Comparison tests for:
#   F. Johansson. Computing hypergeometric functions rigorously. ACM Trans. Math. Software,346 45:30:1–30:26, 2019.
#
# Author: Caleb Jacobs
# DLM: October 2, 2024
=#

using Nemo: hypergeometric_2f1 as nemo_2f1, ComplexField, RealField  # Arb/Flint interface for Julia

CF = ComplexField()                           # Complex box field
convert_to_field(z::Complex) = CF(z.re, z.im) # Convert complex floats to complex field elements
convert_to_field(x::Real)    = CF(x)          # Convert real to complex field elements

function johansson_2f1(a, b, c, z)::ComplexF64
    (A,B,C,Z) = convert_to_field.((a,b,c,z))

    return nemo_2f1(A, B, C, Z)
end
