#=
# Comparison tests for:
#   F. Johansson. Computing hypergeometric functions rigorously. ACM Trans. Math. Software,346 45:30:1–30:26, 2019.
#
# Author: Caleb Jacobs
# DLM: November 4, 2024
=#

# using Nemo: hypergeometric_2f1 as nemo_2f1, ComplexField, RealField  # Arb/Flint interface for Julia

using ArbNumerics: hypergeometric_2F1 as arb_2f1, ArbComplex, ArbFloat

johansson_2f1(a, b, c, z; bits = 2048)::ComplexF64 = arb_2f1(ArbComplex.((a, b, c, z), bits = bits)...)
    
