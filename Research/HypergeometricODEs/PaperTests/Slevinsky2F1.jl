#=
# Comparison tests for:
#   Slevinsky, R., Fast and stable rational approximation of generalized hypergeometric functions, arXiv:2307.06221
#
# Author: Caleb Jacobs
# DLM: October 2, 2024
=#

using HypergeometricFunctions: pFqweniger as weniger_pfq, pFqdrummond as drummond_pfq

"Levin-type factorial"
function weniger_2f1(a, b, c, z::Number) 
    (a,b,c,z) = convert.(ComplexF64, (a,b,c,z))
    return weniger_pfq((a,b), (c,), z)   
end

"Drummond 2F1"
function drummond_2f1(a, b, c, z::Number)
    (a,b,c,z) = convert.(ComplexF64, (a,b,c,z))
    return drummond_pfq((a,b), (c,), z)
end
