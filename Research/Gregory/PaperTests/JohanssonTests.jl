#=
# Comparison tests for:
#   F. Johansson. Computing hypergeometric functions rigorously. ACM Trans. Math. Software,346 45:30:1–30:26, 2019.
#
# Author: Caleb Jacobs
# DLM: May 24, 2024
=#

include("../HypergeometricTests.jl")
using Nemo: hypergeometric_1f1, hypergeometric_2f1, ComplexField # Arb/Flint interface for Julia

CC = ComplexField()                 # Complex box field
cc(z::ComplexF64) = CC(z.re, z.im)  # Convert complex floats to complex field elements

function joTest(a, b; title = "")
    r = 2.49
    n = 50
    np = 5
    Tr = .55
    cr = 8 
    sr = 10

    println("Running test with a = ", a, " and b = ", b, ".")
    (z, f, h, tru, p) = pFqTest(a, b, r = r, n = n, np = np, Tr = Tr, cr = cr, sr = sr, dir = 1, exclude = true, modifyZ1 = true)

    Z = cc.(z)
    A = cc.(complex.(a))
    B = cc.(complex.(b))

    @time F = [hypergeometric_2f1(A[1], A[2], B[1], Z) for Z ∈ Z]

    fj = [convert(ComplexF64, F) for F ∈ F]

    pj = getGraphics(z, fj, f, title = title, exclude = true)

    return (p, pj)
end

function runJohanssonTests()
    a = [[.9,1.9],
         [.9,1.1]]
    b = [[2.9],
         [1.2]]

    p  = Vector{Vector{PlotlyJS.SyncPlot}}(undef, length(a))
    pj = Vector{Vector{PlotlyJS.SyncPlot}}(undef, length(a))
    
    for (c, d, i) ∈ zip(a, b, 1 : length(a))
        (p[i], pj[i]) = joTest(c, d, title = "Nemo/Arb")
    end

    return (p, pj)
end
