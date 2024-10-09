#=
# Comparison tests for:
#   F. Johansson. Computing hypergeometric functions rigorously. ACM Trans. Math. Software,346 45:30:1–30:26, 2019.
#
# Author: Caleb Jacobs
# DLM: October 9, 2024
=#

include("../HypergeometricTests.jl")
using Nemo: hypergeometric_1f1, hypergeometric_2f1, ComplexField # Arb/Flint interface for Julia

CC = ComplexField()                 # Complex box field
cc(z::ComplexF64) = CC(z.re, z.im)  # Convert complex floats to complex field elements

function joTest(a, b; title = "")
    r = 2.49; n = 50; Tr = .55; np = 5
    corrR = .5; innerR = .6; outerR = .8
    z1N = 70; modifyZ1 = true
    a = [.9, 1.9]; b = [1.91]

    test = (a, b, r, n, np, Tr, modifyZ1, corrR, innerR, outerR, z1N)

    println("Running test with a = ", a, " and b = ", b, ".")
    (z, f, h, tru, p) = gridtest(test...)

    Z = cc.(z)
    A = cc.(complex.(a))
    B = cc.(complex.(b))

    @time F = hypergeometric_2f1.(A[1], A[2], B[1], Z)

    fj = convert.(ComplexF64, F)

    pj = generate_graphics(z, fj, tru, title = title)

    return (p, pj)
end

function runJohanssonTests()
    a = [[.9,1.9],
         [.9,1.1]]
    b = [[2.9],
         [1.2]]

    p  = Vector{NTuple{5, PlotlyJS.SyncPlot}}(undef, length(a))
    ps = Vector{NTuple{5, PlotlyJS.SyncPlot}}(undef, length(a))
    
    for (c, d, i) ∈ zip(a, b, 1 : length(a))
        (p[i], ps[i]) = joTest(c, d, title = "Nemo/Arb")
    end

    savefig(p[1][1],  string("Images/Johansson/ECTrap1.png"),       width = 700, height = 700)
    savefig(p[1][3],  string("Images/Johansson/ECTrap1AbsArg.png"), width = 700, height = 700)
    savefig(p[1][4],  string("Images/Johansson/ECTrap1Re.png") ,    width = 700, height = 700)
    savefig(p[1][5],  string("Images/Johansson/ECTrap1Im.png") ,    width = 700, height = 700)
    savefig(ps[1][1], string("Images/Johansson/Jo1.png"),           width = 700, height = 700)

    savefig(p[2][1],  string("Images/Johansson/ECTrap2.png"), width = 700, height = 700)
    savefig(p[2][4],  string("Images/Johansson/ECTrap2AbsArg.png"), width = 700, height = 700)
    savefig(p[2][4],  string("Images/Johansson/ECTrap2Re.png"), width = 700, height = 700)
    savefig(p[2][5],  string("Images/Johansson/ECTrap2Im.png"), width = 700, height = 700)
    savefig(ps[2][1], string("Images/Johansson/Jo2.png"), width = 700, height = 700)
    
    return (p, ps)
end
