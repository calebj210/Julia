#=
# Comparison tests for:
#   Slevinsky, R., Fast and stable rational approximation of generalized hypergeometric functions, arXiv:2307.06221
#
# Author: Caleb Jacobs
# DLM: October 9, 2024
=#

import HypergeometricFunctions: pFqweniger, pFqdrummond
include("../HypergeometricTests.jl")

function runSlevinskyTests()
    r  = 19.99; n  = 81; Tr = 1.4; np = 5
    modifyZ1 = false
    a = [1.25]; b = [1.5]

    test = (a, b, r, n, np, Tr, modifyZ1)
    
    println("Running test with a = ", a, " and b = ", b, ".")
    @time (z, f, h, tru, p) = gridtest(test...)

    testF = Vector{Vector{ComplexF64}}(undef, 2)
    @time testF[1] = [pFqdrummond(Tuple(a), Tuple(b), z) for z ∈ z]             # Drummond
    @time testF[2] = [pFqweniger(Tuple(a), Tuple(b), z) for z ∈ z]              # Levin-type factorial

    titles = ["Drummond", "Factorial Levin-Type"]

    ps = [generate_graphics(z, F, tru, title = title) for (F, title) = zip(testF, titles)]

    for (img, fileName) ∈ zip(ps, titles)
        savefig(img[2], string("Images/Slevinsky/", fileName, ".png"), width = 700, height = 700)
    end
    savefig(p[2], string("Images/Slevinsky/ECTrap.png"),       width = 700, height = 700)
    savefig(p[3], string("Images/Slevinsky/ECTrapAbsArg.png"), width = 700, height = 700)
    savefig(p[4], string("Images/Slevinsky/ECTrapRe.png"),     width = 700, height = 700)
    savefig(p[5], string("Images/Slevinsky/ECTrapIm.png"),     width = 700, height = 700)

    return (p, ps)
end
