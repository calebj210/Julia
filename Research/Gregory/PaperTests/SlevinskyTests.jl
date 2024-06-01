#=
# Comparison tests for:
#   Slevinsky, R., Fast and stable rational approximation of generalized hypergeometric functions, arXiv:2307.06221
#
# Author: Caleb Jacobs
# DLM: April 16, 2024
=#

import HypergeometricFunctions: pFqweniger, pFqdrummond

function runSlevinskyTests()
    r  = 19.99 
    n  = 101
    Tr = 1.3 
    np = 5
    a = [1.25]
    b = [1.5]
    
    println("Running test with a = ", a, " and b = ", b, ".")
    @time (z, f, h, tru, p) = pFqTest(a, b, r = r, n = n, np = np, Tr = Tr, dir = -1, exclude = false)

    testF = Vector{Vector{ComplexF64}}(undef, 2)
    @time testF[1] = [pFqdrummond(Tuple(a), Tuple(b), z) for z ∈ z]             # Drummond
    @time testF[2] = [pFqweniger(Tuple(a), Tuple(b), z) for z ∈ z]              # Levin-type factorial

    titles = ["Drummond", "Factorial Levin-Type"]

    ps = [getGraphics(z, F, tru, title = title, exclude = false) for (F, title) = zip(testF, titles)]

    for (img, fileName) ∈ zip(ps, titles)
        savefig(img[5], string("Images/Slevinsky/", fileName, ".png"), width = 700, height = 700)
    end
    savefig(p[5], string("Images/Slevinsky/ECTrap.png"), width = 700, height = 700)
    savefig(p[2], string("Images/Slevinsky/ECTrapAbsArg.png"), width = 700, height = 700)
    savefig(p[3], string("Images/Slevinsky/ECTrapRe.png"), width = 700, height = 700)
    savefig(p[4], string("Images/Slevinsky/ECTrapIm.png"), width = 700, height = 700)

    return (p, ps)
end
