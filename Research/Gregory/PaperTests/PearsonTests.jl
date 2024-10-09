#=
# Comparison tests for:
#   Pearson, J.W., Olver, S. & Porter, M.A. Numerical methods for the computation of the confluent and Gauss hypergeometric functions. Numer Algor 74, 821–866 (2017).
#
# Author: Caleb Jacobs
# DLM: October 9, 2024
=#

include("PearsonMethods.jl")
include("../HypergeometricTests.jl")

function runPearsonTests()
    r  = 1.99; n  = 41; Tr = .6; np = 5
    corrR = .5; innerR = .6; outerR = .8
    z1N = 70; modifyZ1 = true
    a = [.9, 1.9]; b = [1.91]

    test = (a, b, r, n, np, Tr, modifyZ1, corrR, innerR, outerR, z1N)
    
    println("Running test with a = ", a, " and b = ", b, ".")
    @time (z, f, h, tru, p) = gridtest(test...)

    testF = Vector{Vector{ComplexF64}}(undef, 5)
    @time testF[1] = taylorA.(a..., b..., z, 1e-15)            # Taylor-A test
    @time testF[2] = taylorB.(a..., b..., z, 1e-15)            # Taylor-B test
    @time testF[3] = singleFraction.(a..., b..., z, 1e-15)     # Single fraction test
    @time testF[4] = buhring.(a..., b..., z, 0.0 + 0im, 1e-15) # Buhring test
    @time testF[5] = gjQuad(a..., b..., z, 150)                # Gauss-Jacobi quadrature test

    titles = ["Taylor-A", "Taylor-B", "Single Fraction", "Buhring", "Gauss-Jacobi"]

    ps = [generate_graphics(z, F, tru, title = title) for (F, title) = zip(testF, titles)]

    for (img, fileName) ∈ zip(ps, titles)
        savefig(img[1], string("Images/Pearson/", fileName, ".png"), width = 700, height = 700)
    end
    savefig(p[1], string("Images/Pearson/ECTrap.png"), width = 700, height = 700)
    savefig(p[3], string("Images/Pearson/ECTrapAbsArg.png"), width = 700, height = 700)
    savefig(p[4], string("Images/Pearson/ECTrapRe.png"), width = 700, height = 700)
    savefig(p[5], string("Images/Pearson/ECTrapIm.png"), width = 700, height = 700)

    return (p, ps)
end
