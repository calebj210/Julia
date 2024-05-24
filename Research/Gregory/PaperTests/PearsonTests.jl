#=
# Comparison tests for:
#   Pearson, J.W., Olver, S. & Porter, M.A. Numerical methods for the computation of the confluent and Gauss hypergeometric functions. Numer Algor 74, 821–866 (2017).
#
# Author: Caleb Jacobs
# DLM: May 23, 2024
=#

include("PearsonMethods.jl")
include("../HypergeometricTests.jl")

function runPearsonTests()
    r  = 1.99 
    n  = 40
    Tr = .6
    np = 5
    cr = 9
    sr = 10
    a = [.9, 1.9]
    b = [1.91]
    
    println("Running test with a = ", a, " and b = ", b, ".")
    (z, f, h, tru, p) = pFqTest(a, b, r = r, n = n, np = np, Tr = Tr, cr = cr, sr = sr, dir = 1, exclude = true, modifyZ1 = true)

    testF = Vector{Vector{ComplexF64}}(undef, 5)
    @time testF[1] = [taylorA(a[1], a[2], b[1], z, 1e-15) for z ∈ z]            # Taylor-A test
    @time testF[2] = [taylorB(a[1], a[2], b[1], z, 1e-15) for z ∈ z]            # Taylor-B test
    @time testF[3] = [singleFraction(a[1], a[2], b[1], z, 1e-15) for z ∈ z]     # Single fraction test
    @time testF[4] = [buhring(a[1], a[2], b[1], z, 0.0 + 0im, 1e-15) for z ∈ z] # Buhring test
    @time testF[5] = gjQuad(  a[1], a[2], b[1], z, 150)                         # Gauss-Jacobi quadrature test

    titles = ["Taylor-A", " Taylor-B", "Single Fraction", "Buhring", "Gauss-Jacobi"]

    ps = [getGraphics(z, F, tru, title = title, exclude = true) for (F, title) = zip(testF, titles)]

    for (img, fileName) ∈ zip(ps, titles)
        savefig(img[1], string("Images/Pearson/", fileName, ".png"))
    end

    return (p, ps)
end
