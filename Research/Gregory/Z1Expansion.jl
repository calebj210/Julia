#=
# Linear solver for z = 1 (p + 1)Fp expansion weights
#
# Author: Caleb Jacobs
# Date last modified: March 1 , 2024
=#

include("HypergeometricGrid.jl")
include("Visuals.jl")
include("HypergeometricTests.jl")

"""
    getModifiedVand(a, b, z)
Compute the modified Vanermonde matrix for computing pFq expansion weights about z = 1.
"""
function getModifiedVand(a, b, z)
    Aα = [(z - 1)^(i - 1) 
        for z ∈ z, i ∈ 1 : length(z) / 2] # Regular shifted Vandermonde for α coefficients

    Aβ = [oneMinusZα(z, sum(b) - sum(a)) * (z - 1)^(i - 1) 
        for z ∈ z, i ∈ 1 : length(z) / 2] # Modified Vandermonde for β coefficients

    return [Aα Aβ]
end

"""
    getZ1ExpansionWeights(a, b, z)
Compute z = 1 pFq expansion weights.
"""
function getZ1ExpansionWeights(a, b, z, f)
    A = getModifiedVand(a, b, z)

    ω = A \ f

    ωα = ω[1 : round(Int64, length(ω) / 2)]
    ωβ = ω[round(Int64, length(ω) / 2) + 1 : end]

    return (ωα, ωβ)
end

"""
    pFqZ()
"""
function z1Expansion(a, b, ωα, ωβ, z)
#     γ = (1 - z)^(sum(b) - sum(a))
    γ = oneMinusZα(z, sum(b) - sum(a))

    f = sum((ωα + ωβ * γ) .* (z - 1).^(0 : length(ωα) - 1))

    return f
end

"""
    modified2F1(a, b, c, z; n)
"""
function modified2F1(a, b, c; n = 40, r = 1.99, cr = 10, mr = 5, np = 3)
    (z, f, fh) = pFq([a,b], [c]; n = n, r = r, np = np)  # Initial 2F1

    h = abs(z[1] - z[2])
    zm1 = abs.(z .- 1)
    cIdx = findall((cr - .5) * h .< zm1 .&& zm1 .<= cr * h)
    cIdx = length(cIdx) % 2 == 0 ? cIdx : cIdx[1 : end - 1]

    z1 = z[cIdx]
    f1 = f[cIdx]
    (ωα, ωβ) = getZ1ExpansionWeights([a, b], [c], z1, f1)

    mIdx = findall(zm1 .<= mr * h)
    for i ∈ mIdx
        f[i] = z1Expansion([a, b], [c], ωα, ωβ, z[i])
    end

    g = zeros(ComplexF64, length(f))
    g[mIdx] .= im

    return (z, f, fh)
end

"""
    modifiedPFQTest
"""
function modifiedPFQTest(a,b,c; r = 1.99, n = 40, cr = 10, mr = 5, np = 3)
    generateGrids("grid.csv", n, r)

    println("Press enter after running Mathematica to update values!")
    readline()
    
    (z, fm, hm) = modified2F1(a,b,c, n = n, r = r, np = np, cr = cr, mr = mr)
    (z, f, h) = pFq([a,b],[c], n = n, r = r, np = np)

    tru = getComplexVals("Data/pfq.csv")

    p =  complexAbsPlot(z, (f  - tru) ./ abs.(tru), logscale = true, title = "Default")
    pm = complexAbsPlot(z, (fm - tru) ./ abs.(tru), logscale = true, title = "Expansion")

    display(p)
    display(pm)

    return (z, f, tru)
end
