#=
# Test suite for grid based hypergeometric calculations
#
# Author: Caleb Jacobs
# DLM: October 5, 2023
=#

using CSV, Tables
include("Visuals.jl")
include("HypergeometricGrid.jl")

"Generate CSVs of grid nodes"
function generateGrids(path::String, n::Int64, r)
    g = getGrid(n, r)
    CSV.write("Data/" * path, Tables.table([real(g.z) imag(g.z)]), writeheader = false)
end

"Read CSV to complex vector"
function getComplexVals(path::String)
    reims = CSV.File(path, header = false, delim = ',', types = Float64) |> Tables.matrix

    vals = vec(reims[:, 1] + im * reims[:, 2])

    return vals
end

"Generate pFq test"
function pFqTest(a,b; r = 1, n = 20, np = 3, Tr = 0.5)
    generateGrids("grid.csv", n, r)

    println("Press enter after running Mathematica to update values!")
    readline()

    (z, f) = pFq(a, b, r = r, n = n, np = np, Tr = Tr)
    tru = getComplexVals("Data/pfq.csv")

    title = string("(p + 1)Fp (", a, "; ", b, "; z)") 

    display(complexAbsPlot(z, f - tru, logscale = true, title = title))

    return (z, f, tru)
end

# pFp Tests
pFp12(; r = 1.5, n = 40, np = 3, Tr = 0.4) = pFqTest([.1], [.2], r = r, n = n, np = np, Tr = Tr) 
pFp1122(; r = 1.5, n = 40, np = 3, Tr = 0.4) = pFqTest([.1,.1], [.2,.2], r = r, n = n, np = np, Tr = Tr) 
pFp111222(; r = 1.5, n = 40, np = 3, Tr = 0.4) = pFqTest([.1,.1,.1], [.2,.2,.2], r = r, n = n, np = np, Tr = Tr) 
pFp11112222(; r = 1.5, n = 40, np = 3, Tr = 0.4) = pFqTest([.1,.1,.1,.1], [.2,.2,.2,.2], r = r, n = n, np = np, Tr = Tr) 
pFp1111122222(; r = 1.5, n = 40, np = 3, Tr = 0.4) = pFqTest([.1,.1,.1,.1,.1], [.2,.2,.2,.2,.2], r = r, n = n, np = np, Tr = Tr) 

# pFp+1 Tests
pFpp122(; r = 1.5, n = 40, np = 3, Tr = 0.3) = pFqTest([1], [2,2], r = r, n = n, np = np, Tr = Tr) 
pFpp11222(; r = 1.5, n = 40, np = 3, Tr = 0.3) = pFqTest([1,1], [2,2,2], r = r, n = n, np = np, Tr = Tr) 
pFpp1112222(; r = 1.5, n = 40, np = 3, Tr = 0.3) = pFqTest([1,1,1], [2,2,2,2], r = r, n = n, np = np, Tr = Tr) 
pFpp111122222(; r = 1.5, n = 40, np = 3, Tr = 0.3) = pFqTest([1,1,1,1], [2,2,2,2,2], r = r, n = n, np = np, Tr = Tr) 
pFpp11111222222(; r = 1.5, n = 40, np = 3, Tr = 0.3) = pFqTest([1,1,1,1,1], [2,2,2,2,2,2], r = r, n = n, np = np, Tr = Tr) 

# ppFp Tests
ppFp112(; r = 1.2, n = 40, np = 3, Tr = 0.25) = pFqTest([1,1], [2], r = r, n = n, np = np, Tr = Tr) 
ppFp11122(; r = 1.2, n = 40, np = 3, Tr = 0.25) = pFqTest([1,1,1], [2,2], r = r, n = n, np = np, Tr = Tr) 
ppFp1111222(; r = 1.2, n = 40, np = 3, Tr = 0.25) = pFqTest([1,1,1,1], [2,2,2], r = r, n = n, np = np, Tr = Tr) 
ppFp111112222(; r = 1.2, n = 40, np = 3, Tr = 0.25) = pFqTest([1,1,1,1,1], [2,2,2,2], r = r, n = n, np = np, Tr = Tr) 
ppFp11111122222(; r = 1.2, n = 40, np = 3, Tr = 0.25) = pFqTest([1,1,1,1,1,1], [2,2,2,2,2], r = r, n = n, np = np, Tr = Tr) 
