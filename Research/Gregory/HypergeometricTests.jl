#=
# Test suite for grid based hypergeometric calculations
#
# Author: Caleb Jacobs
# DLM: November 9, 2023
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

    if length(a) != length(b) + 1
        (z, f) = pFq(a, b, r = r, n = n, np = np, Tr = Tr)
        h = []
    else    
        (z, f, h) = pFq(a, b, r = r, n = n, np = np, Tr = Tr)
    end

    tru = getComplexVals("Data/pfq.csv")

    title = string(length(a), "F", length(b), "(", a, "; ", b, "; z)") 

    p1 = complexAbsPlot(z, f - tru, logscale = true, title = title)
    p2 = complexPlot3d(z, f)
    p3 = complexPlot3d(z, f, T = 2)
    p4 = complexPlot3d(z, f, T = 3)

    camera = attr(
        eye=attr(x = -1.25, y = -1.25, z = 1.25)
    )
    relayout!(p2, scene_camera = camera)
    relayout!(p3, scene_camera = camera)
    relayout!(p4, scene_camera = camera)

    display([p1 p2; p3 p4])

    return (z, f, h, tru, p1, p2, p3, p4)
end

# pFp Tests
pFp12(; r = 1.5, n = 40, np = 3, Tr = 0.75) = pFqTest([.1], [.2], r = r, n = n, np = np, Tr = Tr) 
pFp1234(; r = 1.5, n = 40, np = 3, Tr = 0.75) = pFqTest([.1,.2], [.3,.4], r = r, n = n, np = np, Tr = Tr) 
pFp123456(; r = 1.5, n = 40, np = 3, Tr = 0.75) = pFqTest([.1,.2,.3], [.4,.5,.6], r = r, n = n, np = np, Tr = Tr) 
pFp12345678(; r = 1.5, n = 40, np = 3, Tr = 0.75) = pFqTest([.1,.2,.3,.4], [.5,.6,.7,.8], r = r, n = n, np = np, Tr = Tr) 
pFp1234567890(; r = 1.5, n = 40, np = 3, Tr = 0.75) = pFqTest([.1,.2,.3,.4,.5], [.6,.7,.8,.9, 1], r = r, n = n, np = np, Tr = Tr) 

# pFp+1 Tests
pFpp122(; r = 1.5, n = 40, np = 3, Tr = 0.5) = pFqTest([1], [2,2], r = r, n = n, np = np, Tr = Tr) 
pFpp11222(; r = 1.5, n = 40, np = 3, Tr = 0.5) = pFqTest([1,1], [2,2,2], r = r, n = n, np = np, Tr = Tr) 
pFpp1112222(; r = 1.5, n = 40, np = 3, Tr = 0.5) = pFqTest([1,1,1], [2,2,2,2], r = r, n = n, np = np, Tr = Tr) 
pFpp111122222(; r = 1.5, n = 40, np = 3, Tr = 0.5) = pFqTest([1,1,1,1], [2,2,2,2,2], r = r, n = n, np = np, Tr = Tr) 
pFpp56(; r = 1.5, n = 40, np = 3, Tr = 0.5) = pFqTest([1,1,1,1,1], [2,2,2,2,2,2], r = r, n = n, np = np, Tr = Tr) 

# ppFp Tests
ppFp112(; r = 1.2, n = 40, np = 3, Tr = 0.5) = pFqTest([1,1], [2], r = r, n = n, np = np, Tr = Tr) 
ppFp11122(; r = 1.2, n = 40, np = 3, Tr = 0.5) = pFqTest([1,1,1], [2,2], r = r, n = n, np = np, Tr = Tr) 
ppFp1111222(; r = 1.2, n = 40, np = 3, Tr = 0.5) = pFqTest([1,1,1,1], [2,2,2], r = r, n = n, np = np, Tr = Tr) 
ppFp111112222(; r = 1.2, n = 40, np = 3, Tr = 0.5) = pFqTest([1,1,1,1,1], [2,2,2,2], r = r, n = n, np = np, Tr = Tr) 
ppFp65(; r = 1.2, n = 40, np = 3, Tr = 0.5) = pFqTest([.1,.2,.3,.4,.5,.6], [.21,.31,.41,.51,.61], r = r, n = n, np = np, Tr = Tr) 
