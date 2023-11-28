#=
# Test suite for grid based hypergeometric calculations
#
# Author: Caleb Jacobs
# DLM: November 15, 2023
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
function pFqTest(a,b; r = 1, n = 20, np = 3, Tr = 0.5, dir = -1)
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
        eye=attr(x = 1.5dir, y = 1.5dir, z = 1.5)
    )
    relayout!(p2, scene_camera = camera)
    relayout!(p3, scene_camera = camera)
    relayout!(p4, scene_camera = camera)

    return (z, f, h, tru, p1, p2, p3, p4)
end


function runTests(; N = 0)
    test = 
        [
#                 a             b                        r   n np   Tr, cam
            ([1.9,1],        [2.9],                   2.49, 80, 3, 0.6, -1),    # Test 1
            ([1,-1.9],       [2.9],                   2.49, 80, 3, 0.6,  1),    # Test 2
            ([.9,1.1,-2.1],  [1.2,1.3],               2.49, 80, 3, 0.6,  1),    # Test 3
            ([1.1,2.1,.9],   [1.2,1.3],               2.49, 80, 3, 0.6,  1),    # Test 4
            ([.9,1,1.1,1.2], [1.05,1.15,1.25],        2.49, 80, 3, 0.6,  1)     # Test 5
        ]

    if N != 0
        tests = N
    else
        tests = eachindex(a)
    end

    for testN ∈ tests
        (a, b, r, n, np, Tr, dir) = test[testN]

        println("Running test ", testN, " with a = ", a, " and b = ", b, ".")

        (z, f, h, tru, p1, p2, p3, p4) = pFqTest(a, b, r = r, n = n, np = np, Tr = Tr, dir = dir)

        display(p1)
        display(p2)
        display(p3)
        display(p4)

        println("Press enter when ready for next test.")
        readline()
    end

    return nothing
end
