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
function pFqTest(a,b; r = 1, n = 20, np = 3, Tr = 0.5, cr = 9, sr = 10, modifyZ1 = true, dir = -1, exclude = false)
    generateGrids("grid.csv", n, r)

    println("Press enter after running Mathematica to update values!")
    readline()

    if length(a) != length(b) + 1
        (z, f) = pFq(a, b, r = r, n = n, np = np, Tr = Tr)
        h = []
    else    
        (z, f, h) = pFq(a, b, r = r, n = n, np = np, Tr = Tr, cr = cr, sr = sr, modifyZ1 = modifyZ1)
    end

    tru = getComplexVals("Data/pfq.csv")

    title = string(length(a), "F", length(b), "(", a, "; ", b, "; z)") 

    p1 = complexAbsPlot(z, f - tru, logscale = true, title = title)
    p2 = complexPlot3d(z, f, exclude = exclude)
    p3 = complexPlot3d(z, f, T = 2, exclude = exclude, mesh = true)
    p4 = complexPlot3d(z, f, T = 3, exclude = exclude, mesh = true)
    p5 = complexAbsPlot(z, (f - tru) ./ abs.(tru), logscale = true, title = title)

    camera = attr(
        eye=attr(x = 1.5dir, y = 1.5dir, z = 1.5)
    )
    relayout!(p2, scene_camera = camera, template = "plotly_white")
    relayout!(p3, scene_camera = camera, template = "plotly_white")
    relayout!(p4, scene_camera = camera, template = "plotly_white")

    return (z, f, h, tru, p1, p2, p3, p4, p5)
end


function runTests(; N = 0)
    test = 
        [
#                 a                  b                    r   n np   Tr, cam  exclusion
            ([1.9,1],            [2.9],                2.49, 80, 3, 0.6, -1,  true),    # Test 1
            ([1,-1.9],           [2.9],                2.49, 80, 3, 0.6,  1,  true),    # Test 2
            ([1.0,1.1,0.9],      [1.2,1.3],            1.99, 80, 3, 0.6, -1,  true),    # Test 3
            ([1.0,1.1,-0.9],     [1.2,1.3],            1.99, 80, 3, 0.6,  1,  true),    # Test 4
            ([1.0,1.1,1.2,0.9],  [1.3,1.4,1.5],        1.99, 80, 3, 0.6, -1,  true),    # Test 5
            ([1.0,1.1,1.2,-0.9], [1.3,1.4,1.5],        1.99, 80, 3, 0.6,  1,  true),    # Test 6
        ]

    if N != 0
        tests = N
    else
        tests = eachindex(test)
    end

    for testN ∈ tests
        (a, b, r, n, np, Tr, dir, ex) = test[testN]

        println("Running test ", testN, " with a = ", a, " and b = ", b, ".")

        (z, f, h, tru, p1, p2, p3, p4) = pFqTest(a, b, r = r, n = n, np = np, Tr = Tr, dir = dir, exclude = ex)

        display(p1)
        display(p2)
        display(p3)
        display(p4)

        println("Press enter when ready for next test.")
        readline()
    end

    return nothing
end



#             ([.9,1.1,-2.1],  [1.2,1.3],               2.00, 60, 3, 0.6,  1),    # Test 3
#             ([1.1,2.1,.9],   [1.2,1.3],               2.49, 80, 3, 0.6,  1),    # Test 4
