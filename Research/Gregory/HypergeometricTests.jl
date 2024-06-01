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

    return nothing
end

"Read CSV to complex vector"
function getComplexVals(path::String)
    reims = CSV.File(path, header = false, delim = ',', types = Float64) |> Tables.matrix

    vals = vec(reims[:, 1] + im * reims[:, 2])

    return vals
end

"Generate graphics"
function getGraphics(z, f, tru; title = "", dir = -1, exclude = true)
    p1 = complexAbsPlot(z, f - tru, logscale = true, title = title)
    p2 = complexPlot3d(z, f, exclude = exclude)
    p3 = complexPlot3d(z, f, T = 2, exclude = exclude, mesh = true)
    p4 = complexPlot3d(z, f, T = 3, exclude = exclude, mesh = true)
    p5 = complexAbsPlot(z, (f - tru) ./ abs.(tru), logscale = true, title = title)

    camera = attr(
        eye=attr(x = 1.5dir, y = 1.5dir, z = 1.5)
    )
    relayout!(p1, scene_camera = camera, template = "plotly_white")
    relayout!(p2, scene_camera = camera, template = "plotly_white")
    relayout!(p3, scene_camera = camera, template = "plotly_white")
    relayout!(p4, scene_camera = camera, template = "plotly_white")

    return [p1, p2, p3, p4, p5]
end

"Generate pFq test"
function pFqTest(a,b; r = 1, n = 20, np = 3, Tr = 0.5, interpN = 10, circR = .8, circN = 150, corrR = .25, branchN = 150, modifyZ1 = true, dir = -1, exclude = false)
    generateGrids("grid.csv", n, r)

    println("Press enter after running Mathematica to update values!")
    readline()

    if length(a) != length(b) + 1
        (z, f) = pFq(a, b, r = r, n = n, np = np, Tr = Tr)
        h = []
    else    
        (z, f, h) = pFq(a, b, r = r, n = n, np = np, Tr = Tr, interpN = interpN, circR = circR, circN = circN, corrR = corrR, branchN = branchN, modifyZ1 = modifyZ1)
    end

    tru = getComplexVals("Data/pfq.csv")

    title = string(length(a), "F", length(b), "(", a, "; ", b, "; z)") 

    p = getGraphics(z, f, tru, title = title, dir = dir, exclude = exclude)

    return (z, f, h, tru, p)
end


function runTests(; N = 0)
    test = 
        [
#                 a                  b                    r   n np   Tr circR circN corrR interpN branchN cam  exclusion
            ([1.1,1.9],          [2.9],                1.99, 41, 5, 0.5,   .8,  200,   .5,    10,     150, -1,  true),    # Test 1
            ([1.1,-1.9],         [2.9],                1.99, 41, 5, 0.5,   .8,  200,   .5,    10,     150,  1,  true),    # Test 2
            ([1.0,1.1,0.9],      [1.2,1.3],            1.99, 41, 5, 0.6,   .8,  200,   .5,    10,     150, -1,  true),    # Test 3
            ([1.0,1.1,-0.9],     [1.2,1.3],            1.99, 41, 5, 0.6,   .8,  200,   .5,    10,     150,  1,  true),    # Test 4
            ([1.0,1.1,1.2,0.9],  [1.3,1.4,1.5],        1.99, 41, 5, 0.6,   .8,  200,   .5,    10,     150, -1,  true),    # Test 5
            ([1.0,1.1,1.2,-0.9], [1.3,1.4,1.5],        1.99, 41, 5, 0.6,   .8,  200,   .5,    10,     150,  1,  true),    # Test 6
        ]

    if N != 0
        tests = N
    else
        tests = eachindex(test)
    end

    for testN ∈ tests
        (a, b, r, n, np, Tr, circR, circN, corrR, interpN, branchN, dir, ex) = test[testN]

        println("Running test ", testN, " with a = ", a, " and b = ", b, ".")

        (z, f, h, tru, p) = pFqTest(a, b, r = r, n = n, np = np, Tr = Tr, circR = circR, circN = circN, corrR = corrR, interpN = interpN, branchN = branchN, dir = dir, exclude = ex)

        display(p[1])
        display(p[2])
        display(p[3])
        display(p[4])

        println("Press enter when ready for next test.")
        readline()
    end

    return nothing
end
