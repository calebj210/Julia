#=
# Test suite for grid based hypergeometric calculations
#
# Author: Caleb Jacobs
# DLM: October 4, 2023
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
