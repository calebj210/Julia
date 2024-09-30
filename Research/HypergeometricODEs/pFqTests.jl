#=
# Tests for hypergeometric ODE approach
# 
# Author: Caleb Jacobs
# DLM: September 25, 2024
=#

using CSV, Tables
using ComplexVisuals
using BenchmarkTools
include("pFq.jl")

"Generate CSVs of grid nodes"
function generategrids(path::String, n::Int64, r)
    g = getGrid(n, r)
    CSV.write("Data/" * path, Tables.table([real(g.z) imag(g.z)]), writeheader = false)

    return nothing
end

"Read CSV to complex vector"
function getcomplexvals(path::String)
    reims = CSV.File(path, header = false, delim = ',', types = Float64) |> Tables.matrix

    vals = vec(reims[:, 1] + im * reims[:, 2])

    return vals
end

function timetest(a, b, z; order = 20, N = 150, H = .1)
    t = @benchmark fast2f1($(a[1]), $(a[2]), $(b[1]), $(z), order = $(order), N = $(N), H = $(H))
    val = fast2f1(a[1], a[2], b[1], z, order = order, N = N, H = H)
    tru = mathematica_2f1(a[1], a[2], b[1], z)
    
    med = round(median(t.times) * 1e-3, digits = 2)
    dig = round(Int, -log10(abs(val - tru)))

    return "Median = $(med) μs, Correct Digits = $(dig)"
end

