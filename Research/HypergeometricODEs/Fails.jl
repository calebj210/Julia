#=
# Tests for inspected lower accuracy random tests
#
# Author: Caleb Jacobs
# DLM: July 1, 2025
=#

using CairoMakie, LaTeXStrings
using Random
include("pFq.jl")

function random_test(;N = 10000, arng = 25, brng = 25, crng = 25, zrng = 2, seed = 997, complextest = false)
    Random.seed!(seed)

    # Setup random tests
    if complextest
        as = arng * complexrand(N)
        bs = brng * complexrand(N)
        cs = crng * complexrand(N)
        tests = Vector{NTuple{4, ComplexF64}}()
    else
        as = arng * (1 .- 2rand(N))
        bs = brng * (1 .- 2rand(N))
        cs = crng * (1 .- 2rand(N))
        tests = Vector{Tuple{Float64, Float64, Float64, ComplexF64}}()
    end
    zs = zrng * complexrand(N)

    print("Getting tests: ")
    tru = Vector{ComplexF64}()
    for (a,b,c,z) ∈ zip(as, bs, cs, zs)
        val = johansson_2f1(a,b,c,z; bits = 512)
        if isnan(val) || isinf(val)
            continue
        end
        push!(tests, (a,b,c,z))
        push!(tru, convert(ComplexF64, val))
    end
    println("count: $(length(tests))")

    # Generate graphics and results
    set_theme!(theme_latexfonts())
    fig = Figure()
    err = clean_error.([taylor_2f1(test...) for test ∈ tests], tru)

    ax = Axis(fig[1,1], 
        limits = (nothing, (0,1)), 
        xscale = log10, 
        xlabel = "Relative Error", 
    )

    bin = 10.0 .^ (-16:2:2)

    hist!(ax, err, 
        bins = bin, 
        color = :values, 
        normalization = :probability
    )

    sorted_tests = sort_errors(tests, err, 10.0 .^ (-15:1:0))

    return fig, sorted_tests
end

function sort_errors(tests, errs, bins)
    sorted_tests = ntuple(x -> Vector{Tuple{Float64, Float64, Float64, ComplexF64}}(), length(bins))

    for (err, test) ∈ zip(errs, tests)
        for (i, bin) ∈ enumerate(bins)
            if err <= bin
                push!(sorted_tests[end - i + 1], test)
                break
            end
        end
    end

    # return sorted_tests
    return Dict(zip(bins, sorted_tests))
end

function clean_error(f,t)
    err = abs.((f - t) ./ t)

    if err > 1 || isnan(err) || isinf(err)
        err = 1
    elseif err <= 1e-16
        err = 1e-16
    end

    return err
end

function complexrand(N)
    vals = -(1 + im) .+ 2rand(ComplexF64, N)

    return vals
end
