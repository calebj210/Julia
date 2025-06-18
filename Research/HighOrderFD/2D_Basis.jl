#=
#   2D Harmonic basis definitions
#
# Author: Caleb Jacobs
# DLM: June 15, 2025
=#

using LinearAlgebra, GenericLinearAlgebra

function gridnodes(N::T; scale = false) where T <: Integer
    n = N ÷ 2
    if scale
        return nodes = collect(zip(repeat(-n:n, inner = N) / n, repeat(n:-1:-n, outer = N) / n))
    else
        return nodes = collect(zip(repeat(-n:n, inner = N), repeat(n:-1:-n, outer = N)))
    end
    
    return nodes
end

function getlinearlyindependentterms(N)
    nodes = gridnodes(N)
    terms = [1]
    term = 1

    while length(terms) < N^2
        term += 1

        A = ([big(float(harmonic(n, big(float(x)), big(float(y))))) for (x, y) ∈ nodes, n ∈ [terms; term]])
        print(term, ": ")
        if rank(A) > length(terms)
            push!(terms, term)
            println(length(terms))
        else
            println()
        end
    end

    print("terms = [1")
    for i ∈ terms[2:end]
        print(',', i)
    end
    println(']')

    return terms
end

function linearlyindependent(N; reduced = true)
    if reduced
        terms = collect(1:(11-4N+N^2))
    elseif N == 3
        terms = [1,2,3,4,5,6,7,8,16]
    elseif N == 5
        terms = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,18,19,20,21,22,23,24,32,40]
    elseif N == 7
        terms = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,34,35,36,37,38,39,40,42,43,44,45,46,47,48,56,64,72]
    elseif N == 9
        terms = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,58,59,60,61,62,63,64,66,67,68,69,70,71,72,74,75,76,77,78,79,80,88,96,104,112]
    elseif N == 11
        terms = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,90,91,92,93,94,95,96,98,99,100,101,102,103,104,106,107,108,109,110,111,112,114,115,116,117,118,119,120,128,136,144,152,160]
    elseif N == 13
        terms = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,130,131,132,133,134,135,136,138,139,140,141,142,143,144,146,147,148,149,150,151,152,154,155,156,157,158,159,160,162,163,164,165,166,167,168,176,184,192,200,208,216]
    else
        terms = getlinearlyindependentterms(N)
    end

    return terms
end

function harmonic(N::T, x, y) where T <: Integer
    r = norm([x, y])
    θ = atan(y, x)

    if N == 1
        return 1.0
    elseif N % 2 == 0
        n = N ÷ 2

        return r^n * cos(n * θ)
    elseif N % 2 == 1
        n = N ÷ 2

        return r^n * sin(n * θ)
    else
        error("N must be a positive integer")
    end
end

function cardinalweights(N::T)::Matrix{Float64} where T<: Integer
    nodes = gridnodes(N)
    terms = linearlyindependent(N)

    A = big.(float([harmonic(n, big(float(x)), big(float(y))) for (x, y) ∈ nodes, n ∈ terms]))
    B = diagm(ones(N^2))

    return A \ B
end

function cardinal(N::T, w, x, y) where T <: Integer
    terms = linearlyindependent(N)

    return sum(w .* harmonic.(terms, x, y))
end
