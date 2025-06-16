#=
# 2D Harmonic Finite Differences
#
# Author: Caleb Jacobs
# DLM: June 16, 2025
=#

include("2D_Basis.jl")
using LinearAlgebra, GenericLinearAlgebra

function getA(N::T, extendedprecision = true; scale = false) where T <: Integer
    nodes = gridnodes(N, scale = scale)
    terms = linearlyindependent(N)

    if extendedprecision
        A = big.(float([harmonic(n, big(float(x)), big(float(y))) for n ∈ terms, (x, y) ∈ nodes]))
    else
        A = [harmonic(n, x, y) for n ∈ terms, (x, y) ∈ nodes]
    end

    return A
end

function Δweights(N::T; prec = 256, matrixform = true) where T <: Integer
    if prec == 64
        A = getA(N, false)[1:end-1,:]
    else
        setprecision(prec)
        A = getA(N, true)[1:end-1,:]
    end
    scaledA = [zeros(1, N^2); A]
    scaledA[1, N^2 ÷ 2 + 1] = 1

    b = zeros(N^2)
    b[1] = 1

    w = scaledA \ b

    if matrixform
        return reshape(w, N, :)
    else
        return w
    end
end

function ∂xweights(N::T, n = 1; prec = 256, matrixform = true, scale = false, offset = 0) where T <: Integer
    if prec == 64
        A = getA(N, false, scale = scale)[1:end - offset,:]
    else
        setprecision(prec)
        A = getA(N, true, scale = scale)[1:end - offset,:]
    end

    b = zeros(N^2 - offset)
    b[2n] = factorial(n)

    w = A \ b

    if matrixform
        return reshape(w, N, :)
    else
        return w
    end
end
