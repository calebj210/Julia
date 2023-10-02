#=
# Generalized Gregory quadrature for computing high order hypergeometric pFq over a grid
#
# Author: Caleb Jacobs
# DLM: October 2, 2023
=#

function _0F0(g)
    f = exp.(g.z)

    return f
end

function _1F0(a, g)
    f = (1 .- g.z).^(-a)

    return f
end

function _pFq(g::Grid, a, b, r, n, np)
    if length(a) == 1 && length(b) == 1
    elseif length(a) == 2 && length(b) == 1
    else 
        throw(error("Non-supported pFq base case"))
    end
end

function pFq(a, b; r = 1, n = 5, np = 3)
    g = getGrid(r, n, np = np, nl = length(b))  # Initial grid

    if length(a) == length(b)
        f = _0F0(g)                       # Initial function
    elseif length(a) = length(b) - 1
        f = _1F0(a[end], g)
    else
        throw(error("Non-supported pFq size"))
    end

    for nl = length(b) - 1 : -1 : 1
        g = getGrid(r, n, np = np, nl = nl)

        f = _pFq(g)
    end
    
    return (z, f)
end
