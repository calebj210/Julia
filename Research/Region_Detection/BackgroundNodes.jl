#=
# Background node sets
#
# DLM: 27-06-2022
=#

function circ(N::Int)
    t = range(0, 2π * (1 - 1 / N), length = N)

    return hcat((s -> [cos(s), sin(s)]).(t)...)
end

function cartGrid(N, M, xlims, ylims)
    x = range(xlims[1], xlims[2], length = N)
    y = range(ylims[1], ylims[2], length = M)

    return [repeat(x, inner = M)'; 
            repeat(y, outer = N)']
end
