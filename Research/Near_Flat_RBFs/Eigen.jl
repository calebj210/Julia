#=
# Plot consecutive eigenvalues for various surface curvatures
#
# Author: Caleb Jacobs
# DLM: 22-04-2022
=#

include("Conditioning.jl")
using Printf

function eigPlot(ϕ, ε, N)
    κ = 10 .^ range(-2.5, 0, length = 100) 
    n = length(κ)

    X = getNodes(N, 0)
    NNew = size(X, 2)
    @printf "Number of nodes used: %d\n" NNew
    
    λ = zeros(NNew, n)
    for i ∈ 1:n
        X = getNodes(N, κ[i])
        A = colloc(X, ϕ, ε)
        λ[:, i] = eigvals(A)

        if i == 1
            display(λ[:,1])
        end
    end

    κ = repeat(κ', NNew)
    scatter(κ[:], abs.(λ[:]),
            # mz  = log10.(κ[:]),
            mz  = repeat(1:NNew, n, 1),
            ms  = 2,
            msα = 0,
            xscale = :log10,
            yscale = :log10,
            legend = false,
            colorbar = false,
            aspect_ratio = :auto)
end
