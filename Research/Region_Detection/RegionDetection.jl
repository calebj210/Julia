#=
# Approximate normal based region detection algorithm
#
# Author: Caleb Jacobs
# DLM: 27-06-2022
=#

# Needed packages
include("Normals.jl")
using LinearAlgebra
using NearestNeighbors

"""
    classifyNodes(Γ, Ω, n)

Classify each node in `Ω` as interior or exterior to the hypersurface nodes `Γ`.
Use `n` nearest neighbors in normal calculation.
"""
function classifyNodes(Γ, Ω, n)
    d, N = size(Ω)                      # Dimension and number of nodes

    kdtree = KDTree(Γ)                  # Get KD-Tree over Γ
    ΓIdx = knn(kdtree, Γ, n, true)[1]   # Find nearest neighbors over Γ
    ΩIdx = nn(kdtree, Ω)[1]             # Find nearest neighbors to Ω on Γ

    n⃗ = getn⃗(Γ, ΓIdx)                   # Get outward oriented normals to surface

    io = zeros(N)                       # Inside-outside vector
    for i ∈ 1:N
        v⃗ = Ω[:, i] - Γ[:, ΩIdx[i]]     # Vector from Γ to Ω

        io[i] = -sign(v⃗ ⋅ n⃗[:, ΩIdx[i]]) # Determine if node is inside or outside
    end
    
    return io
end
