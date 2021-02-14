using LinearAlgebra
using Plots

#=
# Monte Carlo Simulation for π
# Date Last Modified: January 7, 2021
# Author: Caleb Jacobs
=#

"""
    getPoints(n)

Generate `n` uniformly random points over the unit range [-1,1)×[-1,1)
"""
function getPoints(n)
    data = 2*rand(Float64, (2, n)) .- 1

    return data
end

"""
    innerCount(`points``)

Count the number of data `points` that are in the unit range [-1,1)×[-1,1)
"""
function innerCount(data)
    count = 0                # inner point counter

    # Count the number of points
    for i ∈ 1:size(data,2)
        if (norm(data[:,i]) < 1)
            count += 1
        end
    end

    return count
end

"""
    computeπ(n)

Compute π using the Monte Carlo method using `n` data points
"""
function computeπ(n)
    points = getPoints(n)      # Generate random data
    count = innerCount(points) # number of inner circle points

    πApprox = 4 * count / n

    return πApprox
end
