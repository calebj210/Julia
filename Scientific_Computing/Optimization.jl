#=
  Steepest decent optimization for Group Project #11

  Author: Caleb Jacobs
  Date Last Modified: 12/08/2020
=#

"""
    armijoSearch(f, df, x, p; α₀ = 5, τ = 0.5, β = 0.5, maxIts = 50)

Armijo Backtracking line search algorithm.

# Arguments
- `f::Function`: the function to minimize
- `df::Function`: derivative of function to minimize
- `x::Array{Float64,1}`: Current position
- `p::Array{Float64,1}`: Search direction
- `α₀::Float64`: Initial search parameter
- `τ::Float64`: Armijo step size scaling factor
- `β::Float64`: Stopping criteria
- `maxIts`: Maximum number of iterations to allow for convergence
"""
function armijoSearch(f::Function, df::Function, x::Array{Float64,1}, p::Array{Float64,1};
                      α₀ = 5, τ = 0.5, β = 0.5, maxIts = 50)
    α = α₀
    its = 0
    while f(x + α*p) > (f(x) + β*α*p ⋅ df(x))
        α = its + 1
        its += 1
        if its > maxIts
            println("MAXIMUM ITERATIONS USED: LINE SEARCH DID NOT CONVERGE")
            return
        end
    end

    return (α, its)
end

"""
    minimize()

Steepest decent minimization algorithm using an Armijo line search.
"""
