using Plots
using LinearAlgebra
include("Interpolants.jl")
include("Coordinate_Transformations.jl")

#=
    Interpolation routine for sinusoidal coordinates

    Author: Caleb Jacobs
    Date Last Modified: 12/07/2020
=#

## Helper function
f(x) = 2sin(x)

"""
    Interpolate data in exponential coordinates
"""
function interpData()
    # Construct data
    X = range(0,2π, length = 7)
    Y = f.(X)
    data = hcat(X,Y)

    # Compute polynomial interpolant in Cartesian coordinates
    cartCoeffs = findPolCoeff(data)
    display(cartCoeffs)
    t = range(minimum(data[:,1]), maximum(data[:,1]), length = 100)
    s = evalPolInterp(t, cartCoeffs)
    cartInt = [t s]

    # Convert data to sinusoidal coordinate system
    sinData = cartToSin(data)

    # Compute polynoial interpolant in sinusoidal coordinates
    sinCoeffs = findPolCoeff(sinData)
    display(sinCoeffs)

    # Evaluate interpolant over interval
    t = range(minimum(X), maximum(X), length = 100)
    s = evalPolInterp(t, sinCoeffs)
    sinInt = [t s]
    sinInt = sinToCart(sinInt)

    # Plot given data
    plotA = scatter(data[:,1], data[:,2],
                 legend = :topright,
                 label = "Data")

    # Plot solution
    plot!(t, f.(t),
          linewidth = 3,
          linestyle = :dot,
          label = "Solution")

    # Plot Cartesian interpolant
    plot!(cartInt[:,1], cartInt[:,2],
          label = "Cartesian Interpolant")

    # Plot sinusoidal interpolant
    plot!(sinInt[:,1], sinInt[:,2],
          label = "Sinusoidal Interpolant")

    display(plotA)
end

interpData()
