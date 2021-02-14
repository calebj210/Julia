using Plots
using LinearAlgebra
include("Interpolants.jl")
include("Coordinate_Transformations.jl")

#=
    Interpolation routine for exponential coordinates

    Author: Caleb Jacobs
    Date Last Modified: 12/07/2020
=#

"""
    Interpolate data in exponential coordinates
"""
function interpData()
    # Construct data
    X = range(0,π, length = 3)
    Y = sinc.(X)
    data = hcat(X,Y)

    # Compute polynomial interpolant in Cartesian coordinates
    cartCoeffs = findPolCoeff(data)
    display(cartCoeffs)
    t = range(minimum(data[:,1]), maximum(data[:,1]), length = 100)
    s = evalPolInterp(t, cartCoeffs)
    cartInt = [t s]

    # Convert data to exponential coordinate system
    expData = cartToExp(data)

    # Compute polynoial interpolant in exponential coordinates
    expCoeffs = findPolCoeff(expData)
    display(expCoeffs)

    # Evaluate interpolant over interval
    t = range(minimum(X), maximum(X), length = 100)
    s = evalPolInterp(t, expCoeffs)
    expInt = [t s]
    expInt = expToCart(expInt)

    # Plot given data
    # plotA = scatter(data[:,1], data[:,2],
    #              legend = :topright,
    #              label = "Data")

    # Plot Cartesian interpolant
    plotA = plot(cartInt[:,1], cartInt[:,2] - sinc.(t),
          label = "Cartesian Interpolant")

    # Plot exponential interpolant
    plot!(expInt[:,1], expInt[:,2] - sinc.(t),
          label = "Exponential Interpolant")

    display(plotA)
end

interpData()
