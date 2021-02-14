using Plots
using LinearAlgebra
include("Interpolants.jl")
include("Coordinate_Transformations.jl")

#=
    Interpolation routine for polar coordinates

    Author: Caleb Jacobs
    Date Last Modified: 12/07/2020
=#

"""
    Interpolate data in polar coordinates
"""
function interpData()
    # Create semicircle
    X = range(-1, 1, length = 5)
    Y = sqrt.(1 .- X.^2)

    # Create cartesian data set
    data = [X Y]

    # Find polynomial interpolant under Cartesian coordinates
    cartCoeffs = findPolCoeff(data)
    display(cartCoeffs)
    t = range(minimum(data[:,1]), maximum(data[:,1]), length = 100)
    s = evalPolInterp(t, cartCoeffs)
    cartInt = [t s]

    # Find polynomial interpolant under polar coordinates
    polData = cartToPolar(data)
    coeffs = findPolCoeff(polData)
    display(coeffs)
    t = range(minimum(polData[:,1]), maximum(polData[:,1]), length = 100)
    s = evalPolInterp(t, coeffs)
    polInt = [t s]
    polInt = polarToCart(polInt)

    # Plot given data
    plotA = plot(data[:,1], data[:,2],
                 legend = :topright,
                 aspect_ratio = 1,
                 label = "Data")

    # Plot Cartesian interpolant
    plot!(cartInt[:,1], cartInt[:,2],
          label = "Cartesian Interpolant")

    # Plot polar interpolant
    plot!(polInt[:,1], polInt[:,2],
          label = "Polar Interpolant")

    display(plotA)
end

interpData()
