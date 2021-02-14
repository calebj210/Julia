#=
  Routines to compute and evalute interpolants

  Date last modified: 12-07-2020
  Author: Caleb Jacobs
=#

"""
    Compute the coefficients of an interpolating polynomial

    x is the given data set
"""
function findPolCoeff(x)
    # Number of data points
    n = size(x, 1)

    # Setup up polynomial Interpolation matrix
    A = x[:,1] .^ [0:(n - 1)...]'

    # Compute and return the coefficients
    return A \ x[:,2]
end

"""
    Evaluate interpolant given data and coefficients

    x is the coordinates to evaluate the interpolant at
    coeff are the coefficients of the polynomial interpolant
"""
function evalPolInterp(x, coeff)
    tmp = x .^ [0:(length(coeff) - 1)...]'

    tmp .*= coeff'

    y = sum(tmp, dims = 2)

    return y
end
