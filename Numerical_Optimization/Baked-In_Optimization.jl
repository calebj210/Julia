using Optim

#=
# Optimizer for the Rosenbach Function Using the Optim.jl Package
#
# Date Last Modified: 01/16/2021
# Author: Caleb Jacobs
=#

# Rosenbach function
f(x) = (1 - x[1])^2 + 100(x[2] - x[1]^2)^2

# Optimization function
function optBach(x₀, method)
    # Minimize the Rosenbach function
    res = optimize(f, x₀, method)

    # Display algorithm details
    display(res)

    # Get minimizer
    minimizer = Optim.minimizer(res)

    # Display minimizer
    display(minimizer)

    return minimizer
end

optBach([0.0,0.0], LBFGS())
