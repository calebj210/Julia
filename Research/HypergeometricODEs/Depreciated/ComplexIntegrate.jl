#=
#   Integral techniques for complex variables
#
# Author: Caleb Jacobs
# DLM: September 11, 2024
=#

"""
    gauss_mvt(f, z0; r = 0.1, reltol = 1e-14, maxiter = 10)
Compute `f`(`z0`) by taking a contour integral about a circle of radius `r` to a relative tolerance of `reltol`.
"""
function gauss_mvt(f::Function, z0; r = 0.1, reltol = 1e-14, maxiter = 10)
    z(n) = z0 .+ r * cispi.(2(1:n) / n)     # nth roots of unity
    N = 2                                   # Initial size of stencil

    fold = sum(f, z(N)) / N                 # First function value approximation
    fnew = zero(Complex)                    # Initialize function value
    tol = 0
    for n ∈ 1:maxiter
        N *= 2                              # Double size of stencil

        fnew = sum(f, z(N)) / N             # New function value

        tol = abs(1 - fnew/fold)
        if tol <= reltol                    # Check relative stopping criteria
            return fnew
        end

        fold = fnew                         # Old function value
    end

    @warn "Gauss MVT did not converge to desired tolerance $(reltol) in $(maxiter) iterations. Relative tolerance = $(tol)."
    return fnew
end
