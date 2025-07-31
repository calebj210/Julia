#=
# Compute 2F1 via transformations and a conformal mapping
#
# Author: Caleb Jacobs
# DLM: July 31, 2025
=#

using DSP
using LinearAlgebra
using ComplexVisuals, CairoMakie, LaTeXStrings

include("pFq.jl")

# Generate series weights for z -> 1 - (1 - ζ)^4
function conformalweights(ord)
    a = [0, 4, -6, 4, -1]
    A = zeros(Float64, ord, ord)
    A[1,1] = 1
    for i ∈ 2:ord
        A[:,i] = conv(a,A[:,i - 1])[1:ord]
    end

    return LowerTriangular(A)
end

function conformal_2f1(a, b, c, z; rtol = eps(), mord = 1000, raw = false)
    old = zeros(5)
    old[1] = 1
    for n ∈ 1:4
        old[n+1] = old[n] * (a + n - 1) * (b + n - 1) / n / (c + n - 1)
    end

    α = conformalweights(5) * old

    if !raw
        z = 1 - (1 - z)^.25
    end

    # if abs2(z) >= 1
    #     return NaN + im * NaN
    # end

    zn = z^4
    S = dot(conj(z .^ (0:4)), α)

    A(n) = -10n*(n-1)-4n*(4a+4b+1)-16a*b
    B(n) = 10(n-1)*(n-2)+6(n-1)*(4a+4b+1)+3*16a*b
    C(n) = -5(n-2)*(n-3)-4(n-2)*(4a+4b+1)-3*16a*b
    D(n) = (n-3)*(n-4)+(n-3)*(4a+4b+1)+16a*b
    E(n) = 4(n+1)*(c+n)

    crit_flag = false
    for n ∈ 4:mord
        push!(α, -(A(n) * α[n + 1] + B(n) * α[n] + C(n) * α[n - 1] + D(n) * α[n - 2]) / E(n))

        zn *= z
        val = last(α) * zn
        if abs(val / S) <= rtol
            if crit_flag
                break
            else
                crit_flag = true
            end
        elseif crit_flag
            crit_flag = false
        end
        S += val
    end

    return S
end

function convolutional_maclaurin_2f1(a, b, c, z; max_ord = 25, rtol = 1e-15)
    z = 1 - (1 - z)^.25
    S = zn = coeff = 1.0
    coeffs = [1.0 + 0im]

    # crit_flag = false
    for n ∈ 1:max_ord
        coeff *= (a + n - 1) / (c + n - 1) * (b + n - 1) / n

        # zn *= z
        # val = coeff * zn
        # if abs(val / S) <= rtol
        #     if crit_flag
        #         break
        #     else
        #         crit_flag = true
        #     end
        # elseif crit_flag
        #     crit_flag = false
        # end
        # S += val

        push!(coeffs, coeff)
    end

    return dot(conj(z .^ (0:length(coeffs) - 1)), conformalweights(length(coeffs)), coeffs)
end

function boundarytest()
    test1 = (1.1,1.2,1.3)
    z1 = ComplexGrid(range(-17, 9, 300), range(-13, 13, 300))
    f1 = conformal_2f1.(test1..., z1)
    tru1 = johansson_2f1.(test1..., z1, bits = 106)
    err1 = cleanerror.(f1, tru1)

    test2 = (1,-9/2,-9/4)
    z2 = ComplexGrid(range(-27, 15, 300), range(-21, 21, 300))
    f2 = conformal_2f1.(test2..., z2)
    tru2 = johansson_2f1.(test2..., z2, bits = 106)
    err2 = cleanerror.(f2, tru2)

    fig1 = Figure()
    ax1 = Axis3(fig1[1,1],
       title = L"{_2}F_1(1.1,1.2;1.3;z)",
        titlesize = 20,
        xlabel = L"\mathrm{Re}(z)",
        ylabel = L"\mathrm{Im}(z)",
        zlabel = L"\log_{10}(\text{Error})",
        elevation = .25π,
    )
    fig2 = Figure()
    ax2 = Axis3(fig2[1,1],
        title = L"{_2}F_1(1,-9/2;-9/4;z)",
        titlesize = 20,
        xlabel = L"\mathrm{Re}(z)",
        ylabel = L"\mathrm{Im}(z)",
        zlabel = L"\log_{10}(\text{Error})",
        elevation = .25π,
    )

    p1 = surface!(ax1, reim(z1)..., log10.(err1))
    p2 = surface!(ax2, reim(z2)..., log10.(err2))

    Colorbar(fig1[1,2], p1)
    Colorbar(fig2[1,2], p2)
    resize_to_layout!(fig1)
    resize_to_layout!(fig2)

    return (fig1, fig2)
end

function multiple_sheets()
    ζ = complex_square_grid(sqrt(2) + .1, 100)
    f = conformal_2f1.(1,-9/2,-9/4,ζ; raw = true)

    fig = Figure()
    ax = Axis3(fig[1,1])
    complexsurface!(ax, ζ, sign.(f) .* log10.(abs.(f) .+ 1))

    resize_to_layout!(fig)

    return fig
end

function cleanerror(f, g)
    err = abs(f - g) / (abs(g) + eps())
    if isinf(err) || isnan(err) || err > 1
        err = 1
    elseif err < 1e-16
        err = 1e-16
    end

    return err
end
