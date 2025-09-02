#=
# Tests for conformal approach to 2F1.
#
# Author: Caleb Jacobs
# DLM: August 27, 2025
=#

using ComplexVisuals, CairoMakie, LaTeXStrings
include("pFq.jl")
include("Conformal2F1.jl")

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

function conformal_degree_errors(test)
    a, b, c = test[1:3]
    z = complex_square_grid(3,300)

    println("Getting true solution")
    tru = johansson_2f1.(a, b, c, z, bits = 256)
    println("Getting conformal solutions")
    names = 2:6
    pos = [(1,1), (1,2), (1,3), (2,1), (2,2)]
    vals = [f.(a, b, c, z) for f ∈ (conformal_2_2f1, conformal_3_2f1, conformal_2f1, conformal_5_2f1, conformal_6_2f1)]
    errs = [cleanerror.(val, tru) for val ∈ vals]

    set_theme!(complex_theme)
    fig = Figure()
    for i ∈ 1:5
        ax = Axis3(fig[pos[i]...],
                   title = string(names[i]),
            zlabel = L"\log_{10}|\text{Error}|",
            zlabeloffset = 40,
            viewmode = :fit
        )
        surface!(ax, reim(z)..., log10.(errs[i]))
    end
    ax = Axis3(fig[2,3],
        title = "Phase Portrait",
        zlabel = L"pseudolog$_{10}(f)$",
        zlabeloffset = 40,
    )
    complexsurface!(ax, z, sign.(tru) .* log10.(abs.(tru) .+ 1))

    axiscolorwheel(ax, position = :rt)

    resize_to_layout!(fig)

    return fig
end

function best_of_all(test)
    a, b, c = test[1:3]
    z = complex_square_grid(3,300)

    println("Getting true solution")
    tru = johansson_2f1.(a, b, c, z, bits = 256)

    println("Getting transformed solutions")
    mac = [first.(t.(a, b, c, z, maclaurin_2f1)) for t ∈ transformations]
    con = [first.(t.(a, b, c, z, conformal_2f1)) for t ∈ transformations]

    errm = [cleanerror.(val, tru) for val ∈ mac]
    errc = [cleanerror.(val, tru) for val ∈ con]

    err = min.(errm..., errc...)

    # errc = cat([cleanerror.(val, tru) for val ∈ con]..., dims = 3)

    # errcidx = mapslices(argmin, errc, dims = 3)
    # errcidx = mapslices(minimum, errc, dims = 3)
    # errcidx = mapslices(findmin, errc, dims = 3)
    # err = first.(errs)
    # errcidx = last.(errs)

    fig = Figure()
    ax1 = Axis3(fig[1,1],
        title = "Relative Error",
        zlabel = L"\log_{10}|\text{Error}|",
        zlabeloffset = 40,
        viewmode = :fit
    )
    ax2 = Axis3(fig[1,2],
        title = "Phase Portrait",
        zlabel = L"pseudolog$_{10}(f)$",
        zlabeloffset = 40,
    )

    # heatmap!(ax1, z.real, z.imag, errcidx)
    surface!(ax1, reim(z)..., log10.(err))
    complexsurface!(ax2, z, sign.(tru) .* log10.(abs.(tru) .+ 1))

    axiscolorwheel(ax2, position = :rt)

    rowsize!(fig.layout, 1, Aspect(1, 1))

    resize_to_layout!(fig)

    return fig
end
