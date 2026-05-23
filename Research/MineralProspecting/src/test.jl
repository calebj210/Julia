using LinearAlgebra
using GLMakie
using LaTeXStrings

function single_step_grad(h, λ)
    xs = [(x,y) for x in -h:h:h, y in -h:h:h]

    shift = 2
    F((x,y), z) = (z+shift) / (x^2 + y^2 + (z+shift)^2)^(3/2)
    dF((x,y), z) = (x^2 + y^2 - 2(shift + z)^2) / (x^2 + y^2 + (z+shift)^2)^(5/2)

    u₀ = F.(xs, 0)
    du₀ = dF.(xs, 0)

    D = gradient_stencils(h, λ)

    u_grad = sum(D.Au .* u₀ + D.Adu .* du₀)
    du_grad = sum(D.dAu .* u₀ + D.dAdu .* du₀)
    u_tru = F((0,0), λ * h)
    du_tru = dF((0,0), λ * h)

    err = (
        norm(u_grad - u_tru) / norm(u_tru),
        norm(du_grad - du_tru) / norm(du_tru)
    )

    return err
end

function single_step_layr(h, λ)
    xs = [(x,y,z) for x in -h:h:h, y in -h:h:h, z in [0,λ * h]]

    shift = 2
    F((x,y,z)) = (z+shift) / (x^2 + y^2 + (z+shift)^2)^(3/2)

    u₀ = F.(xs)

    D = two_layer_stencils(h, λ)

    u_grad = sum(D .* u₀)
    u_tru = F((0,0,-λ * h))

    err = norm(u_grad - u_tru) / norm(u_tru)

    return err
end

function single_step_convergence_plot()
    hs = 1 ./ 2 .^(1:14)
    λs = 1 ./ 2.0 .^(-1:3)
    colors = [:purple, :blue, :orange, :green, :red]

    set_theme!(theme_latexfonts())
    fig = Figure()
    ax = Axis(
        fig[1,1],
        xreversed = true,
        xscale = log10,
        yscale = log10,
        xlabel = "Step Size (h)",
        ylabel = "Relative Error",
        xlabelsize = 16,
        ylabelsize = 16,
    )
    for (λ, c) in zip(λs, colors)
        grad_errs = single_step_grad.(hs, -λ)
        layr_errs = single_step_layr.(hs, λ)

        lines!(
            ax, hs, first.(grad_errs), 
            color = c, 
            linewidth = 2,
            label = latexstring("Grad \$\\lambda = \$", λ)
        )
        lines!(
            ax, hs, layr_errs, 
            color = c, 
            linestyle = :dash, 
            linewidth = 2,
            label = latexstring("Layer \$\\lambda = \$", λ)
        )
    end

    lines!(
        ax, hs, 1e2 * hs.^4,
        label = L"\mathcal{O}(h^4)",
        linestyle = :dot,
        linewidth = 3,
    )

    axislegend(ax, labelsize = 12)

    resize_to_layout!(fig)

    return fig
end

function multi_step_grad(h, λ, N, F, dF)
    xs = [(x,y) for x in -100h:h:100h, y in -100h:h:100h]
    # xs = [(x,y) for x in -50:h:50, y in -50:h:50]

    u₀ = F.(xs, 0)
    du₀ = dF.(xs, 0)
    z₀ = 0.0

    D = gradient_stencils(h, λ)

    for _ in 1:N
        utmp = apply_operator(D.Au, u₀) + apply_operator(D.Adu, du₀)
        dutmp = apply_operator(D.dAu, u₀) + apply_operator(D.dAdu, du₀)

        u₀ = utmp
        du₀ = dutmp

        z₀ += λ * h
    end

    u_tru = F.(xs, z₀)
    du_tru = dF.(xs, z₀)

    idx = norm.(xs) .<= 50h
    # idx = norm.(xs) .<= 25

    u_err = norm(u₀[idx] - u_tru[idx]) / norm(u_tru[idx])
    du_err = norm(du₀[idx] - du_tru[idx]) / norm(du_tru[idx])

    return (u_err, du_err)
end

function multi_step_layr(h, λ, N, F)
    xs = [(x,y) for x in -100h:h:100h, y in -100h:h:100h]
    # xs = [(x,y) for x in -50:h:50, y in -50:h:50]

    u₀ = cat(F.(xs, 0), F.(xs, λ * h), dims = 3)
    z₀ = 0.0

    D = two_layer_stencils(h, λ)

    for _ in 1:N
        utmp = apply_operator(D, u₀)

        u₀[:,:,2] = u₀[:,:,1]
        u₀[:,:,1] = utmp

        z₀ -= λ * h
    end

    u_tru = F.(xs, z₀)

    idx = norm.(xs) .<= 50h
    # idx = norm.(xs) .<= 25

    u_err = norm(u₀[idx,1] - u_tru[idx]) / norm(u_tru[idx])

    return u_err
end

function multi_step_convergence_plot()
    hs = 1 ./ 2 .^(1:5)
    λs = 1 ./ 2.0 .^(0:3)
    colors = [:blue, :orange, :green, :red]

    shift = 2
    F((x,y), z) = (z+shift) / (x^2 + y^2 + (z+shift)^2)^(3/2)
    dF((x,y), z) = (x^2 + y^2 - 2(shift + z)^2) / (x^2 + y^2 + (z+shift)^2)^(5/2)

    set_theme!(theme_latexfonts())
    fig = Figure()
    ax = Axis(
        fig[1,1],
        xreversed = true,
        xscale = log2,
        yscale = log10,
        xlabel = "Step Size (h)",
        ylabel = "Relative Error",
        xlabelsize = 16,
        ylabelsize = 16,
        limits = (nothing, (1e-5, 1e-1)),
    )

    for (λ, c) in zip(λs, colors)
        Ns = round.(Int64, 0.5 ./ (λ * hs))

        grad_errs = multi_step_grad.(hs, -λ, Ns, F, dF)
        layr_errs = multi_step_layr.(hs, λ, Ns, F)

        lines!(
            ax, hs, first.(grad_errs), 
            color = c, 
            linewidth = 2,
            label = latexstring("Grad \$\\lambda = \$", λ)
        )

        lines!(
            ax, hs, layr_errs, 
            color = c, 
            linestyle = :dash, 
            linewidth = 2,
            label = latexstring("Layer \$\\lambda = \$", λ)
        )
    end

    lines!(
        ax, hs, 2e-2 * hs.^2,
        label = L"\mathcal{O}(h^2)",
        linestyle = :dot,
        linewidth = 3,
    )

    axislegend(ax, labelsize = 12, position = :lb)

    resize_to_layout!(fig)

    return fig
end
