using LinearAlgebra
using GLMakie

function continuation_compare(h, N, zf)
    x₀ = [(x, y, 0) for x in -1:h:1, y in -1:h:1][1:end-1,1:end-1]
    x₁ = [(x, y, h) for x in -1:h:1, y in -1:h:1][1:end-1,1:end-1]
    xf = [(x, y, -zf) for x in -1:h:1, y in -1:h:1][1:end-1,1:end-1]

    shift = 2
    F((x,y,z)) = (z+shift) / (x^2 + y^2 + (z+shift)^2)^(3/2)
    dF((x,y,z)) = (x^2 + y^2 - 2(shift + z)^2) / (x^2 + y^2 + (z+shift)^2)^(5/2)

    u₀ = F.(x₀)
    u₁ = F.(x₁)
    du₀ = dF.(x₀)

    println(size(x₀))
    println(x₀[end ÷ 2 + 1])

    u_grad = gradient_continue(; h, u₀, du₀, zf, N)
    # u_layr = two_layer_continue(h, Δz, cat(u₀, u₁; dims = 3), -zf)
    u_tru = F.(xf)

    # err_grad = norm(u_grad - u_tru) / norm(u_tru)
    err_grad = abs(u_grad[end ÷ 2 + 1] - u_tru[end ÷ 2 + 1]) / abs(u_tru[end ÷ 2 + 1])
    # err_layr = norm(u_layr - u_tru) / norm(u_tru)

    println("Grad error: ", err_grad)
    # println("Layr error: ", err_layr)
    
    fig = Figure()

    axl = Axis3(fig[1,1])
    axr = Axis3(fig[1,2])
    axc = Axis3(fig[2,1:2])

    surface!(axl, first.(x₀), getindex.(x₀,2), u_tru)
    surface!(axr, first.(x₀), getindex.(x₀,2), u_grad)
    surface!(axc, first.(x₀), getindex.(x₀,2), log10.(abs.(u_tru - u_grad)))

    return fig
end

function single_step_grad(h, λ)
    xs = [(x,y) for x in -100h:h:99h, y in -100h:h:99h]
    # xs = [(x,y) for x in -2h:h:2h, y in -2h:h:2h]
    # xs = [(x,y) for x in -h:h:h, y in -h:h:h]

    shift = 2
    F((x,y), z) = (z+shift) / (x^2 + y^2 + (z+shift)^2)^(3/2)
    dF((x,y), z) = (x^2 + y^2 - 2(shift + z)^2) / (x^2 + y^2 + (z+shift)^2)^(5/2)

    u₀ = F.(xs, 0)
    du₀ = dF.(xs, 0)

    D = gradient_stencils(h, λ)

    # u_grad = sum(D.Au .* u₀ + D.Adu .* du₀)
    # du_grad = sum(D.dAu .* u₀ + D.dAdu .* du₀)
    # u_tru = F((0,0), λ * h)
    # du_tru = dF((0,0), λ * h)
    u_grad = (apply_operator(D.Au, u₀) + apply_operator(D.Adu, du₀))[3:end-2,3:end-2]
    du_grad = (apply_operator(D.dAu, u₀) + apply_operator(D.dAdu, du₀))[3:end-2,3:end-2]
    u_tru = F.(xs, λ * h)[3:end-2,3:end-2]
    du_tru = dF.(xs, λ * h)[3:end-2,3:end-2]

    err = (
        norm(u_grad - u_tru) / norm(u_tru),
        norm(du_grad - du_tru) / norm(du_tru)
    )
    println(err)

    return err
end

function single_step_two(h, λ)
    # xs = [(x,y) for x in -100h:h:99h, y in -100h:h:99h]
    # xs = [(x,y) for x in -2h:h:2h, y in -2h:h:2h]
    xs = [(x,y,z) for x in -h:h:h, y in -h:h:h, z in [0,λ * h]]

    shift = 2
    F((x,y,z)) = (z+shift) / (x^2 + y^2 + (z+shift)^2)^(3/2)

    u₀ = F.(xs)

    D = two_layer_stencils(h, λ)

    u_grad = sum(D .* u₀)
    u_tru = F((0,0,-λ * h))
    # u_grad = (apply_operator(A, u₀))[3:end-2,3:end-2]
    # u_tru = F.(xs, λ * h)[3:end-2,3:end-2]

    err = norm(u_grad - u_tru) / norm(u_tru)
    println(err)

    return err
end

function single_step_grad_convergence_plot()
    hs = 1 ./ 2 .^ (1:14)
    errs = single_step_grad.(hs, .5)

    fig = Figure()
    ax = Axis(
        fig[1,1],
        xreversed = true,
        xscale = log10,
        yscale = log10,
    )
    lines!(ax, hs, first.(errs), label = "u")
    lines!(ax, hs, last.(errs), label = "du")
    lines!(ax, hs, hs.^3, label = "O(h^3)", linestyle = :dash)
    lines!(ax, hs, hs.^4, label = "O(h^4)", linestyle = :dash)

    Legend(fig[1,2], ax)

    return fig
end

function single_step_two_convergence_plot()
    hs = 1 ./ 2 .^ (1:14)
    errs = single_step_two.(hs, .5)

    fig = Figure()
    ax = Axis(
        fig[1,1],
        xreversed = true,
        xscale = log10,
        yscale = log10,
    )
    lines!(ax, hs, first.(errs), label = "u")
    lines!(ax, hs, hs.^4, label = "O(h^4)", linestyle = :dash)

    Legend(fig[1,2], ax)

    return fig
end
