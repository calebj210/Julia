using ComplexVisuals, GLMakie

function taylor_recurrence(a, cn, h, maxn = nothing)
    crit_flag = false
    S = cn
    for n ∈ 1:(isnothing(maxn) ? 1000 : maxn)
        cn *= -a * h / n
        if !isnothing(maxn)
            if abs(cn / S) <= eps()
                if crit_flag
                    break
                else
                    crit_flag = true
                end
            else
                crit_flag = false
            end
        end
        S += cn
    end

    return S
end

function taylor_exp(a, z, mstep = .1, maxn = nothing)
    if abs(z) <= 0.1
        return exp(-a * z)
    end
    z0 = 0.1sign(z)
    f0 = exp(-a * z0)
    for n ∈ 1:1000
        h_sing = abs(z0) * exp(-2)
        h_dir  = abs(z - z0)
        h = sign(z - z0) * min(h_sing, h_dir, mstep)

        f0 = taylor_recurrence(a, f0, h, maxn)
        if h == sign(z - z0) * h_dir
            break
        end
        z0 += h
    end

    return f0
end

function test_exp(a, mstep = .1, maxn = nothing)
    z = complex_square_grid(1, 200)
    f = taylor_exp.(a, z, mstep, maxn)
    g = exp.(-a * z)
    
    err = abs.((f - g) ./ g)

    set_theme!(theme_latexfonts())
    fig = Figure()
    ax = Axis(fig[1,1], xlabel = L"\mathrm{Re}(z)", ylabel = L"\mathrm{Im}(z)", title = "Exponential Error with a = $(a)")
    plt = heatmap!(ax, z, err, colorscale = log10, colorrange = (1e-16, 1e1))
    Colorbar(fig[1,2], plt)

    return fig
end
