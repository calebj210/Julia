#=
# Collection of functions for visualizing properties of complex functions
#
# Author: Caleb Jacobs
# DLM: August 14, 2024
=#

using PlotlyJS
using Colors

function complexAbsPlot(z⃗, f⃗; logscale = false, title = "", smooth = "none")
    x⃗ = real(z⃗[:])
    y⃗ = imag(z⃗[:])
    if logscale
        z⃗ = log10.(abs.(f⃗))
        z⃗[isinf.(z⃗)] .= -17
    else
        z⃗ = abs.(f⃗)
    end

    layout = Layout(
        width = 800, height = 800,
        xaxis = attr(
            title = attr(
                text = "Re(z)",
                font_size = 20
            )
        ),
        yaxis = attr(
            title = attr(
                text = "Im(z)",
                font_size = 20
            )
        ),
        title = attr(
            text = title,
            font_size = 25,
            y = 0.96,
            x = 0.5,
            xanchor = "center",
            yanchor = "top"
        ),
        font=attr(
            size = 18
        )
    )

    if logscale
        plt = plot(heatmap(
                x = x⃗,
                y = y⃗,
                z = z⃗,
                zsmooth = smooth,
                zmin = -16, zmax = 1,
                colorscale = colors.viridis),
               layout)
    else
        plt = plot(heatmap(
                x = x⃗,
                y = y⃗,
                z = z⃗,
                zsmooth = smooth,
                colorscale = colors.viridis),
               layout)
    end
    relayout!(plt, template = "plotly_white")
    
    return plt
end

function complexPlot(z⃗, f⃗)
    x⃗ = real(z⃗[:])
    y⃗ = imag(z⃗[:])
    args = angle.(f⃗[:])
    args[args .< 0] .+= 2π

    plt = plot(heatmap(
            x = x⃗,
            y = y⃗,
            z = args,
            zsmooth = "none",
            zmin = 0, zmax = 2π,
            colorscale = colors.hsv,
            colorbar = attr(
                tickmode = "array",
                tickvals = [0, π, 2π],
                ticktext = ["0", "π", "2π"])),
            Layout(width = 800, height = 800))
    relayout!(plt, template = "plotly_white")
 
    return plt
end

function colorWheel()
    x = range(-1, 1, length = 501)

    z = [x + im * y for x ∈ x, y ∈ x]
    f = z
    f[abs.(z) .> 1] .= NaN

    x⃗ = real(z[:])
    y⃗ = imag(z[:])
    args = angle.(f[:])
    args[args .< 0] .+= 2π

    plt = plot(heatmap(
            x = x⃗,
            y = y⃗,
            z = args,
            zsmooth = "none",
            zmin = 0, zmax = 2π,
            colorscale = colors.hsv,
            showscale = false),
            Layout(
                width  = 1000, 
                height = 1000,
                yaxis_visible = false,
                xaxis_visible = false))

    relayout!(plt, template = "plotly_white")

    return plt
end

function complexPlot3d(z⃗::Matrix, f⃗::Matrix...; T = 1, exclude = false, mesh = false, logscale = false)
    x⃗, y⃗ = reim(z⃗)
    f⃗ = deepcopy(f⃗)
        
    plts = Vector{GenericTrace}()

    for p ∈ f⃗
        if exclude
            p[y⃗ .== 0 .&& x⃗ .> 1] .= NaN + NaN * im     # Exclude common branch cut
        end

        if T == 2
            z⃗ = real(p)
            args = zeros(length(z⃗))
        elseif T == 3
            z⃗ = imag(p)
            args = zeros(length(z⃗))
        else
            z⃗ = abs.(p)
            args = angle.(p)
            args[args .< 0] .+= 2π
        end
        
        if logscale
            z⃗ = sign.(z⃗) .* log10.(abs.(z⃗))
        end

        if !mesh
            push!(plts, surface(
                        x = x⃗, y = y⃗, z = z⃗,
                        connectgaps = false,
                        surfacecolor = args,
                        cmin = 0, cmax = 2π,
                        colorscale = colors.hsv,
                        showscale = false))
        else
            xmin, xmax = extrema(x⃗)
            ymin, ymax = extrema(y⃗)
            hx = (xmax - xmin) / (size(z⃗)[1] + 1)
            hy = (ymax - ymin) / (size(z⃗)[2] + 1)

            push!(plts, surface(
                        x = x⃗, y = y⃗, z = z⃗,
                        connectgaps = false,
                        opacity = 0,
                        showscale = false,
                        contours = attr(
                            x = attr(
                                show = true,
                                start = xmin,
                                size  = hx,
                                color = "black",
                                width = 1),
                            x_end = xmax,
                            y = attr(
                                show = true,
                                start = ymin,
                                size  = hy, 
                                color = "black",
                                width = 1),
                            y_end = ymax)
                        )
            )

        end
    end

    if T == 2
        title = "Re(f)"
    elseif T == 3
        title = "Im(f)"
    else
        title = "Abs-Arg(f)"
    end
    
    layout = Layout(
        width = 800, height = 800,
        title = attr(
            text = title,
            font_size = 28,
            y = 0.86,
            x = 0.5,
            xanchor = "center",
            yanchor = "top",
        ),
        scene = attr(
            xaxis_title = "Re(z)",
            yaxis_title = "Im(z)",
            zaxis_title = "Abs(f)",
            zaxis = attr(
                range = (T == 1) ? [0, 20] : [-10, 10]
            ),
        ),
        font=attr(
            size = 15
        ),
        scene_aspectratio=attr(x = 1, y = 1, z = 1)
    )

    plt = plot(plts, layout)
    relayout!(plt, template = "plotly_white")

    return plt
end

function complexPlot3d(z::Vector, f::Vector...; T = 1, exclude = false, mesh = false, logscale = false)
    Z = reshape(z, round(Int64, sqrt(length(z))), :)

    N = round(Int64, sqrt(length(f[1])))

    F = reshape.(f, N , :)

    return complexPlot3d(Z, F..., T = T, exclude = exclude, mesh = mesh, logscale = logscale)
end
