using PlotlyJS
using Images, Base64
const CS = PlotlyJS.PlotlyBase.ColorSchemes

const plot_template_3d = (;
    template = "plotly_white",
    width  = 1000, height = 1000,

    title = attr(
        font_size = 25,
        y = 0.8,
        x = 0.5
    ),

    font_size = 15,

    scene = attr(
        xaxis_title = attr(
            font_size = 20,
        ),
        yaxis_title = attr(
            font_size = 20,
        ),
        zaxis_title = attr(
            font_size = 20,
        ),
        aspectratio = attr(x = 1, y = 1, z = .6),
    )
)

const plot_template_3d_complex = (;
    plot_template_3d...,

    scene_xaxis_title_text = "Re(z)",
    scene_yaxis_title_text = "Im(z)"
)

function complex_surface_plot(z::T, f::T; kwargs...) where T <: Matrix{<: Number}
    x, y   = vec.(reim(z))
    absval = vec(abs.(f))
    phase  = mod.(angle.(f), 2π)

    i, j, k, c = grid_mesh(phase)

    trace = mesh3d(;
        x = x,
        y = y,
        z = absval,
        i = i,
        j = j,
        k = k,

        intensity = c,
        intensitymode = "cell",
        colorscale = CS.hsv,
        cmin = 0, cmax = 2π,
        showscale = false,
    )

    png_data = read("$(@__DIR__)/AngleWheel.png", String) # Read the PNG file as a string
    png_base64 = base64encode(png_data)

    layout = Layout(;
        plot_template_3d_complex...,

        scene_zaxis_title_text = "abs(f)",

        images = [attr(
            source = "data:image/png;base64,$png_base64",
            xref = "paper",  # Reference the paper coordinates (0 to 1)
            yref = "paper",
            x = 0.8,         # Position horizontally (0.5 is center)
            y = 0.2,         # Position vertically
            sizex = 0.15,      # Size relative to paper width
            sizey = 0.15,      # Size relative to paper height
            opacity = 0.8,    # Adjust transparency as needed
            layer = "above"   # Place the image above the 3D graph
        )],

        kwargs...
    )

    plot(trace, layout)
end

function complex_reim_surface_plot(z::T, f::T; kwargs...) where T <: Matrix{<: Number}
    x, y = reim(z)
    f_re, f_im = reim(f)

    xmin, xmax = extrema(x)
    ymin, ymax = extrema(y)
    hx = abs(x[1,1] - x[1,2])
    hy = abs(y[1,1] - y[2,1])

    trace_real = surface(
        x = x, y = y, z = f_re,
        connectgaps = false,
        opacity = 0,
        showscale = false,
        contours = attr(
            x = attr(
                show = true,
                size = hx,
                color = "black",
                width = 1
            ),
            x_start = xmin, x_end = xmax,
            y = attr(
                show = true,
                size = hy,
                color = "black",
                width = 1
            ),
            y_start = ymin, y_end = ymax,
        )
    )

    trace_imag = surface(
        x = x, y = y, z = f_im,
        connectgaps = false,
        opacity = 0,
        showscale = false,
        contours = attr(
            x = attr(
                show = true,
                size = hx,
                color = "black",
                width = 1
            ),
            x_start = xmin, x_end = xmax,
            y = attr(
                show = true,
                size = hy,
                color = "black",
                width = 1
            ),
            y_start = ymin, y_end = ymax,
        )
    )

    pltre = Plot(trace_real)
    relayout!(pltre; plot_template_3d_complex..., kwargs..., scene_zaxis_title_text = "Re(f)")
    pltim = Plot(trace_imag)
    relayout!(pltim; plot_template_3d_complex..., kwargs..., scene_zaxis_title_text = "Im(f)")

    plot_real = plot(pltre)
    plot_imag = plot(pltim)

    return (plot_real, plot_imag)
end
