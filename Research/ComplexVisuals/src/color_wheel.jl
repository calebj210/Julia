function Colorwheel(fig_or_scene; kwargs...)
    ax = Axis(fig_or_scene;
        backgroundcolor = :transparent,
        aspect = AxisAspect(1),
        kwargs...
    )
    hidespines!(ax)
    hidedecorations!(ax)

    z = complex_square_grid(1, 500)
    f = collect(z)
    f[abs.(z) .>= 1] .= NaN

    phase!(ax, z, f)

    return ax
end

function axiscolorwheel(ax; position = :rb, kwargs...)
    wheel = Colorwheel(
        ax.parent;
        bbox = ax.scene.viewport,
        width = Relative(.2), height = Relative(.2),
        wheel_position_to_aligns(position)...,
        kwargs...
    )

    translate!(wheel.blockscene, 0, 0, 150)

    return wheel
end

function wheel_position_to_aligns(s::Symbol)
    p = string(s)
    length(p) != 2 && throw(ArgumentError("Position symbol must have length == 2"))

    haligns = Dict(
        'l' => .05,
        'r' => .95,
        'c' => .5,
    )
    haskey(haligns, p[1]) || throw(ArgumentError("First letter can be l, r or c, not $(p[1])."))

    valigns = Dict(
        't' => .95,
        'b' => .05,
        'c' => .5,
    )
    haskey(valigns, p[2]) || throw(ArgumentError("Second letter can be b, t or c, not $(p[2])."))

    return (halign = haligns[p[1]], valign = valigns[p[2]])
end
