"Complex phase portrait"
@recipe Phase (z, f) begin
    "Sets whether colors should be interpolated"
    interpolate = false

    MakieCore.mixin_generic_plot_attributes()...
    MakieCore.mixin_colormap_attributes()...

    colormap = :hsv
    colorrange = (0, 2π)
end

function Makie.plot!(ph::Phase{<: Tuple{ComplexGrid, AbstractMatrix{<: Number}}})
    x = to_value(ph.z).real
    y = to_value(ph.z).imag
    f = to_value(ph.f)

    θ = mod.(angle.(f), 2π)

    heatmap!(ph, x, y, θ;
             attributes(ph)...)

    return ph
end

const complex_theme = Theme(
    Axis = (
        xautolimitmargin = (0, 0),
        yautolimitmargin = (0, 0),
        xlabel = L"Re$(z)$",
        ylabel = L"Im$(z)$",
    ),
    
    Axis3 = (
        xlabel = L"Re$(z)$",
        ylabel = L"Im$(z)$",
        zlabel = L"abs$(f)$",
    )
)

Makie.convert_arguments(P::Type{<: Phase}, z::ComplexGrid, f::Function) = convert_arguments(P, z, f.(z))
Makie.convert_arguments(P::Type{<: Phase}, x::AbstractRange{<: Real}, y::AbstractRange{<: Real}, f) = convert_arguments(P, ComplexGrid(x, y), f)
Makie.convert_arguments(P::Type{<: Heatmap}, z::ComplexGrid, f::Matrix{<: Real}) = convert_arguments(P, z.real, z.imag, f)
Makie.convert_arguments(P::Type{<: Heatmap}, z::ComplexGrid, f::Function) = convert_arguments(P, z.real, z.imag, f)
