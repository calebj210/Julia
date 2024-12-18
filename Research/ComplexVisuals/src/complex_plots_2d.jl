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
    @extractvalue ph (z, f)

    θ = mod.(angle.(f), 2π)

    valid_attributes = Makie.shared_attributes(ph, Heatmap)

    heatmap!(ph, z.real, z.imag, θ;
            valid_attributes...
           )

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

function Makie.convert_arguments(P::Type{<: Phase}, z::ComplexGrid, f::Function) 
    return convert_arguments(P, z, f.(z))
end

function Makie.convert_arguments(P::Type{<: Phase}, x::AbstractRange{<: Real}, y::AbstractRange{<: Real}, f) 
    return convert_arguments(P, ComplexGrid(x, y), f)
end

function Makie.convert_arguments(P::Type{<: Heatmap}, z::ComplexGrid, f::Matrix{<: Real}) 
    return convert_arguments(P, z.real, z.imag, f)
end

function Makie.convert_arguments(P::Type{<: Heatmap}, z::ComplexGrid, f::Function) 
    return convert_arguments(P, z.real, z.imag, f)
end
