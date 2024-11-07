@recipe ComplexSurface (z, f) begin
    "[(W)GLMakie only] Specifies whether the surface matrix gets sampled with interpolation."
    interpolate = false

    MakieCore.mixin_generic_plot_attributes()...
    MakieCore.mixin_shading_attributes()...
    MakieCore.mixin_colormap_attributes()...

    colormap = :hsv
    colorrange = (0, 2π)
end

Makie.preferred_axis_type(plot::ComplexSurface) = Axis3

function Makie.plot!(cs::ComplexSurface{<:Tuple{ComplexGrid, AbstractMatrix{<:Number}}})
    @extractvalue cs (z, f)

    height = abs.(f)
    θ = mod.(angle.(f), 2π)

    surface!(cs, z.real, z.imag, height;
             attributes(cs)...,
             color = θ
            )

    return cs
end

Makie.convert_arguments(P::Type{<:ComplexSurface}, z::ComplexGrid, f::Function) = convert_arguments(P, z, f.(z))
Makie.convert_arguments(P::Type{<:ComplexSurface}, x::AbstractRange{<:Real}, y::AbstractRange{<:Real}, f) = convert_arguments(P, ComplexGrid(x, y), f)
