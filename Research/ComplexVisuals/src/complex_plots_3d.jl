"Complex absolute value-argument plot of complex valued functions"
@recipe ComplexSurface (z, f) begin
    "[(W)GLMakie only] Specifies whether the surface matrix gets sampled with interpolation."
    interpolate = false

    Makie.mixin_generic_plot_attributes()...
    Makie.mixin_shading_attributes()...
    Makie.mixin_colormap_attributes()...

    colormap = :hsv
    colorrange = (0, 2π)
end

Makie.preferred_axis_type(plot::ComplexSurface) = Axis3

function Makie.plot!(cs::ComplexSurface{<:Tuple{ComplexGrid, AbstractMatrix{<:Number}}})
    @extractvalue cs (z, f)

    height = abs.(f)
    θ = mod.(angle.(f), 2π)

    valid_attributes = Makie.shared_attributes(cs, Surface)
    surface!(cs, z.real, z.imag, height;
             valid_attributes...,
             color = θ
            )

    return cs
end

function Makie.convert_arguments(P::Type{<:ComplexSurface}, z::ComplexGrid, f::Function) 
    return convert_arguments(P, z, f.(z))
end
function Makie.convert_arguments(P::Type{<:ComplexSurface}, x::T, y::T, f) where {T <: AbstractRange{<:Real}}
    return convert_arguments(P, ComplexGrid(x, y), f)
end
