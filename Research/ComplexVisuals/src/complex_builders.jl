"""
    complex_grid(xminmax, yminmax, Nx, Ny)
Generate a complex grid.
"""
function complex_grid(xminmax::Tuple{Real, Real}, yminmax::Tuple{Real, Real}, Nx::Int, Ny::Int)
    x = range(xminmax[1], xminmax[2], length = Nx)
    y = range(yminmax[1], yminmax[2], length = Ny)

    z = x' .+ y*im

    return z
end
complex_grid(xminmax::Tuple{Real, Real}, yminmax::Tuple{Real, Real}, N::Int) = complex_grid(xminmax, yminmax, N, N)
complex_grid(r::Real, N::Int) = complex_grid((-r, r), (-r, r), N)

"""
    grid_mesh(x)
Convert matrix grid `x` to a mesh triangulation suitable for mesh plots.
"""
function grid_mesh(x::Matrix{T}) where T <: Number
    L = LinearIndices(x) .- 1
    C = CartesianIndices(x)
    Ii = first(C)
    If = last(C)
    D = CartesianIndex(1, 0)
    R = CartesianIndex(0, 1)

    i = Vector{Int64}()
    j = Vector{Int64}()
    k = Vector{Int64}()
    c = Vector{T}()

    for n ∈ Ii : If - oneunit(If)
        push!(i, L[n])
        push!(j, L[n + D])
        push!(k, L[n + D + R])
        push!(c, (x[n]))

        push!(i, L[n])
        push!(k, L[n + R])
        push!(j, L[n + D + R])
        push!(c, (x[n]))
    end

    return (i, j, k, c)
end

function color_wheel()
    z = complex_grid(1, 500)
    x, y = vec.(reim(z))

    f = copy(vec(z))
    f[abs.(f) .> 1] .= NaN
    ang = mod.(angle.(f), 2π)

    trace = heatmap(
        x = x, y = y, z = ang,

        zsmooth = "none",
        colorscale = CS.hsv,
        showscale = false
    )

    layout = Layout(
        width = 1000, height = 1000,
        xaxis_visible = false,
        yaxis_visible = false
    )

    p = Plot(trace, layout)

    savefig(p, "AngleWheel.svg", width = 1000, height = 1000)
end
