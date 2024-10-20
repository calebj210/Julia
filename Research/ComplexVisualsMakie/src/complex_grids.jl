struct ComplexGrid{T} <: AbstractMatrix{Complex{T}}
    "Real range"
    x::AbstractRange{T}

    "Imaginary range"
    y::AbstractRange{T}
end

Base.size(A::ComplexGrid) = (length(A.x), length(A.y))
Base.getindex(A::ComplexGrid, i::Int) = A.x[1 + (i - 1) % length(A.x)] + A.y[1 + (i - 1) ÷ length(A.x)] * im
Base.getindex(A::ComplexGrid, i::Int, j::Int) = A.x[i] + A.y[j] * im
Base.:*(A::ComplexGrid, c::Real) = ComplexGrid(c * A.x, c * A.y)
Base.:*(c::Real, A::ComplexGrid) = ComplexGrid(c * A.x, c * A.y)
Base.:+(A::ComplexGrid, c::Number) = ComplexGrid(A.x .+ real(c), A.y .+ imag(c))
Base.:-(A::ComplexGrid, c::Number) = ComplexGrid(A.x .- real(c), A.y .- imag(c))

"Create grid data in the complex plane"
function complex_grid(xlims::T, ylims::T, Nx::S, Ny::S) where {T <: NTuple{2, <: Real}, S <: Int}
    x = range(xlims..., length = Nx)
    y = range(ylims..., length = Ny)

    return ComplexGrid(x, y)
end
complex_grid(r::Real, N::Int) = complex_grid((-r,r), (-r,r), N, N)
