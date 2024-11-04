struct ComplexGrid{T} <: AbstractMatrix{Complex{T}}
    "Real range"
    real::AbstractRange{T}

    "Imaginary range"
    imag::AbstractRange{T}
end

Base.size(A::ComplexGrid) = (length(A.real), length(A.imag))
Base.getindex(A::ComplexGrid, i::Int) = A.real[1 + (i - 1) % length(A.real)] + 
                                        A.imag[1 + (i - 1) ÷ length(A.real)] * im
Base.getindex(A::ComplexGrid, i::Int, j::Int) = A.real[i] + A.imag[j] * im
Base.:*(A::ComplexGrid, c::Real) = ComplexGrid(c * A.real, c * A.imag)
Base.:*(c::Real, A::ComplexGrid) = ComplexGrid(c * A.real, c * A.imag)
Base.:+(A::ComplexGrid, c::Number) = ComplexGrid(A.real .+ real(c), A.imag .+ imag(c))
Base.:-(A::ComplexGrid, c::Number) = ComplexGrid(A.real .- real(c), A.imag .- imag(c))

"Create a square grid of radius `r` of size `N`×`N`."
complex_square_grid(r::Real, N::Int) = ComplexGrid(range(-r, r, N), range(-r, r, N))
