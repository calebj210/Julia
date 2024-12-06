function fourier_coefficients(f::Vector{<:Real})
    coeffs = rfft(f) / length(f)

    con = real(coeffs[1])
    cos_coeffs, sin_coeffs = 2 .* reim(coeffs[2:end])

    return (con, cos_coeffs, -sin_coeffs)
end
