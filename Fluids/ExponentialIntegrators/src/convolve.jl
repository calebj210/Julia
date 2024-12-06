function convolve(a::Vector, b::Vector)
    N = length(a); M = length(b)
    
    A = [a; zeros(M - 1)]
    B = [b; zeros(N - 1)]

    conv = ifft(fft(A) .* fft(B)) 

    return conv
end
