function cardinalweights(nodes::Vector{<:Tuple{Vararg{T}}}; N = length(nodes)) where T<:Number
    A = vand(nodes; N = N - 1)'

    w = Vector{Vector{T}}(undef, N)
    for n in 1:N
        b = zeros(N)
        b[n] = 1

        w[n] = A \ b
    end

    return w
end

function cardinalfunction(x::Tuple, w)
    return sum(w[n + 1] * harmonic(x, n) for n ∈ 0:length(w) - 1)
end
