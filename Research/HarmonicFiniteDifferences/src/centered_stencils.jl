function centered_nodes(x::T, y::T, z::T; extended_precision = false) where T <: Integer
    if extended_precision
        return vec([BigFloat.((i,j,k)) for i in -x:x, j in -y:y, k in -z:z])
    else
        return vec([(i,j,k) for i in -x:x, j in -y:y, k in -z:z])
    end
end
centered_nodes(N; kwargs...) = centered_nodes(N, N, N; kwargs...)
