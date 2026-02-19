onesided_nodes(x::T, y::T, z::T; type::Type = Float64) where T <: Integer = 
    vec([type.((i,j,k)) for i in -x:x, j in -y:y, k in 0:z])

onesided_nodes(x::T, y::T; type::Type = Float64) where T <: Integer = 
    vec([type.((i,j)) for i in -x:x, j in 0:y])

onesided_nodes(x::T; type::Type = Float64) where T <: Integer = 
    [type.((i,)) for i in 0:x]
