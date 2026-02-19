centered_nodes(x::T, y::T, z::T; type::Type = Float64) where T <: Integer = 
    vec([type.((i,j,k)) for i in -x:x, j in -y:y, k in -z:z])

centered_nodes(x::T, y::T; type::Type = Float64) where T <: Integer = 
    vec([type.((i,j)) for i in -x:x, j in -y:y])

centered_nodes(x::T; type::Type = Float64) where T <: Integer = 
    [type.((i,)) for i in -x:x]
