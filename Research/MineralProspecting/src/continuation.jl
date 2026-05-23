function apply_operator(F::AbstractArray{TF, N}, A::AbstractArray{TA, N}) where {TF, TA, N}
    ny, nx = size(A, 1), size(A, 2)
    h, w = size(F, 1), size(F, 2)
    
    # Assert dimensions match for all depth/additional dimensions beyond 2D
    if N > 2
        for d in 3:N
            @assert size(F, d) == size(A, d) "Dimension $d of stencil F ($(size(F, d))) must match A ($(size(A, d)))."
        end
    end
    
    # Calculate padding based on stencil dimensions (supports both odd and even dimensions)
    r_pad_top = (h + 1) ÷ 2 - 1
    r_pad_bot = h - r_pad_top - 1
    c_pad_left = (w + 1) ÷ 2 - 1
    c_pad_right = w - c_pad_left - 1
    
    # Create periodic ghost-node padded matrix in the first two dimensions
    row_indices = [mod1(i, ny) for i in (1 - r_pad_top):(ny + r_pad_bot)]
    col_indices = [mod1(j, nx) for j in (1 - c_pad_left):(nx + c_pad_right)]
    
    # Dynamically slice the trailing dimensions using static compile-time tuple resolution
    A_padded = A[row_indices, col_indices, ntuple(i -> :, N-2)...]
    
    C = zeros(promote_type(eltype(F), eltype(A)), ny, nx)
    
    @inbounds for I in CartesianIndices(C)
        val = zero(eltype(C))
        for J in CartesianIndices(F)
            # Construct a dynamic CartesianIndex: shift the first two dimensions, splat the rest
            idx = CartesianIndex(I[1] + J[1] - 1, I[2] + J[2] - 1, Tuple(J)[3:end]...)
            val += F[J] * A_padded[idx]
        end
        C[I] = val
    end
    return C
end

function gradient_stencils(h, λ)
    # I = [
    #         0 0 0 0 0;
    #         0 0 0 0 0;
    #         0 0 1 0 0;
    #         0 0 0 0 0;
    #         0 0 0 0 0
    # ]
    # Δ = [
    #      0 -2  -1   -2   0;
    #     -2  16  52   16 -2;
    #     -1  52 -252  52 -1;
    #     -2  16  52   16 -2;
    #      0 -2  -1   -2   0
    # ]
    # Δ2 = [
    #     0  1   1   1  0;
    #     1 -2  -10 -2  1;
    #     1 -10  36 -10 1;
    #     1 -2  -10 -2  1;
    #     0  1   1   1  0
    # ]
    #
    # Au = I - λ^2 / 120 * Δ + λ^4 / 72 * Δ2
    # Adu = λ * h * (I - λ^2 / 360 * Δ + λ^4 / 360 * Δ2)
    #
    # dAu = 1 / (λ * h) * (-λ^2 / 60 * Δ + λ^4 / 18 * Δ2)
    # dAdu = I - λ^2 / 120 * Δ + λ^4 / 72 * Δ2

    I = [
        0 0 0;
        0 1 0;
        0 0 0
    ]
    Δ = [
        1  4  1;
        4 -20 4;
        1  4  1
    ]

    Au = I - (λ^2 / 12) * Δ
    Adu = λ * h * (I - (λ^2 / 36) * Δ)
    dAu = -λ / (6h) * Δ
    dAdu = I - λ^2 / 12 * Δ

    return (; Au, dAu, Adu, dAdu)
end

function two_layer_stencils(_, λ)
    A = zeros(3, 3, 2)
    A[:, :, 1] = [0 0 0; 0 2 0; 0 0 0] + 
        λ^2 * [0 -1 0; -1 4 -1; 0 -1 0]
    A[:, :, 2] = [0 0 0; 0 -1 0; 0 0 0]

    return A
end

function gradient_continue(; h, u₀, du₀, zf = 2, N)
    zs = range(0, zf; length = N)
    λ = Float64(zs.step) / h

    u = copy(u₀)
    du = copy(du₀)

    D = gradient_stencils(h, λ)
    
    z0 = 0

    for _ in 2:N
        z0 += Float64(zs.step)
        utmp = apply_operator(D.Au, u) + apply_operator(D.Adu, du)
        dutmp = apply_operator(D.dAu, u) + apply_operator(D.dAdu, du)

        u = utmp
        du = dutmp
    end

    return u
end

function two_layer_continue(h, Δz, u₀, zf)
    λ = Δz / h

    u = copy(u₀)

    D = two_layer_stencils(h, λ)

    for _ in 0:-Δz:zf
        utmp = apply_operator(D, u)

        u = cat(utmp, u[:,:,1]; dims = 3)
    end

    return u[:,:,1]
end
