### Simple algorithm to compute the inverse of a upper triangular matrix

# Packages to include
using LinearAlgebra

# Inverse solver
function computeInverse(U)
    # Size of U
    n = size(U,1)

    # Constuct n×n identity matrix
    A = Matrix{Float64}(I, (n,n))
    for i∈n:-1:1
        A[i,:] /= U[i,i]
        for j∈1:i-1
            A[j,:] -= U[j,i]*A[i,:]
        end
    end

    return A
end

computeInverse([1 2 3;
                0 1 2;
                0 0 1])
