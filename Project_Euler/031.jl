using Polynomials

"""
    polinv(a, d = 16)
Compute the polynomial inverse of `a` to degree `d`.
"""
function polinv(a::AbstractPolynomial, d = 16)
    b = Polynomial([1 / a[0]])
    
    for m ∈ 1 : d
        nextterm = -sum([a[k] * b[m - k] for k ∈ 1 : m]) / a[0]

        b[m] = nextterm
    end


    return b
end

"Count number of money combinations using generating functions"
function countmoney()
    ms = [1, 2, 5, 10, 20, 50, 100, 200]     # Money values

    num = SparsePolynomial(Dict(sum(ms) => 1))
    den = prod([SparsePolynomial(Dict(0 => 1, m => -1)) for m ∈ ms])

    expansion = num * polinv(den, 200)

    return round(Int, expansion[200 + sum(ms)])
end

countmoney()
