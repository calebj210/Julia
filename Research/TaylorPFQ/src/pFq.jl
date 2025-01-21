function taylor_pFq(a, b, z)
end

function mathematica_pFq(a, b, z)
    val = weval( W`N[HypergeometricPFQ[a,b,z]]`, a = a, b = b, z = z)
    try 
        return Complex(val.args...)
    catch 
        return Real(val)
    end
end
