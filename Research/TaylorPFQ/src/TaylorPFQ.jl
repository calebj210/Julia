module TaylorPFQ

using Polynomials
using SpecialFunctions
using MathLink

include("pFq.jl")
export taylor_pFq,
       mathematica_pFq

end
