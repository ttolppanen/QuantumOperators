using ITensors
using LinearAlgebra
using SparseArrays
#include("PartialTrace.jl")

export entanglement

#cut : cut between the left and right bipartion of the system; cut = 2 cuts the system in to A = {1, 2} and B {3,..., L}

function entanglement(d, L, state, cut)
    rho_p = ptrace(d, L, state, (cut+1):L)
    if issparse(rho_p)
        vals = eigvals(Matrix(rho_p))
    else
        vals = eigvals(rho_p)
    end
    out = 0.0
    for e_val in vals
        out -= (e_val + 1 ≈ 1.0 ? 0.0 : e_val * log(e_val)) 
    end
    return out
end
function entanglement(mps::MPS, cut::Integer)
    ITensors.orthogonalize!(mps, cut)
    if cut == 1
        S = svd(mps[cut], siteind(mps, cut)).S
    else
        S = svd(mps[cut], (linkind(mps, cut - 1), siteind(mps, cut))).S
    end
    out = 0.0
    for i in 1:dim(S, 1)
        p = S[i, i]^2
        out -= (p == 0.0 ? 0.0 : p * log(p)) 
    end
    return out
end