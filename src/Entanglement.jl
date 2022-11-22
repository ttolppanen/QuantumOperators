using ITensors
using LinearAlgebra
using SparseArrays
#include("PartialTrace.jl")

export entanglement

#cut : cut between the left and right bipartion of the system; cut = 2 cuts the system in to A = {1, 2} and B {3,..., L}

function entanglement(state::AbstractVector{<:Number})
    if issparse(state)
        S = svdvals(state * state')
    else
    end
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