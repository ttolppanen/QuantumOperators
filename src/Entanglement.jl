using ITensors
using LinearAlgebra
using SparseArrays
#include("PartialTrace.jl")
#include("Utility/ConvertToReal.jl")

export entanglement

#cut : cut between the left and right bipartion of the system; cut = 2 cuts the system in to A = {1, 2} and B {3,..., L}
#      for entanglement for a density matrix, cut is the sites to trace over

function entanglement(d::Integer, L::Integer, state::AbstractVector{<:Number}, cut::Integer)
    m = schmidt_form(state, d^(cut), d^(L - cut))
    S = svdvals(m)
    out = 0.0
    for s in S
        p = real_with_warning(s^2)
        out -= (p + 1 ≈ 1.0 ? 0.0 : p * log(p)) 
    end
    return out
end
function schmidt_form(state::AbstractVector{<:Number}, d_a, d_b)
    mat_out = zeros(d_a, d_b)
    for b_i in 1:d_b
        for a_i in 1:d_a
            mat_out[a_i, b_i] = state[d_b * (a_i - 1) + b_i]
        end
    end
    return mat_out
end

function entanglement(d::Integer, L::Integer, state::AbstractMatrix{<:Number}, cut)
    rho_p = ptrace(d, L, state, cut)
    if issparse(rho_p)
        vals = eigvals(Matrix(rho_p))
    else
        vals = eigvals(rho_p)
    end
    out = 0.0
    for e_val in vals
        p = real_with_warning(e_val)
        out -= (p + 1 ≈ 1.0 ? 0.0 : p * log(p)) 
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
        p = real_with_warning(S[i, i]^2)
        out -= (p == 0.0 ? 0.0 : p * log(p)) 
    end
    return out
end