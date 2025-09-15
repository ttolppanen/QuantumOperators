# using ITensors
# using LinearAlgebra
# using SparseArrays
# include("PartialTrace.jl")
# include("Utility/ConvertToReal.jl")

export entanglement

# cut : cut between the left and right bipartion of the system; cut = 2 cuts the system in to A = {1, 2} and B {3,..., L}
#      for entanglement for a density matrix, cut is the sites to trace over
function schmidt_form(state::AbstractVector{<:Number}, d_a::Integer, d_b::Integer)
    mat_out = complex(zeros(d_a, d_b))
    for b_i in 1:d_b
        for a_i in 1:d_a
            mat_out[a_i, b_i] = state[d_b * (a_i - 1) + b_i]
        end
    end
    return mat_out
end

# sub_system : list of intergers referring to the sites; i.e. you want the schmidt matrix where A is the first, second and fourth site, you set sub_system = [1,2,4]
function schmidt_form(state::AbstractVector{<:Number}, d::Integer, L::Integer, sub_system::AbstractArray)
    if d^L > size(state)[1]
        throw(ArgumentError("d^L > dim(state), check you have set d and L correctly"))
    end
    sub_a = [L - (i - 1) for i in sub_system] # reverse order
    sort!(sub_a)
    sub_b = setdiff(1:L, sub_a)
    out = complex(zeros(d^length(sub_a), d^length(sub_b)))
    for i in eachindex(state)
        k_i = 1
        for sub_i in eachindex(sub_a)
            k_i += (floor(Int, (i - 1) / d^(sub_a[sub_i] - 1)) % d) * d^(sub_i - 1)
        end
        k_j = 1
        for sub_i in eachindex(sub_b)
            k_j += (floor(Int, (i - 1) / d^(sub_b[sub_i] - 1)) % d) * d^(sub_i - 1)
        end
        out[k_i, k_j] = state[i]
    end
    return out
end

function entanglement(d::Integer, L::Integer, states, cut)
    return [entanglement(d, L, s, cut) for s in states]
end
function entanglement(states::Vector{MPS}, cut::Integer)
    return [entanglement(s, cut) for s in states]
end
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
function entanglement(d::Integer, L::Integer, state::AbstractVector{<:Number}, sub_systems::AbstractArray) # look at schmidt_form for sub_systems
    m = schmidt_form(state, d, L, sub_systems)
    S = svdvals(m)
    out = 0.0
    for s in S
        p = real_with_warning(s^2)
        out -= (p + 1 ≈ 1.0 ? 0.0 : p * log(p)) 
    end
    return out
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
    ITensorMPS.orthogonalize!(mps, cut)
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