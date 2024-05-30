# using SparseArrays
# using LinearAlgebra

export apply_diss_deco!

# operators : list of operators related to dissipation and decoherence, which should be multiplied by
#             the corresponding rates, i.e. sqrt(kappa) * a.

function apply_diss_deco!(state::AbstractVector{<:Number}, operators::Vector{<:AbstractMatrix}; normalize = true)::Int64
    probabilities = []
    n_outcomes = length(operators)
    for i in 1:n_outcomes
        op = operators[i]
        prob = norm(op * state)
        push!(probabilities, prob)
    end
    probabilities .= probabilities ./ sum(probabilities)
    result = findfirst(cumsum(real.(probabilities)) .> rand())
    state .= operators[result] * state
    if normalize 
        normalize!(state)
    end
    return result
end