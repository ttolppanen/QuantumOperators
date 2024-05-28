# using SparseArrays
# using LinearAlgebra

export apply_diss_deco!

function apply_diss_deco!(state::AbstractVector{<:Number}, operators::Vector{AbstractMatrix{<:Number}})::Int64
    probabilities = []
    n_outcomes = length(operators)
    for i in 1:n_outcomes
        op = operators[i]
        prob = calc_msr_probability(op, state) # this should calculate the norm of op * |s>
        push!(probabilities, prob)
    end
    probabilities .= probabilities ./ sum(probabilities)
    result = findfirst(cumsum(real.(probabilities)) .> rand())
    state .= operators[result] * state
    return result
end