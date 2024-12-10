# using SparseArrays
# using LinearAlgebra

export apply_diss_deco!

# operators : list of operators related to dissipation and decoherence, which should be multiplied by
#             the corresponding rates, i.e. sqrt(kappa) * a.

function apply_diss_deco!(state::AbstractVector{<:Number}, operators::Vector{<:AbstractMatrix}; normalize = true)::Int64
    probabilities = []
    n_outcomes = length(operators)
    if n_outcomes == 1
        result = 1
    else
        for i in 1:n_outcomes
            op = operators[i]
            prob = norm(op * state)^2 # norm^2 = <psi|psi>
            push!(probabilities, prob)
        end
        probabilities .= probabilities ./ sum(probabilities)
        result = findfirst(cumsum(real.(probabilities)) .> rand())
    end
    state .= operators[result] * state
    if normalize 
        normalize!(state)
    end
    return result
end
# diss_op and deco_op should be of the form Vector{<:Vector{<:AbstractMatrix}}
function apply_diss_deco!(state::Vector{<:AbstractVector{<:Number}}, subspace_id::Integer, diss_op, deco_op; normalize = true)::Tuple{Int64, Int64}
    probabilities = []
    n_diss = length(diss_op)
    n_deco = length(deco_op)

    if (n_diss + n_deco) == 1
        result = 1
    else
        # calculating probabilities
        for i in 1:n_diss
            op = diss_op[i]
            prob = norm(op[subspace_id] * state[subspace_id])^2 # norm^2 = <psi|psi>
            push!(probabilities, prob)
        end
        for i in 1:n_deco
            op = deco_op[i]
            prob = norm(op[subspace_id] * state[subspace_id])^2 # norm^2 = <psi|psi>
            push!(probabilities, prob)
        end

        # normalizing and drawing the result
        probabilities .= probabilities ./ sum(probabilities)
        result = findfirst(cumsum(real.(probabilities)) .> rand())
    end

    # if dissipation happened, we need to change the subspace
    if result <= length(diss_op)
        new_id = subspace_id - 1
        mul!(state[new_id], diss_op[result][subspace_id], state[subspace_id])
        if normalize
            normalize!(state[new_id])
        end
        return new_id, result
    else
        state[subspace_id] .= deco_op[result - n_diss][subspace_id] * state[subspace_id]
        if normalize 
            normalize!(state[subspace_id])
        end
        return subspace_id, result
    end
end