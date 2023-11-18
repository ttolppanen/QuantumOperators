# using SparseArrays
# using DataStructures
# using LinearAlgebra

export subspace_info
export subspace_indices
export subspace_tools
export find_subspace
export subspace_split
export total_boson_number_subspace_info
export total_boson_number_subspace_indices
export total_boson_number_subspace_tools
export measurement_subspace
export feedback_measurement_subspace

# here find_property takes in the basis vector index as an argument, instead of the state
function subspace_info(dim::Int, find_property::Function)
    out = SortedDict()
    for i in 1:dim
        property = find_property(i)
        if haskey(out, property)
            push!(out[property], i)
        else
            out[property] = [i]
        end
    end
    return out
end

function subspace_indices(dim::Int, find_property::Function)
    dict = subspace_info(dim, find_property)
    return collect(values(dict))
end

function subspace_tools(dim::Int, find_property::Function)
    dict = subspace_info(dim, find_property)
    return subspace_tools(dict)
end
function subspace_tools(subspace_dict::SortedDict)
    dim = sum(length.(collect(values(subspace_dict))))
    permutation_matrix = spzeros(dim, dim)
    subspace_ranges::Vector{UnitRange} = []
    range_end = 0
    for indices in values(subspace_dict)
        range_start = range_end + 1
        range_end = range_start + length(indices) - 1
        range = range_start:range_end
        push!(subspace_ranges, range)
        for (i, val) in pairs(indices)
            permutation_matrix[range_start + i - 1, val] = 1
        end
    end
    return permutation_matrix, subspace_ranges
end

function find_subspace(state, subspace_indices::Dict; kwargs...)
    return find_subspace(state, collect(values(subspace_indices)); kwargs...)
end
function find_subspace(state, subspace_indices; digit_error::Int = 15, id_initial_guess::Int = 1, iterate_order::Int = 1)
    d = length(subspace_indices)
    if iterate_order > 0
        subspace_id = circshift(1:iterate_order:d, -(id_initial_guess - 1))
    else
        subspace_id = circshift(d:iterate_order:1, -(d - id_initial_guess))
    end
    for i in subspace_id
        range = subspace_indices[i]
        n = norm(@view(state[range]))
        if round(n; digits = digit_error) == 1.0
            return i
        end
    end
    throw(ErrorException("Could not find subspace. Maybe digit_error = $digit_error is too large."))
end

function subspace_split(op::AbstractArray, subspace_indices::Array)
    out::Vector{typeof(op)} = []
    rank = length(size(op)) # differentiates between vectors, matrices etc...
    for indices in subspace_indices
        push!(out, op[[indices for _ in 1:rank]...]) # vectors -> op[indices], matrices -> op[indices, indices],
    end
    return out
end

function total_boson_number_subspace_functions(d::Int, L::Int, f::Function)
    op = nall(d, L)
    find_property(i) = real(op[i, i])
    return f(d^L, find_property)
end
total_boson_number_subspace_info(d::Int, L::Int) = total_boson_number_subspace_functions(d, L, subspace_info)
total_boson_number_subspace_indices(d::Int, L::Int) = total_boson_number_subspace_functions(d, L, subspace_indices)
total_boson_number_subspace_tools(d::Int, L::Int) = total_boson_number_subspace_functions(d, L, subspace_tools)


function measurement_subspace(msrop::MsrOpMatrixType, subspace_indices::Array)
    out::Vector{MsrOpMatrixType} = []
    for id in eachindex(subspace_indices)
        push!(out, [])
        for L in eachindex(msrop)
            push!(out[id], [])
            for msr_op in msrop[L]
                indices = subspace_indices[id]
                push!(out[id][L], msr_op[indices, indices])
            end
        end
    end
    return out
end

function feedback_measurement_subspace(feedback::Vector{<:AbstractMatrix}, msrop::MsrOpMatrixType, subspace_indices::Array; digit_error = 15, kwargs...) # kwargs for find_subspace
    d = size(feedback[1])[1]
    out = []
    for id in eachindex(subspace_indices)
        push!(out, [])
        for L in eachindex(msrop)
            push!(out[id], [])
            for msr_outcome in eachindex(msrop[L])
                fb_op = feedback[L]
                indices = subspace_indices[id]
                msr_op = msrop[L][msr_outcome]
                subspace_d = length(indices)

                state = complex(spzeros(d))
                state[indices] .= ones(subspace_d) ./ sqrt(subspace_d)
                state[indices] .= msr_op[indices, indices] * state[indices] # measurements are assumed to stay in the subspace
                if round(norm(state); digits = digit_error) == 0.0
                    push!(out[id][L], (id, spzeros(1,1)))
                    continue
                end
                normalize!(state)
                # measurements should not be able to change the subspace, but since here the measurement probability is not calculated,
                # measurements that are not possible can happen. This shouldn't affect calculations, it just means that
                # out has matrices that won't be used.
                state .= fb_op * state
                normalize!(state)
                
                new_id = find_subspace(state, subspace_indices; digit_error, kwargs...)
                feedback_in_subspace = fb_op[subspace_indices[new_id], indices]
                push!(out[id][L], (new_id, feedback_in_subspace))
            end
        end
    end
    # structure of out is: out[subspace_id][site][msr_outcome]
    # the content is a tuple, which contains the id of the subspace after feedback and that specific measurement,
    # and the feedback operator in that subspace: (id, fb_op).
    # The dimensions of the fb_op in subspace is d(fb_op) = d(new_subspace) * d(old_subspace),
    # so it takes the vector a new subspace when it operates.
    return out, measurement_subspace(msrop, subspace_indices)
end