# using SparseArrays
# using DataStructures

export generate_subspace_info
export find_subspace
export total_boson_number_subspace

function generate_subspace_info(dim::Int, find_property::Function)
    subspace_dict = SortedDict()
    for i in 1:dim
        state = spzeros(dim)
        state[i] = 1.0
        property = find_property(state)
        if haskey(subspace_dict, property)
            push!(subspace_dict[property], i)
        else
            subspace_dict[property] = [i]
        end
    end
    permutation_matrix = spzeros(dim, dim)
    subspace_ranges = []
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
    return subspace_dict, permutation_matrix, subspace_ranges
end

function find_subspace(state, subspace_indices::Dict; kwargs...)
    return find_subspace(state, collect(values(subspace_indices)); kwargs...)
end
function find_subspace(state, subspace_indices; digit_error = 15)
    max_norm, index = findmax(x -> norm(@view state[x]), subspace_indices)
    if round(max_norm ; digits = digit_error) == 1.0
        return subspace_indices[index]
    else
        throw(ErrorException("Could not find subspace. Max norm was $max_norm, with index $index. digit_error = $digit_error is too large."))
    end
end

function total_boson_number_subspace(d::Int, L::Int)
    op = nall(d, L)
    find_property(state) = expval(state, op)
    return generate_subspace_info(d^L, find_property)
end