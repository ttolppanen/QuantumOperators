# using SparseArrays
# using DataStructures

export subspace_info
export subspace_tools
export find_subspace
export split_operator
export total_boson_number_subspace_info
export total_boson_number_subspace_tools

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
function find_subspace(state, subspace_indices; digit_error = 15)
    max_norm, index = findmax(x -> norm(@view state[x]), subspace_indices)
    if round(max_norm ; digits = digit_error) == 1.0
        return index, subspace_indices[index]
    else
        throw(ErrorException("Could not find subspace. Max norm was $max_norm, with index $index. digit_error = $digit_error is too large."))
    end
end

function split_operator(op::AbstractMatrix, proj_mat::AbstractMatrix, ranges::Vector{<:UnitRange})
    proj_op = proj_mat * op * proj_mat'
    out = []
    for range in ranges
        push!(out, proj_op[range, range])
    end
    return out
end

function total_boson_number_subspace_info(d::Int, L::Int)
    op = nall(d, L)
    find_property(i) = real(op[i, i])
    return subspace_info(d^L, find_property)
end
function total_boson_number_subspace_tools(d::Int, L::Int)
    op = nall(d, L)
    find_property(i) = real(op[i, i])
    return subspace_tools(d^L, find_property)
end