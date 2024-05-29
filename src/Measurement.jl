# using ITensors
# using KrylovKit
# using SparseArrays
# using LinearAlgebra

# include("CompleteSpaceOperators.jl")

export measurementoperators
export calc_msr_probability
export measuresite!
export MsrOpMatrixType
export MsrOpITensorType
export unitary_after_measurement!

# op : operator; a complex matrix representing the operators
# indices : physical indices, see ITensors

MsrOpMatrixType = Vector{Vector{SparseMatrixCSC{ComplexF64, Int64}}}
MsrOpITensorType = Vector{Vector{ITensor}}

function measurementoperators(op::AbstractMatrix{<:Number}, L::Integer)::MsrOpMatrixType
    vals, eigenvectors = eigsolve(op) # projection operators are formed from the eigenvectors of the matrix 
    out = []
    for i in 1:L
        push!(out, [])
        for v in eigenvectors # loop over the width of the matrix
            push!(out[i], singlesite(sparse(v * v'), L, i))
        end
    end
    return out # structure is out[site][measurement result]
end
function measurementoperators(op::AbstractMatrix{<:Number}, indices::Vector{Index{Int64}})::MsrOpITensorType
    vals, eigenvectors = eigsolve(op) # projection operators are formed from the eigenvectors of the matrix 
    out = []
    for i in eachindex(indices)
        push!(out, [])
        ind = indices[i]
        for v in eigenvectors # loop over the width of the matrix
            vtensor = ITensor(v * v', ind', ind)
            push!(out[i], vtensor)
        end
    end
    return out # structure is out[site][measurement result]
end

# This is assuming that proj_op is unitary?
function calc_msr_probability(proj_op::AbstractMatrix{<:Number}, state::AbstractVector{<:Number})
    out = 0.0
    for j in axes(proj_op, 2)
        for r in nzrange(proj_op, j)
            i = rowvals(proj_op)[r]
            out += abs2(nonzeros(proj_op)[r]' * state[i])
        end
    end
    return out
end

function measuresite!(state::AbstractVector{<:Number}, msrop::MsrOpMatrixType, siteIndex::Integer)::Int64
    msr_result = rand(Float64)
    sum_of_msr = 0.0
    msr_results = length(msrop[siteIndex])
    for i in 1:(msr_results - 1)
        proj_op = msrop[siteIndex][i]
        sum_of_msr += calc_msr_probability(proj_op, state)
        if msr_result < sum_of_msr
            state .= proj_op * state
            normalize!(state)
            return i
        end
    end
    proj_op = msrop[siteIndex][msr_results] # the last msr happened
    state .= proj_op * state
    normalize!(state)
    return msr_results
end
function measuresite!(mps::MPS, msrop::MsrOpITensorType, siteIndex::Integer; kwargs...) # kwargs for ITensors.apply
    ITensors.orthogonalize!(mps, siteIndex) # Is this necessary?
    msr_result = rand(Float64)
    sum_of_msr = 0.0
    for i in eachindex(msrop[siteIndex])
        proj_op = msrop[siteIndex][i]
        proj_mps = apply(proj_op, mps; kwargs...)
        sum_of_msr += real(inner(proj_mps, proj_mps))
        if msr_result < sum_of_msr
            mps .= proj_mps
            normalize!(mps)
            return i
        end
    end
end


function unitary_after_measurement!(U::AbstractMatrix{<:Complex}, msrop::T) where {T<:MsrOpMatrixType}
    L = length(msrop)
    for site in 1:L
        U_site = singlesite(U, L, site)
        for i in eachindex(msrop[site])
            msrop[site][i] .= U_site * msrop[site][i]
        end
    end
    if !isa(msrop, T) throw(ArgumentError("Type of U changed the type of msrop. Try Matrix(U)")) end
end
function unitary_after_measurement!(U::Vector{<:AbstractMatrix{<:Complex}}, msrop::T) where {T<:MsrOpMatrixType}
    L = length(msrop)
    for site in 1:L
        U_site = singlesite(U[site], L, site)
        for i in eachindex(msrop[site])
            msrop[site][i] .= U_site * msrop[site][i]
        end
    end
    if !isa(msrop, T) throw(ArgumentError("Type of U changed the type of msrop. Try Matrix(U)")) end
end
function unitary_after_measurement!(U::AbstractMatrix{<:Complex}, msrop::T) where {T<:MsrOpITensorType}
    L = length(msrop)
    for site in 1:L
        for i in eachindex(msrop[site])
            msr_res = msrop[site][i]
            mat_inds = inds(msr_res)
            msr_res_mat = U * Matrix(msr_res, mat_inds[1], mat_inds[2])
            msrop[site][i] = ITensor(msr_res_mat, mat_inds[1], mat_inds[2])
        end
    end
    if !isa(msrop, T) throw(ArgumentError("Type of U changed the type of msrop. Try Matrix(U)")) end
end
function unitary_after_measurement!(U::Vector{<:AbstractMatrix{<:Complex}}, msrop::T) where {T<:MsrOpITensorType}
    L = length(msrop)
    for site in 1:L
        for i in eachindex(msrop[site])
            msr_res = msrop[site][i]
            mat_inds = inds(msr_res)
            msr_res_mat = U[site] * Matrix(msr_res, mat_inds[1], mat_inds[2])
            msrop[site][i] = ITensor(msr_res_mat, mat_inds[1], mat_inds[2])
        end
    end
    if !isa(msrop, T) throw(ArgumentError("Type of U changed the type of msrop. Try Matrix(U)")) end
end

# alternatively use:
# opname = "whatever"
# mat = [12312313 123 2 1323]
# eval(Meta.parse("ITensors.op(::OpName\"$opname\", ::SiteType\"Boson\", d::Int) = mat"))