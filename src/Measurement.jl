# using ITensors
# using KrylovKit
# using SparseArrays
# using LinearAlgebra

# include("CompleteSpaceOperators.jl")

export measurementoperators
export measuresite!
export unitary_after_measurement!
export MsrOpMatrixType
export MsrOpITensorType

# op : operator; a complex matrix representing the operators
# indices : physical indices, see ITensors

MsrOpMatrixType = Vector{Vector{AbstractMatrix{<:Number}}}
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

function measuresite!(state::AbstractVector{<:Number}, msrop::MsrOpMatrixType, siteIndex::Integer)
    msr_result = rand(Float64)
    sum_of_msr = 0
    for proj_op in msrop[siteIndex]
        proj_state = proj_op * state
        sum_of_msr += sum(abs2.(proj_state))
        if msr_result < sum_of_msr
            state .= proj_state
            normalize!(state)
            break
        end
    end
end
function measuresite!(mps::MPS, msrop::MsrOpITensorType, siteIndex::Integer; kwargs...) # kwargs for ITensors.apply
    ITensors.orthogonalize!(mps, siteIndex) # Is this necessary?
    msr_result = rand(Float64)
    sum_of_msr = 0
    for proj_op in msrop[siteIndex]
        proj_mps = apply(proj_op, mps; kwargs...)
        sum_of_msr += real(inner(proj_mps, proj_mps))
        if msr_result < sum_of_msr
            mps .= proj_mps
            normalize!(mps)
            break
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