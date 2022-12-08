# using ITensors
# using SparseArrays

export n_bosons_projector

# d : dimension; e.g. with qubits d = 2
# n : number of bosons; how many bosons after projection

function n_bosons_projector(d::Integer, n::Integer)
    if n < 0 throw(ArgumentError("n < 0")) end
    if n >= d throw(ArgumentError("n >= d")) end
    out = spzeros(d, d)
    out[n + 1, :] .= 1
    return complex(out)
end