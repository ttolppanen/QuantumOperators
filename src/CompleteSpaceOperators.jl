# using SparseArrays
# using LinearAlgebra

export aop
export adagop
export nop
export singlesite
export singlesite_n
export singlesite_a
export singlesite_adag
export nall
export bosehubbard

# d : dimension; e.g. with qubits d = 2
# L : number of systems;

⊗(a, b) = kron(a, b)

function aop(d::Integer)
    out = spzeros(d, d)
    for i in 1:(d - 1)
        out[i, i + 1] = sqrt(i)
    end
    return complex(out)
end

adagop(d::Integer) = sparse(aop(d)')

function nop(d::Integer)
    out = spzeros(d, d)
    for i in 1:d
        out[i, i] = i - 1
    end
    return complex(out)
end

function nall(d::Integer, L::Integer)
    n = nop(d)
    out = singlesite(n, L, 1)
    for i in 2:L
        out += singlesite(n, L, i)
    end
    return out
end

function singlesite(op::AbstractMatrix{<:Number}, L::Integer, target::Integer)
    d = size(op)[1]
    return Matrix(I, d^(target-1), d^(target-1)) ⊗ op ⊗ Matrix(I, d^(L-target), d^(L-target))
end

function singlesite_n(d::Integer, L::Integer, target::Integer)
    n = nop(d)
    return singlesite(n, L, target)
end

function singlesite_a(d::Integer, L::Integer, target::Integer)
    a = aop(d)
    return singlesite(a, L, target)
end

function singlesite_adag(d::Integer, L::Integer, target::Integer)
    ad = adagop(d)
    return sparse(singlesite(ad, L, target))
end

function bosehubbard(d::Integer, L::Integer; w=1, U=1, J=1) # Order might be reversed here? The last site could be the first?
    H = spzeros(d^L, d^L)
    for i in 1:d^L
        for site in 1:L
            n = floor((i-1) / d^(site-1)) % d
            n_next = floor((i-1) / d^(site)) % d
            if site != L && n > 0 && n_next + 1 < d
                i_j = i - d^(site - 1) + d^site
                H[i, i_j] = J * sqrt(n * (n_next + 1))
                H[i_j, i] = H[i, i_j]'
            end
            H[i, i] += w * n - U / 2 * n * (n - 1)
        end
    end
    return complex(H)
end

#Not used but needed in tests

function bosehubbard_old(d::Integer, L::Integer; w=1, U=1, J=1) # This should definetly be correct
    a = aop(d)
    n = nop(d)
    H = spzeros(d^L, d^L)
    for i in 1:L
        ñ = singlesite(n, L, i)
        H .+= w .* ñ
        H .+= -U/2 .* ñ * (ñ - I)
        if(i != L)
            hopping = Matrix(I, d^(i-1), d^(i-1)) ⊗ a ⊗ a' ⊗ Matrix(I, d^(L-i - 1), d^(L-i - 1))
            H .+= J .* (hopping + hopping')
        end
    end
    return complex(H)
end