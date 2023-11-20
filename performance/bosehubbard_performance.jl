using QuantumOperators
using SparseArrays

function bosehubbard_old_old(d::Integer, L::Integer; w=1, U=1, J=1) # Order might be reversed here? The last site could be the first?
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

function f()
    d = 2; L = 22;
    
    bosehubbard(d, L)
    @time bosehubbard(d, L)


    # @time bosehubbard_old_old(d, L)
end

f();
