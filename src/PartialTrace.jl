using ITensors

export ptrace

#rho : density matrix; for pure states rho = |Ψ><Ψ|
#d : dimension; e.g. with qubits d = 2
#L : number of systems;
#cut : cut between the left and right bipartion of the system; cut = 2 cuts the system in to A = {1, 2} and B {3,..., L}
function recursive_ptrace_set_index(d::Integer, i::Integer, j::Integer, out, rho, to_sum_over; i_original=i, j_original=j)
    k = to_sum_over[1]
    for k_i in 0:d-1
        if i <= d^(k-1)
            i_rho = i + k_i * d^(k-1)
        else
            j_m = (i-1) % d^(k-1) + 1
            i_rho = (i - j_m) * d + j_m + k_i * d^(k-1)
        end
        if j <= d^(k-1)
            j_rho = j + k_i * d^(k-1)
        else
            j_m = (j-1) % d^(k-1) + 1
            j_rho = (j - j_m) * d + j_m + k_i * d^(k-1)
        end
        if length(to_sum_over) > 1
            recursive_ptrace_set_index(d, i_rho, j_rho, out, rho, to_sum_over[2:end]; i_original, j_original)
        else
            out[i_original, j_original] += rho[i_rho, j_rho]
        end
    end
end
function ptrace(d::Integer, L::Integer, rho::AbstractMatrix{<:Number}, to_sum_over)
    sum_over = sort(union(to_sum_over))
    out = complex(spzeros(d^(L-length(sum_over)), d^(L-length(sum_over))))
    for j in axes(out, 2)
        for i in axes(out, 1)
            recursive_ptrace_set_index(d, i, j, out, rho, sum_over)
        end
    end
    return out
end
function ptrace(d::Integer, L::Integer, state::AbstractVector{<:Number}, to_sum_over)
   return ptrace(d, L, state * state', to_sum_over) 
end
#=
if j <= d^(k-1)
    j_rho = j + k_i * d^(k-1)
else
    j_m = (j-1) % d^(k-1) + 1
    j_rho = (j - j_m) * d + j_m + k_i * d^(k-1)
end
if i <= d^(k-1)
    i_rho = i + k_i * d^(k-1)
else
    j_m = (i-1) % d^(k-1) + 1
    i_rho = (i - j_m) * d + j_m + k_i * d^(k-1)
end
=#