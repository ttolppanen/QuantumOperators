using ITensors

export ptrace

#rho : density matrix; for pure states rho = |Ψ><Ψ|
#d : dimension; e.g. with qubits d = 2
#L : number of systems;
#to_sum_over : sites that will be traced out; e.g. to trace out the first and third site: to_sum_over = [1, 3]
#                                             or to trace out second to fourth site: to_sum_over = 2:4


function recursive_ptrace_set_index(d::Integer, i::Integer, j::Integer, out::AbstractMatrix{<:Number}, rho::AbstractMatrix{<:Number}, to_sum_over::AbstractVector{<:Integer}; i_original=i, j_original=j)
    k = to_sum_over[1]
    dk = d^(k-1)
    for k_i in 0:d-1
        if i <= dk
            i_rho = i + k_i * dk
        else
            j_m = (i-1) % dk + 1
            i_rho = (i - j_m) * d + j_m + k_i * dk
        end
        if j <= dk
            j_rho = j + k_i * dk
        else
            j_m = (j-1) % dk + 1
            j_rho = (j - j_m) * d + j_m + k_i * dk
        end
        if length(to_sum_over) > 1
            @views recursive_ptrace_set_index(d, i_rho, j_rho, out, rho, to_sum_over[2:end]; i_original, j_original)
        else
            out[i_original, j_original] += rho[i_rho, j_rho]
        end
    end
end
function ptrace(d::Integer, L::Integer, rho::AbstractMatrix{<:Number}, to_sum_over)
    sum_over = union(collect(to_sum_over))
    sum_over = sort([L - i + 1 for i in sum_over]) #The sum over is reversed, since with the tensor notation first site is the last site
    out = complex(spzeros(d^(L-length(sum_over)), d^(L-length(sum_over))))
    for j in axes(out, 2)
        for i in axes(out, 1)
            recursive_ptrace_set_index(d, i, j, out, rho, sum_over)
        end
    end
    return out
end