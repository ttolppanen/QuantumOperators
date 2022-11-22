using ITensors

export ptrace

#rho : density matrix; for pure states rho = |Ψ><Ψ|
#d : dimension; e.g. with qubits d = 2
#L : number of systems;
#cut : cut between the left and right bipartion of the system; cut = 2 cuts the system in to A = {1, 2} and B {3,..., L}

function ptrace(d::Integer, L::Integer, rho::AbstractMatrix{<:Number}, k::Integer)
    out = complex(spzeros(d^(L-1), d^(L-1)))
    for j in axes(out, 2)
        for i in axes(out, 1)
            for k_i in 0:d-1
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
                @show i
                @show j
                @show i_rho
                @show j_rho
                out[i, j] += rho[i_rho, j_rho]
            end
        end
    end
    return out
end
function ptrace(d::Integer, L::Integer, state::AbstractVector{<:Number}, cut::Integer)
   return ptrace(d, L, state * state', cut) 
end