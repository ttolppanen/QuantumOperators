# using ITensors
# include("Utility/ConvertToReal.jl")

export expval
export trajmean
export bosonmean

# state : quantum state; either a complex vector or a MPS representing the quantum state
# op : operator; a complex matrix representing the operator, for MPS you can also pass the operator as a string (see ITensors for supported operators as strings)
# f : function; takes the time-evolution as an argument and calculates a quantity e.g. f(states) = expval(states, nall(d, L))
# subspace_indeces : indeces for a generic subspace; A list of integers, e.g. [2, 3] for two qubits and the 1 total boson subspace.

function expval(state::AbstractVector{<:Number}, op::AbstractMatrix{<:Number})
    out = dot(state, op, state)
    return real_with_warning(out)
end
function expval(state::AbstractVector{<:Number}, op::AbstractMatrix{<:Number}, subspace_indeces::AbstractVector{<:Integer})
    @views out = dot(state[subspace_indeces], op[subspace_indeces, subspace_indeces], state[subspace_indeces])
    return real_with_warning(out)
end
function expval(state::MPS, op::Union{Matrix{<:Number}, String}; kwargs...) # keyword arguments for ITensors.expect
    expect(state, op; kwargs...) # kwargs can be {sites} 
end

# expval for time-evolution
function expval(states, op; kwargs...) # kwargs for expval (ITensor.expect)
    return [expval(state, op; kwargs...) for state in states]
end

function trajmean(traj, f::Function)
    traj_num = length(traj)
    out = f.(traj[1])
    for i in 2:traj_num
        out .+= f.(traj[i])
    end
    return out ./ traj_num
end
function trajmean(traj, op; kwargs...) # kwargs for expval (ITensor.expect)
    f(state) = expval(state, op; kwargs...)
    return trajmean(traj, f)
end

function bosonmean(d::Integer, L::Integer, state::AbstractVector{<:Number}, sites::AbstractArray{<:Integer})
    out = 0;
    sites_reversed = [L - (site - 1) for site in sites] # withouth reversing 1 would refer to the last site etc..
    for (i, value) in enumerate(state)
        n_total = 0
        for site in sites_reversed
            n = floor((i-1) / d^(site-1)) % d
            n_total += n
        end
        out += n_total * abs2(value)
    end
    return out
end