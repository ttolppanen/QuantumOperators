using ITensors

export expval
export trajmean

#state : quantum state; either a complex vector or a MPS representing the quantum state
#op : operator; a complex matrix representing the operator, for MPS you can also pass the operator as a string (see ITensors for supported operators as strings)
#f : function; takes the time-evolution as an argument and calculates a quantity e.g. f(states) = expval(states, nall(d, L))

function expval(state::AbstractVector{<:Number}, op::AbstractMatrix{<:Number})
    out = state' * op * state
    if !(isapprox(imag(out), 0; atol=eps(Float64))) #Check if imaginary is not zero
        @warn sprint(showerror, InexactError(:Real, Real, out))
    end
    return real(out)
end
function expval(state::MPS, op::Union{Matrix{<:Number}, String}; kwargs...) #keyword arguments for ITensors.expect
    expect(state, op; kwargs...) #kwargs can be {sites} 
end

#expval for time-evolution
function expval(states, op; kwargs...) #kwargs for expval (ITensor.expect)
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
function trajmean(traj, op; kwargs...) #kwargs for expval (ITensor.expect)
    f(state) = expval(state, op; kwargs...)
    return trajmean(traj, f)
end