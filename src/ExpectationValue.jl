using ITensors

export expval

#state : quantum state; either a complex vector or a MPS representing the quantum state
#op : operator; a complex matrix representing the operator, for MPS you can also pass the operator as a string (see ITensors for supported operators as strings)

function expval(state::AbstractVector{<:Number}, op::AbstractMatrix{<:Number})
    out = state' * op * state
    return isapprox(imag(out), 0; atol=eps(Float64)) ? real(out) : throw(InexactError(:Real, Real, out)) #Check if imaginary is zero
end

function expval(states::AbstractVector{<:AbstractVector{<:Number}}, op::AbstractMatrix{<:Number})
   return [expval(state, op) for state in states]
end
function expval(states::AbstractVector{MPS}, op::Union{Matrix{<:Number}, String}; kwargs...) #keyword arguments for ITensors.expect
    return [expect(state, op; kwargs...) for state in states] #kwargs can be {sites} 
 end