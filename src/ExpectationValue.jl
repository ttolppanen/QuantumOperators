using ITensors

export expval

function expval(state::AbstractVector{<:Number}, op::AbstractMatrix{<:Number})
    out = state' * op * state
    return isapprox(imag(out), 0; atol=eps(Float64)) ? real(out) : throw(InexactError(:Real, Real, out)) #Check if imaginary is zero
end

function expval(states::AbstractVector{<:AbstractVector{<:Number}}, op::AbstractMatrix{<:Number})
   return [expval(state, op) for state in states]
end
function expval(states::AbstractVector{MPS}, op::String; kwargs...) #keyword arguments for ITensors.expect
    return [expect(state, op; kwargs...) for state in states] #kwargs can be {sites} 
 end
 function expval(states::AbstractVector{MPS}, op::Matrix{<:Number}; kwargs...) #keyword arguments for ITensors.expect
    return [expect(state, op; kwargs...) for state in states] #kwargs can be {sites} 
 end