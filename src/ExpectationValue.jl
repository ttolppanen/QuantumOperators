export expval

function expval(state::AbstractVector{<:Number}, op::AbstractMatrix{<:Number})
    out = state' * op * state
    return imag(out) == 0.0 ? real(out) : throw(InexactError(:Real, Real, out))
end

function expval(states::AbstractVector{<:AbstractVector{<:Number}}, op::AbstractMatrix{<:Number})
   return [expval(state, op) for state in states]
end