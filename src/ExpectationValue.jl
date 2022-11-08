export expval

function expval(state, op)
    out = state' * op * state
    return imag(out) == 0.0 ? real(out) : throw(InexactError(:Real, Real, out))
end