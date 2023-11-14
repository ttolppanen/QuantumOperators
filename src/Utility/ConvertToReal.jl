function real_with_warning(x::Number)
    if !(isapprox(imag(x), 0; atol=eps(Float64) * 10)) #Check if imaginary is not zero
        @warn sprint(showerror, InexactError(:Real, Real, x))
    end
    return real(x)
end