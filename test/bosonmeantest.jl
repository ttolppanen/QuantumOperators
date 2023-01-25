# using QuantumStates
# using LinearAlgebra

@testset "Boson Mean" begin
    
d = 3; L = 4
state = allone(d, L)
ntot = nall(d, L)
@test expval(state, ntot) == bosonmean(d, L, state, 1:L)
n = singlesite_n(d, L, 1)
@test expval(state, n) == bosonmean(d, L, state, [1])
a = singlesite_a(d, L, 1)
state = a * state
@test expval(state, n) == bosonmean(d, L, state, [1])
state = a' * a' * a' * state
@test expval(state, n) == bosonmean(d, L, state, [1])
state = normalize(allone(d, L) + zeroone(d, L))
@test bosonmean(d, L, state, [1]) ≈ 0.5
@test bosonmean(d, L, state, [2]) ≈ 1.0
@test bosonmean(d, L, state, [3, 4]) ≈ 1.5
state = normalize(allone(d, L) + zeroone(d, L) + 1im * onezero(d, L))
@test bosonmean(d, L , state, [1]) ≈ 2.0 / 3
state = normalize(a' * state)
@test expval(state, n) ≈ bosonmean(d, L , state, [1])

end # testset