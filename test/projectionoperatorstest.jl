# using QuantumStates
# using LinearAlgebra

@testset "ProjectionOperators" begin
    
d = 3; L = 5
state = zeroone(d, L)
proj_op = n_bosons_projector(d, 1)
projector = singlesite(proj_op, L, 2)
@test projector * state == state
projector = singlesite(proj_op, L, 1)
@test state' * projector * state == 0.0
state = zeroone(d, L) + onezero(d, L)
normalize!(state)
proj_op = n_bosons_projector(d, 2)
projector = singlesite(proj_op, L, 1)
state .= projector * state
@test expval(state, nall(d, L)) ≈ 4

end # testset