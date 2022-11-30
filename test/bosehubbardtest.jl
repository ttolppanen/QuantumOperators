# using QuantumStates
# using Random

@testset "Bosehubbard" begin
    
d = 2; L = 2
zo = zeroone(d, L)
H = bosehubbard(d, L)
@test expval(zo, H) == 1
H = bosehubbard(d, L; w=3)
@test expval(zo, H) == 3
state = bosonstack(3, L, 1)
H = bosehubbard(4, L; U=-3)
@test expval(state, H) == 12
wUJ = 10 * rand(3)
@test bosehubbard(d, L; w=wUJ[1], U=wUJ[2], J=wUJ[3]) == QuantumOperators.bosehubbard_old(d, L; w=wUJ[1], U=wUJ[2], J=wUJ[3])

end # testset
