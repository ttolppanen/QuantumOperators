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
d = 4; L = 5
H = bosehubbard(d, L; w=wUJ[1], U=wUJ[2], J=wUJ[3])
H_correct = QuantumOperators.bosehubbard_old(d, L; w=wUJ[1], U=wUJ[2], J=wUJ[3])
@test norm(H - H_correct) + 1.0 ≈ 1.0

w = [1, 0, 1, 0, 1]
H = bosehubbard(d, L; w, U = 0, J = 0)
@test expval(onezero(d, L), H) == 3
@test expval(zeroone(d, L), H) == 0

d = 3; L = 3
H = bosehubbard(d, L; w = [0, 1, 2], 
                      U = [0, 3, 6],
                      J = [0, 1])

@test expval(onezero(d, L), H) == 2
@test expval(zeroone(d, L), H) == 1
@test expval(productstate(d, [0, 0, 2]), H) == 2 * 2 - 6
@test productstate(d, [0, 1, 1])' * H * productstate(d, [0, 0, 2]) == sqrt(2)

d = 3; L = 4
Ws = rand(4)
Us = rand(4)
Js = rand(3)
@test norm(bosehubbard(d, L; w=Ws, U=Us, J=Js) - bosehubbard(d, L, [(1, 2, Js[3]), (2, 3, Js[2]), (3, 4, Js[1])]; w=Ws, U=Us)) ≈ 0.0
end # testset
