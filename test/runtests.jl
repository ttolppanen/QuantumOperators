using QuantumOperators, Test
using QuantumStates

@test begin
    d = 3; L = 3
    state = zeroone(d, L)
    op = nall(d, L)
    expval(state, op)
end

@testset "Hermitian" begin
    d = 3; L = 4
    @test nop(d) == nop(d)'
    @test bosehubbard(d, L) == bosehubbard(d, L)
end