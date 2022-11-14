using SparseArrays
using LinearAlgebra
using ITensors
using QuantumOperators.QuantumStates

@testset "Measurement Operator Generation" begin
    d = 2; L = 2
    n = nop(d)
    msrop = measurementoperators(n, L)
    @test msrop[1][1] ≈ kron([0 0; 0 1], Matrix(I, 2, 2))
    @test msrop[2][2] ≈ kron(Matrix(I, 2, 2), [1 0; 0 0])
    indices = siteinds("Boson", L; dim = d)
    msrop = measurementoperators(n, indices)
    @test Matrix(msrop[1][1], indices[1]', indices[1]) ≈ [0 0; 0 1]
end

@testset "Measuresite!" begin
    d = 2; L = 2
    state = onezero(d, L) + zeroone(d, L)
    normalize!(state)
    n = nop(d)
    msrop = measurementoperators(n, L)
    msrresults = [onezero(d, L), zeroone(d, L)]
    @test !(real(state' * msrresults[1]) ≈ 1.0 || real(state' * msrresults[2]) ≈ 1.0)
    measuresite!(state, msrop, 1)
    @test real(state' * msrresults[1]) ≈ 1.0 || real(state' * msrresults[2]) ≈ 1.0

    indices = siteinds("Boson", L, dim = d)
    state = onezeromps(indices) + zeroonemps(indices)
    normalize!(state)
    n = nop(d)
    msrop = measurementoperators(n, indices)
    msrresults = [onezeromps(indices), zeroonemps(indices)]
    @test !(real(inner(state', msrresults[1])) ≈ 1.0 || real(inner(state', msrresults[2])) ≈ 1.0)
    measuresite!(state, msrop, 1)
    @test real(inner(state', msrresults[1])) ≈ 1.0 || real(inner(state', msrresults[2])) ≈ 1.0
end