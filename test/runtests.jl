using QuantumOperators, Test
using SparseArrays
using LinearAlgebra
using QuantumStates

include("measurementtest.jl")
include("partialtracetest.jl")

@test begin
    d = 3
    a = aop(d)
    ad = adagop(d)
    n = nop(d)
    return n == ad * a
end

@test begin
    d = 3; L = 3
    state = zeroone(d, L)
    op = nall(d, L)
    expval(state, op) == 1.0
end

function testtype(op)
    @test isa(op, AbstractSparseMatrix)
    @test typeof(op) == typeof(complex(op))
end
@testset "Type" begin
    d = 4
    @testset "a" begin
        testtype(aop(d)) 
    end
    @testset "adag" begin
        testtype(adagop(d)) 
    end
    @testset "n" begin
        testtype(nop(d)) 
    end
    @testset "singlesite_a" begin
        testtype(singlesite_a(d, 3, 1)) 
    end
    @testset "nall" begin
        testtype(nall(d, 3)) 
    end
    @testset "bosehubbard" begin
        testtype(bosehubbard(d, 3)) 
    end
end

@testset "Hermitian" begin
    d = 3; L = 4
    @test nop(d) == nop(d)'
    @test singlesite_n(d, L, 2) == singlesite_n(d, L, 2)'
    @test bosehubbard(d, L) == bosehubbard(d, L)'
end

@testset "ExpectationValue" begin
    d = 3; L = 4
    state = allone(d, L)
    ntot = nall(d, L)
    @test expval(state, ntot) == L
    n = singlesite_n(d, L, 1)
    @test expval(state, n) == 1
    a = singlesite_a(d, L, 1)
    state = a * state
    @test expval(state, n) == 0
    ad = singlesite_adag(d, L, 2)
    state = ad * state
    normalize!(state)
    n = singlesite_n(d, L, 2)
    @test expval(state, n) == 2
    @test expval([state, state], n) == [2, 2]
    @testset "Trajectories" begin
        traj1 = [state, state]
        @test trajmean([traj1, traj1], n) == [2, 2]
        traj2 = [zeroone(d, L), zeroone(d, L)]
        @test trajmean([traj1, traj2, traj1, traj2], n) == [1.5, 1.5] 
    end
    
    @testset "MPS" begin
        state = allonemps(d, L)
        @test expval([state, state], "N"; sites=1) == [1, 1]
        @test expval([state, state], "N"; sites=2:3) == [[1, 1], [1, 1]]
        result = [sum(res) for res in expval([state, state], "N")]
        @test result == [L, L]
    end
end

@testset "Entanglement" begin
    d = 3; L = 2
    mps = zeroonemps(d, L)
    @test entanglement(mps, 1) == 0.0
    mps += onezeromps(siteinds(mps))
    normalize!(mps)
    @test entanglement(mps, 1) ≈ log(2)
end