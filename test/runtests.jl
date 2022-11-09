using QuantumOperators, Test
using SparseArrays
using LinearAlgebra
using QuantumOperators.QuantumStates

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
    @show expval(state, op)
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
end