using QuantumOperators, Test
using QuantumStates
using SparseArrays
using LinearAlgebra
using ITensors

include("measurementtest.jl")
include("partialtracetest.jl")
include("bosehubbardtest.jl")
include("projectionoperatorstest.jl")
include("bosonmeantest.jl")
include("subspacetest.jl")
include("dissipationdecoherencetest.jl")

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
    @testset "anhall" begin
        testtype(anhall(d, 3)) 
    end
    @testset "n_bosons_projector" begin
        testtype(n_bosons_projector(d, 2)) 
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
    state = productstate(d, [2 for _ in 1:L])
    anh = anhall(d, L)
    @test expval(state, anh) == 2 * L
    anh_23 = anhall(d, L; sites = [2, 3])
    @test expval(state, anh_23) == 2 * 2

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

    @testset "Subspace" begin
        state = allone(d, L)
        dict = total_boson_number_subspace_info(d, L)
        ranges, perm_mat = total_boson_number_subspace_tools(d, L)
        n = nall(d, L)
        max_n = L * (d - 1)
        for i in 0:max_n
            if i == L
                @test expval(state, n, dict[i]) == i
                @test expval(perm_mat * state, perm_mat * n * perm_mat', ranges[i + 1]) == i # works since the dict, and ranges, are sorted
            else
                @test expval(state, n, dict[i]) == 0.0
                @test expval(perm_mat * state, perm_mat * n * perm_mat', ranges[i + 1]) == 0.0 # works since the dict, and ranges, are sorted
            end
        end
    end
end

@testset "Entanglement" begin
    d = 3; L = 2
    @testset "Complete" begin
        state = zeroone(d, L) * 1im
        @test entanglement(d, L, state, 1) == 0.0
        state += onezero(d, L)
        normalize!(state)
        @test entanglement(d, L, state, 1) ≈ log(2)

        L = 4
        state = singleone(d, L, 1)
        @test entanglement(d, L, state, 1) == 0.0
        state += singleone(d, L, 2)
        normalize!(state)
        @test entanglement(d, L, state, 1) ≈ log(2)
        @test entanglement(d, L, [state, state], 1) ≈ [log(2), log(2)]
        @test (entanglement(d, L, state, 2) + 1) ≈ 1.0
        
        rho = state * state'
        @test entanglement(d, L, rho, 1) ≈ log(2)
        @test entanglement(d, L, rho, 2) ≈ log(2)
        @test entanglement(d, L, [rho, rho], 2) ≈ [log(2), log(2)]
        @test (entanglement(d, L, rho, 3) + 1) ≈ 1.0
    end
    @testset "MPS" begin
        mps = zeroonemps(d, L)
        @test entanglement(mps, 1) == 0.0
        mps += onezeromps(siteinds(mps))
        normalize!(mps)
        @test entanglement(mps, 1) ≈ log(2)
        @test entanglement([mps, mps], 1) ≈ [log(2), log(2)]
    end
    @testset "Not in a chain" begin
        d = 3; L = 6
        state = normalize(zeroone(3, L) + allone(d, L) + onezero(d, L))
        @test entanglement(d, L, state, 1:3) ≈ entanglement(d, L, state, 3)
        @test entanglement(d, L, state, [1,5]) ≈ entanglement(d, L, state * state', [1, 5]) # the one on the right does this with partial trace
        @test entanglement(d, L, state, 5:6) ≈ entanglement(d, L, state, 1:4)
    end
end