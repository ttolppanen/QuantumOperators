# using SparseArrays
# using LinearAlgebra
# using ITensors
# using QuantumStates

@testset "measurementtest" begin
    
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
    @test !(real(inner(state, msrresults[1])) ≈ 1.0 || real(inner(state, msrresults[2])) ≈ 1.0)
    measuresite!(state, msrop, 1)
    @test real(inner(state, msrresults[1])) ≈ 1.0 || real(inner(state, msrresults[2])) ≈ 1.0
end

@testset "Operate After Measurement" begin
    d = 3; L = 3
    @testset "Total Space" begin
        state = zeroone(d, L)
        msrop = measurementoperators(nop(d), L)
        proj_op = n_bosons_projector(d, 2)
        unitary_after_measurement!(proj_op, msrop)
        measuresite!(state, msrop, 1)
        n = nall(d, L)
        @test expval(state, n) == 3.0
        state .= normalize!(state + onezero(d, L)) # 210 + 101
        measuresite!(state, msrop, 3) # -> 212 or 102
        bosons = expval(state, n)
        @test bosons == 5 || bosons == 3
    end
    @testset "MPS" begin
        state = zeroonemps(d, L)
        msrop = measurementoperators(nop(d), siteinds(state))
        proj_op = n_bosons_projector(d, 2)
        unitary_after_measurement!(proj_op, msrop)
        measuresite!(state, msrop, 1)
        bosons = sum(expval(state, "N"))
        @test bosons == 3.0
        super_pos = normalize!(zeroone(d, L) + onezero(d, L)) # 010 + 101
        state = MPS(Vector(super_pos), siteinds(state))
        measuresite!(state, msrop, 2) # -> 020 or 121
        bosons = sum(expval(state, "N"))
        @test bosons ≈ 4 || bosons ≈ 2 
    end
    d = 4; L = 4
    @testset "Total Space List Of Unitary" begin
        state = zeroone(d, L) # 0101
        n = nall(d, L)
        msrop = measurementoperators(nop(d), L)
        proj_op1 = n_bosons_projector(d, 2)
        proj_op2 = n_bosons_projector(d, 0)
        v_p = [isodd(i) ? proj_op1 : proj_op2 for i in 1:L]
        m = msrop[1][1]
        unitary_after_measurement!(v_p, msrop)
        measuresite!(state, msrop, 1) # 2101
        measuresite!(state, msrop, 2) # 2001
        n = nall(d, L)
        @test expval(state, n) == 3.0
    end
    @testset "MPS List Of Unitary" begin
        state = zeroonemps(d, L)
        msrop = measurementoperators(nop(d), siteinds(state))
        proj_op1 = n_bosons_projector(d, 2)
        proj_op2 = n_bosons_projector(d, 0)
        v_p = [isodd(i) ? proj_op1 : proj_op2 for i in 1:L]
        unitary_after_measurement!(v_p, msrop)
        measuresite!(state, msrop, 1) # 2010
        measuresite!(state, msrop, 3) # 2010
        @test round.(expval(state, "N"), digits = 15) == [2, 0, 2, 0]
    end
end

end # testset
