# using QuantumStates

@testset "Dissipation and Decoherence" begin

@testset "state changed" begin
    d = 2; L = 2;
    state = onezero(d, L)
    operators = [singlesite_a(d, L, 1)]
    apply_diss_deco!(state, operators)
    @test abs(state' * allzero(d, L)) ≈ 1.0
end

@testset "50% outcome" begin
    d = 2; L = 2;
    state = allone(d, L)
    operators = [singlesite_n(d, L, 1), singlesite_n(d, L, 2)]
    outcomes = []
    n = 1000
    for _ in 1:n
        result = apply_diss_deco!(deepcopy(state), operators)
        push!(outcomes, result)
    end
    @test abs(sum(outcomes .- 1) / n - 0.5) < 3 * 0.005 # 3 times standard deviation
end

end # testset