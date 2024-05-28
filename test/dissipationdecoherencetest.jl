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
    state = normalize!(allzero(d, L) + onezero(d, L))
    operators = [singlesite_n(d, L, 1)]
    outcomes = []
    n = 100
    for _ in 1:n
        result = apply_diss_deco!(deepcopy(state), operators)
        push!(outcomes, result)
    end
    @test abs(sum(outcomes) / n - 0.5) < 2 * 0.05 # 2 times standard deviation
end

end # testset