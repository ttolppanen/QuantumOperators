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
    @test abs(sum(outcomes .- 1) / n - 0.5) < 3 * 0.5 / sqrt(n) # 3 times standard deviation
end

@testset "Dissipation in Subspace" begin
    d = 2; L = 2;
    state = onezero(d, L)
    
    indices = total_boson_number_subspace_indices(d, L)
    ranges, perm_mat = total_boson_number_subspace_tools(d, L)
    id = find_subspace(state, indices)

    state = subspace_split(state, indices)
    a_relations = [-1 for _ in eachindex(ranges)]
    a_relations[1] = 0
    diss_op = [subspace_split(singlesite_a(d, L, 1), ranges, perm_mat, a_relations)]
    deco_op = [subspace_split(singlesite_n(d, L, 1), indices)]
    new_id, result = apply_diss_deco!(state, id, diss_op, deco_op)
    if result == 1
        @test new_id == id - 1
        @test state[new_id] == [1.0]
    elseif result == 2
        @test new_id == id
    else
        @test false
    end
end

end # testset