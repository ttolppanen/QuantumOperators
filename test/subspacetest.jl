# using LinearAlgebra

@testset "Subspace" begin

@testset "FindSubspace" begin
    d = 2
    L = 2
    dict = total_boson_number_subspace_info(d, L)
    perm_mat, ranges = total_boson_number_subspace_tools(d, L)
    finder(state) = find_subspace(state, collect(values(dict)))
    finder_range(state) = find_subspace(state, ranges)
    test_states = Dict(
        normalize!([1, 0, 0, 0]) => 1,
        normalize!([0, 1, 0, 0]) => 2,
        normalize!([0, 0, 1, 0]) => 2,
        normalize!([0, 0, 0, 1]) => 3,
        normalize!([0, 1.0, 1.0, 0]) => 2
    )
    for (state, answer) in test_states
        @test finder(state) == answer
        @test finder_range(perm_mat * state) == answer
    end

    @test_throws ErrorException finder(normalize!([1.0, 0, 0, 1.0]))

    finder_id(state) = find_subspace(state, ranges; id_initial_guess = 3, iterate_order = -1)
    for (state, answer) in test_states
        @test finder_id(state) == answer
    end
end

@testset "OP splitting" begin
    d = 2; L = 2;
    indices = total_boson_number_subspace_indices(d, L)

    state = onezero(d, L)
    split_s = subspace_split(state, indices)
    @test split_s[1] == [0.0]
    @test split_s[2] == [0.0, 1.0]
    @test split_s[3] == [0.0] 

    op = bosehubbard(d, L)
    split_H = subspace_split(op, indices)
    @test split_H[1][1, 1] == 0
    @test split_H[3][1, 1] == 2.0
    @test split_H[2][1:2, 1:2] == [1.0 1.0; 1.0 1.0]
end

@testset "Measurement in subspace" begin
    d = 2; L = 2;
    msrop = measurementoperators(nop(d), L) # here msrop[L][1] -> measurement outcome 1 boson, here msrop[L][2] -> 0 boson
    feedback = [singlesite(n_bosons_projector(d, 0), L, i) for i in 1:L] # project to |0>
    indices = total_boson_number_subspace_indices(d, L)
    feedback, msrop = feedback_measurement_subspace(feedback, msrop, indices)
    
    #feedback[subspace_id][site][msr_outcome], and looking at the index
    @test (feedback[3][1][1])[1] == 2
    @test (feedback[2][1][1])[1] == 1
    @test (feedback[1][1][1])[1] == 1

    @test (feedback[3][1][2])[1] == 3 # this shouldn't be possible, but for now is intented.
    @test (feedback[2][1][2])[1] == 2
end

end # testset