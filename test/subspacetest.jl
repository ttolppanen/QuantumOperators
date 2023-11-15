# using LinearAlgebra

@testset "subspace" begin

@testset "FindSubspace" begin
    d = 2
    L = 2
    dict = total_boson_number_subspace_info(d, L)
    perm_mat, ranges = total_boson_number_subspace_tools(d, L)
    finder(state) = find_subspace(state, collect(values(dict)))
    finder_range(state) = find_subspace(state, ranges)
    test_states = Dict(
        normalize!([1, 0, 0, 0]) => ([1], 1:1),
        normalize!([0, 1, 0, 0]) => ([2, 3], 2:3),
        normalize!([0, 0, 1, 0]) => ([2, 3], 2:3),
        normalize!([0, 0, 0, 1]) => ([4], 4:4),
        normalize!([0, 1.0, 1.0, 0]) => ([2, 3], 2:3)
    )
    for (state, answer) in test_states
        @test finder(state)[2] == answer[1]
        @test finder_range(perm_mat * state)[2] == answer[2]
    end

    @test_throws ErrorException finder(normalize!([1.0, 0, 0, 1.0]))
end

@testset "OP splitting" begin
    d = 2; L = 2;
    op = bosehubbard(d, L)
    perm_mat, ranges = total_boson_number_subspace_tools(d, L)
    split_H = split_operator(op, perm_mat, ranges)
    @test split_H[1][1, 1] == 0
    @test split_H[3][1, 1] == 2.0
    @test split_H[2][1:2, 1:2] == [1.0 1.0; 1.0 1.0]
end

end # testset