# using LinearAlgebra

@testset "subspace" begin

@testset "FindSubspace" begin
    d = 2
    L = 2
    finder = generate_total_boson_number_subspace_finder(d, L)
    @test finder(normalize!([1, 0, 0, 0])) == [1]
    @test finder(normalize!([0, 1, 0, 0])) == [2, 3]
    @test finder(normalize!([0, 0, 1, 0])) == [2, 3]
    @test finder(normalize!([0, 0, 0, 1])) == [4]
    @test finder(normalize!([0, 1.0, 1.0, 0])) == [2, 3]

    @test_throws ErrorException finder(normalize!([1.0, 0, 0, 1.0])) == [2, 3]

    bn_ss_indeces = total_boson_number_subspace(d, L)
    @test find_subspace(normalize!([0, 1, 0, 0]), bn_ss_indeces) == [2, 3]
end

end # testset