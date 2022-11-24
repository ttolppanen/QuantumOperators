#using LinearAlgebra
#using QuantumStates

@testset "Partial Trace" begin

function testtrace(d, L, all_rho...)
    @testset "Trace" begin
        for rho in all_rho
            for i in 1:L
                rho_p = ptrace(d, L, rho, i)
                @test tr(rho_p) ≈ 1.0
                rho_p = ptrace(d, L, rho, i:L)
                @test tr(rho_p) ≈ 1.0
            end
            rho_p = ptrace(d, L, rho, 1:2:L)
            @test tr(rho_p) ≈ 1.0
        end 
    end
end

function testpurity(d, L, pure_rho, pure_ent_rho, mixed_ent_rho)
    @testset "Purity" begin
        @testset "Pure" begin
            for i in 1:L
                rho_p = ptrace(d, L, pure_rho, i)
                @test tr(rho_p * rho_p) ≈ 1.0
                rho_p = ptrace(d, L, pure_rho, i:L)
                @test tr(rho_p * rho_p) ≈ 1.0
            end
            rho_p = ptrace(d, L, pure_rho, 1:2:L)
            @test tr(rho_p) ≈ 1.0
        end
        @testset "Entangled" begin
            for i in 1:L
                rho_p = ptrace(d, L, pure_ent_rho, i)
                @test !(tr(rho_p * rho_p) ≈ 1.0)
                rho_p = ptrace(d, L, pure_ent_rho, i:L)
                if i == 1
                    @test tr(rho_p * rho_p) ≈ 1.0
                else
                    @test !(tr(rho_p * rho_p) ≈ 1.0)
                end
            end
            rho_p = ptrace(d, L, pure_ent_rho, 1:2:L)
            @test !(tr(rho_p * rho_p) ≈ 1.0)
            rho_p = ptrace(d, L, pure_ent_rho, 1:Int(floor(L/2)))
            @test tr(rho_p * rho_p) ≈ 0.5
        end
        @testset "Mixed Entangled" begin
            for i in 1:L
                rho_p = ptrace(d, L, mixed_ent_rho, i)
                @test !(tr(rho_p * rho_p) ≈ 1.0)
                rho_p = ptrace(d, L, mixed_ent_rho, i:L)
                if i == 1
                    @test tr(rho_p * rho_p) ≈ 1.0
                else
                    @test !(tr(rho_p * rho_p) ≈ 1.0)
                end
            end
            rho_p = ptrace(d, L, mixed_ent_rho, 1:2:L)
            @test !(tr(rho_p * rho_p) ≈ 1.0)
        end
    end    
end
function teststate(d, L)
    @testset "States" begin
        half = Int(floor(L/2))
        state = onezero(d, L)
        rho = state * state'
        rho_p = ptrace(d, L, rho, 1:half)
        stateshouldbe(rho_p, onezero(d, L - half))
        
        even_sites = 2:2:(half * 2)
        rho_p = ptrace(d, L, rho, even_sites)
        stateshouldbe(rho_p, allone(d, L - length(even_sites)))
        
        state = singleone(d, L , 1) + singleone(d, L, 2)
        normalize!(state)
        rho = state * state'
        rho_p = ptrace(d, L, rho, 3:L)
        should_be = normalize(zeroone(d, 2) + onezero(d, 2))
        stateshouldbe(rho_p, should_be)
        rho_p = ptrace(d, L, rho, (1,2))
        stateshouldbe(rho_p, allzero(d, L-2))
    end
end
function stateshouldbe(rho_p, state)
    @test rho_p ≈ state * state'
end

d = 3; L = 4
state = onezero(d, L)
pure_rho = state * state'
state = normalize(onezero(d, L) + 1im * zeroone(d, L))
pure_ent_rho = state * state'
state = allone(d, L)
mixed_ent_rho = 0.5 * state * state' + 0.5 * pure_ent_rho

testtrace(d, L, pure_rho, pure_ent_rho, mixed_ent_rho)
testpurity(d, L, pure_rho, pure_ent_rho, mixed_ent_rho)
teststate(d, L)

end # testset