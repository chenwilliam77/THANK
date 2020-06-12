include("../../../src/includeall.jl")
using ModelConstructors, DSGE, Test, HDF5, MAT

m = JPT()
Γ0, Γ1, C, Ψ, Π = eqcond(m)
sys = compute_system(m)

Γ0m = h5read("eqcond.h5", "GAM0")
Γ1m = h5read("eqcond.h5", "GAM1")
Cm  = h5read("eqcond.h5", "C")
Ψm  = h5read("eqcond.h5", "PSI")
Πm  = h5read("eqcond.h5", "PPI")
TTTm = h5read("solve.h5", "TTT")
RRRm = h5read("solve.h5", "RRR")
ZZm  = h5read("solve.h5", "ZZ")
DDm  = h5read("solve.h5", "DD")
QQm  = h5read("solve.h5", "QQ")
steady_state = h5read("solve.h5", "steadystate")

@testset "Steady state" begin
    @test @test_matrix_approx_eq steady_state map(x -> x.value, m.steady_state)
end

@testset "JPT equilibrium conditions" begin
    @test @test_matrix_approx_eq Γ0 Γ0m
    @test @test_matrix_approx_eq Γ1 Γ1m
    @test @test_matrix_approx_eq C  Cm
    @test @test_matrix_approx_eq Ψ  Ψm
    @test @test_matrix_approx_eq Π  Πm
end

@testset "JPT reduced form transition matrices" begin
    @test @test_matrix_approx_eq sys[:TTT] TTTm
    @test @test_matrix_approx_eq sys[:RRR] RRRm
    @test @test_matrix_approx_eq sys[:ZZ]  ZZm
    @test @test_matrix_approx_eq sys[:DD]  DDm
    @test @test_matrix_approx_eq sys[:QQ]  QQm
    @test all(sys[:EE]  .== 0.)
    @test all(sys[:CCC] .== 0.)
end
