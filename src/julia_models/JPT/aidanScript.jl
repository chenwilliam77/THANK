include("../../includeall.jl")

using MAT, Test

# instantiate model
m = JPT();

# equilibrium conditions
Γ0, Γ1, C, Ψ, Π = eqcond(m);

GAM0, GAM1, PSI, PPI, ssvec = matopen("../../models/JPT/testmatrices.mat") do file;
    read(file, "GAM0"), read(file, "GAM1"), read(file, "PSI"), read(file, "PPI"), read(file,"ssvec");
end

@testset "Compare eqcond against reference matrices" begin
    @test @test_matrix_approx_eq Γ0 GAM0
    @test @test_matrix_approx_eq Γ1 GAM1
    @test @test_matrix_approx_eq Ψ PSI
    @test @test_matrix_approx_eq Π PPI
end
