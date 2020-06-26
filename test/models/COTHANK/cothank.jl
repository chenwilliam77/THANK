include("../../../src/includeall.jl")

m = COTHANK()
Γ0, Γ1, C, Ψ, Π = eqcond(m)
#=
rows_zero_Γ0 = []
for i in 1:size(Γ0, 1)
    if all(Γ0[i, :] .== 0.)
        push!(rows_zero_Γ0, i)
    end
end

cols_zero_Γ0 = []
for i in 1:size(Γ0, 1)
    if all(Γ0[:, i] .== 0.)
        push!(rows_zero_Γ0, i)
    end
end
=#

eqcondmats = matopen("eqcond.mat")

Γ0m = read(eqcondmats, "GAM0")
Γ1m = read(eqcondmats, "GAM1")
Cm  = read(eqcondmats, "C")
Ψm  = read(eqcondmats, "PSI")
Πm  = read(eqcondmats, "PPI")
ssmat = read(eqcondmats,"ssvec")
