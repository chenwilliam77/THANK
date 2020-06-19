include("../../../src/includeall.jl")

m = COTHANK()
Γ0, Γ1, C, Ψ, Π = eqcond(m)

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
