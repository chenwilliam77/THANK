
"""
```
eqcond(m::JPT)
```
 Expresses the equilibrium conditions in canonical form using Γ0, Γ1, C, Ψ, and Π matrices.
Using the mappings of states/equations to integers defined in jpt.jl, coefficients are
specified in their proper positions.

Γ0 (n_states x n_states) holds coefficients of current time states.
Γ1 (n_states x n_states) holds coefficients of lagged states.
C  (n_states x 1) is a vector of constants
Ψ  (n_states x n_shocks_exogenous) holds coefficients of iid shocks.
Π  (n_states x n_states_expectational) holds coefficients of expectational states.
"""
function eqcond(m::JPT)
    endo = m.endogenous_states
    exo  = m.exogenous_shocks
    ex   = m.expected_shocks
    eq   = m.equilibrium_conditions

    Γ0 = zeros(n_states(m), n_states(m))
    Γ1 = zeros(n_states(m), n_states(m))
    C  = zeros(n_states(m))
    Ψ  = zeros(n_states(m), n_shocks_exogenous(m))
    Π  = zeros(n_states(m), n_shocks_expectational(m))

    ### ENDOGENOUS STATES ###

    ### 1. Production function

    # Sticky prices
    Γ0[eq[:eq_y], endo[:y_t]] = 1.
    Γ0[eq[:eq_y], endo[:k_t]] = -((m[:y_ss] + m[:F_ss]) / m[:y_ss]) * m[:α]
    Γ0[eq[:eq_y], endo[:L_t]] = -((m[:y_ss] + m[:F_ss]) / m[:y_ss]) * (1. - m[:α])

    # Flexible prices
    Γ0[eq[:eq_y_f], endo[:y_f_t]] = 1.
    Γ0[eq[:eq_y_f], endo[:k_f_t]] = -((m[:y_ss] + m[:F_ss]) / m[:y_ss]) * m[:α]
    Γ0[eq[:eq_y_f], endo[:L_f_t]] = -((m[:y_ss] + m[:F_ss]) / m[:y_ss]) * (1. - m[:α])

    ### 2. Cost minimization

    # Sticky prices
    Γ0[eq[:eq_L], endo[:Rk_t]] = 1.
    Γ0[eq[:eq_L], endo[:k_t]]  = 1.
    Γ0[eq[:eq_L], endo[:w_t]]  = -1.
    Γ0[eq[:eq_L], endo[:L_t]]  = -1.

    # Flexible prices
    Γ0[eq[:eq_L_f], endo[:Rk_f_t]] = 1.
    Γ0[eq[:eq_L_f], endo[:k_f_t]]  = 1.
    Γ0[eq[:eq_L_f], endo[:w_f_t]]  = -1.
    Γ0[eq[:eq_L_f], endo[:L_f_t]]  = -1.

    ### 3. Marginal cost

    # Sticky prices
    Γ0[eq[:eq_s], endo[:s_t]]  = 1.
    Γ0[eq[:eq_s], endo[:Rk_t]] = -m[:α]
    Γ0[eq[:eq_s], endo[:w_t]]  = -(1. - m[:α])

    # Flexible prices
    Γ0[eq[:eq_s_f], endo[:s_f_t]]  = 1.
    Γ0[eq[:eq_s_f], endo[:Rk_f_t]] = -m[:α]
    Γ0[eq[:eq_s_f], endo[:w_f_t]]  = -(1. - m[:α])

    ### 4. Phillips Curve


    # Sticky prices
    Γ0[eq[:eq_π], endo[:π_t]]   = 1.
    Γ0[eq[:eq_π], endo[:Eπ_t]]  = -m[:β] / (1 + m[:ι_p] * m[:β])
    Γ1[eq[:eq_π], endo[:π_t]]   = m[:ι_p] / (1 + m[:ι_p] * m[:β])
    Γ0[eq[:eq_π], endo[:s_t]]   = -(1. - m[:β] * m[:ξ_p]) * (1. - m[:ξ_p]) / ((1 + m[:ι_p] * m[:β]) * m[:ξ_p])
    Γ0[eq[:eq_π], endo[:λ_p_t]] = -1.

    # Flexible prices
    Γ0[eq[:eq_R_f], endo[:s_f_t]] = 1.

    ### 5. Consumption FOC

    # Sticky prices
    expγ = exp(m[:γ])
    Γ0[eq[:eq_c], endo[:λ_t]]  = (expγ - m[:h] * m[:β]) * (expγ - m[:h])
    Γ0[eq[:eq_c], endo[:b_t]]  = -(expγ - m[:h] * m[:β] * m[:ρ_b]) * (expγ - m[:h]) /
        ((1. - m[:ρ_b]) * (expγ - m[:h] * m[:β] * m[:ρ_b]) * (expγ - m[:h]) /
         (expγ * m[:h] + expγ ^ 2 + m[:β] * m[:h] ^ 2))
    Γ0[eq[:eq_c], endo[:z_t]]  = -(m[:β] * m[:h] * expγ * m[:ρ_z] - m[:h] * expγ)
    Γ0[eq[:eq_c], endo[:c_t]]  = expγ ^ 2 + m[:β] * m[:h] ^ 2
    Γ0[eq[:eq_c], endo[:Ec_t]] = -m[:β] * m[:h] * expγ
    Γ1[eq[:eq_c], endo[:c_t]]  = m[:h] * expγ

    # Flexible prices
    expγ = exp(m[:γ])
    Γ0[eq[:eq_c_f], endo[:λ_f_t]]  = (expγ - m[:h] * m[:β]) * (expγ - m[:h])
    Γ0[eq[:eq_c_f], endo[:b_t]]  = -(expγ - m[:h] * m[:β] * m[:ρ_b]) * (expγ - m[:h]) /
        ((1. - m[:ρ_b]) * (expγ - m[:h] * m[:β] * m[:ρ_b]) * (expγ - m[:h]) /
         (expγ * m[:h] + expγ ^ 2 + m[:β] * m[:h] ^ 2))
    Γ0[eq[:eq_c_f], endo[:z_t]]  = -(m[:β] * m[:h] * expγ * m[:ρ_z] - m[:h] * expγ)
    Γ0[eq[:eq_c_f], endo[:c_f_t]]  = expγ ^ 2 + m[:β] * m[:h] ^ 2
    Γ0[eq[:eq_c_f], endo[:Ec_f_t]] = -m[:β] * m[:h] * expγ
    Γ1[eq[:eq_c_f], endo[:c_f_t]]  = m[:h] * expγ

    ### 6. Euler equation

    # Sticky prices
    Γ0[eq[:eq_λ], endo[:λ_t]]  = 1.
    Γ0[eq[:eq_λ], endo[:R_t]]  = -1.
    Γ0[eq[:eq_λ], endo[:Eλ_t]] = -1.
    Γ0[eq[:eq_λ], endo[:Eπ_t]] = 1.
    Γ0[eq[:eq_λ], endo[:z_t]]  = m[:ρ_z]

    # Flexible prices
    Γ0[eq[:eq_λ_f], endo[:λ_f_t]]  = 1.
    Γ0[eq[:eq_λ_f], endo[:R_f_t]]  = -1.
    Γ0[eq[:eq_λ_f], endo[:Eλ_f_t]] = -1.
    Γ0[eq[:eq_λ_f], endo[:z_t]]  = m[:ρ_z]

    ### 7. Capital utilization FOC

    # Sticky prices
    Γ0[eq[:eq_Rk], endo[:Rk_t]] = 1.
    Γ0[eq[:eq_Rk], endo[:u_t]]  = -m[:χ]


    # Flexible prices
    Γ0[eq[:eq_Rk_f], endo[:Rk_f_t]] = 1.
    Γ0[eq[:eq_Rk_f], endo[:u_f_t]] = -m[:χ]

    ### 8. Capital FOC

    #Sticky prices
    Γ0[eq[:eq_ϕ], endo[:ϕ_t]] = 1.
    Γ0[eq[:eq_ϕ], endo[:Eϕ_t]] = -m[:β] * exp(-m[:γ]) * (1 - m[:δ])
    Γ0[eq[:eq_ϕ], endo[:z_t]] = m[:ρ_z]
    Γ0[eq[:eq_ϕ], endo[:Eλ_t]] = -1 * (1 - m[:β] * exp(-m[:γ]) * (1 - m[:δ]))
    Γ0[eq[:eq_ϕ], endo[:ERk_t]] = -(1 - m[:β] * exp(-m[:γ]) * (1 - m[:δ]))

    #Flexible prices
    Γ0[eq[:eq_ϕ_f], endo[:ϕ_f_t]] = 1
    Γ0[eq[:eq_ϕ_f], endo[:Eϕ_f_t]] = -m[:β] * exp(-m[:γ]) * (1 - m[:δ])
    Γ0[eq[:eq_ϕ_f], endo[:z_t]] = m[:ρ_z]
    Γ0[eq[:eq_ϕ_f], endo[:Eλ_f_t]] = -(1 - m[:β] * exp(-m[:γ]) * (1 - m[:δ]))
    Γ0[eq[:eq_ϕ_f], endo[:ERk_f_t]] = -(1 - m[:β] * exp(-m[:γ]) *(1 - m[:δ]))

    ### 9. Investment FOC

    #Sticky prices
    Γ0[eq[:eq_i], endo[:λ_t]] = 1 /( m[:S′′] * expγ^2)
    Γ0[eq[:eq_i], endo[:ϕ_t]] = -1 / (m[:S′′] * expγ^2)
    Γ0[eq[:eq_i], endo[:μ_t]] = -1 / (m[:S′′] * expγ^2)
    Γ0[eq[:eq_i], endo[:i_t]] = (1 + m[:β])
    Γ0[eq[:eq_i], endo[:z_t]] = (1 - m[:β] * m[:ρ_z])
    Γ0[eq[:eq_i], endo[:Ei_t]] = -m[:β]
    Γ1[eq[:eq_i], endo[:i_t]] = 1

    #Flexible prices
    Γ0[eq[:eq_i_f], endo[:λ_f_t]] = 1 / (m[:S′′] * expγ^2)
    Γ0[eq[:eq_i_f], endo[:ϕ_f_t]] = -1 / (m[:S′′] * expγ^2)
    Γ0[eq[:eq_i_f], endo[:μ_t]] = -1 / (m[:S′′] * expγ^2)
    Γ0[eq[:eq_i_f], endo[:i_f_t]] = (1 + m[:β])
    Γ0[eq[:eq_i_f], endo[:z_t]] = (1 - m[:β]*m[:ρ_z])
    Γ0[eq[:eq_i_f], endo[:Ei_f_t]] = -m[:β]
    Γ1[eq[:eq_i_f], endo[:i_f_t]] = 1

    ### Capital input

    #Sticky prices
    Γ0[eq[:eq_k], endo[:k_t]] = 1
    Γ0[eq[:eq_k], endo[:u_t]] = -1
    Γ0[eq[:eq_k], endo[:z_t]] = 1
    Γ1[eq[:eq_k], endo[:kbar_t]] =1

    #Flexible prices
    Γ0[eq[:eq_k_f], endo[:k_f_t]] = 1
    Γ0[eq[:eq_k_f], endo[:u_f_t]] = -1
    Γ0[eq[:eq_k_f], endo[:z_t]] = 1
    Γ1[eq[:eq_k_f], endo[:kbar_f_t]] = 1

    ### Capital Accumulation

    #Sticky prices
    Γ0[eq[:eq_kbar], endo[:kbar_t]] = 1
    Γ0[eq[:eq_kbar], endo[:μ_t]] = -(1 - (1 - m[:δ]) * exp(-m[:γ]))
    Γ0[eq[:eq_kbar], endo[:i_t]] = -(1 - (1 - m[:δ]) * exp(-m[:γ]))
    Γ1[eq[:eq_kbar], endo[:kbar_t]] = (1 - m[:δ]) * exp(-m[:γ])
    Γ0[eq[:eq_kbar], endo[:z_t]] = (1 - m[:δ]) * exp(-m[:γ])

    #Flexible Prices
    Γ0[eq[:eq_kbar_f], endo[:kbar_f_t]] = 1
    Γ0[eq[:eq_kbar_f], endo[:μ_t]] = -(1 - (1 - m[:δ]) * exp(-m[:γ]))
    Γ0[eq[:eq_kbar_f], endo[:i_f_t]] = -(1 - (1 - m[:δ]) * exp(-m[:γ]))
    Γ1[eq[:eq_kbar_f], endo[:kbar_f_t]] = (1 - m[:δ]) * exp(-m[:γ])
    Γ0[eq[:eq_kbar_f], endo[:z_t]] = (1 - m[:δ]) * exp(-m[:γ])

    ### wage Phillips curve

    #Sticky prices

    κw = ((1. - m[:ξ_w] * m[:β]) * (1. - m[:ξ_w])) / (m[:ξ_w] * (1. + m[:β]) * (1. + m[:ν] * (1. + 1. / m[:λ_w_ss])))
    Γ0[eq[:eq_w], endo[:w_t]] = 1
    Γ0[eq[:eq_w], endo[:Ew_t]] = -m[:β] / (1 + m[:β])
    Γ0[eq[:eq_w], endo[:wgap_t]] = κw
    Γ0[eq[:eq_w], endo[:π_t]] = (1 + m[:β] * m[:ι_w]) / (1 + m[:β])
    Γ0[eq[:eq_w], endo[:Eπ_t]] = -m[:β] / (1 + m[:β])
    Γ0[eq[:eq_w], endo[:z_t]] = (1 + m[:β] * m[:ι_w] - m[:β] * m[:ρ_z]) / (1 + m[:β])
    Γ0[eq[:eq_w], endo[:λ_w_t]] = -1
    Γ1[eq[:eq_w], endo[:w_t]] = 1 / (1 + m[:β])
    Γ1[eq[:eq_w], endo[:π_t]] = m[:ι_w] / (1 + m[:β])
    Γ1[eq[:eq_w], endo[:z_t]] = m[:ι_w] / (1 + m[:β])
    #Flexible prices
    Γ0[eq[:eq_w_f], endo[:wgap_f_t]] = 1

    ### Wage gap

    #Flexible prices
    Γ0[eq[:eq_wgap], endo[:wgap_t]] = 1
    Γ0[eq[:eq_wgap], endo[:w_t]] = -1
    Γ0[eq[:eq_wgap], endo[:b_t]] = 1 / ((1 - m[:ρ_b]) * (expγ - m[:h] * m[:β] * m[:ρ_b]) /
         *(expγ - m[:h]) / (expγ * m[:h] + expγ^2 + m[:β] * m[:h]^2))
    Γ0[eq[:eq_wgap], endo[:L_t]] = m[:ν]
    Γ0[eq[:eq_wgap], endo[:λ_t]] = -1

    #Sticky prices
    Γ0[eq[:eq_wgap_f], endo[:wgap_f_t]] = 1
    Γ0[eq[:eq_wgap_f], endo[:w_f_t]] = -1
    Γ0[eq[:eq_wgap_f], endo[:b_t]] = 1 / ((1 - m[:ρ_b]) * (expγ - m[:h] * m[:β] * m[:ρ_b]) /
         *(expγ - m[:h]) / (expγ * m[:h] + expγ^2 + m[:β] * m[:h]^2))
    Γ0[eq[:eq_wgap_f], endo[:L_f_t]] = m[:ν]
    Γ0[eq[:eq_wgap_f], endo[:λ_f_t]] = -1

    ### Monetary policy rule

    #Flexible prices
    Γ0[eq[:eq_R], endo[:R_t]] = 1
    Γ1[eq[:eq_R], endo[:R_t]] = m[:ρ_R]
    Γ0[eq[:eq_R], endo[:π_t]] = -(1 - m[:ρ_R]) *m[:ψ1]
    Γ0[eq[:eq_R], endo[:gdp_t]] = -(1 - m[:ρ_R]) * m[:ψ2] - m[:ψ3]
    Γ0[eq[:eq_R], endo[:gdp_f_t]] = (1 - m[:ρ_R]) * m[:ψ2] + m[:ψ3]
    Γ1[eq[:eq_R], endo[:gdp_t]] = -m[:ψ3]
    Γ1[eq[:eq_R], endo[:gdp_f_t]] = m[:ψ3]
    Γ0[eq[:eq_R], endo[:mp_t]] = -1

    ### definition of gdp

    #Sticky prices
    Γ0[eq[:eq_gdp], endo[:gdp_t]] = -1
    Γ0[eq[:eq_gdp], endo[:y_t]] = 1
    Γ0[eq[:eq_gdp], endo[:u_t]] = -m[:k_ss] * m[:Rk_ss] / m[:y_ss]

    #Flexible prices
    Γ0[eq[:eq_gdp_f], endo[:gdp_f_t]] = -1
    Γ0[eq[:eq_gdp_f], endo[:y_f_t]] = 1
    Γ0[eq[:eq_gdp_f], endo[:u_f_t]] = -m[:k_ss] * m[:Rk_ss] / m[:y_ss]

    ### market clearing

    #Sticky prices
    Γ0[eq[:eq_u], endo[:c_t]] = m[:c_ss] / m[:y_ss]
    Γ0[eq[:eq_u], endo[:i_t]] = m[:i_ss] / m[:y_ss]
    Γ0[eq[:eq_u], endo[:y_t]] = -1 / m[:g_ss]
    Γ0[eq[:eq_u], endo[:g_t]] = 1 / m[:g_ss]
    Γ0[eq[:eq_u], endo[:u_t]] = m[:k_ss] * m[:Rk_ss] / m[:y_ss]

    #Flexible Prices
    Γ0[eq[:eq_u_f], endo[:c_f_t]] = m[:c_ss] / m[:y_ss]
    Γ0[eq[:eq_u_f], endo[:i_f_t]] = m[:i_ss] / m[:y_ss]
    Γ0[eq[:eq_u_f], endo[:y_f_t]] = -1 / m[:g_ss]
    Γ0[eq[:eq_u_f], endo[:g_t]] = 1 / m[:g_ss]
    Γ0[eq[:eq_u_f], endo[:u_f_t]] = m[:k_ss] * m[:Rk_ss] / m[:y_ss]

    ### Exogenous shocks

    # z
    Γ0[eq[:eq_z], endo[:z_t]] = 1
    Γ1[eq[:eq_z], endo[:z_t]] = m[:ρ_z]
    Ψ[eq[:eq_z], exo[:z_sh]] = 1

    # g
    Γ0[eq[:eq_g], endo[:g_t]] = 1
    Γ1[eq[:eq_g], endo[:g_t]] = m[:ρ_g]
    Ψ[eq[:eq_g], exo[:g_sh]] = 1

    # λ_w
    Γ0[eq[:eq_λ_w], endo[:λ_w_t]] = 1
    Γ1[eq[:eq_λ], endo[:λ_w_t]] = m[:ρ_λ_w]
    Γ0[eq[:eq_λ_w], endo[:λ_w_t1]] = -1
    Γ1[eq[:eq_λ_w], endo[:λ_w_t1]] = -m[:η_λ_w]
    Γ0[eq[:eq_λ_w1], endo[:λ_w_t1]] = 1
    Ψ[eq[:eq_λ_w1], exo[:λ_w_sh]] = 1

    # μ
    Γ0[eq[:eq_μ], endo[:μ_t]] = 1
    Γ1[eq[:eq_μ], endo[:μ_t]] = m[:ρ_μ]
    Ψ[eq[:eq_μ], exo[:μ_sh]] = 1

    # λ_p
    Γ0[eq[:eq_λ_p], endo[:λ_p_t]] = 1
    Γ1[eq[:eq_λ_p], endo[:λ_p_t]] = m[:ρ_λ_p]
    Γ0[eq[:eq_λ_p], endo[:λ_p_t1]] = -1
    Γ1[eq[:eq_λ_p], endo[:λ_p_t1]] = -m[:η_λ_p]
    Γ0[eq[:eq_λ_p1], endo[:λ_p_t1]] = 1
    Ψ[eq[:eq_λ_p1], exo[:λ_p_sh]] = 1

    ### Expectational terms

    # π sticky
    Γ0[eq[:eq_Eπ], endo[:π_t]] = 1
    Γ1[eq[:eq_Eπ], endo[:Eπ_t]] = 1
    Π[eq[:eq_Eπ], ex[:Eπ_sh]] = 1

    # c sticky
    Γ0[eq[:eq_Ec], endo[:c_t]] = 1
    Γ1[eq[:eq_Ec], endo[:Ec_t]] = 1
    Π[eq[:eq_Ec], ex[:Ec_sh]] = 1

    # c flexible
    Γ0[eq[:eq_Ec_f], endo[:c_f_t]] = 1
    Γ1[eq[:eq_Ec_f], endo[:Ec_f_t]] = 1
    Π[eq[:eq_Ec_f], ex[:Ec_f_sh]] = 1

    # λ sticky
    Γ0[eq[:eq_Eλ], endo[:λ_t]] = 1
    Γ1[eq[:eq_Eλ], endo[:Eλ_t]] = 1
    Π[eq[:eq_Eλ], ex[:Eλ_sh]] = 1

    # λ flexible
    Γ0[eq[:eq_Eλ_f], endo[:λ_f_t]] = 1
    Γ1[eq[:eq_Eλ_f], endo[:Eλ_f_t]] = 1
    Π[eq[:eq_Eλ_f], ex[:Eλ_f_sh]] = 1

    # ϕ sticky
    Γ0[eq[:eq_Eϕ], endo[:ϕ_t]] = 1
    Γ1[eq[:eq_Eϕ], endo[:Eϕ_t]] = 1
    Π[eq[:eq_Eϕ], ex[:Eϕ_sh]] = 1

    # ϕ flexible
    Γ0[eq[:eq_Eϕ_f], endo[:ϕ_f_t]] = 1
    Γ1[eq[:eq_Eϕ_f], endo[:Eϕ_f_t]] = 1
    Π[eq[:eq_Eϕ_f], ex[:Eϕ_f_sh]] = 1

    # Rk sticky
    Γ0[eq[:eq_ERk], endo[:Rk_t]] = 1
    Γ1[eq[:eq_ERk], endo[:ERk_t]] = 1
    Π[eq[:eq_ERk], ex[:ERk_sh]] = 1

    # Rk flexible
    Γ0[eq[:eq_ERk_f], endo[:Rk_f_t]] = 1
    Γ1[eq[:eq_ERk_f], endo[:ERk_f_t]] = 1
    Π[eq[:eq_ERk_f], ex[:ERk_f_sh]] = 1

    # i sticky
    Γ0[eq[:eq_Ei], endo[:i_t]] = 1
    Γ1[eq[:eq_Ei], endo[:Ei_t]] = 1
    Π[eq[:eq_Ei], ex[:Ei_sh]] = 1

    #  i flexible
    Γ0[eq[:eq_Ei_f], endo[:i_f_t]] = 1
    Γ1[eq[:eq_Ei_f], endo[:Ei_f_t]] = 1
    Π[eq[:eq_Ei_f], ex[:Ei_f_sh]] = 1

    # w sticky
    Γ0[eq[:eq_Ew], endo[:w_t]] = 1
    Γ1[eq[:eq_Ew], endo[:Ew_t]] = 1
    Π[eq[:eq_Ew], ex[:Ew_sh]] = 1

    ### Lagged variables

    # gdp
    Γ0[eq[:eq_gdp1], endo[:gdp_t1]] = 1
    Γ1[eq[:eq_gdp1], endo[:gdp_t]] = 1

    # c
    Γ0[eq[:eq_c1], endo[:c_t1]] = 1
    Γ1[eq[:eq_c1], endo[:c_t]] = 1

    # i
    Γ0[eq[:eq_i1], endo[:i_t1]] = 1
    Γ1[eq[:eq_i1], endo[:i_t]] = 1

    # w
    Γ0[eq[:eq_w1], endo[:w_t1]] = 1
    Γ1[eq[:eq_w1], endo[:w_t]] = 1

    return Γ0, Γ1, C, Ψ, Π

end
