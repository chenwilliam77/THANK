"""
```
eqcond(m::COTHANK)
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
function eqcond(m::COTHANK)
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

    ### 1. Intermediate goods producers

    # Sticky prices
    Γ0[eq[:eq_KL_ratio_good1], endo[:klR1_t]] = -1.
    Γ0[eq[:eq_KL_ratio_good1], endo[:w1_t]]   = 1.
    Γ0[eq[:eq_KL_ratio_good1], endo[:ρ_t]]    = -1.

    Γ0[eq[:eq_mc_good1], endo[:mc1_t]]    = -1.
    Γ0[eq[:eq_mc_good1], endo[:ρ_t]]      = m[:α]
    Γ0[eq[:eq_mc_good1], endo[:w1_t]]     = (1. - m[:α])

    Γ0[eq[:eq_KL_ratio_good2], endo[:klR2_t]] = -1.
    Γ0[eq[:eq_KL_ratio_good2], endo[:w2_t]]   = 1.
    Γ0[eq[:eq_KL_ratio_good2], endo[:ρ_t]]    = -1.

    Γ0[eq[:eq_mc_good2], endo[:mc2_t]]   = -1.
    Γ0[eq[:eq_mc_good2], endo[:ρ_t]]     = m[:α]
    Γ0[eq[:eq_mc_good2], endo[:w2_t]]    = (1. - m[:α])
    Γ0[eq[:eq_mc_good2], endo[:a_t]]     = -(1. - m[:α])
#=
    # Flexible prices
    Γ0[eq[:eq_klR1_f], endo[:klR1_t_f]] = -1.
    Γ0[eq[:eq_klR1_f], endo[:w1_t_f]]   = 1.
    Γ0[eq[:eq_klR1_f], endo[:ρ_t_f]]    = -1.

    Γ0[eq[:eq_mc1_f], endo[:mc1_f_t]]     = -1
    Γ0[eq[:eq_mc1_f], endo[:ρ_t_f]]      = m[:α]
    Γ0[eq[:eq_mc1_f], endo[:w1_t_f]]     = (1 - m[:α])

    Γ0[eq[:eq_klR2_f], endo[:klR2_f_t]] = -1
    Γ0[eq[:eq_klR2_f], endo[:w2_f_t]]   = 1
    Γ0[eq[:eq_klR2_f], endo[:ρ_f_t]]    = -1

    Γ0[eq[:eq_mc2_f], endo[:mc2_f_t]]   = -1
    Γ0[eq[:eq_mc2_f], endo[:ρ_f_t]]   = m[:α]
    Γ0[eq[:eq_mc2_f], endo[:w2_f_t]]   = (1 - m[:α])
    Γ0[eq[:eq_mc2_f], endo[:a_f_t]]   = (1 - m[:α])
=#
    ### 2. Production of goods 1 and 2

    # Sticky prices
    Γ0[eq[:eq_production_good1], endo[:y1_t]]   = -1.
    Γ0[eq[:eq_production_good1], endo[:klR1_t]] = (1. + m[:λ_p_ss]) * m[:α]
    Γ0[eq[:eq_production_good1], endo[:L1_t]]   = (1. + m[:λ_p_ss])

    Γ0[eq[:eq_production_good2], endo[:y2_t]]   = -1
    Γ0[eq[:eq_production_good2], endo[:klR2_t]] = (1. + m[:λ_p_ss]) * m[:α]
    Γ0[eq[:eq_production_good2], endo[:L2_t]]   = (.1 + m[:λ_p_ss])
    Γ0[eq[:eq_production_good2], endo[:a_t]]    = -m[:λ_p_ss]

#=
    # Flexible prices
    Γ0[eq[:eq_y1_f], endo[:y1_f_t]] = -1
    Γ0[eq[:eq_y1_f], endo[:klR1_f_t]] = (1 + m[:λ_p]) * m[:α]
    Γ0[eq[:eq_y1_f], endo[:L1_f_t]] = (1 + m[:λ_p])
    Γ0[eq[:eq_y1_f], endo[:a_f_t]] = -m[:λ_p]
=#
    ### 3. Production of final gods

    # Sticky prices
    Γ0[eq[:eq_production_finalgood], endo[:y_t]] = -1.
    Γ0[eq[:eq_production_finalgood], endo[:y1_t]] = m[:ν]
    Γ0[eq[:eq_production_finalgood], endo[:y2_t]] = (1 - m[:ν])
    Γ0[eq[:eq_production_finalgood], endo[:d_t]] = -(1 - m[:ν]) / (1 - m[:ζ])

    Γ0[eq[:eq_price_index], endo[:π_t]] = -1.
    Γ0[eq[:eq_price_index], endo[:d_t]] = (1. - m[:ν]) / (1. - m[:ζ])
    Γ0[eq[:eq_price_index], endo[:π1_t]] = m[:ν]
    Γ0[eq[:eq_price_index], endo[:π2_t]] = (1. - m[:ν])
    Γ1[eq[:eq_price_index], endo[:d_t]] = (1. - m[:ν]) / (1. - m[:ζ])

    ### 4. Price Phillips Curve

    # Sticky prices
    κ = (1. - m[:ξ_p] * m[:β]) * (1. - m[:ξ_p]) / ( m[:ξ_p] * (1. + m[:ι_p] * m[:β]))
    Γ0[eq[:eq_price_phillips_good1], endo[:π1_t]] = -1.
    Γ0[eq[:eq_price_phillips_good1], endo[:Eπ1_t]] = m[:β] / (1. + m[:ι_p] * m[:β])
    Γ0[eq[:eq_price_phillips_good1], endo[:mc1_t]] = κ
    Γ0[eq[:eq_price_phillips_good1], endo[:λ_p_t]] = 1. # IN THE NOTES, IT SAYS WE REDEFINE λ̂_p_t = κ λ̂_p_t, DOES THIS APPLY THE REST OF THE DOCUMENT OR JUST FOR THE PRICE_PHILLIPS CURVE EQUATIONS?
    Γ1[eq[:eq_price_phillips_good1], endo[:π1_t]]  = -m[:ι_p] / (1. + m[:ι_p] * m[:β])

    Γ0[eq[:eq_price_phillips_good2], endo[:π2_t]]  = -1.
    Γ0[eq[:eq_price_phillips_good2], endo[:Eπ2_t]] = m[:β] / (1. + m[:ι_p] * m[:β])
    Γ0[eq[:eq_price_phillips_good2], endo[:mc2_t]] = κ
    Γ0[eq[:eq_price_phillips_good2], endo[:λ_p_t]] = 1.
    Γ1[eq[:eq_price_phillips_good2], endo[:π2_t]]  = -m[:ι_p] / (1. + m[:ι_p] * m[:β])

    # Flexible prices
#    Γ0[eq[:eq_R_f], endo[:mc_f_t]] = 1 # THIS CAME FROM THANK, CHECK IF IT'S RIGHT

    ### 5. Consumption
    expγ = exp(m[:γ])

    # Sticky prices
    Γ0[eq[:eq_marg_utility_S], endo[:λ_S_t]] = -1.
    Γ0[eq[:eq_marg_utility_S], endo[:b_t]]   = 1.
    Γ0[eq[:eq_marg_utility_S], endo[:z_t]]   = -m[:h] / (expγ - m[:h])
    Γ0[eq[:eq_marg_utility_S], endo[:c_S_t]] = -expγ / (expγ - m[:h])
    Γ1[eq[:eq_marg_utility_S], endo[:c_S_t]] = m[:h] / (expγ - m[:h])
#=
    # Flexible
    Γ0[eq[:eq_marg_utility_S_f], endo[:λ_S_f_t]] = -1.
    Γ0[eq[:eq_marg_utility_S_f], endo[:b_f_t]]   = 1.
    Γ0[eq[:eq_marg_utility_S_f], endo[:z_f_t]]   = - m[:h] / (expγ - m[:h])
    Γ0[eq[:eq_marg_utility_S_f], endo[:c_S_f_t]] = -expγ / (expγ - m[:h])
    Γ1[eq[:eq_marg_utility_S_f], endo[:C_S_f_t]] = m[:h] / (expγ - m[:h])
=#
    # Sticky
    Γ0[eq[:eq_marg_utility_H1], endo[:λ_H1_t]] = -1.
    Γ0[eq[:eq_marg_utility_H1], endo[:b_t]]    = 1.
    Γ0[eq[:eq_marg_utility_H1], endo[:z_t]]    = -m[:h] / (expγ - m[:h])
    Γ0[eq[:eq_marg_utility_H1], endo[:c_H1_t]] = -expγ / (expγ - m[:h])
    Γ1[eq[:eq_marg_utility_H1], endo[:c_H1_t]] = m[:h] / (expγ - m[:h])
#=
    # Flexible
    Γ0[eq[:eq_marg_utility_H1_f], endo[:λ_H1_f_t]] = -1.
    Γ0[eq[:eq_marg_utility_H1_f], endo[:b_f_t]]   = 1.
    Γ0[eq[:eq_marg_utility_H1_f], endo[:z_f_t]]   = - m[:h] / (expγ - m[:h])
    Γ0[eq[:eq_marg_utility_H1_f], endo[:c_H1_f_t]] = -expγ / (expγ - m[:h])
    Γ1[eq[:eq_marg_utility_H1_f], endo[:c_H1_f_t]] = m[:h] / (expγ - m[:h])
=#
    # Sticky
    Γ0[eq[:eq_marg_utility_H2], endo[:λ_H2_t]] = -1.
    Γ0[eq[:eq_marg_utility_H2], endo[:b_t]]    = 1.
    Γ0[eq[:eq_marg_utility_H2], endo[:z_t]]    = -m[:h] / (expγ - m[:h])
    Γ0[eq[:eq_marg_utility_H2], endo[:c_H2_t]] = -expγ / (expγ - m[:h])
    Γ1[eq[:eq_marg_utility_H2], endo[:c_H2_t]] = m[:h] / (expγ - m[:h])
#=
    # Flexible
    Γ0[eq[:eq_marg_utility_H2_f], endo[:λ_H2_f_t]] = -1
    Γ0[eq[:eq_marg_utility_H2_f], endo[:b_f_t]]   = 1
    Γ0[eq[:eq_marg_utility_H2_f], endo[:z_f_t]]   = - m[:h] / (expγ - m[:h])
    Γ0[eq[:eq_marg_utility_H2_f], endo[:c_H2_f_t]] = -expγ / (expγ - m[:h])
    Γ1[eq[:eq_marg_utility_H2_f], endo[:c_H2_f_t]] = m[:h] / (expγ - m[:h])
=#
    # Sticky

    Γ0[eq[:eq_average_c], endo[:c_t]]   = -1
    Γ0[eq[:eq_average_c], endo[:c_S_t]] = (1 - m[:θ]) * m[:c_S_ss] / m[:c_ss]
    Γ0[eq[:eq_average_c], endo[:c_H1_t]] = m[:θ] * m[:f_H1] * m[:c_H1_ss] / m[:c_ss]
    Γ0[eq[:eq_average_c], endo[:c_H2_t]] = m[:θ] * (1 - m[:f_H1]) * m[:c_H2_ss] / m[:c_ss]
#=
    # Flexible
    Γ0[eq[:eq_average_c_f], endo[:c_f_t]]   = -1
    Γ0[eq[:eq_average_c_f], endo[:c_S_t]] = (1 - m[:θ]) * m[:c_S_ss] / m[:c_ss]
    Γ0[eq[:eq_average_c_f], endo[:c_H1_t]] = m[:θ] * m[:f_H1] * m[:c_H1_ss] / m[:c_ss]
    Γ0[eq[:eq_average_c_f], endo[:c_H2_t]] = m[:θ] * (1 - m[:f_H1]) * m[:c_H2_ss] / m[:c_ss]
=#
    # Sticky
    Γ0[eq[:eq_euler], endo[:λ_S_t]]   = -1. # Euler equation also pins down the bonds demand/self-insurance motive
    Γ0[eq[:eq_euler], endo[:R_t]]     = 1.
    Γ0[eq[:eq_euler], endo[:Ez_t]]    = -1.
    Γ0[eq[:eq_euler], endo[:Eπ_t]]    = -1.
    Γ0[eq[:eq_euler], endo[:Eλ_S_t]]  = m[:β] * m[:R_ss] / (expγ * m[:π_ss]) * m[:s]
    Γ0[eq[:eq_euler], endo[:Eλ_H1_t]] =  m[:β] * m[:R_ss] / (expγ * m[:π_ss]) *
                                         m[:f_S1] * m[:λ_H1_ss] * (1. - m[:s]) / m[:λ_S_ss]
    Γ0[eq[:eq_euler], endo[:Eλ_H2_t]] =  m[:β] * m[:R_ss] / (expγ * m[:π_ss]) * (1. - m[:f_S1]) *
                                         m[:λ_H2_ss] * (1. - m[:s]) / m[:λ_S_ss]
    Γ0[eq[:eq_euler], endo[:Ex_t]]    = m[:β] * m[:R_ss] / (expγ * m[:π_ss]) * m[:s] *
                                        (1. - m[:f_S1] * m[:λ_H1_ss] / m[:λ_S_ss] -
                                        (1. - m[:f_S1]) * m[:λ_H2_ss] / m[:λ_S_ss]) *
                                        m[:s_x] / (m[:s] * m[:x_ss])
#=
    # Flexible
    Γ0[eq[:eq_euler_f], endo[:λ_S_f_t]] = -1.
    Γ0[eq[:eq_euler_f], endo[:R_f_t]] = 1.
    Γ0[eq[:eq_euler_f], endo[:Ez_f_t]] = -1.
    Γ0[eq[:eq_euler_f], endo[:Eπ_f_t]] = -1.
    Γ0[eq[:eq_euler_f], endo[:Eλ_S_f_t]] = m[:β] * m[:R_ss] / (expγ * m[:π_ss]) * m[:s]
    Γ0[eq[:eq_euler_f], endo[:Eλ_H1_f_t]] =  m[:β] * m[:R_ss] / (expγ * m[:π_ss]) *
                                           m[:f_S1] * m[:λ_H1_ss] * (1. - m[:s]) / m[:λ_S_ss]
    Γ0[eq[:eq_euler_f], endo[:Eλ_H2_f_t]] =  m[:β] * m[:R_ss] / (expγ * m[:π_ss]) * (1. - m[:f_S1]) *
                                           m[:λ_H2_ss] * (1. - m[:s]) / m[:λ_S_ss]
    Γ0[eq[:eq_euler_f], endo[:Ex_f_t]] = m[:β] * m[:R_ss] / (expγ * m[:π_ss]) * m[:s] *
                                       (1. - m[:f_S1] * m[:λ_H1_ss] / m[:λ_S_ss] -
                                       (1. - m[:f_S1]) * m[:λ_H2_ss] / m[:λ_S_ss]) *
                                       m[:s_x] / (m[:s] * m[:x_ss])
=#
    # Sticky
    Γ0[eq[:eq_c_H1], endo[:c_H1_t]] = -1.
    Γ0[eq[:eq_c_H1], endo[:w1_t]]   = m[:w1_ss] * m[:H1] / m[:c_H1_ss]
    Γ0[eq[:eq_c_H1], endo[:L1_t]]   = m[:w1_ss] * m[:H1] / m[:c_H1_ss]
    Γ0[eq[:eq_c_H1], endo[:t_H1_t]] = -m[:x_ss] / m[:c_H1_ss]
    Γ0[eq[:eq_c_H1], endo[:τ_H1_t]] = m[:x_ss] / m[:c_H1_ss]
    Γ0[eq[:eq_c_H1], endo[:z_t]]    = -m[:f_S1] * (1. - m[:s]) / (m[:θ] * m[:f_H1]) *
                                       m[:R_ss] / (expγ * m[:π_ss]) * m[:bᴿ_ss] / m[:c_H1_ss]
    Γ0[eq[:eq_c_H1], endo[:π_t]]    = -m[:f_S1] * (1. - m[:s]) / (m[:θ] * m[:f_H1]) *
                                       m[:R_ss] / (expγ * m[:π_ss]) * m[:bᴿ_ss] / m[:c_H1_ss]
    Γ0[eq[:eq_c_H1], endo[:x_t]]    = -m[:f_S1] * m[:s] / (m[:θ] * m[:f_H1]) *
                                       m[:R_ss] / (expγ * m[:π_ss]) * m[:bᴿ_ss] / m[:c_H1_ss] *
                                       m[:s_x]/(m[:s]* m[:x_ss])
    Γ1[eq[:eq_c_H1], endo[:R_t]]    = -m[:f_S1] * (1. - m[:s]) / (m[:θ] * m[:f_H1]) *
                                       m[:R_ss] / (expγ * m[:π_ss]) * m[:bᴿ_ss] / m[:c_H1_ss]
    Γ1[eq[:eq_c_H1], endo[:bᴿ_t]]   = -m[:f_S1] * (1. - m[:s]) / (m[:θ] * m[:f_H1]) *
                                       m[:R_ss] / (expγ * m[:π_ss]) * m[:bᴿ_ss] / m[:c_H1_ss] *
                                       m[:x_ss] / m[:bᴿ_ss]
#=
    # Flexible
    Γ0[eq[:eq_c_H1_f], endo[:c_f_H1_t]] = -1
    Γ0[eq[:eq_c_H1_f], endo[:w1_f_t]] = m[:w1_ss] * m[:H1] / m[:c_H1_ss]
    Γ0[eq[:eq_c_H1_f], endo[:H1_f_t]] = m[:w1_ss] * m[:H1] / m[:c_H1_ss]
    Γ0[eq[:eq_c_H1_f], endo[:t_H1_f_t]] = -m[:x_ss] / m[:c_H1_ss]
    Γ0[eq[:eq_c_H1_f], endo[:τ_H1_f_t]] = m[:x_ss] / m[:c_H1_ss]
    Γ0[eq[:eq_c_H1_f], endo[:z_f_t]] = -m[:f_S1] * (1 - m[:s]) / (m[:θ] * m[:f_H1]) *
                                        m[:R_ss] / (expγ * m[:π_ss]) * m[:bᴿ_ss] / m[:c_H1_ss]
    Γ0[eq[:eq_c_H1_f], endo[:π_f_t]] = -m[:f_S1] * (1 - m[:s]) / (m[:θ] * m[:f_H1]) *
                                        m[:R_ss] / (expγ * m[:π_ss]) * m[:bᴿ_ss] / m[:c_H1_ss]
    Γ0[eq[:eq_c_H1_f], endo[:x_f_t]] = -m[:f_S1] * m[:s] / (m[:θ] * m[:f_H1]) *
                                        m[:R_ss] / (expγ * m[:π_ss]) * m[:bᴿ_ss] / m[:c_H1_ss] *
                                        m[:s_x]/(m[:s]*m[:x_ss])
    Γ1[eq[:eq_c_H1_f], endo[:R_f_t]] = -m[:f_S1] * (1 - m[:s]) / (m[:θ] * m[:f_H1]) *
                                        m[:R_ss] / (expγ * m[:π_ss]) * m[:bᴿ_ss] / m[:c_H1_ss]
    Γ1[eq[:eq_c_H1_f], endo[:bᴿ_f_t]] = -m[:f_S1] * (1 - m[:s]) / (m[:θ] * m[:f_H1]) *
                                        m[:R_ss] / (expγ * m[:π_ss]) * m[:bᴿ_ss] / m[:c_H1_ss] *
                                        m[:x_ss] / m[:bᴿ_ss]
=#
    # Sticky
    Γ0[eq[:eq_c_H2], endo[:c_H2_t]] = -1
    Γ0[eq[:eq_c_H2], endo[:w2_t]]   = m[:w2_ss] * m[:H2] / m[:c_H2_ss]
    Γ0[eq[:eq_c_H2], endo[:L2_t]]   = m[:w2_ss] * m[:H2] / m[:c_H2_ss]
    Γ0[eq[:eq_c_H2], endo[:t_H2_t]] = -m[:x_ss] / m[:c_H2_ss]
    Γ0[eq[:eq_c_H2], endo[:τ_H2_t]] = m[:x_ss] / m[:c_H2_ss]
    Γ0[eq[:eq_c_H2], endo[:z_t]]    = -(1. - m[:f_S1]) * (1 - m[:s]) / (m[:θ] * (1. - m[:f_H1])) *
                                       m[:R_ss] / (expγ * m[:π_ss]) * m[:bᴿ_ss] / m[:c_H2_ss]
    Γ0[eq[:eq_c_H2], endo[:π_t]]    = -(1. - m[:f_S1]) * (1 - m[:s]) / (m[:θ] * (1. - m[:f_H1])) *
                                       m[:R_ss] / (expγ * m[:π_ss]) * m[:bᴿ_ss] / m[:c_H2_ss]
    Γ0[eq[:eq_c_H2], endo[:x_t]]    = -(1. - m[:f_S1]) * m[:s] / (m[:θ] * (1. - m[:f_H1])) *
                                       m[:R_ss] / (expγ * m[:π_ss]) * m[:bᴿ_ss] / m[:c_H2_ss] *
                                       m[:s_x] / (m[:s] * m[:x_ss])
    Γ1[eq[:eq_c_H2], endo[:R_t]]    = -(1. - m[:f_S1]) * (1 - m[:s]) / (m[:θ] * (1. - m[:f_H1])) *
                                       m[:R_ss] / (expγ * m[:π_ss]) * m[:bᴿ_ss] / m[:c_H2_ss]
    Γ1[eq[:eq_c_H2], endo[:bᴿ_t]]    = -(1. - m[:f_S1]) * (1 - m[:s]) / (m[:θ] * (1. - m[:f_H1])) *
                                        m[:R_ss] / (expγ * m[:π_ss]) * m[:bᴿ_ss] / m[:c_H2_ss] *
                                        m[:x_ss] / m[:bᴿ_ss]
#=
    # Flexible
    Γ0[eq[:eq_c_H2_f], endo[:c_f_H2_t]] = -1
    Γ0[eq[:eq_c_H2_f], endo[:w2_f_t]] = m[:w2_ss] * m[:H2] / m[:c_H2_ss]
    Γ0[eq[:eq_c_H2_f], endo[:H2_f_t]] = m[:w2_ss] * m[:H2] / m[:c_H2_ss]
    Γ0[eq[:eq_c_H2_f], endo[:t_H2_f_t]] = -m[:x_ss] / m[:c_H2_ss]
    Γ0[eq[:eq_c_H2_f], endo[:τ_H2_f_t]] = m[:x_ss] / m[:c_H2_ss]
    Γ0[eq[:eq_c_H2_f], endo[:z_f_t]] = -(1. - m[:f_S1]) * (1 - m[:s]) / (m[:θ] * (1. - m[:f_H1])) *
                                        m[:R_ss] / (expγ * m[:π_ss]) * m[:bᴿ_ss] / m[:c_H2_ss]
    Γ0[eq[:eq_c_H2_f], endo[:π_f_t]] = -(1. - m[:f_S1]) * (1 - m[:s]) / (m[:θ] * (1. - m[:f_H1])) *
                                        m[:R_ss] / (expγ * m[:π_ss]) * m[:bᴿ_ss] / m[:c_H2_ss]
    Γ0[eq[:eq_c_H2_f], endo[:x_f_t]] = -(1. - m[:f_S1]) * m[:s] / (m[:θ] * (1. - m[:f_H1])) *
                                        m[:R_ss] / (expγ * m[:π_ss]) * m[:bᴿ_ss] / m[:c_H2_ss] *
                                        m[:s_x] / (m[:s] * m[:x_ss])
    Γ1[eq[:eq_c_H2_f], endo[:R_f_t]] = -(1. - m[:f_S1]) * (1 - m[:s]) / (m[:θ] * (1. - m[:f_H1])) *
                                        m[:R_ss] / (expγ * m[:π_ss]) * m[:bᴿ_ss] / m[:c_H2_ss]
    Γ1[eq[:eq_c_H2_f], endo[:bᴿ_f_t]] = -(1. - m[:f_S1]) * (1 - m[:s]) / (m[:θ] * (1. - m[:f_H1])) *
                                        m[:R_ss] / (expγ * m[:π_ss]) * m[:bᴿ_ss] / m[:c_H2_ss] *
                                        m[:x_ss] / m[:bᴿ_ss]
=#


    ### 6. Investment and capital accumulation

    # Sticky prices
    Γ0[eq[:eq_cap_util], endo[:ρ_t]] = 1.
    Γ0[eq[:eq_cap_util], endo[:u_t]] = -m[:χ]

    Γ0[eq[:eq_ϕ], endo[:ϕ_t]]    = -1.
    Γ0[eq[:eq_ϕ], endo[:Eϕ_t]]   = (1. - m[:δ]) * m[:β] * exp(-m[:γ])
    Γ0[eq[:eq_ϕ], endo[:Ez_t]]   = -1.
    Γ0[eq[:eq_ϕ], endo[:Eλ_S_t]] = (1. - (1. - m[:δ]) * m[:β] * exp(-m[:γ]))
    Γ0[eq[:eq_ϕ], endo[:Eρ_t]]   = (1. - (1. - m[:δ]) * m[:β] * exp(-m[:γ]))

    Γ0[eq[:eq_i], endo[:λ_S_t]]  = -1.
    Γ0[eq[:eq_i], endo[:ϕ_t]]    = 1.
    Γ0[eq[:eq_i], endo[:μ_t]]    = 1.
    Γ0[eq[:eq_i], endo[:ι_S_t]]  = -exp(2. * m[:γ]) * m[:S′′] * (1. + m[:β])
    Γ0[eq[:eq_i], endo[:z_t]]    = -exp(2. * m[:γ]) * m[:S′′]
    Γ0[eq[:eq_i], endo[:Eι_S_t]] = m[:β] * exp(2. * m[:γ]) * m[:S′′]
    Γ0[eq[:eq_i], endo[:Ez_t]]   = m[:β] * exp(2. * m[:γ]) * m[:S′′]
    Γ1[eq[:eq_i], endo[:ι_S_t]]  = -exp(2. * m[:γ]) * m[:S′′]

    Γ0[eq[:eq_cap_input], endo[:k_t]]   = -1.
    Γ0[eq[:eq_cap_input], endo[:u_t]]   = 1.
    Γ0[eq[:eq_cap_input], endo[:z_t]]   = -1.
    Γ1[eq[:eq_cap_input], endo[:k_S_t]] = -1.

    Γ0[eq[:eq_cap_accum], endo[:k_S_t]] = 1.
    Γ0[eq[:eq_cap_accum], endo[:z_t]]   = (1. - m[:δ]) * exp(-m[:γ])
    Γ0[eq[:eq_cap_accum], endo[:μ_t]]   = (1. - (1. - m[:δ]) * exp(-m[:γ]))
    Γ0[eq[:eq_cap_accum], endo[:ι_S_t]] = (1. - (1. - m[:δ]) * exp(-m[:γ]))
    Γ1[eq[:eq_cap_accum], endo[:k_S_t]] = (1. - m[:δ]) * exp(-m[:γ])
#=
    # Flexible prices
    Γ0[eq[:eq_k_S_f], endo[:k_S_f_t]] = 1
    Γ0[eq[:eq_k_S_f], endo[:z_t]] = (1 - m[:δ]) * exp(-m[:γ])
    Γ0[eq[:eq_k_S_f], endo[:μ_t]] = (1 - (1 - m[:δ]) * exp(-m[:γ]))
    Γ0[eq[:eq_k_S_f], endo[:ι_S_f_t]] = (1 - (1 - m[:δ]) * exp(-m[:γ]))
    Γ1[eq[:eq_k_S_f], endo[:k_S_f_t]] = (1 - m[:δ]) * exp(-m[:γ])
=#
    ### 7. Wage Phillips Curve
# wage phillips curve, marg sub evolution, average marginal utility labor
    # Sticky prices
    κ_w = (1. - m[:ξ_w] * m[:β]) * (1. - m[:ξ_w]) / (m[:ξ_w] * (1. + m[:β]) * (1. +
              m[:ν] * (1. + 1. / m[:λ_w_ss])))
    Γ0[eq[:eq_wage_phillips_good1], endo[:w1_t]]   = -1.
    Γ0[eq[:eq_wage_phillips_good1], endo[:Ew1_t]]  = m[:β] / (1. + m[:β])
    Γ0[eq[:eq_wage_phillips_good1], endo[:g_w1_t]] = -κ_w
    Γ0[eq[:eq_wage_phillips_good1], endo[:π_t]]    = -(1. + m[:β] * m[:ι_w]) / (1. + m[:β])
    Γ0[eq[:eq_wage_phillips_good1], endo[:Eπ_t]]   = m[:β] / (1. + m[:β])
    Γ0[eq[:eq_wage_phillips_good1], endo[:z_t]]    = -(1. + m[:β] * m[:ι_w] - m[:ρ_z] * m[:β]) / (1. + m[:β])
    Γ0[eq[:eq_wage_phillips_good1], endo[:λ_w_t]]  = 1.
    Γ1[eq[:eq_wage_phillips_good1], endo[:w1_t]]   = -1. / (1. + m[:β])
    Γ1[eq[:eq_wage_phillips_good1], endo[:π_t]]    = -m[:ι_w] / (1. + m[:β])
    Γ1[eq[:eq_wage_phillips_good1], endo[:z_t]]    = -m[:ι_w] / (1. + m[:β])

    Γ0[eq[:eq_marg_sub_good1], endo[:g_w1_t]]  = 1.
    Γ0[eq[:eq_marg_sub_good1], endo[:w1_t]]    = -1.
    Γ0[eq[:eq_marg_sub_good1], endo[:L1_t]]    = m[:ν]
    Γ0[eq[:eq_marg_sub_good1], endo[:b_t]]     = 1.
    Γ0[eq[:eq_marg_sub_good1], endo[:λ_SH1_t]] = -1.

    Γ0[eq[:eq_avg_marg_util_SH1], endo[:λ_SH1_t]] = -1.
    Γ0[eq[:eq_avg_marg_util_SH1], endo[:λ_S_t]]   = (1. - m[:θ]) * m[:f_S1] / ((1. - m[:θ]) * m[:f_S1] + m[:θ] * m[:f_H1]) *
                                                    m[:λ_S_ss] / m[:λ_SH1_ss]
    Γ0[eq[:eq_avg_marg_util_SH1], endo[:λ_H1_t]]  = m[:θ] * m[:f_H1] / ((1. - m[:θ]) * m[:f_S1] + m[:θ] * m[:f_H1]) *
                                                    m[:λ_H1_ss] / m[:λ_SH1_ss]

    Γ0[eq[:eq_wage_phillips_good2], endo[:w2_t]]   = -1.
    Γ0[eq[:eq_wage_phillips_good2], endo[:Ew2_t]]  = m[:β] / (1. + m[:β])
    Γ0[eq[:eq_wage_phillips_good2], endo[:g_w2_t]] = -κ_w
    Γ0[eq[:eq_wage_phillips_good2], endo[:π_t]]    = -(1. + m[:β] * m[:ι_w]) / (1. + m[:β])
    Γ0[eq[:eq_wage_phillips_good2], endo[:Eπ_t]]   = m[:β] / (1. + m[:β])
    Γ0[eq[:eq_wage_phillips_good2], endo[:z_t]]    = -(1. + m[:β] * m[:ι_w] - m[:ρ_z] * m[:β]) / (1. + m[:β])
    Γ0[eq[:eq_wage_phillips_good2], endo[:λ_w_t]]  = 1.
    Γ1[eq[:eq_wage_phillips_good2], endo[:w2_t]]   = -1. / (1. + m[:β])
    Γ1[eq[:eq_wage_phillips_good2], endo[:π_t]]    = -m[:ι_w] / (1. + m[:β])
    Γ1[eq[:eq_wage_phillips_good2], endo[:z_t]]    = -m[:ι_w] / (1. + m[:β])

    Γ0[eq[:eq_marg_sub_good2], endo[:g_w2_t]]  = 1.
    Γ0[eq[:eq_marg_sub_good2], endo[:w2_t]]    = -1.
    Γ0[eq[:eq_marg_sub_good2], endo[:L2_t]]    = m[:ν]
    Γ0[eq[:eq_marg_sub_good2], endo[:b_t]]     = 1.
    Γ0[eq[:eq_marg_sub_good2], endo[:λ_SH2_t]] = -1.

    Γ0[eq[:eq_avg_marg_util_SH2], endo[:λ_SH2_t]] = -1.
    Γ0[eq[:eq_avg_marg_util_SH2], endo[:λ_S_t]]   = (1. - m[:θ]) * (1. - m[:f_S1]) / ((1. - m[:θ]) *
                                                    (1. - m[:f_S1]) + m[:θ] * (1 - m[:f_H1])) *
                                                    m[:λ_S_ss] / m[:λ_SH2_ss]
    Γ0[eq[:eq_avg_marg_util_SH2], endo[:λ_H2_t]]  = m[:θ] * (1. - m[:f_H1]) / ((1. - m[:θ]) *
                                                    (1. - m[:f_S1]) + m[:θ] * (1. - m[:f_H1])) *
                                                    m[:λ_H2_ss] / m[:λ_SH2_ss]

    ### 8. Monetary and fiscal policy
    Γ0[eq[:eq_mp], endo[:R_t]]    = -1.
    Γ1[eq[:eq_mp], endo[:R_t]]    = -m[:ρ_R]
    Γ0[eq[:eq_mp], endo[:π_t]]    = (1. - m[:ρ_R]) * m[:ψ_1]
    Γ0[eq[:eq_mp], endo[:x_t]]    = (1. - m[:ρ_R]) * m[:ψ_3]
    Γ0[eq[:eq_mp], endo[:η_mp_t]] = 1.
    Γ1[eq[:eq_mp], endo[:x_t]]    = (1. - m[:ρ_R]) * m[:ψ_3]

    Γ0[eq[:eq_revenue], endo[:t_t]]    = -1.
    Γ0[eq[:eq_revenue], endo[:t_S_t]]  = (1. - m[:θ])
    Γ0[eq[:eq_revenue], endo[:t_H1_t]] = m[:θ] * m[:f_H1]
    Γ0[eq[:eq_revenue], endo[:t_H2_t]] = m[:θ] * (1. - m[:f_H1])

    Γ0[eq[:eq_transfer], endo[:τ_t]]    = -1.
    Γ0[eq[:eq_transfer], endo[:τ_S_t]]  = (1. - m[:θ])
    Γ0[eq[:eq_transfer], endo[:τ_H1_t]] = m[:θ] * m[:f_H1]
    Γ0[eq[:eq_transfer], endo[:τ_H2_t]] = m[:θ] * (1. - m[:f_H1])

    Γ0[eq[:eq_realdebt], endo[:bᴿ_t]] = m[:x_ss] / m[:bᴿ_ss]
    Γ0[eq[:eq_realdebt], endo[:π_t]]  = m[:R_ss] / (expγ * m[:π_ss])
    Γ0[eq[:eq_realdebt], endo[:z_t]]  = m[:R_ss] / (expγ * m[:π_ss])
    Γ0[eq[:eq_realdebt], endo[:g_t]]  = -(1. - m[:R_ss] / (expγ * m[:π_ss])) * m[:x_ss] / (m[:g_ss] - m[:t_ss] + m[:τ_ss])
    Γ0[eq[:eq_realdebt], endo[:t_t]]  = (1. - m[:R_ss] / (expγ * m[:π_ss])) * m[:x_ss] / (m[:g_ss] - m[:t_ss] + m[:τ_ss])
    Γ0[eq[:eq_realdebt], endo[:τ_t]]  = -(1. - m[:R_ss] / (expγ * m[:π_ss])) * m[:x_ss] / (m[:g_ss] - m[:t_ss] + m[:τ_ss])
    Γ1[eq[:eq_realdebt], endo[:R_t]]  = m[:R_ss] / (expγ * m[:π_ss])
    Γ1[eq[:eq_realdebt], endo[:bᴿ_t]] = m[:R_ss] / (expγ * m[:π_ss]) * m[:x_ss] / m[:bᴿ_ss]

    Γ0[eq[:eq_fiscal_rule], endo[:t_t]] = -m[:x_ss] / m[:t_ss]
    Γ0[eq[:eq_fiscal_rule], endo[:x_t]] = (1. + (1. - m[:ρ_T]) * m[:ζ_X] + (1. - m[:ρ_T]) * m[:ζ_G])
    Γ0[eq[:eq_fiscal_rule], endo[:g_t]] = (1. - m[:ρ_T]) * m[:ζ_G] * m[:x_ss] / m[:g_ss]
    Γ1[eq[:eq_fiscal_rule], endo[:t_t]] = -m[:ρ_T] * m[:x_ss] / m[:t_ss]
    Γ1[eq[:eq_fiscal_rule], endo[:bᴿ_t]] = -(1. - m[:ρ_T]) * m[:ζ_B]  * m[:x_ss] / m[:bᴿ_ss]
    Γ1[eq[:eq_fiscal_rule], endo[:x_t]] = m[:ρ_T] + (1. - m[:ρ_T]) * m[:ζ_B]

    Γ0[eq[:eq_redistribute_H1], endo[:t_H1_t]] = m[:θ] * m[:f_H1] # IN THE TEXT, THIS IS f_H, but I think it should be f_H1
    Γ0[eq[:eq_redistribute_H1], endo[:t_t]]    = -m[:ψ_H1]

    Γ0[eq[:eq_redistribute_H2], endo[:t_H2_t]] = m[:θ] * (1. - m[:f_H1]) # IN THE TEXT, THIS IS f_H, but I think it should be f_H1
    Γ0[eq[:eq_redistribute_H2], endo[:t_t]]    = -m[:ψ_H2]

    ### 9. Market clearing and aggregate resource constraint
    Γ0[eq[:eq_real_profits], endo[:π_S_t]]  = -1.
    Γ0[eq[:eq_real_profits], endo[:y_t]]    = m[:y_ss] / ((1. - m[:θ]) * m[:π_S_ss])
    Γ0[eq[:eq_real_profits], endo[:mc1_t]]  = -m[:mc1_ss] * m[:klR1_ss] * m[:L1] / ((1. - m[:θ]) * m[:π_S_ss])
    Γ0[eq[:eq_real_profits], endo[:klR1_t]] = -m[:mc1_ss] * m[:klR1_ss] * m[:L1] / ((1. - m[:θ]) * m[:π_S_ss])
    Γ0[eq[:eq_real_profits], endo[:L1_t]]   = -m[:mc1_ss] * m[:klR1_ss] * m[:L1] / ((1. - m[:θ]) * m[:π_S_ss])
    Γ0[eq[:eq_real_profits], endo[:mc2_t]]  = -m[:A2] ^ (1. - m[:α]) * m[:mc2_ss] * m[:klR2_ss] * m[:L2] /
                                               ((1. - m[:θ]) * m[:π_S_ss])
    Γ0[eq[:eq_real_profits], endo[:klR2_t]] = -m[:A2] ^ (1. - m[:α]) * m[:mc2_ss] * m[:klR2_ss] * m[:L2] /
                                               ((1. - m[:θ]) * m[:π_S_ss])
    Γ0[eq[:eq_real_profits], endo[:L2_t]]   = -m[:A2] ^ (1. - m[:α]) * m[:mc2_ss] * m[:klR2_ss] * m[:L2] /
                                               ((1. - m[:θ]) * m[:π_S_ss])

    Γ0[eq[:eq_cap_market_clear], endo[:klR1_t]] = m[:klR1_ss] * m[:L1] / m[:k_ss]
    Γ0[eq[:eq_cap_market_clear], endo[:L1_t]]   = m[:klR1_ss] * m[:L1] / m[:k_ss]
    Γ0[eq[:eq_cap_market_clear], endo[:klR2_t]] = m[:klR2_ss] * m[:L2] / m[:k_ss]
    Γ0[eq[:eq_cap_market_clear], endo[:L2_t]]   = m[:klR2_ss] * m[:L2] / m[:k_ss]
    Γ0[eq[:eq_cap_market_clear], endo[:k_t]]    = -1.

    Γ0[eq[:eq_goods_market_clear], endo[:c_t]]   = m[:c_ss] / m[:y_ss]
    Γ0[eq[:eq_goods_market_clear], endo[:i_S_t]] = (1 - m[:θ]) * m[:i_S_ss] / m[:y_ss]
    Γ0[eq[:eq_goods_market_clear], endo[:g_t]]   = m[:x_ss] / m[:y_ss]
    Γ0[eq[:eq_goods_market_clear], endo[:u_t]]   = m[:ρ_ss] * m[:k_ss] / m[:y_ss]
    Γ0[eq[:eq_goods_market_clear], endo[:y_t]]   = -1

    Γ0[eq[:eq_gdp], endo[:x_t]] = -1.
    Γ0[eq[:eq_gdp], endo[:y_t]] = 1.
    Γ0[eq[:eq_gdp], endo[:u_t]] = -m[:ρ_ss] * m[:k_ss] / m[:y_ss]

    ### 10. Exogenous processes

    # z
    Γ0[eq[:eq_z], endo[:z_t]] = 1
    Γ1[eq[:eq_z], endo[:z_t]] = m[:ρ_z]
    Ψ[eq[:eq_z], exo[:z_sh]]  = 1

    # a
    Γ0[eq[:eq_a], endo[:a_t]] = 1.
    Γ1[eq[:eq_a], endo[:a_t]] = m[:ρ_a]
    Ψ[eq[:eq_a],  exo[:a_sh]] = 1.

    # d
    Γ0[eq[:eq_d], endo[:d_t]] = 1.
    Γ1[eq[:eq_d], endo[:d_t]] = m[:ρ_d]
    Ψ[eq[:eq_d],  exo[:d_sh]] = 1.

    # φ
    Γ0[eq[:eq_φ], endo[:φ_t]] = 1.
    Γ1[eq[:eq_φ], endo[:φ_t]] = m[:ρ_φ]
    Ψ[eq[:eq_φ],  exo[:φ_sh]] = 1.

    # μ
    Γ0[eq[:eq_μ], endo[:μ_t]] = 1.
    Γ1[eq[:eq_μ], endo[:μ_t]] = m[:ρ_μ]
    Ψ[eq[:eq_μ],  exo[:μ_sh]] = 1.

    # b
    Γ0[eq[:eq_b], endo[:b_t]] = 1.
    Γ1[eq[:eq_b], endo[:b_t]] = m[:ρ_b]
    Ψ[eq[:eq_b],  exo[:b_sh]] = 1.

    # λ_p
    Γ0[eq[:eq_λ_p], endo[:λ_p_t]]   = 1.
    Γ1[eq[:eq_λ_p], endo[:λ_p_t]]   = m[:ρ_λ_p]
#=  Γ0[eq[:eq_λ_p], endo[:λ_p_t1]]  = -1
    Γ1[eq[:eq_λ_p], endo[:λ_p_t1]]  = -m[:η_λ_p]
    Γ0[eq[:eq_λ_p1], endo[:λ_p_t1]] = 1 =#
    Ψ[eq[:eq_λ_p], exo[:λ_p_sh]]   = 1.

    # λ_w
    Γ0[eq[:eq_λ_w], endo[:λ_w_t]] = 1.
    Γ1[eq[:eq_λ_w], endo[:λ_w_t]] = m[:ρ_λ_w]
#=  Γ0[eq[:eq_λ_w], endo[:λ_w_t1]]  = -1
    Γ1[eq[:eq_λ_w], endo[:λ_w_t1]]  = -m[:η_λ_w]
    Γ0[eq[:eq_λ_w1], endo[:λ_w_t1]] = 1 =#
    Ψ[eq[:eq_λ_w],  exo[:λ_w_sh]] = 1.

    # η_mp
    Γ0[eq[:eq_η_mp], endo[:η_mp_t]] = 1.
    Γ1[eq[:eq_η_mp], endo[:η_mp_t]] = m[:ρ_η_mp]
    Ψ[eq[:eq_η_mp],  exo[:η_mp_sh]] = 1.

    # g
    Γ0[eq[:eq_g], endo[:g_t]] = 1.
    Γ1[eq[:eq_g], endo[:g_t]] = m[:ρ_g]
    Ψ[eq[:eq_g],  exo[:g_sh]] = 1.

    # τ_S
    Γ0[eq[:eq_τ_S], endo[:τ_S_t]] = 1.
    Γ1[eq[:eq_τ_S], endo[:τ_S_t]] = m[:ρ_τ_S]
    Ψ[eq[:eq_τ_S],  exo[:τ_S_sh]] = 1.

    # τ_H1
    Γ0[eq[:eq_τ_H1], endo[:τ_H1_t]] = 1.
    Γ1[eq[:eq_τ_H1], endo[:τ_H1_t]] = m[:ρ_τ_H1]
    Ψ[eq[:eq_τ_H1],  exo[:τ_H1_sh]] = 1.

    # τ_H2
    Γ0[eq[:eq_τ_H2], endo[:τ_H2_t]] = 1.
    Γ1[eq[:eq_τ_H2], endo[:τ_H2_t]] = m[:ρ_τ_H2]
    Ψ[eq[:eq_τ_H2],  exo[:τ_H2_sh]] = 1.

    ### Expectational terms

    # π
    Γ0[eq[:eq_Eπ], endo[:π_t]]  = 1.
    Γ1[eq[:eq_Eπ], endo[:Eπ_t]] = 1.
    Π[eq[:eq_Eπ],  ex[:Eπ_sh]]  = 1.

    # π1
    Γ0[eq[:eq_Eπ1], endo[:π1_t]]  = 1.
    Γ1[eq[:eq_Eπ1], endo[:Eπ1_t]] = 1.
    Π[eq[:eq_Eπ1],  ex[:Eπ1_sh]]  = 1.

    # π2
    Γ0[eq[:eq_Eπ2], endo[:π2_t]]  = 1.
    Γ1[eq[:eq_Eπ2], endo[:Eπ2_t]] = 1.
    Π[eq[:eq_Eπ2],  ex[:Eπ2_sh]]  = 1.

    # z
    Γ0[eq[:eq_Ez], endo[:z_t]]  = 1.
    Γ1[eq[:eq_Ez], endo[:Ez_t]] = 1.
    Π[eq[:eq_Ez],  ex[:Ez_sh]]  = 1.

    # λ_S
    Γ0[eq[:eq_Eλ_S], endo[:λ_S_t]]  = 1.
    Γ1[eq[:eq_Eλ_S], endo[:Eλ_S_t]] = 1.
    Π[eq[:eq_Eλ_S],  ex[:Eλ_S_sh]]  = 1.

    # λ_H1
    Γ0[eq[:eq_Eλ_H1], endo[:λ_H1_t]]  = 1.
    Γ1[eq[:eq_Eλ_H1], endo[:Eλ_H1_t]] = 1.
    Π[eq[:eq_Eλ_H1],  ex[:Eλ_H1_sh]]  = 1.

    # λ_H2
    Γ0[eq[:eq_Eλ_H2], endo[:λ_H2_t]]  = 1.
    Γ1[eq[:eq_Eλ_H2], endo[:Eλ_H2_t]] = 1.
    Π[eq[:eq_Eλ_H2],  ex[:Eλ_H2_sh]]  = 1.

    # x
    Γ0[eq[:eq_Ex], endo[:x_t]]  = 1
    Γ1[eq[:eq_Ex], endo[:Ex_t]] = 1
    Π[eq[:eq_Ex],  ex[:Ex_sh]]  = 1

    # ϕ
    Γ0[eq[:eq_Eϕ], endo[:ϕ_t]]  = 1.
    Γ1[eq[:eq_Eϕ], endo[:Eϕ_t]] = 1.
    Π[eq[:eq_Eϕ],  ex[:Eϕ_sh]]  = 1.

    # ρ
    Γ0[eq[:eq_Eρ], endo[:ρ_t]]  = 1.
    Γ1[eq[:eq_Eρ], endo[:Eρ_t]] = 1.
    Π[eq[:eq_Eρ],  ex[:Eρ_sh]]  = 1.

    # ι_S
    Γ0[eq[:eq_Eι_S], endo[:ι_S_t]]  = 1.
    Γ1[eq[:eq_Eι_S], endo[:Eι_S_t]] = 1.
    Π[eq[:eq_Eι_S],  ex[:Eι_S_sh]]  = 1.

    # w1
    Γ0[eq[:eq_Ew1], endo[:w1_t]]  = 1.
    Γ1[eq[:eq_Ew1], endo[:Ew1_t]] = 1.
    Π[eq[:eq_Ew1],  ex[:Ew1_sh]]  = 1.

    # w2
    Γ0[eq[:eq_Ew2], endo[:w2_t]]  = 1.
    Γ1[eq[:eq_Ew2], endo[:Ew2_t]] = 1.
    Π[eq[:eq_Ew2],  ex[:Ew2_sh]]  = 1.

#=    # c
    Γ0[eq[:eq_Ec], endo[:c_t]] = 1
    Γ1[eq[:eq_Ec], endo[:Ec_t]] = 1
    Π[eq[:eq_Ec], endo[:Ec_t]] = 1

    # c_H1
    Γ0[eq[:eq_Ec_H1], endo[:c_H1_t]] = 1
    Γ1[eq[:eq_Ec_H1], endo[:Ec_H1_t]] = 1
    Π[eq[:eq_Ec_H1], endo[:Ec_H1_t]] = 1

    # c_H2
    Γ0[eq[:eq_Ec_H2], endo[:c_H2_t]] = 1
    Γ1[eq[:eq_Ec_H2], endo[:Ec_H2_t]] = 1
    Π[eq[:eq_Ec_H2], endo[:Ec_H2_t]] = 1

    # c_H2
    Γ0[eq[:eq_Ec_H2], endo[:c_H2_t]] = 1
    Γ1[eq[:eq_Ec_H2], endo[:Ec_H2_t]] = 1
    Π[eq[:eq_Ec_H2], endo[:Ec_H2_t]] = 1
=#

    return Γ0, Γ1, C, Ψ, Π
end
