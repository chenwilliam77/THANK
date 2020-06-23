"""
```
COTHANK{T} <: AbstractRepModel{T}
```

The `THANK` type defines the structure of THANK, an implementation
of the DSGE model in Bilbiie, Primiceri, Tambalotti (2019)

### Fields

#### Parameters and Steady-States
* `parameters::Vector{AbstractParameter}`: Vector of all time-invariant model
  parameters.

* `steady_state::Vector{AbstractParameter}`: Model steady-state values, computed
  as a function of elements of `parameters`.

* `keys::OrderedDict{Symbol,Int}`: Maps human-readable names for all model
  parameters and steady-states to their indices in `parameters` and
  `steady_state`.

#### Inputs to Measurement and Equilibrium Condition Equations

The following fields are dictionaries that map human-readable names to row and
column indices in the matrix representations of of the measurement equation and
equilibrium conditions.

* `endogenous_states::OrderedDict{Symbol,Int}`: Maps each state to a column in
  the measurement and equilibrium condition matrices.

* `exogenous_shocks::OrderedDict{Symbol,Int}`: Maps each shock to a column in
  the measurement and equilibrium condition matrices.

* `expected_shocks::OrderedDict{Symbol,Int}`: Maps each expected shock to a
  column in the measurement and equilibrium condition matrices.

* `equilibrium_conditions::OrderedDict{Symbol,Int}`: Maps each equlibrium
  condition to a row in the model's equilibrium condition matrices.

* `endogenous_states_augmented::OrderedDict{Symbol,Int}`: Maps lagged states to
  their columns in the measurement and equilibrium condition equations. These
  are added after `gensys` solves the model.

* `observables::OrderedDict{Symbol,Int}`: Maps each observable to a row in the
  model's measurement equation matrices.

* `pseudo_observables::OrderedDict{Symbol,Int}`: Maps each pseudo-observable to
  a row in the model's pseudo-measurement equation matrices.

#### Model Specifications and Settings

* `spec::String`: The model specification identifier, \"m1002\", cached here for
  filepath computation.

* `subspec::String`: The model subspecification number, indicating that
  some parameters from the original model spec (\"ss10\") are initialized
  differently. Cached here for filepath computation.

* `settings::Dict{Symbol,Setting}`: Settings/flags that affect computation
  without changing the economic or mathematical setup of the model.

* `test_settings::Dict{Symbol,Setting}`: Settings/flags for testing mode

#### Other Fields

* `rng::MersenneTwister`: Random number generator. Can be is seeded to ensure
  reproducibility in algorithms that involve randomness (such as
  Metropolis-Hastings).

* `testing::Bool`: Indicates whether the model is in testing mode. If `true`,
  settings from `m.test_settings` are used in place of those in `m.settings`.

* `observable_mappings::OrderedDict{Symbol,Observable}`: A dictionary that
  stores data sources, series mnemonics, and transformations to/from model
  units. DSGE.jl will fetch data from the Federal Reserve Bank of St. Louis's
  FRED database; all other data must be downloaded by the user. See `load_data`
  and `Observable` for further details.

* `pseudo_observable_mappings::OrderedDict{Symbol,PseudoObservable}`: A
  dictionary that stores names and transformations to/from model units. See
  `PseudoObservable` for further details.
"""
mutable struct COTHANK{T} <: AbstractRepModel{T}
    parameters::ParameterVector{T}                         # vector of all time-invariant model parameters
    steady_state::ParameterVector{T}                       # model steady-state values
    keys::OrderedDict{Symbol,Int}                          # human-readable names for all the model
                                                           # parameters and steady-states

    endogenous_states::OrderedDict{Symbol,Int}             # these fields used to create matrices in the
    exogenous_shocks::OrderedDict{Symbol,Int}              # measurement and equilibrium condition equations.
    expected_shocks::OrderedDict{Symbol,Int}               #
    equilibrium_conditions::OrderedDict{Symbol,Int}        #
    endogenous_states_augmented::OrderedDict{Symbol,Int}   #
    observables::OrderedDict{Symbol,Int}                   #
    pseudo_observables::OrderedDict{Symbol,Int}            #

    spec::String                                           # Model specification number (eg "m990")
    subspec::String                                        # Model subspecification (eg "ss0")
    settings::Dict{Symbol,Setting}                         # Settings/flags for computation
    test_settings::Dict{Symbol,Setting}                    # Settings/flags for testing mode
    rng::MersenneTwister                                   # Random number generator
    testing::Bool                                          # Whether we are in testing mode or not

    observable_mappings::OrderedDict{Symbol, Observable}
    pseudo_observable_mappings::OrderedDict{Symbol, PseudoObservable}
end

description(m::COTHANK) = "Bilbiie, Primiceri, Tambalotti (2019), $(m.subspec)."

"""
`init_model_indices!(m::COTHANK)`

Arguments:
`m:: COTHANK`: a model object

Description:
Initializes indices for all of `m`'s states, shocks, and equilibrium conditions.
"""
function init_model_indices!(m::COTHANK)
    # Endogenous states
    endogenous_states = [:y_t, :y1_t, :y2_t, :k_t, :k_S_t, :L1_t,  :L2_t, :klR1_t, :klR2_t,   # check if missing ι terms
                         :w1_t, :w2_t, :π_t, :π_S_t,  :π1_t, :π2_t, :mc1_t, :mc2_t, :λ_p_t, :λ_w_t,      # sticky prices
                         :λ_S_t, :λ_H1_t, :λ_H2_t, :λ_SH1_t, :λ_SH2_t, :c_t,  :c_H1_t, :c_H2_t, :c_S_t, :i_S_t,
                         :R_t, :u_t, :ϕ_t, :z_t, :ρ_t, :a_t, :d_t, :x_t, :b_t, :t_t, :φ_t,
                         :g_t, :g_w1_t, :g_w2_t, :μ_t, :bᴿ_t, :η_mp_t, :t_S_t, :t_H1_t, :t_H2_t,
                         :τ_t, :τ_S_t, :τ_H1_t, :τ_H2_t, :Eπ_t, :Eπ1_t, :Eπ2_t, :Ez_t, :Eλ_S_t,
                         :Eλ_H1_t, :Eλ_H2_t, :Ex_t, :Eϕ_t, :Eρ_t, :Ei_S_t, :Ew1_t, :Ew2_t]
                                                                                # CHECK BACK IF ANY CAN BE AUGMENTED VARIABLES
#=                         :y_f_t, :k_f_t, :L_f_t, :Rk_f_t, :w_f_t, :mc_f_t,       # flexible prices
                         :λ_f_t, :c_f_t, :R_f_t, :u_f_t, :ϕ_f_t,
                         :i_f_t, :kbar_f_t, :wgap_f_t, :gdp_f_t, :Ec_f_t, :Eλ_f_t,
                         :Eϕ_f_t, :ERk_f_t, :Ei_f_t, :λ_h_f_t, :λ_s_f_t, :c_h_f_t, :c_s_f_t,
                         :t_h_f_t, :Ey_f_t, :Eλ_h_f_t, :Eλ_s_f_t]
=#
    # Exogenous shocks
    exogenous_shocks = [:z_sh, :a_sh, :d_sh, :φ_sh, :μ_sh, :b_sh, :λ_p_sh, :λ_w_sh, :η_mp_sh, :g_sh,
                        :τ_S_sh, :τ_H1_sh, :τ_H2_sh]

    # Expectations shocks
    expected_shocks = [:Eπ_sh, :Eπ1_sh, :Eπ2_sh, :Ez_sh, :Eλ_S_sh, :Eλ_H1_sh, :Eλ_H2_sh,
                       :Ex_sh, :Eϕ_sh, :Eρ_sh, :Ei_S_sh, :Ew1_sh, :Ew2_sh] # sticky prices

    # Equilibrium conditions
    equilibrium_conditions = [:eq_KL_ratio_good1, :eq_mc_good1, :eq_KL_ratio_good2, :eq_mc_good2, :eq_production_good1, :eq_production_good2,
                              :eq_production_finalgood, :eq_demand_good1, :eq_price_index, :eq_price_phillips_good1,
                              :eq_price_phillips_good2, :eq_marg_utility_S, :eq_marg_utility_H1, :eq_marg_utility_H2, :eq_average_c,
                              :eq_euler, :eq_c_H1, :eq_c_H2, :eq_cap_util, :eq_ϕ, :eq_i, :eq_cap_input, :eq_cap_accum, :eq_wage_phillips_good1,
                              :eq_marg_sub_good1, :eq_avg_marg_util_SH1, :eq_wage_phillips_good2, :eq_marg_sub_good2, :eq_avg_marg_util_SH2,
                              :eq_mp, :eq_revenue, :eq_transfer, :eq_realdebt, :eq_fiscal_rule, :eq_redistribute_H1, :eq_redistribute_H2,
                              :eq_real_profits, :eq_cap_market_clear, :eq_goods_market_clear, :eq_gdp,
                              :eq_z, :eq_a, :eq_d, :eq_φ, :eq_μ, :eq_b, :eq_λ_p, :eq_λ_w, :eq_η_mp, :eq_g, :eq_τ_S, :eq_τ_H1, :eq_τ_H2,
                              :eq_Eπ, :eq_Eπ1, :eq_Eπ2, :eq_Ez, :eq_Eλ_S, :eq_Eλ_H1, :eq_Eλ_H2, :eq_Ex, :eq_Eϕ,
                              :eq_Eρ, :eq_Ei_S, :eq_Ew1, :eq_Ew2]
#=
                              :eq_y_f, :eq_k_f, :eq_L_f, :eq_Rk_f, :eq_w_f, :eq_mc_f,       # flexible prices
                              :eq_λ_f, :eq_c_f, :eq_R_f, :eq_u_f, :eq_ϕ_f,
                              :eq_i_f, :eq_kbar_f, :eq_wgap_f, :eq_gdp_f, :eq_Ec_f, :eq_Eλ_f,
                              :eq_Eϕ_f, :eq_ERk_f, :eq_Ei_f, :eq_λ_h_f, :eq_λ_s_f, :eq_c_h_f,
                              :eq_c_s_f, :eq_t_h_f, :eq_Ey_f, :eq_Eλ_h_f, :eq_Eλ_s_f]
=#
    # Additional states added after solving model
    # Lagged states and observables measurement error
    endogenous_states_augmented = []

    # Observables
    observables = keys(m.observable_mappings)

    # Pseudo-observables
    pseudo_observables = keys(m.pseudo_observable_mappings)

    for (i,k) in enumerate(endogenous_states);           m.endogenous_states[k]           = i end
    for (i,k) in enumerate(exogenous_shocks);            m.exogenous_shocks[k]            = i end
    for (i,k) in enumerate(expected_shocks);             m.expected_shocks[k]             = i end
    for (i,k) in enumerate(equilibrium_conditions);      m.equilibrium_conditions[k]      = i end
    for (i,k) in enumerate(endogenous_states);           m.endogenous_states[k]           = i end
    for (i,k) in enumerate(endogenous_states_augmented); m.endogenous_states_augmented[k] = i+length(endogenous_states) end
    for (i,k) in enumerate(observables);                 m.observables[k]                 = i end
    for (i,k) in enumerate(pseudo_observables);          m.pseudo_observables[k]          = i end
end

function COTHANK(subspec::String="ss1";
                   custom_settings::Dict{Symbol, Setting} = Dict{Symbol, Setting}(),
                   testing = false)

    # Model-specific specifications
    spec               = split(basename(@__FILE__),'.')[1]
    subspec            = subspec
    settings           = Dict{Symbol,Setting}()
    test_settings      = Dict{Symbol,Setting}()
    rng                = MersenneTwister(0)

    # Initialize empty model
    m = COTHANK{Float64}(
            # model parameters and steady state values
            Vector{AbstractParameter{Float64}}(), Vector{Float64}(), OrderedDict{Symbol,Int}(),

            # model indices
            OrderedDict{Symbol,Int}(), OrderedDict{Symbol,Int}(), OrderedDict{Symbol,Int}(), OrderedDict{Symbol,Int}(), OrderedDict{Symbol,Int}(), OrderedDict{Symbol,Int}(), OrderedDict{Symbol,Int}(),

            spec,
            subspec,
            settings,
            test_settings,
            rng,
            testing,
            OrderedDict{Symbol,Observable}(),
            OrderedDict{Symbol,PseudoObservable}())

    # Set settings
    model_settings!(m)
    DSGE.default_test_settings!(m)
    for custom_setting in values(custom_settings)
        m <= custom_setting
    end

    # Set observable and pseudo-observable transformations
    init_observable_mappings!(m)
    #init_pseudo_observable_mappings!(m)

    # Initialize parameters
    init_parameters!(m)

    init_model_indices!(m)
    init_subspec!(m)
    steadystate!(m)

    return m
end

"""
```
init_parameters!(m::COTHANK)
```

Initializes the model's parameters, as well as empty values for the steady-state
parameters (in preparation for `steadystate!(m)` being called to initialize
those).
"""
function init_parameters!(m::COTHANK)
    # Calibrated parameters#
    m <= parameter(:δ, 0.025, fixed=true,
                   description="δ: The capital depreciation rate.",
                   tex_label = "\\delta" ) # from Giorgio's matlab script

    # Non standard devation parameters
    m <= parameter(:α, 1-0.6/(1.15), (0., 1.), (1e-5, 0.999), ModelConstructors.SquareRoot(), Normal(0.30, 0.05), fixed = false,
                   description="α: Capital elasticity in the intermediate goods sector's production function (also known as the capital share).",
                   tex_label = "\\alpha") # from Giorgio's matlab script

    m <= parameter(:ι_p, 0.21, (0., 1.), (1e-5, 0.999), ModelConstructors.SquareRoot(), Beta(0.5, 0.15), fixed = false,
                   description = "ι_p: The weight attributed to last period's inflation in price indexation. " *
                   "(1-ι_p) is the weight attributed to steady-state inflation.",
                   tex_label = "\\iota_p") # THANK

    m <= parameter(:ι_w, 0.15, (0., 1.), (1e-5, 0.999), ModelConstructors.SquareRoot(), Beta(0.5, 0.15), fixed = false,
                   description="ι_w: The weight attributed to last period's wage in wage indexation. "
                    * "(1-ι_w) is the weight attributed to steady-state wages.",
                   tex_label = "\\iota_w") # THANK

    m <= parameter(:h, 0.5, (-5.0, 5.0), (-5., 5.), ModelConstructors.Untransformed(), Beta(0.5, 0.1), fixed = false,
                   description="h: habit formation parameter.",
                   tex_label = "h") # from Giorgio's matlab script

    m <= parameter(:λ_p_ss, 0.15,(-1000.,1000.),(-1000., 1000.), ModelConstructors.Untransformed(), Normal(0.15, 0.05),  fixed = false,
                   description="λ_p_ss: The steady state net price markup.",
                   tex_label = "\\lambda_p") # from Giorgio's matlab script

    m <= parameter(:λ_w_ss, 0.12,(-1000.,1000.),(-1000., 1000.), ModelConstructors.Untransformed(), Normal(0.15, 0.05), fixed = false,
                   description = "λ_w_ss: The steady state net wage markup, which affects the elasticity of substitution between differentiated labor services.",
                   tex_label = "\\lambda_w") # from Giorgio's matlab script

    m <= parameter(:L_ss, 0.5, (-1000., 1000.), (-1e3, 1e3), ModelConstructors.Untransformed(), Normal(0.0,0.5), fixed=false,
                   description="L_ss: The steady state for log hours.",
                   tex_label = "\\log(L)^{ss}") # from Giogio's matlab script

    m <= parameter(:Fβ, 0.1, (1e-5, 10.), (1e-5, 10.), ModelConstructors.Exponential(), Gamma(0.25, 0.1), fixed=false,
                   description = "Fβ: Discount rate transformed.",

                   tex_label = "100(\\beta^{-1} - 1)") # from THANK

    m <= parameter(:ν, 4., (1e-5, 10.), (1e-5, 10.), ModelConstructors.Exponential(), Gamma(2, 0.75), fixed = false,
                   description="ν_l: The inverse Frisch elasticity.",
                   tex_label = "\\nu") # from thank.jl

    m <= parameter(:ξ_p, 0.85, (1e-5, 0.999), (1e-5, 0.999), ModelConstructors.SquareRoot(), Beta(0.66, 0.1), fixed=false,
                   description="ξ_p: The Calvo parameter. In every period, intermediate goods producers optimize prices with probability (1-ξ_p). "
                   * "With probability ξ_p, prices are adjusted according to a weighted average of the previous period's inflation "
                   * "(π_t1) and steady-state inflation (π_ss).",
                   tex_label = "\\xi_p") # from THANK

    m <= parameter(:ξ_w, 0.75, (1e-5, 0.999), (1e-5, 0.999), ModelConstructors.SquareRoot(), Beta(0.66, 0.1), fixed=false,
                   description="ξ_w: (1-ξ_w) is the probability with which households can freely choose wages in each period. "
                   * "With probability ξ_w, wages increase at a geometrically weighted average of the steady state rate of "
                   * "wage increases and last period's productivity times last period's inflation.",
                   tex_label = "\\xi_w") # from THANK

    m <= parameter(:χ, 5., (0., 10.), (1e-5, 0.), ModelConstructors.Exponential(), Gamma(5., 1.), fixed = false,
                   description="χ: The elasticity of the capital utilization cost function.",
                   tex_label = "\\chi") # from THANK

    m <= parameter(:S′′, 2.5, (-15., 15.), (-15., 15.), ModelConstructors.Untransformed(), Gamma(4., 1), fixed = false,
                   description="S′′: The investment adjust cost.",
                   tex_label = "S^{\\prime\\prime}") # from THANK

    m <= parameter(:ψ_1, 2., (1e-5, 10.), (1e-5, 10.00), ModelConstructors.Exponential(), Normal(1.7, 0.3), fixed = false,
                   description="ψ₁: Weight on inflation gap in monetary policy rule.",
                   tex_label = "\\psi_1") # from THANK

    m <= parameter(:ψ_2, 0.1, (-0.5, 0.5), (-0.5, 0.5), ModelConstructors.Untransformed(), Normal(0.125, 0.05), fixed = false,
                   description="ψ₂: Weight on output gap in monetary policy rule.",
                   tex_label = "\\psi_2") # from THANK

    m <= parameter(:ψ_3, 0.25, (-0.5, 0.5), (-0.5, 0.5), ModelConstructors.Untransformed(), Normal(0.125, 0.05), fixed = false,

                   description="ψ₃: Weight on ouptut growth in the monetary policy rule.",
                   tex_label = "\\psi_3") # from THANK

    # new parameters (from θ to v), all of these are inputs into the steady state calculations given by Giorgio's matlab script
    m <= parameter(:θ, 0.25, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Untransformed(), RootInverseGamma(0.5, 1), fixed = true) # from Giorgio's MATLAB script
    m <= parameter(:τ_x, 0.0, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Untransformed(), RootInverseGamma(0.5, 1), fixed = true) # from Giorgio's MATLAB script
    m <= parameter(:t_x, 0.2, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Untransformed(), RootInverseGamma(0.5, 1), fixed = true) # from Giorgio's MATLAB script

    m <= parameter(:f_H1, 1/4, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Untransformed(), RootInverseGamma(0.5, 1), fixed = true) #from Giorgio's MATLAB script
    m <= parameter(:f_S1, 3/4, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Untransformed(), RootInverseGamma(0.5, 1), fixed = true) # from Giorgio's MATLAB script
    m <= parameter(:H1, 1.0, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Untransformed(), RootInverseGamma(0.5, 1), fixed = true) # from Giorgio's MATLAB script
    m <= parameter(:H2, 1.0, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Untransformed(), RootInverseGamma(0.5, 1), fixed = true) # from Giorgio's MATLAB script
    m <= parameter(:A2, 0.7, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Untransformed(), RootInverseGamma(0.5, 1), fixed = true) # from Giorgio's MATLAB script
    m <= parameter(:s, 0.995, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Untransformed(), RootInverseGamma(0.5, 1), fixed = true) # from Giorgio's MATLAB script
    m <= parameter(:g_x, 0.18, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Untransformed(), RootInverseGamma(0.5, 1), fixed = true) # from Giorgio's MATLAB script
    m <= parameter(:v, 0.7, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Untransformed(), RootInverseGamma(0.5, 1), fixed = true)

    # the following parameters (from ζ to sx) are new but had no given value, so current value is random placeholder
    m <= parameter(:ζ, 0.75, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Untransformed(), RootInverseGamma(0.5, 1), fixed = true) # ???
         # the above parameter appears as \varzeta in the final good producers equilibrium conditions
         # but no \varzeta unicode characters exists in Julia so we use \zeta instead
    m <= parameter(:ζ_X, 0.75, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Untransformed(), RootInverseGamma(0.5, 1), fixed = true) # ???
    m <= parameter(:ζ_G, 0.75, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Untransformed(), RootInverseGamma(0.5, 1), fixed = true) # ???
    m <= parameter(:ζ_B, 0.75, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Untransformed(), RootInverseGamma(0.5, 1), fixed = true) # ???
    m <= parameter(:s_x, 0.18, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Untransformed(), RootInverseGamma(0.5, 1), fixed = true) # ???
    m <= parameter(:sx, 0.18, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Untransformed(), RootInverseGamma(0.5, 1), fixed = true) # ???

    # exogenous processes - autocorrelation
    m <= parameter(:ρ_R, 0.8, (0., 1.), (1e-5, 0.999), ModelConstructors.SquareRoot(), Beta(0.6, 0.2), fixed = false,
                   description="ρ_R: persistence in the monetary policy rule.",
                   tex_label = "\\rho_{R}") # from THANK

    m <= parameter(:ρ_η_mp, 0.8, (0., 1.), (1e-5, 0.999), ModelConstructors.SquareRoot(), Beta(0.6, 0.2), fixed = false,
                   description="ρ_R: persistence in the monetary policy rule.",
                   tex_label = "\\rho_{R}") # from THANK

    m <= parameter(:ρ_z, 0.25, (1e-5, 0.999), (1e-5, 0.999), ModelConstructors.SquareRoot(), Beta(0.6, 0.2), fixed = false,
                   description="ρ_z: AR(1) coefficient in the technology process.",
                   tex_label = "\\rho_{z}") # from THANK

    m <= parameter(:ρ_a, 0.15, (0., 1.), (1e-5, 0.999), ModelConstructors.SquareRoot(), Beta(0.4, 0.2), fixed = false,
                   description="ρ_a: AR(1) coefficient", # fix description
                   tex_label = "\\rho_{a}") # ???

    m <= parameter(:ρ_d, 0.15, (0., 1.), (1e-5, 0.999), ModelConstructors.SquareRoot(), Beta(0.4, 0.2), fixed = false,
                   description="ρ_d: AR(1) coefficient", # fix description
                   tex_label = "\\rho_{d}") # ???

    m <= parameter(:ρ_φ, 0.15, (0., 1.), (1e-5, 0.999), ModelConstructors.SquareRoot(), Beta(0.4, 0.2), fixed = false,
                   description="ρ_φ: AR(1) coefficient.", # fix description
                   tex_label = "\\rho_{\\varphi}") # ???

    m <= parameter(:ρ_μ, 0.7, (1e-5, 0.999), (1e-5, 0.999), ModelConstructors.SquareRoot(), Beta(0.6, 0.2), fixed = false,
                   description="ρ_μ: AR(1) coefficient in capital adjustment cost process.",
                   tex_label = "\\rho_{\\mu}") # from THANK

    m <= parameter(:ρ_b, 0.7, (1e-5, 0.999), (1e-5, 0.999), ModelConstructors.SquareRoot(), Beta(0.6, 0.2), fixed = false,
                   description="ρ_b: AR(1) coefficient in the intertemporal preference shifter process.",
                   tex_label = "\\rho_b") # from THANK

    m <= parameter(:ρ_λ_p, 0.95, (1e-5, 0.999), (1e-5, 0.999), ModelConstructors.SquareRoot(), Beta(0.6, 0.2), fixed = false,
                   description="ρ_λ_p: AR(1) coefficient in the price mark-up shock process.",
                   tex_label = "\\rho_{\\lambda_p}") # from THANK

    m <= parameter(:ρ_λ_w, 0.98, (1e-5, 0.999), (1e-5, 0.999), ModelConstructors.SquareRoot(), Beta(0.6, 0.2), fixed = false,
                   description="ρ_λ_w: AR(1) coefficient in the wage mark-up shock process.",
                   tex_label = "\\rho_{\\lambda_w}") # from THANK

    m <= parameter(:ρ_g, 0.98, (1e-5, 0.999), (1e-5, 0.999), ModelConstructors.SquareRoot(), Beta(0.6, 0.2), fixed = false,
                   description="ρ_g: AR(1) coefficient in the government spending process.",
                   tex_label = "\\rho_g") # from THANK

    m <= parameter(:ρ_τ_S, 0.98, (1e-5, 0.999), (1e-5, 0.999), ModelConstructors.SquareRoot(), Beta(0.6, 0.2), fixed = false,
                   description="ρ_τ_S: AR(1) coefficient.", # fix description
                   tex_label = "\\rho_\\tau_S") # ???

    m <= parameter(:ρ_τ_H1, 0.98, (1e-5, 0.999), (1e-5, 0.999), ModelConstructors.SquareRoot(), Beta(0.6, 0.2), fixed = false,
                   description="ρ_τ_H1: AR(1) coefficient.", # fix description
                   tex_label = "\\rho_\\tau_H1") # ???

    m <= parameter(:ρ_τ_H2, 0.98, (1e-5, 0.999), (1e-5, 0.999), ModelConstructors.SquareRoot(), Beta(0.6, 0.2), fixed = false,
                   description="ρ_τ_H2: AR(1) coefficient.", # fix description
                   tex_label = "\\rho_g") # ???

    m <= parameter(:ρ_T, 0.8, (0., 1.), (1e-5, 0.999), ModelConstructors.SquareRoot(), Beta(0.6, 0.2), fixed = false,
                   description="ρ_T: parameter governing fiscal rule.", # Update this description later
                   tex_label = "\\rho_{R}") # ???

    m <= parameter(:η_λ_p, 0.75, (1e-5, 0.999), (1e-5, 0.999), ModelConstructors.SquareRoot(), Beta(0.50, 0.20), fixed = false,
                   description="η_λ_p: Moving average component in the price markup shock.",
                   tex_label = "\\eta_{\\lambda_p}") # from THANK

    m <= parameter(:η_λ_w, 0.95, (1e-5, 0.999), (1e-5, 0.999), ModelConstructors.SquareRoot(), Beta(0.50, 0.20), fixed = false,
                   description="η_λ_w: Moving average component in the wage markup shock.",
                   tex_label = "\\eta_{\\lambda_w}") # from THANK

    # exogenous processes - standard deviation
    m <= parameter(:σ_z, 0.9, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(0.5, 1), fixed = false,
                   description="σ_z: The standard deviation of the technology process.",
                   tex_label = "\\sigma_{z}") # from THANK

    m <= parameter(:σ_a, 0.9, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(0.5, 1), fixed = false,
                   description="σ_a: The standard deviation of the COVID-19 shock",
                   tex_label = "\\sigma_{a}") # ???

    m <= parameter(:σ_d, 0.9, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(0.5, 1), fixed = false,
                   description="σ_d: ", # fix description
                   tex_label = "\\sigma_{d}") # ???

    m <= parameter(:σ_φ, 4., (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(0.5, 1), fixed = false,
                   description="σ_φ: The standard deviation of the shock to disutility from labor supply", # fix description
                   tex_label = "\\sigma_{\\mu}") # ???, set to 4, following some values used in literature (e.g. Del Negro, Smets, and Schorfheide (2007) (DSGE-VAR paper))

    m <= parameter(:σ_μ, 5., (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(0.5, 1), fixed = false,
                   description="σ_μ: The standard deviation of the exogenous marginal efficiency of investment shock process.",
                   tex_label = "\\sigma_{\\mu}") # from THANK

    m <= parameter(:σ_b, 0.05, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(0.1, 1), fixed = false,
                   description="σ_b: The standard deviation of the intertemporal preference shifter process.",
                   tex_label = "\\sigma_{b}") # from THANK

    m <= parameter(:σ_λ_p, 0.15, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(0.1, 1), fixed = false,
                   description="σ_λ_p: The standard deviation of the price markup shock",
                   tex_label = "\\sigma_{\\lambda_p}") # from THANK

    m <= parameter(:σ_λ_w, 0.2, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(0.1, 1), fixed = false,
                   description="σ_λ_w: The standard deviation of the wage markup shock",
                   tex_label = "\\sigma_{\\lambda_w}") # from THANK

    m <= parameter(:σ_η_mp, 0.238, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(2, 0.10), fixed = false,
                   description="σ_η_mp: The standard deviation of the monetary policy shock.", # fix description
                   tex_label = "\\sigma_{\\eta_mp}") # ???, set to a reasonable calibration from the FRBNY DSGE model

    m <= parameter(:σ_g, 0.35, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(0.5, 1), fixed = false,
                   description="σ_g: The standard deviation of the government spending process.",
                   tex_label = "\\sigma_{g}") # from THANK

    m <= parameter(:σ_τ_S, 0.0, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(0.5, 1), fixed = false,
                   description="σ_τ_S.", # fix description
                   tex_label = "\\sigma_{\\tau_S}") # ???, zeroed out for now

    m <= parameter(:σ_τ_H1, 0.0, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(0.5, 1), fixed = false,
                   description="σ_τ_H1:.", # fix description
                   tex_label = "\\sigma_{\\tau_{H1}}") # ???, zeroed out for now

    m <= parameter(:σ_τ_H2, 0.0, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(0.5, 1), fixed = false,
                   description="σ_τ_H2: .", # fix description
                   tex_label = "\\sigma_{\\tau_{H2}}") # ???, zeroed out for now




    # steady states
    m <= SteadyStateParameter(:γ, NaN, tex_label = "\\gamma")
    m <= SteadyStateParameter(:β, NaN, tex_label = "\\beta")
    m <= SteadyStateParameter(:r_ss, NaN, tex_label = "r_ss")
    m <= SteadyStateParameter(:π_ss, NaN, tex_label = "\\pi_ss")
     m <= SteadyStateParameter(:expγ, NaN, tex_label = "exp^\\gamma")
    m <= SteadyStateParameter(:ψ_H1, NaN, tex_label = "\\psi_H1")
    m <= SteadyStateParameter(:ψ_H2, NaN, tex_label = "\\psi_H2")
    m <= SteadyStateParameter(:L1, NaN, tex_label = "L1")
    m <= SteadyStateParameter(:L2, NaN, tex_label = "L2")
    m <= SteadyStateParameter(:ρ_ss, NaN, tex_label = "\\rho")
    m <= SteadyStateParameter(:w1_ss, NaN, tex_label = "w1")
    m <= SteadyStateParameter(:w2_ss, NaN, tex_label = "w2")
    m <= SteadyStateParameter(:klR1_ss, NaN, tex_label = "klr1_{ss}")
    m <= SteadyStateParameter(:klR2_ss, NaN, tex_label = "klr2_{ss}")
    m <= SteadyStateParameter(:k_ss, NaN, tex_label = "k_{ss}")
    m <= SteadyStateParameter(:F1_ss, NaN, tex_label = "F1_{ss}")
    m <= SteadyStateParameter(:y1_ss, NaN, tex_label = "y1_ss")
    m <= SteadyStateParameter(:y_ss, NaN, tex_label = "y_{ss}")
    m <= SteadyStateParameter(:y2_ss, NaN, tex_label = "y2_{ss}")
    m <= SteadyStateParameter(:F2_ss, NaN, tex_label = "F2_{ss}")
    m <= SteadyStateParameter(:x_ss, NaN, tex_label = "x_{ss}")
    m <= SteadyStateParameter(:i_S_ss, NaN, tex_label = "i_s_{ss}")
    m <= SteadyStateParameter(:g_ss, NaN, tex_label = "g_{ss}")
    m <= SteadyStateParameter(:c_ss, NaN, tex_label = "c_{ss}")
    m <= SteadyStateParameter(:τ_ss, NaN, tex_label = "\\tau_{ss}")
    m <= SteadyStateParameter(:τ_H1_ss, NaN, tex_label = "\\tau_h1_{ss}")
    m <= SteadyStateParameter(:τ_H2_ss, NaN, tex_label = "\\tau_H2_{ss}")
    m <= SteadyStateParameter(:τ_s_ss, NaN, tex_label = "\\tau_s_{ss}")

    m <= SteadyStateParameter(:t_ss, NaN, tex_label = "t_{ss}")
    m <= SteadyStateParameter(:t_H1_ss, NaN, tex_label = "t_H1_{ss}")
    m <= SteadyStateParameter(:t_H2_ss, NaN, tex_label = "t_H2_{ss}")
    m <= SteadyStateParameter(:t_S_ss, NaN, tex_label = "t_s_{ss}")

    m <= SteadyStateParameter(:c_H1_ss, NaN, tex_label = "c_h10_{ss}")
    m <= SteadyStateParameter(:c_H2_ss, NaN, tex_label = "c_h20_{ss}")
    m <= SteadyStateParameter(:c_S_ss, NaN, tex_label = "c_s0_{ss}")
    m <= SteadyStateParameter(:λ_H1_ss, NaN, tex_label = "\\lambda_h10_{ss}")
    m <= SteadyStateParameter(:λ_H2_ss, NaN, tex_label = "\\lambda_h20_{ss}")
    m <= SteadyStateParameter(:λ_S_ss, NaN, tex_label = "\\lambda_s0_{ss}")
    m <= SteadyStateParameter(:λ_SH1_ss, NaN, tex_label = "\\lambda_sh1_{ss}")
    m <= SteadyStateParameter(:λ_SH2_ss, NaN, tex_label = "\\lambda_sh1_{ss}")
    m <= SteadyStateParameter(:bᴿ_ss, NaN, tex_label = "bR0_{ss}")
    m <= SteadyStateParameter(:R_ss, NaN, tex_label = "R0_{ss}")
    m <= SteadyStateParameter(:mc1_ss, NaN, tex_label = "mc1_{ss}")
    m <= SteadyStateParameter(:mc2_ss, NaN, tex_label = "mc2_{ss}")

end

"""
```
steadystate!(m::COTHANK)
```

Calculates the model's steady-state values. `steadystate!(m)` must be called whenever
the parameters of `m` are updated.
"""
function steadystate!(m::COTHANK)

    # These steady state values are inputs given by Giorgio's matlab script
    m[:γ]        = log(1.02)/4
    m[:expγ]     = exp(m[:γ])
    m[:β]        = 0.99
    m[:π_ss]     = 1.02^(1/4)
    m[:r_ss]     = exp(m[:γ]) / m[:β] - 1.

    # Lines 553-595 are taken from Giorgio's MATLAB script for solving steady state, steady states that we
    # are unsure about are in lines 598-604
    m[:ψ_H1]     = m[:θ] * m[:f_H1] / 1.1 # from Giorgio's MATLAB script for solving steady state
    m[:ψ_H2]     = m[:θ] * (1 - m[:f_H1]) / 1.1 # from Giorgio's MATLAB script for solving steady state

    m[:L1]   = ((1-m[:θ]) * m[:f_S1] + m[:θ] * m[:f_H1]) * m[:H1] # from Giorgio's MATLAB script for solving steady state
    m[:L2]   = ((1-m[:θ]) * (1 - m[:f_S1]) + m[:θ] * (1 - m[:f_H1])) *m[:H2] # from Giorgio's MATLAB script for solving steady state

    m[:ρ_ss]    = exp(m[:γ]) / m[:β] - (1 - m[:δ]) # from Giorgio's MATLAB script for solving steady state

    m[:w1_ss]   = (m[:α] ^ m[:α] * (1 - m[:α]) ^ (1 - m[:α]) / (1 + m[:λ_p_ss]) / (m[:ρ_ss] ^ m[:α])) ^ (1 / (1 - m[:α])) # from Giorgio's MATLAB script for solving steady state
    m[:w2_ss]   = m[:A2] * m[:w1_ss] # from Giorgio's MATLAB script for solving steady state

    m[:klR1_ss] = m[:w1_ss] / m[:ρ_ss] * m[:α] / (1 - m[:α]) # from Giorgio's MATLAB script for solving steady state
    m[:klR2_ss] = m[:A2] * m[:klR1_ss] # from Giorgio's MATLAB script for solving steady state
    m[:k_ss]    = m[:klR1_ss] * m[:L1] # from Giorgio's MATLAB script for solving steady state

    m[:F1_ss]   = m[:λ_p_ss] / (1 + m[:λ_p_ss]) * m[:L1] * m[:klR1_ss] ^ m[:α] # from Giorgio's MATLAB script for solving steady state
    m[:y1_ss]   = m[:L1] * m[:klR1_ss] ^ m[:α] - m[:F1_ss] # from Giorgio's MATLAB script for solving steady state
    m[:y_ss]    = m[:y1_ss] /m[:v] # from Giorgio's MATLAB script for solving steady state
    m[:y2_ss]   = (1 - m[:v]) * m[:y_ss] # from Giorgio's MATLAB script for solving steady state
    m[:F2_ss]   = m[:λ_p_ss] * m[:y2_ss] / m[:A2] # from Giorgio's MATLAB script for solving steady state
    m[:x_ss]    = m[:y_ss] # from Giorgio's MATLAB script for solving steady state

    m[:i_S_ss]   = (1 - (1 - m[:δ]) * exp(-m[:γ])) * exp(m[:γ]) * m[:k_ss] / (1 - m[:θ]) # from Giorgio's MATLAB script for solving steady state
    m[:g_ss]    = m[:g_x] * m[:x_ss] # from Giorgio's MATLAB script for solving steady state
    m[:c_ss]    = m[:y_ss] - m[:g_ss] - (1 - m[:θ]) * m[:i_S_ss] # from Giorgio's MATLAB script for solving steady state

    m[:τ_ss]    = m[:τ_x] * m[:x_ss] # from Giorgio's MATLAB script for solving steady state
    m[:τ_H1_ss] = m[:τ_ss] * m[:ψ_H1] / m[:θ] / (1 - m[:f_H1]) # from Giorgio's MATLAB script for solving steady state
    m[:τ_H2_ss] = m[:τ_ss] *m[:ψ_H2] / m[:θ] / (1 - m[:f_H1]) # from Giorgio's MATLAB script for solving steady state
    m[:τ_s_ss]  = m[:τ_ss] - m[:θ] * (m[:f_H1] * m[:τ_H1_ss] + (1 - m[:f_H1]) * m[:τ_H2_ss]) # from Giorgio's MATLAB script for solving steady state

    m[:t_ss]    = m[:t_x] * m[:x_ss] # from Giorgio's MATLAB script for solving steady state
    m[:t_H1_ss] = m[:t_ss] * m[:ψ_H1] / m[:θ] / m[:f_H1] # from Giorgio's MATLAB script for solving steady state
    m[:t_H2_ss] = m[:t_ss] * m[:ψ_H2] / m[:θ] / (1 - m[:f_H1]) # from Giorgio's MATLAB script for solving steady state
    m[:t_S_ss]  = (m[:t_ss] - m[:θ] * ( m[:f_H1] * m[:t_H1_ss] + (1 - m[:f_H1]) * m[:t_H2_ss])) / (1 - m[:θ]) # from Giorgio's MATLAB script for solving steady state

    m[:c_H1_ss] = m[:c_ss] / m[:θ] / m[:f_H1] # from Giorgio's MATLAB script for solving steady state
    m[:c_H2_ss] = m[:c_ss] / m[:θ] / (1 - m[:f_H1]) # from Giorgio's MATLAB script for solving steady state
    m[:c_S_ss]  = m[:c_ss] / (1 - m[:θ]) # from Giorgio's MATLAB script for solving steady state
    m[:λ_H1_ss] = exp(m[:γ]) / (exp(m[:γ]) * m[:c_H1_ss] - m[:h] * m[:c_ss]) # from Giorgio's MATLAB script for solving steady state
    m[:λ_H2_ss] = exp(m[:γ]) / (exp(m[:γ]) * m[:c_H2_ss] - m[:h] * m[:c_ss]) # from Giorgio's MATLAB script for solving steady state
    m[:λ_S_ss]  = exp(m[:γ]) / (exp(m[:γ]) *m[:c_S_ss] - m[:h] * m[:c_ss]) # from Giorgio's MATLAB script for solving steady state
    m[:R_ss]    = 1.025 / 4.
    m[:bᴿ_ss]   = (m[:g_ss] - m[:t_ss] + m[:τ_ss]) / (1 - m[:R_ss] / exp(m[:γ]) / m[:π_ss]) # from Giorgio's MATLAB script for solving steady state

    # PROBLEM STEADY STATES
    m[:λ_SH1_ss]= (1 - m[:θ]) * m[:f_S1] / ((1 - m[:θ]) * m[:f_S1] + m[:θ] * m[:f_H1]) * m[:λ_S_ss] +
                       m[:θ] * m[:f_H1] / ((1 - m[:θ]) * m[:f_S1] + m[:θ] * m[:f_H1]) * m[:λ_H1_ss] # ???
    m[:λ_SH2_ss]= (1 - m[:θ]) * (1 - m[:f_S1]) / ((1 - m[:θ]) * (1 - m[:f_S1]) + m[:θ] * (1 - m[:f_H1])) * m[:λ_S_ss] +
                       m[:θ] * (1 - m[:f_H1]) / ((1 - m[:θ]) * (1 - m[:f_S1]) + m[:θ] * (1 - m[:f_H1])) * m[:λ_H2_ss] # ???
    m[:mc1_ss]  = 1 / (m[:α] ^ m[:α] * (1 - m[:α]) ^ (1 - m[:α])) * m[:ρ_ss] ^ m[:α] * m[:w1_ss] ^ (1 - m[:α]) # ???
    m[:mc1_ss]  = 1 / (m[:α] ^ m[:α] * (1 - m[:α]) ^ (1 - m[:α])) * m[:ρ_ss] ^ m[:α] * (m[:w1_ss]/ (m[:A2])) ^ (1 - m[:α]) # ???


    return m
end

function model_settings!(m::COTHANK)

    DSGE.default_settings!(m)

    # Anticipated shocks
    m <= Setting(:n_anticipated_shocks, 0,
                 "Number of anticipated policy shocks")
    m <= Setting(:n_anticipated_shocks_padding, 0,
                 "Padding for anticipated policy shocks")

    # Data
    m <= Setting(:data_id, 3, "Dataset identifier")
    if get_setting(m, :cond_id) in collect(1:5)
        m <= Setting(:cond_full_names, [:obs_gdp, :obs_corepce, :obs_spread, :obs_nominalrate, :obs_longrate],
                     "Observables used in conditional forecasts")
    elseif get_setting(m, :cond_id) == 6
        m <= Setting(:cond_full_names, [:obs_gdp, :obs_corepce, :obs_spread, :obs_nominalrate, :obs_longrate,
                                        :obs_gdpdeflator], "Observables used in conditional forecasts")
    end
    m <= Setting(:cond_semi_names, [:obs_spread, :obs_nominalrate, :obs_longrate],
                 "Observables used in semiconditional forecasts")

    # Forecast
    m <= Setting(:use_population_forecast, true,
                 "Whether to use population forecasts as data")
    m <= Setting(:shockdec_startdate, Nullable(quartertodate("2007-Q1")),
                 "Date of start of shock decomposition output period. If null, then shockdec starts at date_mainsample_start")
end
