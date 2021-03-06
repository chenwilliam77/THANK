"""
```
THANK{T} <: AbstractRepModel{T}
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
mutable struct THANK{T} <: AbstractRepModel{T}
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

description(m::THANK) = "Bilbiie, Primiceri, Tambalotti (2019), $(m.subspec)."

"""
`init_model_indices!(m::THANK)`

Arguments:
`m:: THANK`: a model object

Description:
Initializes indices for all of `m`'s states, shocks, and equilibrium conditions.
"""
function init_model_indices!(m::THANK)
    # Endogenous states
    endogenous_states = [:y_t, :k_t, :L_t, :Rk_t, :w_t, :π_t, :mc_t, :λ_t,       # sticky prices
                         :c_t, :R_t, :u_t, :ϕ_t, :i_t, :kbar_t, :wgap_t,
                         :gdp_t, :z_t, :g_t, :μ_t, :λ_p_t, :λ_w_t,
                         :b_t, :mp_t, :Eπ_t, :Ec_t, :Eλ_t, :Eϕ_t, :ERk_t, :Ei_t,
                         :Ew_t, :gdp_t1, :c_t1, :i_t1, :w_t1, :λ_p_t1, :λ_w_t1,
                         :λ_h_t, :λ_s_t, :c_h_t, :c_s_t, :t_h_t, :Ey_t, :Eλ_s_t, :Eλ_h_t,
                                                                                # CHECK BACK IF ANY CAN BE AUGMENTED VARIABLES
                         :y_f_t, :k_f_t, :L_f_t, :Rk_f_t, :w_f_t, :mc_f_t,       # flexible prices
                         :λ_f_t, :c_f_t, :R_f_t, :u_f_t, :ϕ_f_t,
                         :i_f_t, :kbar_f_t, :wgap_f_t, :gdp_f_t, :Ec_f_t, :Eλ_f_t,
                         :Eϕ_f_t, :ERk_f_t, :Ei_f_t, :λ_h_f_t, :λ_s_f_t, :c_h_f_t, :c_s_f_t,
                         :t_h_f_t, :Ey_f_t, :Eλ_h_f_t, :Eλ_s_f_t]

    # Exogenous shocks
    exogenous_shocks = [:R_sh, :z_sh, :g_sh, :μ_sh, :λ_p_sh, :λ_w_sh, :b_sh]

    # Expectations shocks
    expected_shocks = [:Eπ_sh, :Ec_sh, :Eλ_sh, :Eϕ_sh, :ERk_sh, :Ei_sh, :Ew_sh, # sticky prices
                       :Ec_f_sh, :Eλ_f_sh, :Eϕ_f_sh, :ERk_f_sh, :Ei_f_sh,       # flexible prices
                       :Ey_sh, :Eλ_s_sh, :Eλ_h_sh, :Ey_f_sh, :Eλ_s_f_sh, :Eλ_h_f_sh] # new parameters

    # Equilibrium conditions
    equilibrium_conditions = [:eq_y, :eq_k, :eq_L, :eq_Rk, :eq_w, :eq_π, :eq_mc, :eq_λ,     # sticky prices
                              :eq_c, :eq_R, :eq_u, :eq_ϕ, :eq_i, :eq_kbar, :eq_wgap,
                              :eq_gdp, :eq_z, :eq_g, :eq_μ, :eq_λ_p, :eq_λ_w,
                              :eq_b, :eq_mp, :eq_Eπ,  :eq_Ec, :eq_Eλ, :eq_Eϕ,
                              :eq_ERk, :eq_Ei, :eq_Ew, :eq_gdp1, :eq_c1, :eq_i1, :eq_w1,    # CHECK BACK IF ANY CAN BE AUGMENTED VARIABLES
                              :eq_λ_p1, :eq_λ_w1, :eq_λ_h, :eq_λ_s, :eq_c_h, :eq_c_s,
                              :eq_t_h, :eq_Ey, :eq_Eλ_s, :eq_Eλ_h,

                              :eq_y_f, :eq_k_f, :eq_L_f, :eq_Rk_f, :eq_w_f, :eq_mc_f,       # flexible prices
                              :eq_λ_f, :eq_c_f, :eq_R_f, :eq_u_f, :eq_ϕ_f,
                              :eq_i_f, :eq_kbar_f, :eq_wgap_f, :eq_gdp_f, :eq_Ec_f, :eq_Eλ_f,
                              :eq_Eϕ_f, :eq_ERk_f, :eq_Ei_f, :eq_λ_h_f, :eq_λ_s_f, :eq_c_h_f,
                              :eq_c_s_f, :eq_t_h_f, :eq_Ey_f, :eq_Eλ_h_f, :eq_Eλ_s_f]

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

function THANK(subspec::String="ss1";
                   custom_settings::Dict{Symbol, Setting} = Dict{Symbol, Setting}(),
                   testing = false)

    # Model-specific specifications
    spec               = split(basename(@__FILE__),'.')[1]
    subspec            = subspec
    settings           = Dict{Symbol,Setting}()
    test_settings      = Dict{Symbol,Setting}()
    rng                = MersenneTwister(0)

    # Initialize empty model
    m = THANK{Float64}(
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
init_parameters!(m::THANK)
```

Initializes the model's parameters, as well as empty values for the steady-state
parameters (in preparation for `steadystate!(m)` being called to initialize
those).
"""
function init_parameters!(m::THANK)
    # Calibrated parameters
    m <= parameter(:gy_ss, 0.2, fixed=true,
                   description="g_ss_ratio: steady state government spending to GDP ratio.",
                   tex_label = "g^{ss}" )
    m <= parameter(:δ, 0.025, fixed=true,
                   description="δ: The capital depreciation rate.",
                   tex_label = "\\delta" )
    m <= parameter(:t_h_0_L_ss, 0.0, fixed=true,
                   description="New parameter! Might estimate", # FIX THIS DESCRIPTION
                   tex_label = "g^{ss}" )
    # Non standard devation parameters
    m <= parameter(:α, 0.25, (0., 1.), (1e-5, 0.999), ModelConstructors.SquareRoot(), Normal(0.30, 0.05), fixed = false,
                   description="α: Capital elasticity in the intermediate goods sector's production function (also known as the capital share).",
                   tex_label = "\\alpha")

    m <= parameter(:ι_p, 0.21, (0., 1.), (1e-5, 0.999), ModelConstructors.SquareRoot(), Beta(0.5, 0.15), fixed = false,
                   description = "ι_p: The weight attributed to last period's inflation in price indexation. " *
                   "(1-ι_p) is the weight attributed to steady-state inflation.",
                   tex_label = "\\iota_p")

    m <= parameter(:ι_w, 0.15, (0., 1.), (1e-5, 0.999), ModelConstructors.SquareRoot(), Beta(0.5, 0.15), fixed = false,
                   description="ι_w: The weight attributed to last period's wage in wage indexation. "
                    * "(1-ι_w) is the weight attributed to steady-state wages.",
                   tex_label = "\\iota_w")

    m <= parameter(:γ100, 0.53, (-5.0, 5.0), (-5., 5.), ModelConstructors.Untransformed(), Normal(0.5, 0.025), fixed = false,

                   description="γ: The log of the steady-state growth rate of technology.",
                   tex_label = "100\\gamma")

    m <= parameter(:h, 0.85, (-5.0, 5.0), (-5., 5.), ModelConstructors.Untransformed(), Beta(0.5, 0.1), fixed = false,
                   description="h: habit formation parameter.",
                   tex_label = "h")

    m <= parameter(:λ_p_ss, 0.25,(-1000.,1000.),(-1000., 1000.), ModelConstructors.Untransformed(), Normal(0.15, 0.05),  fixed = false,
                   description="λ_p_ss: The steady state net price markup.",
                   tex_label = "\\lambda_p")

    m <= parameter(:λ_w_ss, 0.12,(-1000.,1000.),(-1000., 1000.), ModelConstructors.Untransformed(), Normal(0.15, 0.05), fixed = false,
                   description = "λ_w_ss: The steady state net wage markup, which affects the elasticity of substitution between differentiated labor services.",
                   tex_label = "\\lambda_w")

    m <= parameter(:L_ss, 0.5, (-1000., 1000.), (-1e3, 1e3), ModelConstructors.Untransformed(), Normal(0.0,0.5), fixed=false,
                   description="L_ss: The steady state for log hours.",
                   tex_label = "\\log(L)^{ss}")

    m <= parameter(:π_ss100, 0.75, (1e-5, 10.), (1e-5, 10.), ModelConstructors.Exponential(), Normal(0.5, 0.1), fixed = false,

                   description="π_ss100: The steady-state rate of net inflation multiplied by 100.",
                   tex_label = "100 \\pi^{ss}")

    m <= parameter(:Fβ, 0.1, (1e-5, 10.), (1e-5, 10.), ModelConstructors.Exponential(), Gamma(0.25, 0.1), fixed=false,
                   description = "Fβ: Discount rate transformed.",
                   tex_label = "100(\\beta^{-1} - 1)")

    m <= parameter(:ν, 4., (1e-5, 10.), (1e-5, 10.), ModelConstructors.Exponential(), Gamma(2, 0.75), fixed = false,
                   description="ν_l: The inverse Frisch elasticity.",
                   tex_label = "\\nu")

    m <= parameter(:ξ_p, 0.85, (1e-5, 0.999), (1e-5, 0.999), ModelConstructors.SquareRoot(), Beta(0.66, 0.1), fixed=false,
                   description="ξ_p: The Calvo parameter. In every period, intermediate goods producers optimize prices with probability (1-ξ_p). "
                   * "With probability ξ_p, prices are adjusted according to a weighted average of the previous period's inflation "
                   * "(π_t1) and steady-state inflation (π_ss).",
                   tex_label = "\\xi_p")

    m <= parameter(:ξ_w, 0.75, (1e-5, 0.999), (1e-5, 0.999), ModelConstructors.SquareRoot(), Beta(0.66, 0.1), fixed=false,
                   description="ξ_w: (1-ξ_w) is the probability with which households can freely choose wages in each period. "
                   * "With probability ξ_w, wages increase at a geometrically weighted average of the steady state rate of "
                   * "wage increases and last period's productivity times last period's inflation.",
                   tex_label = "\\xi_w")

    m <= parameter(:χ, 5., (0., 10.), (1e-5, 0.), ModelConstructors.Exponential(), Gamma(5., 1.), fixed = false,
                   description="χ: The elasticity of the capital utilization cost function.",
                   tex_label = "\\chi")

    m <= parameter(:S′′, 2.5, (-15., 15.), (-15., 15.), ModelConstructors.Untransformed(), Gamma(4., 1), fixed = false,
                   description="S′′: The investment adjust cost.",
                   tex_label = "S^{\\prime\\prime}")

    m <= parameter(:ψ1, 2., (1e-5, 10.), (1e-5, 10.00), ModelConstructors.Exponential(), Normal(1.7, 0.3), fixed = false,
                   description="ψ₁: Weight on inflation gap in monetary policy rule.",
                   tex_label = "\\psi_1")

    m <= parameter(:ψ2, 0.1, (-0.5, 0.5), (-0.5, 0.5), ModelConstructors.Untransformed(), Normal(0.125, 0.05), fixed = false,
                   description="ψ₂: Weight on output gap in monetary policy rule.",
                   tex_label = "\\psi_2")

    m <= parameter(:ψ3, 0.25, (-0.5, 0.5), (-0.5, 0.5), ModelConstructors.Untransformed(), Normal(0.125, 0.05), fixed = false,
                   description="ψ₃: Weight on rate of change of output gap in the monetary policy rule.",
                   tex_label = "\\psi_3")

    m <= parameter(:η_λ_p, 0.75, (1e-5, 0.999), (1e-5, 0.999), ModelConstructors.SquareRoot(), Beta(0.50, 0.20), fixed = false,
                   description="η_λ_p: Moving average component in the price markup shock.",
                   tex_label = "\\eta_{\\lambda_p}")

    m <= parameter(:η_λ_w, 0.95, (1e-5, 0.999), (1e-5, 0.999), ModelConstructors.SquareRoot(), Beta(0.50, 0.20), fixed = false,
                   description="η_λ_w: Moving average component in the wage markup shock.",
                   tex_label = "\\eta_{\\lambda_w}")

    # exogenous processes - autocorrelation
    m <= parameter(:ρ_R, 0.8, (0., 1.), (1e-5, 0.999), ModelConstructors.SquareRoot(), Beta(0.6, 0.2), fixed = false,
                   description="ρ_R: persistence in the monetary policy rule.",
                   tex_label = "\\rho_{R}")

    m <= parameter(:ρ_z, 0.25, (1e-5, 0.999), (1e-5, 0.999), ModelConstructors.SquareRoot(), Beta(0.6, 0.2), fixed = false,
                   description="ρ_z: AR(1) coefficient in the technology process.",
                   tex_label = "\\rho_{z}")

    m <= parameter(:ρ_g, 0.98, (1e-5, 0.999), (1e-5, 0.999), ModelConstructors.SquareRoot(), Beta(0.6, 0.2), fixed = false,
                   description="ρ_g: AR(1) coefficient in the government spending process.",
                   tex_label = "\\rho_g")

    m <= parameter(:ρ_μ, 0.7, (1e-5, 0.999), (1e-5, 0.999), ModelConstructors.SquareRoot(), Beta(0.6, 0.2), fixed = false,
                   description="ρ_μ: AR(1) coefficient in capital adjustment cost process.",
                   tex_label = "\\rho_{\\mu}")

    m <= parameter(:ρ_λ_p, 0.95, (1e-5, 0.999), (1e-5, 0.999), ModelConstructors.SquareRoot(), Beta(0.6, 0.2), fixed = false,
                   description="ρ_λ_p: AR(1) coefficient in the price mark-up shock process.",
                   tex_label = "\\rho_{\\lambda_p}")

    m <= parameter(:ρ_λ_w, 0.98, (1e-5, 0.999), (1e-5, 0.999), ModelConstructors.SquareRoot(), Beta(0.6, 0.2), fixed = false,
                   description="ρ_λ_w: AR(1) coefficient in the wage mark-up shock process.",
                   tex_label = "\\rho_{\\lambda_w}")

    m <= parameter(:ρ_b, 0.7, (1e-5, 0.999), (1e-5, 0.999), ModelConstructors.SquareRoot(), Beta(0.6, 0.2), fixed = false,
                   description="ρ_b: AR(1) coefficient in the intertemporal preference shifter process.",
                   tex_label = "\\rho_b")

    m <= parameter(:ρ_mp, 0.15, (0., 1.), (1e-5, 0.999), ModelConstructors.SquareRoot(), Beta(0.4, 0.2), fixed = false,
                   description="ρ_mp: AR(1) coefficient in the monetary policy shock process.",
                   tex_label = "\\rho_{mp}")

    # exogenous processes - standard deviation
    m <= parameter(:σ_R, 0.2, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(2, 0.10), fixed = false,
                   description="σ_mp: The standard deviation of the monetary policy shock.",
                   tex_label = "\\sigma_{mp}")

    m <= parameter(:σ_z, 0.9, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(0.5, 1), fixed = false,
                   description="σ_z: The standard deviation of the technology process.",
                   tex_label = "\\sigma_{z}")

    m <= parameter(:σ_g, 0.35, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(0.5, 1), fixed = false,
                   description="σ_g: The standard deviation of the government spending process.",
                   tex_label = "\\sigma_{g}")

    m <= parameter(:σ_μ, 5., (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(0.5, 1), fixed = false,
                   description="σ_μ: The standard deviation of the exogenous marginal efficiency of investment shock process.",
                   tex_label = "\\sigma_{\\mu}")

    m <= parameter(:σ_λ_p, 0.15, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(0.1, 1), fixed = false,
                   description="σ_λ_p: The mean of the process that generates the price elasticity of the composite good. " *
                   "Specifically, the elasticity is (1+λ_{f,t})/(λ_{p_t}).",
                   tex_label = "\\sigma_{\\lambda_p}")

    m <= parameter(:σ_λ_w, 0.2, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(0.1, 1), fixed = false,
                   tex_label = "\\sigma_{\\lambda_w}")

    m <= parameter(:σ_b, 0.05, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Exponential(), RootInverseGamma(0.1, 1), fixed = false,
                   description="σ_b: The standard deviation of the intertemporal preference shifter process.",
                   tex_label = "\\sigma_{b}")
    # new parameters
    m <= parameter(:θ, 0.25, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Untransformed(), RootInverseGamma(0.5, 1), fixed = true)
    m <= parameter(:σ, 0.98, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Untransformed(), RootInverseGamma(0.5, 1), fixed = true)
    m <= parameter(:σ′, 0.0, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Untransformed(), RootInverseGamma(0.5, 1), fixed = true)
    m <= parameter(:τ_d, 0.1, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Untransformed(), RootInverseGamma(0.5, 1), fixed = true)
    m <= parameter(:τ_k, 0.0, (1e-8, 5.), (1e-8, 5.), ModelConstructors.Untransformed(), RootInverseGamma(0.5, 1), fixed = true)
    # steady states
    m <= SteadyStateParameter(:γ, NaN, tex_label = "\\gamma")
    m <= SteadyStateParameter(:expγ, NaN, tex_label = "")
    m <= SteadyStateParameter(:β, NaN, tex_label = "\\beta")
    m <= SteadyStateParameter(:r_ss, NaN, tex_label = "r^{ss}")
    m <= SteadyStateParameter(:r_ss100, NaN, tex_label = "100 r^{ss}")
    m <= SteadyStateParameter(:expL_ss, NaN, tex_label = "L^{ss}")
    m <= SteadyStateParameter(:Rk_ss, NaN, description = "Steady-state short-term rate of return on capital.", tex_label = "R^{k, ss}")
    m <= SteadyStateParameter(:mc_ss, NaN, description = "Steady-state marginal cost", tex_label = "s^{ss}")
    m <= SteadyStateParameter(:w_ss, NaN, tex_label = "w^{ss}")
    m <= SteadyStateParameter(:kL_ss, NaN, tex_label = "k^{ss}/L^{ss}")
    m <= SteadyStateParameter(:FL_ss, NaN, description = "Steady-state ratio of fixed costs to labor", tex_label = "F^{ss}/L^{ss}")
    m <= SteadyStateParameter(:yL_ss, NaN, description = "Steady-state ratio of output to labor", tex_label = "y^{ss}/L^{ss}")
    m <= SteadyStateParameter(:gL_ss, NaN, tex_label = "k^{ss}")
    m <= SteadyStateParameter(:i_s_L_ss, NaN, tex_label = "k^{ss}")
    m <= SteadyStateParameter(:cL_ss, NaN, tex_label = "k^{ss}")
    m <= SteadyStateParameter(:c_h_L_ss, NaN, tex_label = "k^{ss}")
    m <= SteadyStateParameter(:c_s_L_ss, NaN, tex_label = "k^{ss}")
    m <= SteadyStateParameter(:λ_h_L_ss, NaN, tex_label = "k^{ss}")
    m <= SteadyStateParameter(:λ_s_L_ss, NaN, tex_label = "k^{ss}")
    m <= SteadyStateParameter(:k_ss, NaN, tex_label = "k^{ss}")
    m <= SteadyStateParameter(:F_ss, NaN, tex_label = "F^{ss}")
    m <= SteadyStateParameter(:y_ss, NaN, tex_label = "y^{ss}")
    m <= SteadyStateParameter(:c_ss, NaN, tex_label = "c^{ss}")

    # new equations
    m <= SteadyStateParameter(:c_h_ss, NaN, tex_label = "k^{ss}")
    m <= SteadyStateParameter(:c_s_ss, NaN, tex_label = "")
    m <= SteadyStateParameter(:i_s_ss, NaN, tex_label = "k^{ss}")
    m <= SteadyStateParameter(:λ_h_ss, NaN, tex_label = "k^{ss}")
    m <= SteadyStateParameter(:λ_s_ss, NaN, tex_label = "k^{ss}")
    m <= SteadyStateParameter(:λ_ss, NaN, tex_label = "k^{ss}")
end

"""
```
steadystate!(m::THANK)
```

Calculates the model's steady-state values. `steadystate!(m)` must be called whenever
the parameters of `m` are updated.
"""
function steadystate!(m::THANK)

    m[:γ]        = m[:γ100] / 100.
    m[:expγ]     = exp(m[:γ])
    m[:β]        = 100. / (m[:Fβ] + 100.)
    m[:r_ss]     = exp(m[:γ]) / m[:β] - 1.

    m[:r_ss100]  = 100. * m[:r_ss]
    m[:expL_ss]  = exp(m[:L_ss])
    m[:Rk_ss]    = exp(m[:γ]) / m[:β] - 1. + m[:δ]
    m[:mc_ss]    = 1. / (1. + m[:λ_p_ss])

    m[:w_ss]     = (m[:mc_ss] * ((1. - m[:α]) ^ (1. - m[:α])) / ((m[:α] ^ (-m[:α])) * m[:Rk_ss] ^ m[:α])) ^ (1. / (1. - m[:α]))

    m[:kL_ss]    = (m[:w_ss] / m[:Rk_ss]) * (m[:α] / (1. - m[:α]))
    m[:FL_ss]    = m[:kL_ss] ^ m[:α] - m[:Rk_ss] * m[:kL_ss] - m[:w_ss]
    m[:yL_ss]    = m[:kL_ss] ^ m[:α] - m[:FL_ss]
    m[:gL_ss]    = m[:gy_ss] * m[:yL_ss]
    m[:i_s_L_ss] = (1 - (1 - m[:δ]) * exp(-m[:γ])) * exp(m[:γ]) * m[:kL_ss] / (1 - m[:θ]) # check in on this later
    m[:cL_ss]    = m[:yL_ss] - (1 - m[:θ]) * m[:i_s_L_ss] - m[:gL_ss]
    m[:c_h_L_ss] = m[:w_ss] + m[:t_h_0_L_ss] + m[:τ_k] * m[:Rk_ss] * m[:kL_ss] / m[:θ] - m[:gL_ss]
    m[:c_s_L_ss] = (1 / (1 - m[:θ])) * m[:cL_ss] - (m[:θ]) / (1 - m[:θ]) * m[:c_h_L_ss]
    m[:λ_h_L_ss] = exp(m[:γ]) / (exp(m[:γ]) * m[:c_h_L_ss] - m[:h] * m[:cL_ss])
    m[:λ_s_L_ss] = exp(m[:γ]) / (exp(m[:γ]) * m[:c_s_L_ss] - m[:h] * m[:cL_ss])
    m[:k_ss]     = m[:kL_ss] * m[:expL_ss]

    m[:F_ss]     = m[:FL_ss] * m[:expL_ss]
    m[:y_ss]     = m[:yL_ss] * m[:expL_ss]
    m[:c_ss]     = m[:cL_ss] * m[:expL_ss]

    # New equations
    m[:c_h_ss]   = m[:c_h_L_ss] * exp(m[:L_ss])
    m[:c_s_ss]   = m[:c_s_L_ss] * exp(m[:L_ss])
    m[:i_s_ss]   = m[:i_s_L_ss] * exp(m[:L_ss])

    m[:λ_h_ss]   = m[:λ_h_L_ss] / exp(m[:L_ss])
    m[:λ_s_ss]   = m[:λ_s_L_ss] / exp(m[:L_ss])

    m[:λ_ss]     = m[:θ] * m[:λ_h_ss] + (1 - m[:θ]) * m[:λ_s_ss]
    return m
end

function model_settings!(m::THANK)

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
