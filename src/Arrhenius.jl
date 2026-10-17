module Arrhenius
    using YAML
    using NPZ
    using SHA
    using LinearAlgebra
    using SparseArrays

    include("Constants.jl")
    include("BlowersMasel.jl")
    include("DataStructure.jl")
    include("GasMechanismKineticsIO.jl")
    include("GasMechanismThermoIO.jl")
    include("GasMechanismTransportIO.jl")
    include("GasMechanismAudit.jl")
    export audit_mechanism, refit_nasa7_continuity, refit_nasa7_common_tmid, refit_plog_to_chebyshev
    include("Solution.jl")
    include("Magic.jl")
    include("Thermo.jl")
    include("Kinetics.jl")
    include("Transport.jl")
    include("IonTransport.jl")
    include("MulticomponentTransport.jl")
    include("DustyGasTransport.jl")
    include("Equilibrium.jl")
include("SignedEquilibrium.jl")
    include("IdealGasStates.jl")
    include("IdealGasMixing.jl")
    include("Reactors.jl")
    include("ReactorNetworks.jl")
    include("MovingWallReactors.jl")
    include("SurfaceKinetics.jl")
    include("CondensedEquilibrium.jl")
    include("CoverageDependentThermo.jl")

    module Plasma
        using YAML
        using LinearAlgebra
        using SparseArrays
        import ForwardDiff
        using ..Arrhenius: R, one_atm, IdealGasThermo, KineticsWorkspace, PlogData, BlowersMaselData, SurfaceMechanism,
            SurfaceWorkspace, surface_rates!, surface_rates, surface_rhs, surface_problem,
            solve_surface, CoverageThermoModel, coverage_thermo!, coverage_thermo,
            ReactorSurface, CatalyticNetwork, catalytic_state, catalytic_rhs, catalytic_problem,
            solve_catalytic, catalytic_diagnostics
        import ..Arrhenius: reactor_state, reactor_rhs, reactor_jacobian!, reactor_problem,
            reactor_properties, solve_reactor, _ReversibleRatePlan, _nasa7_coefficients,
            _extended_thermo, _reactor_composition, _rate_factors!, _mass_action!, Reaction

        include("PlasmaEEDFTypes.jl")
        include("PlasmaEEDFIO.jl")
        include("PlasmaEEDF.jl")
        include("PlasmaKinetics.jl")
        include("PlasmaIO.jl")
        include("PlasmaRateIO.jl")
        include("PlasmaThermochemistry.jl")
        include("PlasmaEnergyReactor.jl")

        export EEDFModel, EEDFState, TwoTermOptions, EEDFResult, read_eedf_model, solve_eedf
        export PlasmaMechanism, PlasmaState, PlasmaReactor, PlasmaRHS
        export plasma_properties, plasma_rates, set_plasma_state!, set_mean_electron_energy!
        export set_reduced_electric_field!, update_eedf!, plasma_thermodynamics, set_plasma_enthalpy!
        export PlasmaEnergyReactor, PlasmaEnergyRHS
        export SurfaceMechanism, SurfaceWorkspace, surface_rates!, surface_rates, surface_rhs
        export surface_problem, solve_surface, CoverageThermoModel, coverage_thermo!, coverage_thermo
        export ReactorSurface, CatalyticNetwork, catalytic_state, catalytic_rhs, catalytic_problem
        export solve_catalytic, catalytic_diagnostics
    end
    using .Plasma: EEDFModel, EEDFState, TwoTermOptions, EEDFResult, read_eedf_model, solve_eedf,
        PlasmaMechanism, PlasmaState, PlasmaReactor, PlasmaRHS, plasma_properties, plasma_rates,
        set_plasma_state!, set_mean_electron_energy!, set_reduced_electric_field!, update_eedf!,
        plasma_thermodynamics, set_plasma_enthalpy!, PlasmaEnergyReactor, PlasmaEnergyRHS,
        ElectronCollisionKind, ElectronCollision, EffectiveCollision, ElasticCollision,
        ExcitationCollision, IonizationCollision, AttachmentCollision, EEDFWorkspace,
        _eedf_simpson, _integral_pq, _sg_coefficients, _add_collision_pq!, _eedf_norm,
        _plasma_import, _plasma_species, _EEDF_BOLTZMANN, _assemble_operator!, _maxwellian,
        _plasma_chebyshev_rate
    export EEDFModel, EEDFState, TwoTermOptions, EEDFResult, read_eedf_model, solve_eedf
    export PlasmaMechanism, PlasmaState, PlasmaReactor, PlasmaRHS
    export plasma_properties, plasma_rates, set_plasma_state!, set_mean_electron_energy!
    export set_reduced_electric_field!, update_eedf!, plasma_thermodynamics, set_plasma_enthalpy!
    export PlasmaEnergyReactor, PlasmaEnergyRHS
    # Preserve the legacy Arrhenius.<name> bindings for every name formerly
    # defined by direct includes, including private plasma helper names.
    for name in names(Plasma; all=true, imported=false)
        isdefined(@__MODULE__, name) && continue
        @eval const $(name) = Plasma.$(name)
    end
    include("SurfaceFlowReactors.jl")
    include("PureWater.jl")
    include("CriticalProperties.jl")
    include("RealGasThermo.jl")
    include("RealGasReactors.jl")
    include("Flames.jl")
    include("FlameSensitivity.jl")
    include("IonizedFlames.jl")
    include("FlameIO.jl")
    include("CounterflowFlames.jl")
    include("PremixedCounterflowFlames.jl")
    include("CatalyticFlames.jl")
    include("CounterflowContinuation.jl")
    include("Precision.jl")
end
