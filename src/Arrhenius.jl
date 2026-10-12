module Arrhenius
    using YAML
    using NPZ
    using SHA
    using LinearAlgebra
    using SparseArrays

    include("Constants.jl")
    include("PlasmaEEDFTypes.jl")
    include("PlasmaEEDFIO.jl")
    include("PlasmaEEDF.jl")
    export EEDFModel, EEDFState, TwoTermOptions, EEDFResult, read_eedf_model, solve_eedf

    include("BlowersMasel.jl")
    include("DataStructure.jl")
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
    include("PlasmaKinetics.jl")
    include("PlasmaIO.jl")
    include("PlasmaRateIO.jl")
    include("PlasmaThermochemistry.jl")
    include("PlasmaEnergyReactor.jl")
    include("ReactorNetworks.jl")
    include("MovingWallReactors.jl")
    include("SurfaceKinetics.jl")
    include("CondensedEquilibrium.jl")
    include("CoverageDependentThermo.jl")
    include("SurfaceFlowReactors.jl")
    include("PureWater.jl")
    include("CriticalProperties.jl")
    include("RealGasThermo.jl")
    include("RealGasReactors.jl")
    include("Flames.jl")
    include("IonizedFlames.jl")
    include("FlameIO.jl")
    include("CounterflowFlames.jl")
    include("PremixedCounterflowFlames.jl")
    include("CatalyticFlames.jl")
    include("Precision.jl")
end
