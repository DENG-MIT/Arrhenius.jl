# Four coverage-dependent enthalpy/entropy laws and a cross-interaction map.
# First export covdep_lin, covdep_pwlin, covdep_poly, covdep_int and covdep_cross
# from Cantera's example_data/covdepsurf.yaml using export_coverage_thermo.py.
using Arrhenius

const COVERAGE_PHASE_NAMES = ("covdep_lin", "covdep_pwlin", "covdep_poly", "covdep_int", "covdep_cross")

function prepare_coverage_sweep(directory)
    models = [CoverageThermoModel(joinpath(directory, name*".coverage.json")) for name in COVERAGE_PHASE_NAMES]
    workspaces = CoverageThermoWorkspace.(models)
    return (; models, workspaces)
end

function coverage_dependent_surf_calculation(prepared)
    models, workspaces = prepared.models, prepared.workspaces
    length(models) == length(workspaces) == 5 || throw(DimensionMismatch("five surface phases required"))
    coverages = collect(range(0., 1.; length=101))
    curves = Array{Float64}(undef, 2, 101, 4)
    theta = zeros(2)
    for phase in 1:4
        for (point, co) in enumerate(coverages)
            theta[1], theta[2] = 1-co, co
            properties = coverage_thermo!(workspaces[phase], models[phase], 300., theta)
            curves[1, point, phase] = properties.enthalpy_RT[2]
            curves[2, point, phase] = properties.entropy_R[2]
        end
    end
    cross = zeros(101, 101)
    theta_cross = zeros(3)
    for (i, co) in enumerate(coverages)
        for j in 1:102-i
            oxygen = coverages[j]
            # Remove only subtraction roundoff at the occupied-site boundary.
            theta_cross[1] = max(0., 1-co-oxygen)
            theta_cross[2], theta_cross[3] = co, oxygen
            properties = coverage_thermo!(workspaces[5], models[5], 300., theta_cross)
            cross[i, j] = properties.enthalpy_RT[2]
        end
    end
    return (; coverages, curves, cross)
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) == 1 || error("supply the directory containing the five .coverage.json parameter files")
    result = coverage_dependent_surf_calculation(prepare_coverage_sweep(ARGS[1]))
    println("Four 101-point enthalpy/entropy curves; ", sum(!iszero, result.cross), " cross-interaction states")
    println("CO enthalpy h/RT at zero/full coverage: ", result.curves[1, [1, end], :])
end
