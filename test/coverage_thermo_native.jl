using Test
using Arrhenius

const _RGAS = Arrhenius.R

# NASA7 coefficients from the covdepsurf.yaml example (Blondal et al. 2019),
# laid out in the SurfaceMechanism row format: [T_mid, high(7), low(7)].
function nasa7_row(tmid, low, high)
    row = zeros(15)
    row[1] = tmid
    row[2:8] .= high
    row[9:15] .= low
    return row
end

const PT_NASA7 = nasa7_row(1554.81,
    [7.10129478e-03, -4.25609798e-05, 8.98507278e-08, -7.80169595e-11,
     2.32458299e-14, -0.876096726, -0.0311207473],
    [0.16030291, -2.52239722e-04, 1.14183461e-07, -1.21476333e-11,
     3.85825979e-16, -70.8116648, -0.909545048])
const CO_NASA7 = nasa7_row(891.33,
    [-1.38214121, 0.0375305409, -8.29758476e-05, 8.09701555e-08, -2.85470829e-11,
     -3.45176032e+04, 4.3544767],
    [1.3809066, 8.0571901e-03, -4.6430896e-06, 8.91170699e-10, -5.90048361e-14,
     -3.43319289e+04, -4.85318015])
const O_NASA7 = nasa7_row(888.26,
    [-0.759013067, 0.0189868498, -3.82473745e-05, 3.43558395e-08, -1.13974372e-11,
     -1.72389494e+04, 1.76017396],
    [1.89893619, 2.03295425e-03, -1.19976574e-06, 2.32680659e-10, -1.53508282e-14,
     -1.75144954e+04, -9.6410408])

# Independent NASA7 evaluation for cross-checking the shared base-state layout.
function nasa7_hs(row, T)
    a = T <= row[1] ? row[9:15] : row[2:8]
    h = _RGAS*(T*(a[1] + T*(a[2]/2 + T*(a[3]/3 + T*(a[4]/4 + T*a[5]/5)))) + a[6])
    s = _RGAS*(a[1]*log(T) + T*(a[2] + T*(a[3]/2 + T*(a[4]/3 + T*a[5]/4))) + a[7])
    cp = _RGAS*(a[1] + T*(a[2] + T*(a[3] + T*(a[4] + T*a[5]))))
    return h, s, cp
end

function pt_co_model(; dependencies=Arrhenius.AbstractCoverageDependency[],
                     reference_coverage=1.0)
    return CoverageThermoModel(species_names=["Pt(s)", "CO(s)"],
        thermo_type=[1, 1], thermo_coefficients=vcat(PT_NASA7', CO_NASA7'),
        reference_coverage=reference_coverage, dependencies=dependencies)
end

@testset "Coverage-dependent surface thermo" begin

@testset "Low-coverage limit reduces to NASA7 base state" begin
    m = pt_co_model(dependencies=[LinearDependency(2, 2; enthalpy=8.1e7, entropy=-4.9e3),
        HeatCapacityDependency(2, 2; a=1.9e3, b=-1.5e4)])
    # theta_CO = 0: all self-interaction corrections vanish.
    out = coverage_thermo(m, 300.0, [1.0, 0.0])
    h, s, cp = nasa7_hs(CO_NASA7, 300.0)
    @test out.enthalpy_RT[2] ≈ h/(_RGAS*300.0) rtol=1e-14
    @test out.entropy_R[2] ≈ s/_RGAS rtol=1e-14
    @test out.cp_R[2] ≈ cp/_RGAS rtol=1e-14
    hp, sp, cpp = nasa7_hs(PT_NASA7, 300.0)
    @test out.enthalpy_RT[1] ≈ hp/(_RGAS*300.0) rtol=1e-14
    # Both branches of the NASA7 range switch agree with the independent form.
    for T in (200.0, 891.33, 900.0, 2500.0)
        out = coverage_thermo(m, T, [0.4, 0.6])
        h, s, cp = nasa7_hs(CO_NASA7, T)
        eh = 8.1e7*0.6
        es = -4.9e3*0.6
        q2 = 0.36
        ecp = (1.9e3*log(T) - 1.5e4)*q2
        eh += (T*(1.9e3*log(T) - 1.9e3 - 1.5e4) -
               298.15*(1.9e3*log(298.15) - 1.9e3 - 1.5e4))*q2
        es += 0.5*(log(T)*(1.9e3*log(T) - 3.0e4) -
                   log(298.15)*(1.9e3*log(298.15) - 3.0e4))*q2
        @test out.enthalpy_RT[2] ≈ (h+eh)/(_RGAS*T) rtol=1e-13
        @test out.entropy_R[2] ≈ (s+es)/_RGAS rtol=1e-13
        @test out.cp_R[2] ≈ (cp+ecp)/_RGAS rtol=1e-13
        @test out.gibbs_RT[2] ≈ out.enthalpy_RT[2] - out.entropy_R[2] rtol=1e-14
    end
end

@testset "Thermodynamic identities dh/dT = cp and ds/dT = cp/T" begin
    m = pt_co_model(dependencies=[PolynomialDependency(2, 2;
            enthalpy_coefficients=[1.2e7, 6.1e7, 3.1e7, 0.0],
            entropy_coefficients=[6.1e3, -1.5e4, 7.6e3, 0.0]),
        HeatCapacityDependency(2, 2; a=1.9e3, b=-1.5e4)])
    theta = [0.35, 0.65]
    # Fourth-order differences avoid the second-order truncation error that
    # dominates the very small platinum heat capacity at these temperatures.
    derivative(f, T, h) = (8*(f(T+h)-f(T-h))-(f(T+2h)-f(T-2h)))/(12h)
    for T in (250.0, 400.0, 850.0, 1200.0, 3000.0)
        hstep = 1e-4*T
        mid = coverage_thermo(m, T, theta)
        hdim(t) = coverage_thermo(m, t, theta).enthalpy_RT[2]*_RGAS*t
        sdim(t) = coverage_thermo(m, t, theta).entropy_R[2]*_RGAS
        cpdim = mid.cp_R[2]*_RGAS
        @test derivative(hdim, T, hstep) ≈ cpdim rtol=1e-7
        @test derivative(sdim, T, hstep) ≈ cpdim/T rtol=1e-7
        # Target species without corrections still obeys the identities.
        hp(t) = coverage_thermo(m, t, theta).enthalpy_RT[1]*_RGAS*t
        @test derivative(hp, T, hstep) ≈ mid.cp_R[1]*_RGAS rtol=1e-7
    end
end

@testset "Linear and polynomial laws" begin
    m = pt_co_model(dependencies=[LinearDependency(2, 2; enthalpy=8.0e7, entropy=-5.0e3)])
    T = 300.0
    lo = coverage_thermo(m, T, [0.8, 0.2])
    hi = coverage_thermo(m, T, [0.4, 0.6])
    @test (hi.enthalpy_RT[2] - lo.enthalpy_RT[2])*_RGAS*T ≈ 8.0e7*0.4 rtol=1e-12
    @test (hi.entropy_R[2] - lo.entropy_R[2])*_RGAS ≈ -5.0e3*0.4 rtol=1e-12
    @test lo.cp_R[2] ≈ hi.cp_R[2] rtol=1e-14

    coeffs = [0.5e7, -1.0e7, 2.0e7, -0.25e7]
    scoeffs = [1.0e3, 2.0e3, -3.0e3, 0.5e3]
    m = pt_co_model(dependencies=[PolynomialDependency(2, 2;
        enthalpy_coefficients=coeffs, entropy_coefficients=scoeffs)])
    for q in (0.13, 0.47, 1.0)
        out = coverage_thermo(m, T, [1-q, q])
        h, s, _ = nasa7_hs(CO_NASA7, T)
        @test out.enthalpy_RT[2]*_RGAS*T - h ≈ sum(coeffs .* q.^(1:4)) rtol=1e-12
        @test out.entropy_R[2]*_RGAS - s ≈ sum(scoeffs .* q.^(1:4)) rtol=1e-12
    end
end

@testset "Piecewise-linear law: slopes, knot continuity, boundaries" begin
    elow, ehigh, echange = 4.0e7, 1.5e8, 0.44
    slow, shigh, schange = 4.0e3, -3.0e3, 0.22
    m = pt_co_model(dependencies=[PiecewiseLinearDependency(2, 2;
        enthalpy_low=elow, enthalpy_high=ehigh, enthalpy_change=echange,
        entropy_low=slow, entropy_high=shigh, entropy_change=schange)])
    T = 300.0
    htheta(q) = coverage_thermo(m, T, [1-q, q]).enthalpy_RT[2]*_RGAS*T
    stheta(q) = coverage_thermo(m, T, [1-q, q]).entropy_R[2]*_RGAS
    h0, s0, _ = nasa7_hs(CO_NASA7, T)
    # Low/high region slopes.
    @test htheta(0.30) - htheta(0.20) ≈ elow*0.10 rtol=1e-12
    @test htheta(0.70) - htheta(0.60) ≈ ehigh*0.10 rtol=1e-12
    @test stheta(0.15) - stheta(0.05) ≈ slow*0.10 rtol=1e-12
    @test stheta(0.60) - stheta(0.50) ≈ shigh*0.10 rtol=1e-12
    # Continuity across the knots (finite-difference one-sided values agree).
    eps_q = 1e-9
    @test htheta(echange - eps_q) ≈ htheta(echange + eps_q) rtol=1e-7
    @test stheta(schange - eps_q) ≈ stheta(schange + eps_q) rtol=1e-7
    # Boundary values.
    @test htheta(0.0) ≈ h0 rtol=1e-14
    @test htheta(1.0) - h0 ≈ elow*echange + ehigh*(1-echange) rtol=1e-12
    @test stheta(1.0) - s0 ≈ slow*schange + shigh*(1-schange) rtol=1e-12
end

@testset "Interpolative law: knots, segments, boundaries" begin
    covs = collect(0.0:0.25:1.0)
    enthalpies = [0.0, 1.0e7, 3.0e7, 2.0e7, 8.0e7]
    entropies = [0.0, -1.0e3, 4.0e3, 1.0e3, -2.0e3]
    m = pt_co_model(dependencies=[InterpolativeDependency(2, 2;
        enthalpy_coverages=covs, enthalpies=enthalpies,
        entropy_coverages=covs, entropies=entropies)])
    T = 300.0
    h0, s0, _ = nasa7_hs(CO_NASA7, T)
    for (i, q) in enumerate(covs)
        out = coverage_thermo(m, T, [1-q, q])
        @test out.enthalpy_RT[2]*_RGAS*T - h0 ≈ enthalpies[i] rtol=1e-12
        @test out.entropy_R[2]*_RGAS - s0 ≈ entropies[i] rtol=1e-12
    end
    # Mid-segment values are linear between knots.
    out = coverage_thermo(m, T, [0.875, 0.125])
    @test out.enthalpy_RT[2]*_RGAS*T - h0 ≈ 0.5*(enthalpies[1] + enthalpies[2]) rtol=1e-12
    @test out.entropy_R[2]*_RGAS - s0 ≈ 0.5*(entropies[1] + entropies[2]) rtol=1e-12
end

@testset "Heat-capacity law anchored at 298.15 K" begin
    a, b = 1.9e3, -1.5e4
    dep = HeatCapacityDependency(2, 2; a=a, b=b)
    m = pt_co_model(dependencies=[dep])
    base = pt_co_model()
    theta = [0.5, 0.5]
    # Enthalpy and entropy corrections vanish exactly at the anchor temperature.
    anchored = coverage_thermo(m, 298.15, theta)
    reference = coverage_thermo(base, 298.15, theta)
    @test anchored.enthalpy_RT ≈ reference.enthalpy_RT rtol=1e-14
    @test anchored.entropy_R ≈ reference.entropy_R rtol=1e-14
    # Heat capacity itself is not anchored; it is quadratic in coverage.
    @test anchored.cp_R[2] - reference.cp_R[2] ≈ (a*log(298.15) + b)*0.25/_RGAS rtol=1e-12
    quarter = coverage_thermo(m, 800.0, [0.75, 0.25]).cp_R[2]
    full = coverage_thermo(m, 800.0, [0.0, 1.0]).cp_R[2]
    basecp = coverage_thermo(base, 800.0, theta).cp_R[2]
    @test (quarter - basecp)/(full - basecp) ≈ 0.25^2 rtol=1e-12
    # Away from the anchor, enthalpy/entropy pick up the integrated terms.
    out = coverage_thermo(m, 800.0, theta)
    @test out.enthalpy_RT[2] != coverage_thermo(base, 800.0, theta).enthalpy_RT[2]
end

@testset "Reference-coverage standard-state correction" begin
    deps = [LinearDependency(2, 2; enthalpy=8.0e7, entropy=-5.0e3),
        HeatCapacityDependency(2, 2; a=1.9e3, b=-1.5e4)]
    plain = pt_co_model(dependencies=deps)
    shifted = pt_co_model(dependencies=deps, reference_coverage=0.11)
    for T in (300.0, 1000.0)
        theta = [0.61, 0.39]
        a = coverage_thermo(plain, T, theta)
        b = coverage_thermo(shifted, T, theta)
        # Only the standard-state entropy and Gibbs energy shift, uniformly.
        @test a.enthalpy_RT ≈ b.enthalpy_RT rtol=1e-14
        @test a.cp_R ≈ b.cp_R rtol=1e-14
        @test b.entropy_R .- a.entropy_R ≈ fill(log(0.11), 2) rtol=1e-12
        @test b.gibbs_RT .- a.gibbs_RT ≈ fill(-log(0.11), 2) rtol=1e-12
        # The shift is a standard-state reference correction, not a mixing
        # term: it does not depend on the species coverage.
        c = coverage_thermo(shifted, T, [0.05, 0.95])
        d = coverage_thermo(plain, T, [0.05, 0.95])
        @test c.entropy_R .- d.entropy_R ≈ fill(log(0.11), 2) rtol=1e-12
    end
end

@testset "Self- and cross-interactions on a triangular coverage map" begin
    e_self, e_cross = 8.0e7, 4.0e7
    s_self = -5.0e3
    m = CoverageThermoModel(species_names=["Pt(s)", "CO(s)", "O(s)"],
        thermo_type=[1, 1, 1],
        thermo_coefficients=vcat(PT_NASA7', CO_NASA7', O_NASA7'),
        reference_coverage=0.11,
        dependencies=[LinearDependency(2, 2; enthalpy=e_self, entropy=s_self),
                      LinearDependency(2, 3; enthalpy=e_cross)])
    T = 300.0
    w = CoverageThermoWorkspace(m)
    h0, s0, _ = nasa7_hs(CO_NASA7, T)
    covs = range(0.0, 1.0; length=11)
    for qco in covs, qo in covs
        qco + qo > 1 + 1e-12 && continue
        theta = [max(0., 1-qco-qo), qco, qo]
        out = coverage_thermo!(w, m, T, theta)
        @test out.enthalpy_RT[2]*_RGAS*T ≈ h0 + e_self*qco + e_cross*qo rtol=1e-12
        @test out.entropy_R[2]*_RGAS ≈ s0 + s_self*qco + _RGAS*log(0.11) rtol=1e-12
        # O(s) and Pt(s) carry no dependencies.
        @test out.enthalpy_RT[3]*_RGAS*T ≈ nasa7_hs(O_NASA7, T)[1] rtol=1e-13
        # Workspace arrays alias the returned buffers across states.
        @test out.enthalpy_RT === w.enthalpy_RT
    end
end

@testset "Constant-cp base species" begin
    coeffs = zeros(2, 15)
    coeffs[1, 1:4] .= [298.15, 1.0e6, 2.0e5, 2.5e4]
    coeffs[2, 1:4] .= [298.15, -1.1e8, 3.0e5, 2.0e4]
    m = CoverageThermoModel(species_names=["PT(s)", "CO(s)"], thermo_type=[2, 2],
        thermo_coefficients=coeffs,
        dependencies=[LinearDependency(2, 2; enthalpy=8.0e7)])
    T = 900.0
    out = coverage_thermo(m, T, [0.7, 0.3])
    h = -1.1e8 + 2.0e4*(T - 298.15) + 8.0e7*0.3
    s = 3.0e5 + 2.0e4*log(T/298.15)
    @test out.enthalpy_RT[2] ≈ h/(_RGAS*T) rtol=1e-13
    @test out.entropy_R[2] ≈ s/_RGAS rtol=1e-13
    @test out.cp_R[2] ≈ 2.0e4/_RGAS rtol=1e-13
    @test out.enthalpy_RT[1] ≈ (1.0e6 + 2.5e4*(T-298.15))/(_RGAS*T) rtol=1e-13
end

@testset "SI parameter archive" begin
    text = raw"""{"format":"arrhenius-coverage-thermo-v1", "units":"K, J/kmol, J/kmol/K",
      "species_names":["A"], "thermo_type":[2], "thermo_coefficients":[[300,6000,300,20]],
      "reference_state_coverage":0.5,
      "dependencies":[{"target":1,"influencing":1,"kind":"polynomial",
        "enthalpy":[100,0,0,0],"entropy":[0,0,0,0],"heat_capacity_a":0,"heat_capacity_b":0}]}"""
    mktemp() do path, io
        write(io, text)
        close(io)
        model = CoverageThermoModel(path)
        out = coverage_thermo(model, 300., [1.])
        @test out.enthalpy_RT[1] ≈ 6100/(R*300)
        @test out.entropy_R[1] ≈ 300/R+log(.5)
        @test out.cp_R[1] ≈ 20/R
        for invalid in (replace(text, "thermo-v1"=>"thermo-v9"),
                        replace(text, "J/kmol/K"=>"cal/mol/K"),
                        replace(text, "polynomial"=>"unknown"))
            write(path, invalid)
            @test_throws ArgumentError CoverageThermoModel(path)
        end
    end
end

@testset "Input validation" begin
    good_names = ["Pt(s)", "CO(s)"]
    good_types = [1, 1]
    good_coeffs = vcat(PT_NASA7', CO_NASA7')
    @test_throws ArgumentError CoverageThermoModel(species_names=["A", "A"],
        thermo_type=good_types, thermo_coefficients=good_coeffs)
    @test_throws ArgumentError CoverageThermoModel(species_names=String[],
        thermo_type=Int[], thermo_coefficients=zeros(0, 15))
    @test_throws ArgumentError CoverageThermoModel(species_names=good_names,
        thermo_type=[1, 3], thermo_coefficients=good_coeffs)
    @test_throws DimensionMismatch CoverageThermoModel(species_names=good_names,
        thermo_type=good_types, thermo_coefficients=zeros(3, 15))
    @test_throws DimensionMismatch CoverageThermoModel(species_names=good_names,
        thermo_type=good_types, thermo_coefficients=zeros(2, 14))
    @test_throws ArgumentError CoverageThermoModel(species_names=good_names,
        thermo_type=good_types, thermo_coefficients=fill(NaN, 2, 15))
    for bad in (0.0, -0.1, 1.1, NaN)
        @test_throws ArgumentError CoverageThermoModel(species_names=good_names,
            thermo_type=good_types, thermo_coefficients=good_coeffs,
            reference_coverage=bad)
    end
    @test_throws ArgumentError CoverageThermoModel(species_names=good_names,
        thermo_type=good_types, thermo_coefficients=good_coeffs,
        dependencies=[LinearDependency(3, 2; enthalpy=1.0)])
    @test_throws ArgumentError CoverageThermoModel(species_names=good_names,
        thermo_type=good_types, thermo_coefficients=good_coeffs,
        dependencies=[LinearDependency(2, 0; enthalpy=1.0)])
    @test_throws ArgumentError LinearDependency(2, 2; enthalpy=Inf)
    @test_throws ArgumentError pt_co_model(dependencies=[LinearDependency(2, 2, Inf, 0.)])
    bad_coeffs = copy(good_coeffs)
    bad_coeffs[1, 1] = 0
    @test_throws ArgumentError CoverageThermoModel(species_names=good_names,
        thermo_type=good_types, thermo_coefficients=bad_coeffs)
    @test_throws ArgumentError PolynomialDependency(2, 2;
        enthalpy_coefficients=[1.0, 2.0])
    @test_throws ArgumentError PiecewiseLinearDependency(2, 2;
        enthalpy_low=1.0, enthalpy_high=2.0, enthalpy_change=1.5,
        entropy_low=1.0, entropy_high=2.0, entropy_change=0.5)
    @test_throws ArgumentError InterpolativeDependency(2, 2;
        enthalpy_coverages=[0.1, 1.0], enthalpies=[0.0, 1.0],
        entropy_coverages=[0.0, 1.0], entropies=[0.0, 1.0])
    @test_throws ArgumentError InterpolativeDependency(2, 2;
        enthalpy_coverages=[0.0, 0.9], enthalpies=[0.0, 1.0],
        entropy_coverages=[0.0, 1.0], entropies=[0.0, 1.0])
    @test_throws ArgumentError InterpolativeDependency(2, 2;
        enthalpy_coverages=[0.0, 0.5, 0.4, 1.0], enthalpies=zeros(4),
        entropy_coverages=[0.0, 1.0], entropies=[0.0, 1.0])
    @test_throws ArgumentError InterpolativeDependency(2, 2;
        enthalpy_coverages=[0.0, 0.5, 1.0], enthalpies=zeros(2),
        entropy_coverages=[0.0, 1.0], entropies=[0.0, 1.0])

    m = pt_co_model(dependencies=[LinearDependency(2, 2; enthalpy=1.0)])
    w = CoverageThermoWorkspace(m)
    @test_throws DomainError coverage_thermo!(w, m, 0.0, [0.5, 0.5])
    @test_throws DomainError coverage_thermo!(w, m, -300.0, [0.5, 0.5])
    @test_throws DomainError coverage_thermo!(w, m, NaN, [0.5, 0.5])
    @test_throws DimensionMismatch coverage_thermo!(w, m, 300.0, [1.0])
    @test_throws DomainError coverage_thermo!(w, m, 300.0, [1.2, -0.2])
    @test_throws DomainError coverage_thermo!(w, m, 300.0, [NaN, 0.0])
    @test_throws ArgumentError coverage_thermo!(w, m, 300.0, [0.6, 0.6])
    mismatched = CoverageThermoModel(species_names=["A"], thermo_type=[1],
        thermo_coefficients=reshape(PT_NASA7, 1, :))
    @test_throws DimensionMismatch coverage_thermo!(CoverageThermoWorkspace(mismatched), m, 300., [.5, .5])
end

end
