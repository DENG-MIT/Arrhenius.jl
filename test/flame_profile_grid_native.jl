using Arrhenius, LinearAlgebra, Test

# A normalized inert fixture exercises grid/IO interfaces without a flame solve.
function profile_grid_test_flame(gas;grid=[0.,.2,.6,.8,1.],discretization=:conservative)
    n=gas.n_species
    Y=zeros(n);Y[findfirst(==("AR"),gas.species_names)]=1
    state=zeros(n+2,length(grid))
    state[1,:].=.3
    state[2:n+1,:].=Y
    state[end,:].=.04
    return BurnerFlame(gas,copy(grid),one_atm,300.,Y,state,2,300.,argmax(Y),false,.04,
        Float64[],Float64[],Float64[],:full_knots,:mixture_averaged,false,nothing,:mole,discretization)
end

function profile_grid_test_snapshot(f)
    names=fieldnames(typeof(f))
    return NamedTuple{names}(map(names) do name
        value=getfield(f,name)
        value isa AbstractArray ? copy(value) : value
    end)
end

function profile_grid_test_constant(actual,expected,label)
    @test size(actual) == size(expected)
    @test all(iszero(actual[i]) for i in eachindex(expected) if iszero(expected[i]))
    relative=maximum((abs(actual[i]-expected[i])/abs(expected[i])
        for i in eachindex(expected) if !iszero(expected[i]));init=0.)
    # Up to eight successive convex interpolations of an analytic constant.
    # Preserve exact zeros; allow only accumulated Float64 rounding elsewhere.
    @test relative <= 16eps(Float64)
    println("constant-field roundoff ",label,": max_relative_error=",relative,
        ", limit=",16eps(Float64))
end

@testset "prescribed-temperature grid policy and state transfer" begin
    mechanism=joinpath(@__DIR__,"..","mechanism","h2o2.yaml")
    gas=CreateSolution(mechanism)
    data=MultiTransportData(mechanism*".multicomponent.npz",gas)
    @test BurnerFlame(gas;X="H2:1.1,O2:1,AR:5",mdot=.03).profile_grid_policy == :full_knots
    positions=[0.,.3,.4,1.];temperatures=[300.,500.,300.,300.]
    f=profile_grid_test_flame(gas)
    before=copy(f.grid);anchor_z=f.grid[f.anchor]
    set_temperature_profile!(f,positions,temperatures;relative=false)
    @test f.profile_grid_policy == :full_knots
    @test f.grid == sort(unique(vcat(before,positions)))
    @test f.grid[f.anchor] == anchor_z
    profile_grid_test_constant(f.state[end,:],fill(.04,length(f.grid)),"full-knot mass flux")
    profile_grid_test_constant(mass_fractions(f),repeat(f.inlet_Y,1,length(f.grid)),"full-knot composition")

    f=profile_grid_test_flame(gas)
    set_temperature_profile!(f,positions,temperatures;relative=false,grid_policy=:adaptive)
    @test f.grid == before
    @test f.profile_grid_policy == :adaptive
    @test f.profile_positions == positions && f.profile_temperatures == temperatures
    @test !(positions[2] in f.grid) # hidden peak data survives off-grid
    positions[2]=.31;temperatures[2]=550.
    @test f.profile_positions == [0.,.3,.4,1.] && f.profile_temperatures == [300.,500.,300.,300.]
    @test Arrhenius._profile_chord_defects([.2,.6],f.profile_positions,f.profile_temperatures) ≈ [100.]
    split=Arrhenius._profile_chord_defects([.2,.4,.6],f.profile_positions,f.profile_temperatures)
    @test split ≈ [400/3,0.]
    # Chord error need not decrease each bisection; the independent Lipschitz
    # bound L*h/2 does. This peak is resolved under the unchanged 1% budget.
    @test maximum(split) <= 2000*.2/2
    for _ in 1:8
        Arrhenius._refine_flame!(f;slope=1.,curve=1.,ratio=100.) || break
    end
    @test maximum(Arrhenius._profile_chord_defects(f.grid,f.profile_positions,f.profile_temperatures)) <= 5.
    @test f.profile_positions == [0.,.3,.4,1.] && f.profile_temperatures == [300.,500.,300.,300.]
    @test temperature(f) ≈ [Arrhenius._interpolate_profile(f.profile_positions,f.profile_temperatures,z) for z in f.grid] atol=1e-12
    @test f.grid[f.anchor] == anchor_z
    profile_grid_test_constant(f.state[end,:],fill(.04,length(f.grid)),"refined mass flux")
    profile_grid_test_constant(mass_fractions(f),repeat(f.inlet_Y,1,length(f.grid)),"refined composition")

    linear=profile_grid_test_flame(gas)
    redundant=[0.,.13,.35,.7,1.]
    set_temperature_profile!(linear,redundant,300 .+ 200 .* redundant;relative=false,grid_policy=:adaptive)
    @test maximum(Arrhenius._profile_chord_defects(linear.grid,linear.profile_positions,linear.profile_temperatures)) < 1e-12
    @test !Arrhenius._refine_flame!(linear;slope=1.,curve=1.,ratio=100.)
    @test linear.grid == before && linear.profile_positions == redundant
    shifted=profile_grid_test_flame(gas;grid=2 .+ before)
    set_temperature_profile!(shifted,[0.,.3,1.],[300.,600.,900.];grid_policy=:adaptive)
    @test shifted.profile_positions == [2.,2.3,3.]
    @test shifted.grid == 2 .+ before

    # All invalid calls leave every field of the existing flame unchanged.
    saved=profile_grid_test_snapshot(f)
    @test_throws ArgumentError set_temperature_profile!(f,[0.,1.],[300.,500.];grid_policy=:unknown)
    @test isequal(profile_grid_test_snapshot(f),saved)
    @test_throws ArgumentError set_temperature_profile!(f,[0.,1.],[300.,500.];relative=1)
    @test isequal(profile_grid_test_snapshot(f),saved)
    for (z,T,exception) in (([0.,1.],[300.],DimensionMismatch),
            ([0.,.3,.3,1.],[300.,400.,500.,600.],ArgumentError),
            ([0.,1.],[300.,NaN],ArgumentError),([0.,1.],[300.,7000.],ArgumentError),
            ([0.,.9],[300.,500.],ArgumentError),([0.,1.],[301.,500.],ArgumentError))
        @test_throws exception set_temperature_profile!(f,z,T;relative=false,grid_policy=:adaptive)
        @test isequal(profile_grid_test_snapshot(f),saved)
    end
    legacy=profile_grid_test_flame(gas;discretization=:finite_difference)
    saved_legacy=profile_grid_test_snapshot(legacy)
    @test_throws ArgumentError set_temperature_profile!(legacy,[0.,1.],[300.,500.];grid_policy=:adaptive)
    @test isequal(profile_grid_test_snapshot(legacy),saved_legacy)
    set_temperature_profile!(legacy,[0.,.3,1.],[300.,500.,600.])
    @test legacy.grid == before # historical finite-difference sampling

    # Invalid Soret arguments must not insert knots as a side effect.
    f=profile_grid_test_flame(gas)
    set_temperature_profile!(f,[0.,.3,.4,1.],[300.,500.,300.,300.];grid_policy=:adaptive)
    saved=profile_grid_test_snapshot(f)
    @test_throws ArgumentError set_transport!(f,:mixture_averaged;soret=true)
    @test isequal(profile_grid_test_snapshot(f),saved)
    @test_throws ArgumentError set_transport!(f,:multicomponent;data,soret=true,flux_gradient_basis=:invalid)
    @test isequal(profile_grid_test_snapshot(f),saved)
    @test_throws ArgumentError set_transport!(f,:mixture_averaged;data=1)
    @test isequal(profile_grid_test_snapshot(f),saved)
    bad=deepcopy(data);bad.species_names[1]="wrong"
    @test_throws ArgumentError set_transport!(f,:multicomponent;data=bad,soret=true)
    @test isequal(profile_grid_test_snapshot(f),saved)
    expected=profile_grid_test_flame(gas)
    set_temperature_profile!(expected,f.profile_positions,f.profile_temperatures;relative=false)
    set_transport!(f,:multicomponent;data,soret=true)
    @test f.soret_enabled && f.profile_grid_policy == :full_knots && !f.converged
    @test f.grid == expected.grid && f.state == expected.state && f.anchor == expected.anchor
    @test f.profile_positions == expected.profile_positions && f.profile_temperatures == expected.profile_temperatures
    saved=profile_grid_test_snapshot(f)
    @test_throws ArgumentError set_temperature_profile!(f,[0.,1.],[300.,500.];grid_policy=:adaptive)
    @test isequal(profile_grid_test_snapshot(f),saved)
end

@testset "profile policy restart and invalid-archive atomicity" begin
    mechanism=joinpath(@__DIR__,"..","mechanism","h2o2.yaml")
    gas=CreateSolution(mechanism)
    data=MultiTransportData(mechanism*".multicomponent.npz",gas)
    mktempdir() do directory
        path=joinpath(directory,"profile.npz")
        f=profile_grid_test_flame(gas)
        set_temperature_profile!(f,[0.,.3,.4,1.],[300.,500.,300.,300.];grid_policy=:adaptive)
        save_flame(path,f)
        arrays=Arrhenius.npzread(path)
        @test String(copy(vec(arrays["profile_grid_policy_utf8"]))) == "adaptive"
        restored=profile_grid_test_flame(gas)
        restore_flame!(restored,path)
        @test restored.profile_grid_policy == :adaptive
        @test restored.grid == f.grid && restored.state == f.state
        @test restored.profile_positions == f.profile_positions && restored.profile_temperatures == f.profile_temperatures

        # Old conservative archives use full knots, including when the
        # destination currently holds an adaptive profile.
        old=deepcopy(arrays);delete!(old,"profile_grid_policy_utf8")
        Arrhenius.npzwrite(path,old)
        restore_flame!(restored,path)
        @test restored.profile_grid_policy == :full_knots
        @test all(z->z in restored.grid,f.profile_positions)
        profile_grid_test_constant(restored.state[end,:],fill(.04,length(restored.grid)),"historical restart mass flux")
        delete!(old,"discretization_utf8")
        Arrhenius.npzwrite(path,old)
        restore_flame!(restored,path)
        @test restored.discretization == :finite_difference
        @test restored.profile_grid_policy == :full_knots && restored.grid == f.grid

        # A Soret roundtrip retains the restored knots and full-knot policy.
        set_transport!(f,:multicomponent;data,soret=true)
        save_flame(path,f;overwrite=true)
        set_transport!(restored,:multicomponent;data)
        restore_flame!(restored,path)
        @test restored.soret_enabled && restored.profile_grid_policy == :full_knots
        @test restored.grid == f.grid && restored.state == f.state
        @test restored.discretization == :conservative

        # Even failures late in profile validation must leave the target intact.
        saved=profile_grid_test_snapshot(restored)
        for bad in (merge(deepcopy(arrays),Dict("profile_grid_policy_utf8"=>collect(codeunits("bad")))),
                merge(deepcopy(arrays),Dict("profile_temperatures"=>[300.,NaN,300.,300.])),
                merge(deepcopy(arrays),Dict("mass_flux"=>[-.04])),
                merge(deepcopy(arrays),Dict("soret"=>[1])),
                merge(deepcopy(arrays),Dict("discretization_utf8"=>collect(codeunits("finite_difference")))))
            Arrhenius.npzwrite(path,bad)
            @test_throws ArgumentError restore_flame!(restored,path)
            @test isequal(profile_grid_test_snapshot(restored),saved)
        end
        incomplete=deepcopy(arrays);delete!(incomplete,"profile_temperatures")
        Arrhenius.npzwrite(path,incomplete)
        @test_throws ArgumentError restore_flame!(restored,path)
        @test isequal(profile_grid_test_snapshot(restored),saved)
    end
end
