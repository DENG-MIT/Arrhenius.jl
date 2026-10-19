using Arrhenius, LinearAlgebra, Test

@testset "premixed Newton correction allocation" begin
    B = 2
    n = 4
    bandwidth = 2B - 1
    band = zeros(6B - 2,n)
    band[4B-1,:] .= 1
    _,pivots = LinearAlgebra.LAPACK.gbtrf!(bandwidth,bandwidth,n,band)
    residual = reshape([1.0,-2.0,3.0,-4.0],B,2)
    correction = zeros(n)
    step = reshape(Arrhenius._flame_solve_correction!(correction,band,pivots,residual,bandwidth),B,2)
    @test step == -residual
    Arrhenius._flame_solve_correction!(correction,band,pivots,residual,bandwidth)
    @test (@allocated Arrhenius._flame_solve_correction!(correction,band,pivots,residual,bandwidth)) == 0
end

@testset "premixed residual and Jacobian allocation" begin
    mechanism = joinpath(@__DIR__,"..","mechanism","h2o2.yaml")
    gas = CreateSolution(mechanism)
    X = mole_fractions(gas,"H2:1.1,O2:1,AR:5")
    flame = solve!(FreeFlame(gas;X,discretization=:finite_difference))
    workspace = Arrhenius.FlameWorkspace(flame)
    residual = similar(flame.state)
    Arrhenius.flame_residual!(residual,flame,flame.state,workspace)
    @test (@allocated Arrhenius.flame_residual!(residual,flame,flame.state,workspace)) == 0
    Arrhenius._flame_jacobian(flame,flame.state,workspace,residual;analytic=false)
    @test (@allocated Arrhenius._flame_jacobian(flame,flame.state,workspace,residual;analytic=false)) == 0
end


@testset "burner and fixed flame residual and Jacobian allocation" begin
    mechanism = joinpath(@__DIR__,"..","mechanism","h2o2.yaml")
    gas = CreateSolution(mechanism)
    X = "H2:1.1,O2:1,AR:5"
    grid = [0.0,0.0025,0.005,0.0075,0.01]
    for fixed in (false,true)
        flame = BurnerFlame(gas;mdot=0.1,T=300.0,X,width=0.01,grid)
        if fixed
            set_temperature_profile!(flame,grid,fill(300.0,length(grid));relative=false)
        end
        workspace = Arrhenius.FlameWorkspace(flame)
        residual = similar(flame.state)
        Arrhenius.flame_residual!(residual,flame,flame.state,workspace)
        @test (@allocated Arrhenius.flame_residual!(residual,flame,flame.state,workspace)) == 0
        Arrhenius._flame_jacobian(flame,flame.state,workspace,residual;analytic=false)
        @test (@allocated Arrhenius._flame_jacobian(flame,flame.state,workspace,residual;analytic=false)) == 0
    end
end

@testset "flame-owned workspace isolates transport configurations" begin
    mechanism = joinpath(@__DIR__,"..","mechanism","h2o2.yaml")
    gas = CreateSolution(mechanism)
    data = MultiTransportData(mechanism*".multicomponent.npz",gas;mechanism)
    X = "H2:1.1,O2:1,AR:5"
    makeflame(model,basis) = FreeFlame(gas;X,transport_model=model,
        multicomponent_data=model==:multicomponent ? data : nothing,
        flux_gradient_basis=basis,discretization=:finite_difference)
    mixed_mole = makeflame(:mixture_averaged,:mole)
    mixed_mass = makeflame(:mixture_averaged,:mass)
    multi_mole = makeflame(:multicomponent,:mole)
    fresh_mole = solve!(makeflame(:mixture_averaged,:mole))
    fresh_mass = solve!(makeflame(:mixture_averaged,:mass))
    fresh_multi = solve!(makeflame(:multicomponent,:mole))
    solve!(mixed_mole)
    first_workspace = Arrhenius._take_flame_workspace(mixed_mole)
    @test length(mixed_mole.workspace) > 1 # adaptive grids keep flame-local workspaces
    solve!(mixed_mass)
    solve!(multi_mole)
    solve!(mixed_mole)
    @test any(w -> w === first_workspace,mixed_mole.workspace)
    @test mixed_mole.state == fresh_mole.state
    @test mixed_mass.state == fresh_mass.state
    @test multi_mole.state == fresh_multi.state
    @test mixed_mole.workspace !== mixed_mass.workspace
    @test mixed_mass.workspace !== multi_mole.workspace
    old_workspace = Arrhenius._take_flame_workspace(mixed_mass)
    set_transport!(mixed_mass,:mixture_averaged;flux_gradient_basis=:mole)
    @test Arrhenius._take_flame_workspace(mixed_mass) !== old_workspace
    @test !any(w -> w === old_workspace,mixed_mass.workspace)
end
