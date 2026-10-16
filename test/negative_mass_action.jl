using Arrhenius, LinearAlgebra, SparseArrays, Test

@testset "Cantera negative-concentration mass-action rules" begin
    # Cantera C2 repeats a species factor for 2A. It preserves positive and
    # single-negative first-order factors, but a two-negative product is zero.
    c2 = sparse([1], [1], [2.0], 1, 1)
    o2 = copy(c2)
    inds2, orders2 = rowvals(o2), nonzeros(o2)
    c1=sparse([1],[1],[1.0],1,1);o1=copy(c1)
    @test Arrhenius._mass_action_product([-1e-14],rowvals(o1),nonzeros(o1),c1,1,
        nzrange(o1,1)) == -1e-14
    for y in (0.0, 1e-25, -1e-25, -1e-14, -1e-7)
        C = [y]
        got = Arrhenius._mass_action_product(C, inds2, orders2, c2, 1, nzrange(o2, 1))
        expected = y > 0 ? y^2 : 0.0
        @test got == expected
    end

    # For C2 with two distinct species, a single negative factor stays signed;
    # two negative factors are suppressed, matching Cantera StoichManager.h.
    c2distinct = sparse([1, 2], [1, 1], [1.0, 1.0], 2, 1)
    o2distinct = copy(c2distinct)
    inds, orders = rowvals(o2distinct), nonzeros(o2distinct)
    @test Arrhenius._mass_action_product([-1.0, 2.0], inds, orders,
        c2distinct, 1, nzrange(o2distinct, 1)) == -2.0
    @test Arrhenius._mass_action_product([-1.0, -2.0], inds, orders,
        c2distinct, 1, nzrange(o2distinct, 1)) == 0.0
    c3 = sparse([1, 2, 3], [1, 1, 1], ones(3), 3, 1)
    o3 = copy(c3)
    @test Arrhenius._mass_action_product([-1.0, 2.0, 3.0], rowvals(o3),
        nonzeros(o3), c3, 1, nzrange(o3, 1)) == -6.0
    @test Arrhenius._mass_action_product([-1.0, -2.0, 3.0], rowvals(o3),
        nonzeros(o3), c3, 1, nzrange(o3, 1)) == 0.0
    # The analytic mass-action derivative follows the same active product branch.
    for (stoich, orders_matrix, values) in ((c2distinct,o2distinct,[-1.0,2.0]),
            (c3,o3,[-1.0,2.0,3.0]))
        gradient=zeros(length(values))
        Arrhenius._mass_action_product_gradient!(gradient,values,rowvals(orders_matrix),
            nonzeros(orders_matrix),stoich,1,nzrange(orders_matrix,1),1.0)
        h=1e-6
        for k in eachindex(values)
            plus=copy(values);minus=copy(values);plus[k]+=h;minus[k]-=h
            fd=(Arrhenius._mass_action_product(plus,rowvals(orders_matrix),
                    nonzeros(orders_matrix),stoich,1,nzrange(orders_matrix,1)) -
                Arrhenius._mass_action_product(minus,rowvals(orders_matrix),
                    nonzeros(orders_matrix),stoich,1,nzrange(orders_matrix,1)))/(2h)
            @test gradient[k] ≈ fd rtol=2e-8 atol=2e-8
        end
    end
    suppressed=zeros(3)
    Arrhenius._mass_action_product_gradient!(suppressed,[-1.0,-2.0,3.0],rowvals(o3),
        nonzeros(o3),c3,1,nzrange(o3,1),1.0)
    @test suppressed == zeros(3)


    # Cantera's general stoichiometry path suppresses any nonpositive
    # participating concentration, including integer custom orders.
    c4 = sparse(collect(1:4), ones(Int, 4), ones(4), 4, 1)
    o4 = copy(c4)
    @test Arrhenius._mass_action_product([1.0, 2.0, 3.0, 4.0], rowvals(o4),
        nonzeros(o4), c4, 1, nzrange(o4, 1)) == 24.0
    @test Arrhenius._mass_action_product([1.0, 2.0, -3.0, 4.0], rowvals(o4),
        nonzeros(o4), c4, 1, nzrange(o4, 1)) == 0.0

    custom_orders = sparse([1], [1], [1.0], 1, 1)
    custom_stoich = spzeros(1, 1)
    fractional_orders=sparse([1],[1],[0.5],1,1)
    for y in (0.0,1e-25,-1e-25,-1e-14)
        value=Arrhenius._mass_action_product([y],rowvals(fractional_orders),
            nonzeros(fractional_orders),custom_stoich,1,nzrange(fractional_orders,1))
        @test value == (y>0 ? sqrt(y) : 0.0)
    end
    fractional_gradient=zeros(1)
    Arrhenius._mass_action_product_gradient!(fractional_gradient,[1e-25],
        rowvals(fractional_orders),nonzeros(fractional_orders),custom_stoich,1,
        nzrange(fractional_orders,1),1.0)
    @test fractional_gradient[1] ≈ 0.5/sqrt(1e-25) rtol=1e-15

    @test Arrhenius._mass_action_product([-1e-14], rowvals(custom_orders),
        nonzeros(custom_orders), custom_stoich, 1, nzrange(custom_orders, 1)) == 0.0
    # Exercise the requested radical Y edges through ideal-gas Y -> X -> C
    # conversion, with all absent species zero.
    gas = CreateSolution(joinpath(@__DIR__, "data", "nakamura_2017.yaml"))
    nh = findfirst(==("NH"), gas.species_names)
    n2h2 = findfirst(==("N2H2"), gas.species_names)
    n2 = findfirst(==("N2"), gas.species_names)
    reaction = gas.reaction
    reaction115 = findfirst(i -> reaction.reactant_stoich_coeffs[n2h2,i] == 1 &&
        reaction.product_stoich_coeffs[nh,i] == 2, 1:reaction.n_reactions)
    @test reaction115 !== nothing
    Cprobe = zeros(gas.n_species); Cprobe[nh] = -1.0
    legacy = KineticsWorkspace(reaction); fill!(legacy.kf,0); fill!(legacy.kr,0)
    legacy.kr[reaction115] = 1
    Arrhenius._mass_action!(legacy,reaction,Cprobe)
    @test legacy.rates_of_progress[reaction115] == 0.0
    cantera = KineticsWorkspace(reaction); fill!(cantera.kf,0); fill!(cantera.kr,0)
    cantera.kr[reaction115] = 1
    Arrhenius._mass_action!(cantera,reaction,Cprobe)
    @test cantera.rates_of_progress[reaction115] == 0.0
    P, T = 101325.0, 300.0
    for y in (0.0, 1e-25, -1e-25, -1e-14, -1e-7)
        Y = zeros(gas.n_species)
        Y[nh], Y[n2] = y, 1-y
        inverseMW = dot(Y, 1 ./ gas.MW)
        xnh = Y[nh] / gas.MW[nh] / inverseMW
        cnh = P / (Arrhenius.R*T) * xnh
        got = Arrhenius._mass_action_product([cnh], inds2, orders2, c2, 1, nzrange(o2, 1))
        @test got == (y > 0 ? cnh^2 : 0.0)
    end

    # Cantera's NoNorm density uses 1/sum(Y/MW), including when sum(Y) is not
    # one. Transport keeps signed midpoint Y; only X is floored to Tiny.
    MW = [gas.MW[nh], gas.MW[n2]]
    probe = (gas=(MW=MW,), pressure=P)
    Yraw = [0.7, -0.1]
    density = Arrhenius._flame_density(probe, T, Yraw).density
    @test density ≈ P/(dot(Yraw, 1 ./ MW)*Arrhenius.R*T) rtol=2e-15
    density_nonorm = Arrhenius._flame_density(probe, T, Yraw).density
    @test density_nonorm ≈ P / (dot(Yraw, 1 ./ MW) * Arrhenius.R*T) rtol=2e-15
    density_clipped = Arrhenius._flame_density(probe, T, Yraw; clip_species=true).density
    Yclipped = max.(Yraw, 0.0)
    @test density_clipped ≈ P*sum(Yclipped) / (dot(Yclipped, 1 ./ MW) * Arrhenius.R*T) rtol=2e-15

    for y in (0.0, 1e-25, -1e-25, -1e-14, -1e-7)
        Yedge = [y, 1-y]
        ymid, xmid = zeros(2), zeros(2)
        meanMW = Arrhenius._flame_midpoint_composition!(ymid, xmid,
            Yedge, Yedge, MW, false)
        inverseMW = dot(Yedge, 1 ./ MW)
        @test ymid == Yedge
        @test meanMW ≈ 1/inverseMW
        @test xmid[1] ≈ (y/MW[1])/inverseMW rtol=1e-14 atol=0.0
        @test sum(xmid) ≈ 1.0 rtol=1e-14
        @test max(xmid[1], 1e-20) == 1e-20
    end
    ymid, xmid = zeros(2), zeros(2)
    Arrhenius._flame_midpoint_composition!(ymid, xmid, [-1e-7, 1+1e-7],
        [-1e-7, 1+1e-7], MW, true)
    @test ymid == [0.0, 1.0] # missing radical after positive-species fallback

end
