module PlasmaEEDFLoaderTests
using Test
using Arrhenius
using Arrhenius: ElectronCollisionKind, EffectiveCollision, ElasticCollision,
    ExcitationCollision, IonizationCollision, AttachmentCollision

function write_tmp(content)
    path = joinpath(mktempdir(), "case.yaml")
    write(path, content)
    return path
end

const BASE = """
phases:
- name: gas
  thermo: plasma
  electron-energy-distribution:
    type: Boltzmann-two-term
    energy-levels: [0.0, 1.0, 2.0, 3.0]
"""

@testset "plasma EEDF loader" begin

    @testset "root entries and grid" begin
        path = write_tmp(BASE * """
electron-collisions:
- target: N2
  kind: effective
  threshold: 0.0
  energy-levels: [0.0, 0.5, 1.0]
  cross-sections: [1.0e-20, 2.0e-20, 3.0e-20]
- target: O2
  kind: excitation
  threshold: 0.977
  energy-levels: [0.0, 1.0, 2.0]
  cross-sections: [0.0, 1.0e-22, 2.0e-22]
""")
        m = read_eedf_model(path)
        @test m.energy_edges == [0.0, 1.0, 2.0, 3.0]
        @test length(m.collisions) == 2
        @test m.target_names == ["N2", "O2"]
        c1 = m.collisions[1]
        @test c1.kind === EffectiveCollision
        @test c1.threshold == 0.0  # root zero threshold stays zero
        @test c1.origin === :root
        @test c1.energy == [0.0, 0.5, 1.0]
        @test c1.cross_section == [1.0e-20, 2.0e-20, 3.0e-20]
        @test m.collisions[2].kind === ExcitationCollision
        @test m.collisions[2].threshold == 0.977
        # phase keyword selects by name
        @test read_eedf_model(path; phase="gas").energy_edges == m.energy_edges
        @test_throws ArgumentError read_eedf_model(path; phase="nope")
    end

    @testset "plasma reaction kinds and threshold inference" begin
        path = write_tmp(BASE * """
reactions:
- equation: N2 + Electron => N2+ + 2 Electron
  type: electron-collision-plasma
  energy-levels: [0.0, 15.6, 16.0]
  cross-sections: [0.0, 0.0, 1.0e-22]
  threshold: 0.0
- equation: O2 + Electron => O2-
  type: electron-collision-plasma
  energy-levels: [0.0, 0.06, 0.10]
  cross-sections: [0.0, 1.0e-41, 0.0]
  threshold: 0.0
- equation: O2 + Electron => O2 + Electron
  type: electron-collision-plasma
  energy-levels: [0.0, 1.0, 2.0]
  cross-sections: [1.0e-20, 1.0e-20, 1.0e-20]
  threshold: 0.0
- equation: O2 + Electron => O2(a1) + Electron
  type: electron-collision-plasma
  energy-levels: [0.0, 0.977, 2.0]
  cross-sections: [0.0, 1.0e-22, 2.0e-22]
  threshold: 0.977
- equation: O2 + M => O + O + M
  type: three-body
  rate-constant: {A: 1.0, b: 0.0, Ea: 0.0}
""")
        m = read_eedf_model(path)
        @test length(m.collisions) == 4
        kinds = [c.kind for c in m.collisions]
        @test kinds == [IonizationCollision, AttachmentCollision,
                        EffectiveCollision, ExcitationCollision]
        # reaction zero threshold -> first strictly positive energy level
        @test m.collisions[1].threshold == 15.6
        @test m.collisions[2].threshold == 0.06
        @test m.collisions[3].threshold == 0.0
        @test m.collisions[4].threshold == 0.977
        @test all(c.origin === :reaction for c in m.collisions)
        # ionic "+" preserved by spaced-" + " splitting
        @test m.target_names == ["N2", "O2"]
        # duplicate background records for O2 must be rejected
        bad = write_tmp(BASE * """
reactions:
- equation: O2 + Electron => O2 + Electron
  type: electron-collision-plasma
  energy-levels: [0.0, 1.0, 2.0]
  cross-sections: [1.0e-20, 1.0e-20, 1.0e-20]
- equation: O2 + Electron => O2 + Electron
  type: electron-collision-plasma
  energy-levels: [0.0, 1.0, 2.0]
  cross-sections: [1.0e-20, 1.0e-20, 1.0e-20]
""")
        @test_throws ArgumentError read_eedf_model(bad)
    end

    @testset "phase sections, imports, and renamed electron" begin
        mechanism_dir = mktempdir()
        data_dir = mktempdir()
        imported = joinpath(data_dir, "catalog.yaml")
        write(imported, """
species:
- {name: carrier, composition: {E: 1}}
- {name: N2, composition: {N: 2}}
- {name: O2, composition: {O: 2}}
- {name: N2+, composition: {N: 2, E: -1}}
- {name: O2-, composition: {O: 2, E: 1}}
reactions:
- id: imported-attachment
  equation: O2 + carrier => O2-
  type: electron-collision-plasma
  energy-levels: [0.0, 0.2, 1.0]
  cross-sections: [0.0, 1.0e-22, 0.0]
- id: imported-unused
  equation: N2 + carrier => N2+ + 2 carrier
  type: electron-collision-plasma
  energy-levels: [0.0, 9.0, 10.0]
  cross-sections: [0.0, 0.0, 1.0e-22]
""")
        path = joinpath(mechanism_dir, "case.yaml")
        write(path, """
phases:
- name: selected
  species:
  - catalog.yaml/species: all
  reactions:
  - collisions: [elastic, ion]
  - catalog.yaml/reactions: [imported-attachment]
  electron-energy-distribution:
    type: Boltzmann-two-term
    energy-levels: [0.0, 1.0, 2.0, 3.0]
- name: other
  species:
  - catalog.yaml/species: all
  reactions:
  - unused: all
  electron-energy-distribution:
    type: Boltzmann-two-term
    energy-levels: [0.0, 2.0, 4.0]
collisions:
- id: elastic
  equation: N2 + carrier => carrier + N2
  type: electron-collision-plasma
  kind: elastic
  threshold: 0.0
  energy-levels: [0.0, 1.0, 2.0]
  cross-sections: [1.0e-20, 1.0e-20, 1.0e-20]
- id: ion
  equation: N2 + carrier => N2+ + 2 carrier
  type: electron-collision-plasma
  energy-levels: [0.0, 12.0, 13.0]
  cross-sections: [0.0, 0.0, 1.0e-22]
unused:
- id: unused-local
  equation: N2 + carrier => N2 + carrier
  type: electron-collision-plasma
  energy-levels: [0.0, 1.0]
  cross-sections: [1.0e-20, 1.0e-20]
electron-collisions:
- target: O2
  kind: effective
  energy-levels: [0.0, 1.0]
  cross-sections: [1.0e-20, 1.0e-20]
""")

        @test_throws ArgumentError read_eedf_model(path; phase="selected")
        model = read_eedf_model(path; phase="selected", data_paths=[data_dir])
        @test model.target_names == ["O2", "N2"]
        @test [c.kind for c in model.collisions] ==
            [EffectiveCollision, ElasticCollision, IonizationCollision, AttachmentCollision]
        @test [c.threshold for c in model.collisions] == [0.0, 0.0, 12.0, 0.2]
        @test [c.origin for c in model.collisions] == [:root, :reaction, :reaction, :reaction]
        @test length(read_eedf_model(path; phase="other",
                                     data_paths=[data_dir]).collisions) == 2

        all_path = joinpath(mechanism_dir, "all.yaml")
        write(all_path, replace(read(path, String),
            "- collisions: [elastic, ion]\n  - catalog.yaml/reactions: [imported-attachment]" =>
            "- collisions: all"))
        @test length(read_eedf_model(all_path; phase="selected",
                                     data_paths=[data_dir]).collisions) == 3

        bad_selector = joinpath(mechanism_dir, "bad-selector.yaml")
        write(bad_selector, replace(read(path, String),
            "collisions: [elastic, ion]" => "collisions: elastic"))
        @test_throws ArgumentError read_eedf_model(bad_selector; phase="selected",
                                                   data_paths=[data_dir])

        reversible = joinpath(mechanism_dir, "reversible.yaml")
        write(reversible, replace(read(path, String),
            "N2 + carrier => carrier + N2" => "N2 + carrier <=> carrier + N2"))
        @test_throws ArgumentError read_eedf_model(reversible; phase="selected",
                                                   data_paths=[data_dir])

        write(imported, replace(read(imported, String), "species:\n" =>
            "species:\n- {name: second-electron, composition: {E: 1}}\n"))
        @test_throws ArgumentError read_eedf_model(path; phase="selected",
                                                   data_paths=[data_dir])
    end

    @testset "validation failures" begin
        @test_throws ArgumentError read_eedf_model(write_tmp("""
phases:
- name: gas
  electron-energy-distribution:
    type: isotropic
    energy-levels: [0.0, 1.0, 2.0]
electron-collisions:
- {target: N2, kind: effective, energy-levels: [0.0, 1.0], cross-sections: [1.0, 2.0]}
"""))  # unsupported EEDF type

        @test_throws ArgumentError read_eedf_model(write_tmp("""
phases:
- name: gas
  electron-energy-distribution: {type: Boltzmann-two-term, energy-levels: [0.0, 2.0, 1.0]}
electron-collisions:
- {target: N2, kind: effective, energy-levels: [0.0, 1.0], cross-sections: [1.0, 2.0]}
"""))  # grid not increasing

        @test_throws ArgumentError read_eedf_model(write_tmp(BASE * """
electron-collisions:
- {target: N2, kind: banana, energy-levels: [0.0, 1.0], cross-sections: [1.0, 2.0]}
"""))  # unknown kind

        @test_throws ArgumentError read_eedf_model(write_tmp(BASE * """
electron-collisions:
- {target: N2, kind: effective, energy-levels: [0.0, 0.5, 1.0], cross-sections: [1.0, 2.0]}
"""))  # length mismatch

        @test_throws ArgumentError read_eedf_model(write_tmp(BASE * """
electron-collisions:
- {target: N2, kind: effective, energy-levels: [0.0, 0.5, 0.5], cross-sections: [1.0, 2.0, 3.0]}
"""))  # energies not strictly increasing

        @test_throws ArgumentError read_eedf_model(write_tmp(BASE * """
electron-collisions:
- {target: N2, kind: effective, energy-levels: [0.0, 1.0], cross-sections: [1.0, -2.0]}
"""))  # negative cross section

        @test_throws ArgumentError read_eedf_model(write_tmp(BASE * """
electron-collisions:
- {target: N2, kind: effective, energy-levels: [0.0, 1.0], cross-sections: [1.0, 2.0]}
- {target: N2, kind: effective, energy-levels: [0.0, 1.0], cross-sections: [1.0, 2.0]}
"""))  # duplicate effective for one target

        @test_throws ArgumentError read_eedf_model(write_tmp(BASE * """
reactions:
- equation: N2 + O2 + Electron => N2+ + O2+ + 2 Electron
  type: electron-collision-plasma
  energy-levels: [0.0, 1.0, 2.0]
  cross-sections: [0.0, 1.0, 2.0]
"""))  # more than one neutral target

        @test_throws ArgumentError read_eedf_model(write_tmp(BASE))  # no targets at all
    end

    @testset "strict input boundaries" begin
        collision = "electron-collisions:\n- {target: N2, kind: effective, energy-levels: [0.0, 1.0], cross-sections: [1.0, 2.0]}\n"
        @test_throws ArgumentError read_eedf_model(write_tmp(replace(BASE, "[0.0, 1.0, 2.0, 3.0]"=>"[0.0, 1.0, 1.0, 3.0]") * collision))
        @test_throws ArgumentError read_eedf_model(write_tmp(BASE * replace(collision,"[0.0, 1.0]"=>"[-1.0, 1.0]")))
        @test_throws ArgumentError read_eedf_model(write_tmp(BASE * "electron-collisions: [4]\n"))
        @test_throws ArgumentError read_eedf_model(write_tmp(BASE * collision * "- {target: N2, kind: elastic, energy-levels: [0.0, 1.0], cross-sections: [1.0, 2.0]}\n"))
    end

    @testset "reference model" begin
        path = get(ENV, "ARRHENIUS_EEDF_MODEL", "")
        if isempty(path)
            @info "ARRHENIUS_EEDF_MODEL not set; skipping reference model test"
            @test_skip false
        else
            m = read_eedf_model(path)
            counts = Dict(k => count(c -> c.kind === k, m.collisions)
                          for k in instances(ElectronCollisionKind))
            @test length(m.collisions) == 43
            @test counts[EffectiveCollision] == 2
            @test counts[ExcitationCollision] == 36
            @test counts[IonizationCollision] == 3
            @test counts[AttachmentCollision] == 2
        end
    end
end

end # module
