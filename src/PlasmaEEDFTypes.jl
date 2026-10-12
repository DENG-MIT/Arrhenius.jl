# Collision data shared by the independent loader and numerical solver.
@enum ElectronCollisionKind::UInt8 begin
    EffectiveCollision
    ElasticCollision
    ExcitationCollision
    IonizationCollision
    AttachmentCollision
end

struct ElectronCollision
    target::String
    kind::ElectronCollisionKind
    threshold::Float64
    energy::Vector{Float64}
    cross_section::Vector{Float64}
    origin::Symbol
end

"""
    EEDFModel

Fixed electron-energy grid and collision cross-section tables. Construct a model
with [`read_eedf_model`](@ref); energies are eV and cross sections are m².
"""
struct EEDFModel
    energy_edges::Vector{Float64}
    collisions::Vector{ElectronCollision}
    target_names::Vector{String}
end
