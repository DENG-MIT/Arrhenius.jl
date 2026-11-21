"""Supertype for whole-mechanism reaction data.

A concrete subtype stores every reaction's stoichiometry, rate-constant and
special-case (third-body/falloff/Troe/PLOG/Blowers-Masel) data for an entire
mechanism as one struct-of-arrays. `Reactions{T}` stores the mechanism in the
YAML/Cantera file's own reaction order; this is currently the only concrete
subtype.
"""
abstract type AbstractReactions{T<:AbstractFloat} end
