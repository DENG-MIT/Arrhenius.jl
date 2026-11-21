"""Supertype for whole-mechanism reaction data.

A concrete subtype stores every reaction's stoichiometry, rate-constant and
special-case (third-body/falloff/Troe/PLOG/Blowers-Masel) data for an entire
mechanism as one struct-of-arrays, in whatever reaction order that subtype
uses. `Reactions{T}` stores the mechanism in the YAML/Cantera file's own
reaction order; `ReorderedReactions{T}` stores the same data bucketed into
contiguous ranges by reaction kind for faster scalar evaluation at scale,
while still answering every original-order query (`rate_multipliers`,
`reaction_rate_constants`, sensitivities, ...) in the mechanism's own order
via `original_reaction_index`/`reordered_reaction_index`.
"""
abstract type AbstractReactions{T<:AbstractFloat} end
