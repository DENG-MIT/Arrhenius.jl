"""Parse Cantera's level-2 Newton log into per-timestep Jacobian counts."""
function cantera_stat_total(stats, key)
    values = Int.(get(stats, key, Int[]))
    isempty(values) && return 0
    # SteadyStateSystem::saveStats stores m_jac->nEvals() cumulatively, while
    # m_nsteps is reset after each save and must be summed across snapshots.
    key == "jacobian_evals" && return last(values)
    sum(values)
end

function cantera_attempt_trace(path)
    rows = Dict{String,Any}[]
    current = nothing
    last_jac = 0
    pending_refreshes = 0
    for line in eachline(path)
        if occursin(r"Timestep \(\d+\) dt=", line)
            dt = first(split(strip(split(line, "dt=")[2])))
            current = Dict{String,Any}("dt" => parse(Float64, dt),
                "jacobians" => 0, "jacobians_cumulative" => last_jac,
                "jacobian_age_start" => nothing)
            push!(rows, current)
            continue
        end
        if current !== nothing && (occursin("Re-evaluating Jacobian", line) ||
                occursin("Maximum Jacobian age reached", line))
            current["jacobians"] += 1
            pending_refreshes += 1
            current["jacobian_age_start"] === nothing &&
                occursin("Maximum Jacobian age reached", line) &&
                (current["jacobian_age_start"] = 0)
        end
        if current !== nothing && occursin(r"Timestep \(\d+\) failed", line)
            current["accepted"] = false
            current = nothing
        elseif current !== nothing && occursin(r"Timestep \(\d+\) succeeded", line)
            current["accepted"] = true
            current = nothing
        end
        cols = split(strip(line))
        if length(cols) == 8 && tryparse(Int, cols[1]) !== nothing &&
                tryparse(Float64, cols[2]) !== nothing && occursin("/", cols[end])
            njac = tryparse(Int, cols[7])
            age = tryparse(Int, first(split(cols[8], "/")))
            if njac !== nothing && age !== nothing
                delta = max(0, njac - last_jac)
                if current !== nothing
                    newly_observed = max(0, delta - pending_refreshes)
                    current["jacobians"] += newly_observed
                    current["jacobians_cumulative"] = njac
                    if current["jacobian_age_start"] === nothing
                        current["jacobian_age_start"] = newly_observed > 0 ? 0 : age - 1
                    end
                end
                pending_refreshes = max(0, pending_refreshes - delta)
                last_jac = njac
            end
        end
    end
    return rows
end
