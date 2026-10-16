using Arrhenius
length(ARGS) in (1,2) || error("usage: julia --project scripts/audit_nh3.jl <mechanism-dir> [output.md]")
root=abspath(ARGS[1]); output=length(ARGS)==2 ? abspath(ARGS[2]) : nothing
files=("alzueta_2023.yaml","creck_nh3.yaml","glarborg_2018.yaml","h2o2.yaml","nakamura_2017.yaml","shrestha_2020.yaml","stagni_2023.yaml")
reports=[audit_mechanism(joinpath(root,f)) for f in files]
lines=["# NH3 mechanism audit", "", "Audit checks elemental and electron-count balance, duplicate-reaction candidates, unused and missing-transport species, NASA7/NASA9 continuity, negative A declarations, hard-sphere collision screens, and broad transport parameter ranges.", "", "Collision-limit and transport-range findings are screening heuristics, not Cantera validation verdicts. Duplicate groups are stoichiometric candidates; exact Cantera duplicate semantics also depend on reaction and rate type and third-body efficiencies.", "", "| Mechanism | Species | Reactions | Error | Warning | Info |", "|---|---:|---:|---:|---:|---:|"]
for (f,r) in zip(files,reports)
    level(l)=count(x->x.level==l,r.issues)
    push!(lines,"| $f | $(r.species) | $(r.reactions) | $(level(:error)) | $(level(:warning)) | $(level(:info)) |")
end
push!(lines,"", "## Findings", "")
for (f,r) in zip(files,reports)
    push!(lines,"### $f","")
    if isempty(r.issues)
        push!(lines,"No findings.","")
    else
        for x in r.issues
            details=join((string(k,"=",v) for (k,v) in pairs(x) if k ∉ (:level,:check)), ", ")
            push!(lines,"- **$(x.level), $(x.check):** $details.")
        end
        push!(lines,"")
    end
end
report=join(lines,"\n")*"\n"
if output===nothing
    print(report)
else
    write(output,report)
    println("wrote ",output)
end
