# YAML species selection, element metadata and thermodynamic model reader.
function read_species_basics(yaml)
    n_species=length(yaml["phases"][1]["species"])
    n_reactions=length(get(yaml,"reactions",Any[]))
    species_names=yaml["phases"][1]["species"]
    elements=yaml["phases"][1]["elements"]
    n_elements=length(elements)
    ele_matrix=zeros(n_elements,n_species)
    _species_names=[yaml["species"][i]["name"] for i in eachindex(yaml["species"])]
    for (i,species) in enumerate(species_names)
        spec=yaml["species"][findfirst(==(species),_species_names)]
        for j in 1:n_elements
            haskey(spec["composition"],elements[j]) && (ele_matrix[j,i]=spec["composition"][elements[j]])
        end
    end
    n_species,n_reactions,species_names,elements,n_elements,ele_matrix
end

"""Read selected-phase species metadata and construct the YAML thermo model."""
function _read_gas_thermo(yaml)
    n_species,n_reactions,species_names,elements,n_elements,ele_matrix=read_species_basics(yaml)
    model=lowercase(String(yaml["phases"][1]["thermo"]))
    thermo = model=="ideal-gas" ? IdealGasThermo(yaml) :
        (model=="redlich-kwong" ? RedlichKwongThermo(yaml) :
         getfield(@__MODULE__,Symbol(replace(model,"-"=>"_"),:Thermo))(yaml))
    thermo,n_species,n_reactions,species_names,elements,n_elements,ele_matrix
end
