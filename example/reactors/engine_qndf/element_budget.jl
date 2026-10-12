"""Apply a relative inventory limit plus the actual species-atol induced budget."""
function engine_element_budget(element_matrix, molecular_weights, component_mass_atols,
        inventory_reference, balance_defect, quadrature_change;relative_limit=2e-7)
    absolute_budget=abs.(element_matrix)*(component_mass_atols./molecular_weights)
    reference_scale=abs.(inventory_reference)
    total_budget=relative_limit.*reference_scale.+absolute_budget
    all(isfinite,total_budget) && all(total_budget.>=0) || error("invalid elemental error budget")
    ratio(error,budget)=budget==0 ? (error==0 ? 0.0 : Inf) : abs(error)/budget
    balance_ratio=ratio.(balance_defect,total_budget)
    quadrature_ratio=ratio.(quadrature_change,total_budget)
    return (;absolute_budget,reference_scale,total_budget,balance_ratio,quadrature_ratio,
        balance_pass=all(balance_ratio.<1),quadrature_pass=all(quadrature_ratio.<1),relative_limit)
end
