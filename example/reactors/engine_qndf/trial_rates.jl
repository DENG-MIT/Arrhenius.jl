# Exact helper prefix from realgas-signed-model.jl; original SHA256 8636e69785806363013611d63b757dd6b090b04b1959f59450a07415a561085a
# PRIVATE numerical trial extension, validated only for integer-order dodecane.
# The public constructors and all accepted-state gates remain unchanged.
using Arrhenius, LinearAlgebra, SparseArrays, ForwardDiff

function check_trial_reactions(reaction)
    isempty(reaction.blowers_masel.reaction_indices) || error("BM trial rates are unsupported")
    reaction.reactant_orders==reaction.reactant_stoich_coeffs || error("custom orders unsupported")
    all(isinteger,nonzeros(reaction.reactant_orders)) || error("integer reactant orders required")
    all(isinteger,nonzeros(reaction.product_stoich_coeffs)) || error("integer product orders required")
    all(vec(sum(reaction.reactant_orders;dims=1)).<=3) || error("reactant orders above three unsupported")
    all(vec(sum(reaction.product_stoich_coeffs;dims=1))[reaction.is_reversible].<=3) ||
        error("active reverse orders above three unsupported")
    nothing
end

function trial_wdot!(output,reaction,T,C,S0,h,work;activity_concentrations=nothing,kwargs...)
    # Reuse native Arrhenius / equilibrium / collider / falloff / PLOG factors.
    wdot!(output,reaction,T,C,S0,h,work;activity_concentrations,get_rate_constants=true,kwargs...)
    activity=isnothing(activity_concentrations) ? C : activity_concentrations
    reactants,products=reaction.reactant_orders,reaction.product_stoich_coeffs
    ri,rv=rowvals(reactants),nonzeros(reactants)
    pi,pv=rowvals(products),nonzeros(products)
    @inbounds for i in 1:reaction.n_reactions
        forward,reverse=work.kf[i],work.kr[i]
        negative_factors=0.
        for j in nzrange(reactants,i)
            concentration=activity[ri[j]]
            forward*=Arrhenius._concentration_power(concentration,rv[j])
            negative_factors+=concentration<0 ? rv[j] : 0.
        end
        negative_factors>=2 && (forward=zero(forward))
        if reaction.is_reversible[i]
            negative_factors=0.
            for j in nzrange(products,i)
                concentration=activity[pi[j]]
                reverse*=Arrhenius._concentration_power(concentration,pv[j])
                negative_factors+=concentration<0 ? pv[j] : 0.
            end
            negative_factors>=2 && (reverse=zero(reverse))
        end
        work.kf[i],work.kr[i]=forward,reverse
        work.rates_of_progress[i]=forward-reverse
    end
    mul!(output,reaction.vk,work.rates_of_progress)
    output
end


