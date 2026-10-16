struct SignedProductPlan
    indices::Vector{Int}
    orders::Vector{Float64}
    repeated::Vector{Int}
    kind::Int
end

function signed_product_plans(stoich, orders)
    size(stoich)==size(orders)||error("order/stoichiometry dimensions")
    plans=SignedProductPlan[]
    for j in axes(stoich,2)
        si=rowvals(stoich)[nzrange(stoich,j)]; oi=rowvals(orders)[nzrange(orders,j)]
        indices=sort!(union(si,oi)); s=[stoich[k,j] for k in indices]; o=[orders[k,j] for k in indices]
        all(isfinite,s)&&all(isfinite,o)&&all(>=(0),s)||error("invalid product plan")
        general=length(indices)>3||any(k->!isinteger(s[k])||s[k]!=o[k],eachindex(s))
        repeated=Int[]
        if !general
            for (k,value) in zip(indices,s); append!(repeated,fill(k,Int(value))); end
        end
        kind=!general&&1<=length(repeated)<=3 ? length(repeated) : 4
        push!(plans,SignedProductPlan(indices,o,repeated,kind))
    end
    return plans
end

function signed_multiply(plan::SignedProductPlan,C,rate)
    r=plan.repeated
    if plan.kind==1
        return rate*C[r[1]]
    elseif plan.kind==2
        a,b=C[r[1]],C[r[2]]
        return a<0&&b<0 ? zero(rate) : rate*(a*b)
    elseif plan.kind==3
        a,b,c=C[r[1]],C[r[2]],C[r[3]]
        return (a<0)+(b<0)+(c<0)>=2 ? zero(rate) : rate*(a*b*c)
    end
    output=rate
    for (k,order) in zip(plan.indices,plan.orders)
        if order!=0
            output=C[k]>0 ? output*C[k]^order : zero(output)
        end
    end
    return output
end
