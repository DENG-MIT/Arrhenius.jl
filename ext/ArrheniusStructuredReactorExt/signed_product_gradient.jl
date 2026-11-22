function signed_product_gradient!(g,plan::SignedProductPlan,C,rate)
    n=length(plan.indices)
    length(g)==n||throw(DimensionMismatch("g has length $(length(g)), expected $n"))
    r=plan.repeated
    if plan.kind<=3
        if plan.kind==2
            (C[r[1]]<0&&C[r[2]]<0)&&return fill!(g,zero(rate))
        elseif plan.kind==3
            ((C[r[1]]<0)+(C[r[2]]<0)+(C[r[3]]<0)>=2)&&return fill!(g,zero(rate))
        end
        for i in 1:n
            k=plan.indices[i]
            acc=zero(rate)
            for p in eachindex(r)
                r[p]==k||continue
                term=rate
                for q in eachindex(r)
                    q==p&&continue
                    term*=C[r[q]]
                end
                acc+=term
            end
            g[i]=acc
        end
        return g
    end
    indices=plan.indices; orders=plan.orders
    for j in 1:n
        if orders[j]!=0&&C[indices[j]]<=0
            return fill!(g,zero(rate))
        end
    end
    for i in 1:n
        oi=orders[i]
        if oi==0
            g[i]=zero(rate)
        else
            term=rate*oi*C[indices[i]]^(oi-1)
            for j in 1:n
                (j==i||orders[j]==0)&&continue
                term*=C[indices[j]]^orders[j]
            end
            g[i]=term
        end
    end
    return g
end
