struct DenseShiftSolver{C}
    cache::C
    original::Matrix{Float64}
    A::Matrix{Float64}
    b::Vector{Float64}
    factorizations::Base.RefValue{Int}
    solves::Base.RefValue{Int}
end
function new_dense_solver(d)
    J=Matrix{Float64}(d["J"]); A=copy(J); b=zeros(size(J,1))
    alg=Sys.isapple() ? LS.AppleAccelerateLUFactorization() : LS.MKLLUFactorization()
    cache=LS.init(SciMLBase.LinearProblem(A,b),alg;
        alias=LS.LinearAliasSpecifier(alias_A=true,alias_b=true))
    cache.A===A || error("matrix alias contract");cache.b===b || error("RHS alias contract")
    DenseShiftSolver(cache,J,A,b,Ref(0),Ref(0))
end
function factor_shift!(s::DenseShiftSolver,sigma)
    copyto!(s.A,s.original)
    for i in axes(s.A,1);s.A[i,i]-=sigma;end
    LS.reinit!(s.cache;A=s.A,b=s.b)
    s.cache.isfresh || error("freshness trigger failed")
    return nothing
end
function solve_factored!(x,s::DenseShiftSolver,b)
    copyto!(s.b,b);wasfresh=s.cache.isfresh
    LS.reinit!(s.cache;b=s.b,reuse_precs=true)
    s.cache.isfresh==wasfresh || error("RHS update changed matrix freshness")
    sol=LS.solve!(s.cache)
    SciMLBase.successful_retcode(sol) || error("dense linear solve failed")
    !s.cache.isfresh || error("factorization was not reused")
    copyto!(x,sol.u);s.factorizations[]+=wasfresh;s.solves[]+=1
    return nothing
end
solver_report(s::DenseShiftSolver)=Dict("backend"=>string(typeof(s.cache.alg)),
    "factorizations"=>s.factorizations[],"solves"=>s.solves[],"fallbacks"=>0)
