module KLUCoreStructuredSolver

using KLU, LinearAlgebra, SparseArrays, SciMLBase
import LinearSolve
const LS = LinearSolve
include(joinpath(@__DIR__, "dense_solver.jl"))

export new_structured_solver, factor_shift!, solve_factored!, solver_report

struct KLUNumericalFailure <: Exception
    message::String
end
Base.showerror(io::IO, err::KLUNumericalFailure) = print(io, err.message)
_numerical_failure(err) = err isa KLUNumericalFailure || err isa SingularException

const KLUFactor = KLU.KLUFactorization{Float64,Int64}

mutable struct KLUStructuredShiftSolver{D}
    ns::Int
    base_nzval::Vector{Float64}
    q::Vector{Float64}
    u::Vector{Float64}
    t::Vector{Float64}
    e::Vector{Float64}
    a::Float64
    shifted::SparseMatrixCSC{Float64,Int64}
    colptr0::Vector{Int64}
    rowval0::Vector{Int64}
    nzval0::Vector{Float64}
    q0::Vector{Float64}
    u0::Vector{Float64}
    t0::Vector{Float64}
    e0::Vector{Float64}
    diagonal_slots::Vector{Int64}
    factor::Union{Nothing,KLUFactor}
    rhs_species::Vector{Float64}
    y::Vector{Float64}
    r::Vector{Float64}
    z::Vector{Float64}
    small_factor::Matrix{Float64}
    small_pivots::Vector{LinearAlgebra.BlasInt}
    small_rhs::Vector{Float64}
    dense::D
    sigma::Float64
    fallback_active::Bool
    factorization_calls::Int
    solve_calls::Int
    fresh_factors::Int
    refactor_attempts::Int
    successful_refactors::Int
    discarded_refactors::Int
    recovery_rebuilds::Int
    auxiliary_solves::Int
    fallback_shift_calls::Int
    fallback_solve_calls::Int
    fallback_reasons::Dict{String,Int}
    rcond::Float64
    rgrowth::Float64
    minimum_rcond::Float64
    minimum_rgrowth::Float64
    fill_ratio::Float64
end

function _vector(d,key,n)
    x=Vector{Float64}(vec(d[key])); length(x)==n || error("$key length mismatch")
    all(isfinite,x) || error("$key contains nonfinite values"); x
end

function new_structured_solver(d)
    q=Vector{Float64}(vec(d["q"])); ns=length(q); ns>0 || error("empty system")
    all(isfinite,q) || error("q contains nonfinite values")
    u=_vector(d,"rank_vector",ns); e=_vector(d,"energy_row",ns)
    tc=_vector(d,"temperature_column",ns+1); t=tc[1:ns]; a=tc[end]
    colptr=Vector{Int64}(vec(d["core_colptr"])); rowval=Vector{Int64}(vec(d["core_rowval"]))
    nzval=Vector{Float64}(vec(d["core_nzval"])); length(colptr)==ns+1 || error("colptr mismatch")
    colptr[1]==1 && colptr[end]==length(nzval)+1 || error("invalid CSC endpoints")
    length(rowval)==length(nzval) && all(isfinite,nzval) || error("invalid CSC values")
    diagonal=Vector{Int64}(undef,ns)
    for j in 1:ns
        prior=0; slot=0
        for p in colptr[j]:colptr[j+1]-1
            i=rowval[p]; 1<=i<=ns && i>prior || error("invalid CSC row order")
            prior=i; i==j && (slot=p)
        end
        slot>0 || error("missing diagonal in column $j"); diagonal[j]=slot
    end
    shifted=SparseMatrixCSC(ns,ns,copy(colptr),copy(rowval),copy(nzval))
    dense=new_dense_solver(d)
    KLUStructuredShiftSolver(ns,copy(nzval),q,u,t,e,a,shifted,
        copy(colptr),copy(rowval),copy(nzval),copy(q),copy(u),copy(t),copy(e),diagonal,
        nothing,zeros(ns),zeros(ns),zeros(ns),zeros(ns),zeros(2,2),
        zeros(LinearAlgebra.BlasInt,2),zeros(2),dense,NaN,false,
        0,0,0,0,0,0,0,0,0,0,Dict{String,Int}(),NaN,NaN,Inf,Inf,NaN)
end

function _assert_immutable(s)
    s.shifted.colptr==s.colptr0 || error("sparse column pattern changed")
    s.shifted.rowval==s.rowval0 || error("sparse row pattern changed")
    isequal(s.base_nzval,s.nzval0) || error("base sparse values changed")
    isequal(s.q,s.q0) && isequal(s.u,s.u0) && isequal(s.t,s.t0) && isequal(s.e,s.e0) ||
        error("border data changed")
end

function _diagnose!(s)
    f=s.factor; f===nothing && throw(KLUNumericalFailure("missing KLU factor"))
    LinearAlgebra.issuccess(f) || throw(KLUNumericalFailure("KLU status is not successful"))
    f.common.status==KLU.KLU_OK || throw(KLUNumericalFailure("KLU status $(f.common.status)"))
    rc=KLU.rcond(f); rg=KLU.rgrowth(f)
    isfinite(rc) && rc>0 || throw(KLUNumericalFailure("nonpositive KLU rcond"))
    isfinite(rg) && rg>0 || throw(KLUNumericalFailure("nonpositive KLU rgrowth"))
    s.rcond=rc; s.rgrowth=rg; s.minimum_rcond=min(s.minimum_rcond,rc)
    s.minimum_rgrowth=min(s.minimum_rgrowth,rg)
end

function _fresh_factor!(s)
    s.factor=KLU.klu(s.shifted)
    s.fresh_factors+=1
    _diagnose!(s)
end

function _refactor_or_rebuild!(s)
    s.refactor_attempts+=1
    try
        KLU.klu!(s.factor,s.shifted.nzval)
        _diagnose!(s)
        s.successful_refactors+=1
    catch err
        _numerical_failure(err) || rethrow()
        s.factor=nothing
        s.discarded_refactors+=1
        s.recovery_rebuilds+=1
        _fresh_factor!(s)
    end
end

function _fallback!(s,sigma,err)
    reason=string(nameof(typeof(err)),": ",sprint(showerror,err))
    s.fallback_reasons[reason]=get(s.fallback_reasons,reason,0)+1
    factor_shift!(s.dense,sigma); s.fallback_active=true; s.fallback_shift_calls+=1
end

function factor_shift!(s::KLUStructuredShiftSolver,sigma)
    isfinite(sigma) && sigma>0 || error("sigma must be finite and positive")
    _assert_immutable(s); s.factorization_calls+=1; s.sigma=Float64(sigma); s.fallback_active=false
    copyto!(s.shifted.nzval,s.base_nzval)
    @inbounds for p in s.diagonal_slots; s.shifted.nzval[p]-=sigma; end
    try
        s.factor===nothing ? _fresh_factor!(s) : _refactor_or_rebuild!(s)
        copyto!(s.r,s.u); KLU.solve!(s.factor,s.r)
        copyto!(s.z,s.t); KLU.solve!(s.factor,s.z); s.auxiliary_solves+=2
        all(isfinite,s.r) && all(isfinite,s.z) || throw(KLUNumericalFailure("nonfinite border solves"))
        s.small_factor[1,1]=1+dot(s.q,s.r); s.small_factor[1,2]=dot(s.q,s.z)
        s.small_factor[2,1]=-dot(s.e,s.r); s.small_factor[2,2]=s.a-sigma-dot(s.e,s.z)
        all(isfinite,s.small_factor) || throw(KLUNumericalFailure("nonfinite reduced system"))
        _,piv,info=LAPACK.getrf!(s.small_factor); info==0 || throw(KLUNumericalFailure("singular reduced system"))
        copyto!(s.small_pivots,piv)
        isnan(s.fill_ratio) && (s.fill_ratio=nnz(s.factor)/nnz(s.shifted))
    catch err
        err isa InterruptException && rethrow(); _numerical_failure(err) || rethrow()
        s.factor=nothing; _fallback!(s,sigma,err)
    end
    nothing
end

function solve_factored!(x,s::KLUStructuredShiftSolver,b)
    length(x)==s.ns+1 && length(b)==s.ns+1 || error("RHS/solution length mismatch")
    isfinite(s.sigma) || error("factor_shift! required"); s.solve_calls+=1
    if s.fallback_active
        solve_factored!(x,s.dense,b); s.fallback_solve_calls+=1; return nothing
    end
    try
        copyto!(s.rhs_species,@view(b[1:s.ns])); KLU.solve!(s.factor,s.rhs_species)
        all(isfinite,s.rhs_species) || throw(KLUNumericalFailure("nonfinite KLU solve"))
        s.small_rhs[1]=dot(s.q,s.rhs_species); s.small_rhs[2]=b[end]-dot(s.e,s.rhs_species)
        LAPACK.getrs!('N',s.small_factor,s.small_pivots,s.small_rhs)
        alpha,theta=s.small_rhs; isfinite(alpha)&&isfinite(theta) || throw(KLUNumericalFailure("nonfinite reduced solution"))
        @inbounds for i in 1:s.ns; x[i]=s.rhs_species[i]-s.r[i]*alpha-s.z[i]*theta; end
        x[end]=theta; all(isfinite,x) || throw(KLUNumericalFailure("nonfinite reconstruction"))
    catch err
        err isa InterruptException && rethrow(); _numerical_failure(err) || rethrow()
        s.factor=nothing; _fallback!(s,s.sigma,err); solve_factored!(x,s.dense,b); s.fallback_solve_calls+=1
    end
    nothing
end

solver_report(s::KLUStructuredShiftSolver)=Dict(
    "backend"=>"KLU.jl fixed-pivot refactor bordered 2x2","factorizations"=>s.factorization_calls,
    "solves"=>s.solve_calls,"fresh_factors"=>s.fresh_factors,"refactor_attempts"=>s.refactor_attempts,
    "successful_refactors"=>s.successful_refactors,"discarded_refactors"=>s.discarded_refactors,
    "recovery_rebuilds"=>s.recovery_rebuilds,"auxiliary_solves"=>s.auxiliary_solves,
    "fallback_shift_calls"=>s.fallback_shift_calls,"fallback_solves"=>s.fallback_solve_calls,
    "fallbacks"=>s.fallback_solve_calls,"fallback_reasons"=>copy(s.fallback_reasons),
    "rcond"=>s.rcond,"rgrowth"=>s.rgrowth,"minimum_rcond"=>s.minimum_rcond,
    "minimum_rgrowth"=>s.minimum_rgrowth,"core_stored_nnz"=>nnz(s.shifted),
    "factor_fill_ratio"=>s.fill_ratio,"dense_fallback_backend"=>string(typeof(s.dense.cache.alg)))

end
