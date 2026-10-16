using LinearAlgebra, Libdl

"Set or verify the actual numerical backends used by a sequential Julia benchmark."
function benchmark_julia_thread_settings(;enforce=true)
    Threads.nthreads()==1 || error("launch Julia with JULIA_NUM_THREADS=1")
    enforce && BLAS.set_num_threads(1)
    result = Dict("julia_threads"=>Threads.nthreads(),"blas_threads"=>BLAS.get_num_threads())
    result["blas_threads"]==1 || error("Julia BLAS must use one thread")
    if Sys.isapple()
        if enforce
            status=ccall((:BLASSetThreading,"/System/Library/Frameworks/Accelerate.framework/Accelerate"),Cint,(Cuint,),1)
            status==0 || error("Accelerate single-thread setting failed")
        end
        result["accelerate_threading_mode"]=Int(ccall((:BLASGetThreading,"/System/Library/Frameworks/Accelerate.framework/Accelerate"),Cuint,()))
        result["accelerate_threading_mode"]==1 || error("Accelerate must use one thread")
    end
    # LinearSolve can load MKL lazily, independently of Julia's default BLAS.
    for path in unique(filter(path->occursin("mkl_rt",lowercase(basename(path))),Libdl.dllist()))
        handle=Libdl.dlopen(path)
        try
            setter=Libdl.dlsym(handle,:MKL_Set_Num_Threads)
            getter=Libdl.dlsym(handle,:MKL_Get_Max_Threads)
            enforce && ccall(setter,Cvoid,(Cint,),1)
            count=Int(ccall(getter,Cint,()))
            result["mkl_threads:"*basename(path)]=count
            count==1 || error("LinearSolve MKL must use one thread")
        finally
            Libdl.dlclose(handle)
        end
    end
    return result
end
