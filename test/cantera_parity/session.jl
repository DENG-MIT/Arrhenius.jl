# Execute diagnostic/test jobs in one Julia process, retaining compiled code.
# Usage: julia --project=... session.jl QUEUE_DIRECTORY
using Arrhenius, LinearAlgebra
BLAS.set_num_threads(1)
queue=abspath(ARGS[1])
empty!(ARGS)
while !isfile(joinpath(queue,"stop"))
    jobs=sort(filter(x->endswith(x,".jl"),readdir(queue)))
    for job in jobs
        path=joinpath(queue,job)
        isfile(path*".done") && continue
        open(path*".log","w") do io
            redirect_stdout(io) do
                redirect_stderr(io) do
                    try
                        Base.include(Main,path)
                        println("JOB_PASS")
                    catch err
                        showerror(io,err,catch_backtrace()); println(io,"\nJOB_FAIL")
                    end
                end
            end
        end
        touch(path*".done")
    end
    sleep(.25)
end
