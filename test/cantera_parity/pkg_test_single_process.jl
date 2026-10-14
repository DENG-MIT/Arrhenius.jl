# Run Pkg's complete test suite without creating a second Julia process.
# Pkg prepares its test sandbox; the public test_fn hook evaluates its normal
# test entry point, then signals completion before the subprocess launcher.
using Pkg
struct SingleProcessPkgTestsPassed <: Exception end
root=normpath(joinpath(@__DIR__,"..",".."))
Pkg.activate(root)
try
    withenv("JULIA_PKG_OFFLINE"=>"true","JULIA_PKG_PRECOMPILE_AUTO"=>"0") do
        Pkg.test(;allow_reresolve=false,test_fn=()->begin
            Base.include_string(Main,Pkg.Operations.gen_test_code(root;test_args=``),"pkg_single_process_test")
            throw(SingleProcessPkgTestsPassed())
        end)
    end
catch err
    err isa SingleProcessPkgTestsPassed || rethrow()
    println("PKG_FULL_SUITE_PASS (in-process test_fn; no test subprocess)")
end
