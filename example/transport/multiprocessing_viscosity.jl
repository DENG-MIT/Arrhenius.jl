# Native counterpart of Cantera's transport/multiprocessing_viscosity.py.
# julia --threads=4 --project=. example/transport/multiprocessing_viscosity.jl GRI30.yaml TRANSPORT.npz
using Arrhenius

function _transport_property_chunk!(out, indices, mechanism, transport, temperatures, property)
    # Each task owns its phase and mutable storage.
    gas = CreateSolution(mechanism)
    data = MultiTransportData(transport,gas;mechanism)
    multi = MultiTransportWorkspace(data)
    mixture = TransportWorkspace(gas)
    X = zeros(gas.n_species)
    for (name,amount) in (("CH4",1.0),("O2",1.0),("N2",3.76))
        index = findfirst(==(name),gas.species_names)
        isnothing(index) && throw(ArgumentError("mechanism must contain $name"))
        X[index] = amount
    end
    X ./= sum(X)
    for i in indices
        out[i] = if property == :conductivity
            multicomponent_thermal_conductivity!(multi,data,gas,one_atm,temperatures[i],X)
        elseif property == :viscosity
            first(mixture_transport!(mixture,gas,one_atm,temperatures[i],X))
        else
            throw(ArgumentError("unknown transport property"))
        end
    end
    return nothing
end

function _transport_property_sweep(mechanism,transport,temperatures,property;workers=1)
    output = Vector{Float64}(undef,length(temperatures))
    if workers == 1
        _transport_property_chunk!(output,eachindex(temperatures),mechanism,transport,temperatures,property)
    else
        @sync for worker in 1:workers
            first_index = fld((worker-1)*length(temperatures),workers)+1
            last_index = fld(worker*length(temperatures),workers)
            Threads.@spawn _transport_property_chunk!(output,first_index:last_index,mechanism,transport,temperatures,property)
        end
    end
    return output
end

"""
    parallel_transport_calculation(mechanism, transport; npoints=5000, workers=4)

Compute multicomponent thermal conductivity [W/(m K)] and viscosity [Pa s]
for CH4:1, O2:1, N2:3.76 at one atmosphere over 300–900 K. Each property is
calculated in parallel and serially, returning all four ordered arrays and
the temperature grid. Each parallel task owns its phase and transport storage;
the serial and parallel calculations initialize their own models.

`mechanism` needs its kinetic/transport sidecar from `export_sidecar.py`;
`transport` is produced by `export_multicomponent.py`. Launch Julia with at
least `workers` threads, for example `julia --threads=4 --project=.`.
"""
function parallel_transport_calculation(mechanism,transport;npoints::Integer=5000,workers::Integer=4)
    npoints > 0 || throw(ArgumentError("positive temperature count required"))
    workers > 0 || throw(ArgumentError("positive worker count required"))
    workers <= Threads.nthreads() || throw(ArgumentError("launch Julia with at least $workers threads"))
    T = npoints == 1 ? [300.0] : collect(range(300.0,900.0;length=npoints))
    conductivity_parallel = _transport_property_sweep(mechanism,transport,T,:conductivity;workers)
    conductivity_serial = _transport_property_sweep(mechanism,transport,T,:conductivity)
    viscosity_parallel = _transport_property_sweep(mechanism,transport,T,:viscosity;workers)
    viscosity_serial = _transport_property_sweep(mechanism,transport,T,:viscosity)
    return (;T,conductivity_parallel,conductivity_serial,viscosity_parallel,viscosity_serial)
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    length(ARGS) == 2 || error("supply GRI30.yaml and its multicomponent transport archive")
    result = parallel_transport_calculation(ARGS[1],ARGS[2])
    println("Computed ",length(result.T)," temperatures for both properties in parallel and serially.")
    println("Final conductivity [W/(m K)]: ",last(result.conductivity_serial))
    println("Final viscosity [Pa s]: ",last(result.viscosity_serial))
end
