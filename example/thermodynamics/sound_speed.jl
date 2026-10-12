using Arrhenius

"Native equilibrium/frozen sound speeds (m/s); rows follow the supplied K grid."
function sound_speed_calculation(gas;temperatures=collect(300.:200.:4900.))
    x = mole_fractions(gas,"CH4:1,O2:2,N2:7.52")
    data = zeros(length(temperatures),4)
    final_states = zeros(2+gas.n_species,length(temperatures))
    for (j,T) in enumerate(temperatures)
        result = equilibrium_sound_speeds(gas;T,P=one_atm,X=x)
        data[j,:] = [T,result.equilibrium,result.frozen,result.frozen_at_equilibrium]
        # The published example carries the previous equilibrium composition
        # into the next TP state; it preserves the original elemental inventory.
        x = result.perturbed_state.X
        final_states[:,j] = [result.perturbed_state.T,result.perturbed_state.P,x...]
    end
    return (;data,final_states)
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) == 1 || error("provide gri30_highT.yaml and its .yaml.npz sidecar")
    result = sound_speed_calculation(CreateSolution(ARGS[1]))
    println("T [K], equilibrium [m/s], frozen [m/s], frozen at perturbed equilibrium [m/s]")
    display(result.data)
end
