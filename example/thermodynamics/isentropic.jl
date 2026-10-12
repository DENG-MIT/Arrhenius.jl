using Arrhenius, LinearAlgebra

"Native adiabatic nozzle calculation; pressure, temperature and area use SI units."
function isentropic_calculation(gas; points=200)
    T0,P0 = 1200.,10one_atm
    x = mole_fractions(gas,"H2:1,N2:0.1")
    mw = dot(gas.MW,x)
    h0 = cal_hmass_mean(gas,T0,P0,x)
    data = zeros(points,4)
    states = zeros(points,5) # T, P, rho, h, entropy
    for (j,p) in enumerate(range(.01P0,.99P0;length=points))
        state = isentropic_state(gas;T=T0,P=P0,X=x,pressure=p)
        h = cal_hmass_mean(gas,state.T,p,x)
        rho = p*mw/(R*state.T)
        velocity = sqrt(2(h0-h))
        area = 1/(rho*velocity)
        data[j,:] = [area,velocity/frozen_sound_speed(gas;T=state.T,P=p,X=x),state.T/T0,p/P0]
        states[j,:] = [state.T,p,rho,h,cal_smass_mean(gas,state.T,p,x)]
    end
    throat_area = minimum(data[:,1])
    data[:,1] ./= throat_area
    return (;data,states,throat_area)
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) == 1 || error("provide the stock Cantera h2o2.yaml with N2 and its .yaml.npz sidecar")
    result = isentropic_calculation(CreateSolution(ARGS[1]))
    println("area ratio, Mach number, temperature ratio, pressure ratio")
    display(result.data)
end
