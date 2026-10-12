module IonizedFlameTests
using Arrhenius,Test,LinearAlgebra
using ..IonTransportTests: fixture
const A=Arrhenius
@testset "ionized flame state and Jacobian contracts" begin
    mktempdir() do dir
        gas,_,_=fixture(dir)
        f=FreeFlame(gas;T=900.,X="Ar:0.7,He:0.3",grid=[0,.001,.003,.007,.01])
        n=gas.n_species;B=n+3;N=length(f.grid)
        @test f isa IonizedFlame
        @test !f.field_enabled
        u=f.state
        for j in 1:N
            u[:,j].=0
            u[1,j]=.85+.025*j
            u[2,j]=.7+.001*j
            u[3,j]=.3-.001*j
            u[4,j]=(-1)^j*1e-11
            u[5,j]=(-1)^j*1e-15
            u[end-1,j]=.1*j
            u[end,j]=.4+.01*j
        end
        set_electric_field!(f,true)
        @test electric_field(f)==1000 .* vec(u[end-1,:])
        @test mass_fractions(f)==u[2:n+1,:]
        @test velocity(f)==vec(u[end,:])
        @test A._flame_bounds(f,5,B)==(-1e-14,1.)
        @test A._flame_bounds(f,4,B)==(-1e-10,1.)
        for transient in (false,true)
            weights=A._flame_correction_weights(f,u,transient)
            @test weights[5]≈1e-5*sum(abs,u[5,:])/N+1e-20
            @test weights[4]≈1e-5*sum(abs,u[4,:])/N+1e-16
            @test weights[end-1]≈1e-4*sum(abs,u[end-1,:])/N+(transient ? 1e-11 : 1e-9)/1000
        end
        w=A.FlameWorkspace(f);res=similar(u);flame_residual!(res,f,u,w)
        @test mole_fractions(f)≈w.X
        @test density(f)≈w.rho
        @test heat_release_rate(f)≈-vec(sum(w.h.*w.source;dims=1))
        @test set_transport!(f,:ionized_gas)===f
        @test_throws ArgumentError set_transport!(f,:mixture_averaged)
        @test_throws ArgumentError set_transport!(f,:ionized_gas;soret=true)
        old=copy(u);band=copy(A._flame_jacobian(f,u,w,res))
        @test u==old
        bw=2B-1
        J=zeros(B*N,B*N)
        for col in 1:B*N,row in max(1,col-bw):min(B*N,col+bw)
            J[row,col]=band[2bw+1+row-col,col]
        end
        @test all(isfinite,J)
        for j in 2:N
            S=sum(u[k+1,j]/gas.MW[k] for k in 1:n)
            Q=sum(f.ion_data.charges[k]*u[k+1,j]/gas.MW[k] for k in 1:n)
            fac=(f.grid[j]-f.grid[j-1])*A._ION_FLAME_FARADAY*w.rho[j]/(1000*A._ION_FLAME_EPS0)
            row=(j-1)*B+B-1
            @test J[row,(j-1)*B+1]≈fac*Q/u[1,j] rtol=3e-5 atol=1e-6
            for k in 1:n
                exact=-fac/gas.MW[k]*(f.ion_data.charges[k]-Q/S)
                @test J[row,(j-1)*B+k+1]≈exact rtol=3e-5 atol=1e-5
            end
            @test J[row,(j-1)*B+B-1]≈1 rtol=1e-8
            @test J[row,(j-2)*B+B-1]≈-1 rtol=1e-8
        end
        d=copy(u)
        d[1,:].=.1;d[2,:].=.02;d[3,:].=-.01
        d[4,:].=1e-9;d[5,:].=1e-13;d[end-1,:].=range(-.2,.3,length=N);d[end,:].=.05
        rp=similar(u);rm=similar(u);h=1e-4
        flame_residual!(rp,f,u+h*d,w;update_transport=false)
        flame_residual!(rm,f,u-h*d,w;update_transport=false)
        fd=(rp-rm)/(2h);jv=reshape(J*vec(d),B,N)
        # The Jacobian uses finite one-sided secants. Scale their error by the
        # gross directional contributions, which may cancel in the net product.
        gross=reshape(abs.(J)*abs.(vec(d)),B,N)
        for k in 1:B
            @test norm(jv[k,:]-fd[k,:],Inf)<=1e-4*max(norm(gross[k,:],Inf),1e-20)
        end
        @test u==old
        f.converged=true;set_electric_field!(f,false)
        @test !f.converged
        flame_residual!(res,f,u,w)
        @test all(iszero,w.flux[3:4,:])
        @test res[end-1,:]==u[end-1,:]
        @test_throws ArgumentError FreeFlame(gas;X="Ar:1",soret=true)
        @test_throws ArgumentError FreeFlame(gas;X="Ar:1",discretization=:unsupported)
    end
end
end
