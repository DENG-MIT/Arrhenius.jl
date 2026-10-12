module IonConservativeTests
using Arrhenius,Test,LinearAlgebra
using ..IonTransportTests: fixture
const A=Arrhenius
@testset "conservative ionized balances" begin
 mktempdir() do dir
  gas,_,_=fixture(dir)
  f=FreeFlame(gas;T=900.,X="Ar:.7,He:.3",grid=[0,.001,.003,.007,.01],discretization=:conservative)
  @test f.discretization==:conservative
  @test !f.initialized
  @test_throws ArgumentError FreeFlame(gas;X="Ar:1",flux_gradient_basis=:mass)
  @test_throws ArgumentError BurnerFlame(gas;X="Ar:1",mdot=.1,multicomponent_data=:invalid)
  u=f.state;n=gas.n_species;B,N=size(u)
  for j in 1:N
   u[1,j]=.85+.025*j;u[2,j]=.7+.001*j
   u[4,j]=(-1)^j*1e-11;u[5,j]=(-1)^j*1e-15
   u[3,j]=1-u[2,j]-u[4,j]-u[5,j]
   u[end-1,j]=.1*j;u[end,j]=.4+.01*j
  end
  w=A.FlameWorkspace(f);c=w.conservative;r=similar(u)
  for enabled in (false,true)
   set_electric_field!(f,enabled);old=copy(u);flame_residual!(r,f,u,w)
   @test u==old
   @test all(isfinite,r)
   cells=[(f.grid[min(j+1,N)]-f.grid[j-1])/2 for j in 2:N]
   integrated=zeros(n);energy=0.
   for (i,j) in enumerate(2:N)
    integrated .+= r[2:n+1,j].*w.rho[j]*cells[i]/A._flame_timescale
    energy+=r[1,j]*1000*c.cp[j]*w.rho[j]*cells[i]/A._flame_timescale
   end
   @test integrated≈c.species_flux[:,1]-c.species_flux[:,N] rtol=1e-12 atol=1e-18
   @test energy≈c.enthalpy_flux[1]-c.enthalpy_flux[N] rtol=1e-12 atol=1e-8
   @test A._flame_node_mass_flux(f,u,w,N)==w.rho[N]*u[end,N]
   @test c.species_flux[:,N]≈w.rho[N]*u[end,N]*u[2:n+1,N]+c.outlet_flux rtol=1e-14
   if enabled
    expected=-w.rho[N]*u[5,N]*1000*u[end-1,N]*.4
    @test c.outlet_flux[4]≈expected rtol=1e-14
    @test abs(sum(c.outlet_flux))<=1e-14*sum(abs,c.outlet_flux)
   else
    @test all(iszero,c.outlet_flux)
   end
   previous=copy(u);previous[1,:].-=1e-3;previous[2,:].-=1e-4;previous[3,:].+=1e-4
   hp=A._flame_previous_enthalpy(f,w,previous);rt=similar(u)
   flame_residual!(rt,f,u,w;previous,dt=2e-5,previous_enthalpy=hp)
   @test rt[end-1:end,:]≈r[end-1:end,:] rtol=0 atol=1e-12
   @test rt[:,1]≈r[:,1] rtol=0 atol=1e-12
   @test rt[2:n+1,2:end]-r[2:n+1,2:end]≈-5*(u[2:n+1,2:end]-previous[2:n+1,2:end]) rtol=1e-12 atol=1e-15
   @test rt[1,2:end]-r[1,2:end]≈-5*(c.enthalpy[2:end]-hp[2:end])./(1000*c.cp[2:end]) rtol=1e-12 atol=1e-15
  end
  flame_residual!(r,f,u,w);band=copy(A._flame_jacobian(f,u,w,r));bw=2B-1;J=zeros(B*N,B*N)
  for col in 1:B*N,row in max(1,col-bw):min(B*N,col+bw)
   J[row,col]=band[2bw+1+row-col,col]
  end
  d=copy(u);d[1,:].=.1;d[2,:].=.02;d[3,:].=-.01;d[4,:].=1e-9;d[5,:].=1e-13;d[end-1,:].=.2;d[end,:].=.05
  rp=similar(u);rm=similar(u);h=1e-4
  flame_residual!(rp,f,u+h*d,w;update_transport=false);flame_residual!(rm,f,u-h*d,w;update_transport=false)
  fd=(rp-rm)/(2h);jv=reshape(J*vec(d),B,N);gross=reshape(abs.(J)*abs.(vec(d)),B,N)
  for k in 1:B
   @test norm(jv[k,:]-fd[k,:],Inf)<=1e-4*max(norm(gross[k,:],Inf),1e-20)
  end
 end
end
@testset "conservative ionized initialization" begin
 mktempdir() do dir
  gas,_,_=fixture(dir)
  f=BurnerFlame(gas;T=900.,X="Ar:.7,He:.3",mdot=.1,
      grid=[0,.001,.003,.007,.01],discretization=:conservative)
  set_electric_field!(f,true)
  @test_throws ArgumentError solve!(f;ratio=1.)
  @test f.discretization==:conservative && f.field_enabled
  @test !f.initialized && !f.converged
  set_electric_field!(f,false)
  @test solve!(f;refine_grid=false)===f
  @test f.initialized && f.converged && !f.field_enabled
  @test f.discretization==:conservative
  @test temperature(f)≈fill(900.,5) rtol=1e-10
  @test maximum(abs,vec(sum(mass_fractions(f);dims=1)).-1)<1e-10
  @test all(iszero,electric_field(f))
 end
end
end
