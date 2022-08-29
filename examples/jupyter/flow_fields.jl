using MeshArrays, OceanStateEstimation, MITgcmTools
import IndividualDisplacements.NetCDF as NetCDF

"""
    read_velocities(γ::gcmgrid,t::Int,pth::String)

Read velocity components `u,v` from files in `pth`for time `t`
"""
function read_velocities(γ::gcmgrid,t::Int,pth::String)
    u=read_nctiles("$pth"*"UVELMASS/UVELMASS","UVELMASS",γ,I=(:,:,:,t))
    v=read_nctiles("$pth"*"VVELMASS/VVELMASS","VVELMASS",γ,I=(:,:,:,t))
    return u,v
end

"""
    solid_body_rotation(np,nz)

Set up an idealized flow field which consists of 
[rigid body rotation](https://en.wikipedia.org/wiki/Rigid_body), 
plus a convergent term, plus a sinking term.

```
u,v,w=solid_body_rotation(12,4)
```
"""
function solid_body_rotation(np,nz)
    Γ=simple_periodic_domain(np)
    Γ = UnitGrid(Γ.XC.grid;option="full")
    γ=Γ.XC.grid;
    
    #Solid-body rotation around central location ...
    i=Int(np/2+1)
    u=-(Γ.YG.-Γ.YG[1][i,i])
    v=(Γ.XG.-Γ.XG[1][i,i])
    
    #... plus a convergent term to / from central location
    d=-0.01
    u=u+d*(Γ.XG.-Γ.XG[1][i,i])
    v=v+d*(Γ.YG.-Γ.YG[1][i,i])
    
    #Replicate u,v in vertical dimension
    uu=MeshArray(γ,γ.ioPrec,nz)
    [uu[k]=u[1] for k=1:nz]
    vv=MeshArray(γ,γ.ioPrec,nz)
    [vv[k]=v[1] for k=1:nz]
    
    #Vertical velocity component w    
    w=fill(-0.01,MeshArray(γ,γ.ioPrec,nz+1))
    
    return write(uu),write(vv),write(w)
end

"""
    global_ocean_circulation(;k=10,ny=2)

Set up Global Ocean particle simulation in 2D with seasonally varying flow field.

```
𝑃,𝐷=global_ocean_circulation(k=10,ny=2);
```
"""
function global_ocean_circulation(;k=1,ny=2)

  #k=10 #choice of vertical level
  #ny=2 #number of simulated years (20 for k>20)
  r_reset = 0.01 #fraction of the particles reset per month (0.05 for k<=10)

  #read grid and set up connections between subdomains
  p=dirname(pathof(IndividualDisplacements))
  γ=GridSpec("LatLonCap",MeshArrays.GRID_LLC90)
  Γ=GridLoad(γ;option="full")
  Γ=merge(Γ,MeshArrays.NeighborTileIndices_cs(Γ))
  func=(u -> MeshArrays.update_location_llc!(u,Γ))

  #initialize u0,u1 etc
  𝑃,𝐷=set_up_FlowFields(k,Γ,func,ECCOclim_path);

  #add parameters for use in reset!
  tmp=(frac=r_reset, Γ=Γ)
  𝐷=merge(𝐷,tmp)

  return 𝑃,𝐷

end

"""
    OCCA_FlowFields(;backward_in_time::Bool=false,nmax=Inf)

Define gridded variables and return result as NamedTuple
"""
function OCCA_FlowFields(;backward_in_time::Bool=false,nmax=Inf)

   γ=GridSpec("PeriodicChannel",MeshArrays.GRID_LL360)
   Γ=GridLoad(γ;option="full")
   n=length(Γ.RC)
   isfinite(nmax) ? n=min(n,Int(nmax)) : nothing

   g=Γ.XC.grid
   func=(u -> MeshArrays.update_location_dpdo!(u,g))

   jj=[:hFacC, :hFacW, :hFacS, :DXG, :DYG, :RAC, :RAZ, :RAS]
   ii=findall([!in(i,jj) for i in keys(Γ)])
   Γ=(; zip(Symbol.(keys(Γ)[ii]), values(Γ)[ii])...)

   backward_in_time ? s=-1.0 : s=1.0
   s=Float32(s)

   function rd(filename, varname,n)
   fil = NetCDF.open(filename, varname)
   siz = size(fil)
   tmp = zeros(siz[1:2]...,n)
   [tmp .+= fil[:,:,1:n,t] for t=1:12]
   tmp ./= 12.0
   tmp[findall(tmp.<-1e22)] .= 0.0
   return tmp
   end

   fileIn=OCCAclim_path*"DDuvel.0406clim.nc"
   u=s*read(rd(fileIn,"u",n),MeshArray(γ,Float32,n))

   fileIn=OCCAclim_path*"DDvvel.0406clim.nc"
   v=s*read(rd(fileIn,"v",n),MeshArray(γ,Float32,n))

   fileIn=OCCAclim_path*"DDwvel.0406clim.nc"
   w=s*rd(fileIn,"w",n)
   w=-cat(w,zeros(360, 160),dims=3)
   w[:,:,1] .=0.0
   w=read(w,MeshArray(γ,Float32,n+1))

   fileIn=OCCAclim_path*"DDtheta.0406clim.nc"
   θ=read(rd(fileIn,"theta",n),MeshArray(γ,Float32,n))

#   fileIn=OCCAclim_path*"DDsalt.0406clim.nc"
#   𝑆=read(rd(fileIn,"salt",n),MeshArray(γ,Float64,n))

   for i in eachindex(u)
      u[i]=u[i]./Γ.DXC[1]
      v[i]=v[i]./Γ.DYC[1]
   end

   for i in eachindex(u)
      u[i]=circshift(u[i],[-180 0])
      v[i]=circshift(v[i],[-180 0])
      θ[i]=circshift(θ[i],[-180 0])
#      𝑆[i]=circshift(𝑆[i],[-180 0])
   end

   for i in eachindex(w)
      w[i]=w[i]./Γ.DRC[min(i[2]+1,n)]
      w[i]=circshift(w[i],[-180 0])
   end

   tmpx=circshift(Γ.XC[1],[-180 0])
   tmpx[1:180,:]=tmpx[1:180,:] .- 360.0
   Γ.XC[1]=tmpx

   tmpx=circshift(Γ.XG[1],[-180 0])
   tmpx[1:180,:]=tmpx[1:180,:] .- 360.0
   Γ.XG[1]=tmpx
   Γ.Depth[1]=circshift(Γ.Depth[1],[-180 0])

   t0=0.0; t1=86400*366*2.0;

   for k=1:n
    (tmpu,tmpv)=exchange(u[:,k],v[:,k],1)
    u[:,k]=tmpu
    v[:,k]=tmpv
   end
   for k=1:n+1
    tmpw=exchange(w[:,k],1)
    w[:,k]=tmpw
   end

   𝑃=FlowFields(u,u,v,v,w,w,[t0,t1],func)

   𝐷 = (θ0=θ, θ1=θ, XC=exchange(Γ.XC), YC=exchange(Γ.YC), 
   RF=Γ.RF, RC=Γ.RC,ioSize=(360,160,n))

   return 𝑃,𝐷,Γ

end

"""
    test1_setup()

Call `gcmgrid`, initialize a single point,
rely on `dxdt!`, and just output `sol` at the end.

```
using IndividualDisplacements, MeshArrays
import IndividualDisplacements.OrdinaryDiffEq: ODEProblem, solve, Tsit5
𝑃,sol=test1_setup()
```
"""
function test1_setup()

    mygrid=gcmgrid("flt_example/","ll",1,[(80,42)], [80 42], Float32, read, write)
    XC=MeshArray(mygrid,Float32); XC[1]=vec(2500.:5000.:397500.0)*ones(1,42);
    XG=MeshArray(mygrid,Float32); XG[1]=vec(0.:5000.:395000.0)*ones(1,42);
    YC=MeshArray(mygrid,Float32); YC[1]=ones(80,1)*transpose(vec(2500.:5000.:207500.0));
    YG=MeshArray(mygrid,Float32); YG[1]=ones(80,1)*transpose(vec(0.:5000.:205000.0));

    dx=5000.0
    t0=0.0; t1=18001.0*3600.0
    u=-(YG.-YC[1][40,21])/2000000.
    v=(XG.-XC[1][40,21])/2000000.
    u0=u./dx; u1=u./dx
    v0=v./dx; v1=v./dx

    𝑃=FlowFields(u0[1], u1[1], v0[1], v1[1], [t0,t1])
    
    u0=[200000.0;0.0]./dx
    du=fill(0.0,2);
    prob = ODEProblem(dxdt!,u0,[0.0,2998*3600.0],𝑃)
    sol = solve(prob,Tsit5(),reltol=1e-8,abstol=1e-8)

    return 𝑃,sol
end

"""
    test2_periodic_domain(np = 12, nq = 12)

Call `simple_periodic_domain`, initialize 6x6 point cloud,
rely on `dxdt!`, and call `postprocess_xy` at the end.

```
using IndividualDisplacements, MeshArrays
import IndividualDisplacements.OrdinaryDiffEq: ODEProblem, solve, Euler

df,𝑃=test2_periodic_domain()

using Plots
@gif for t in 𝑃.t0:1.0:𝑃.t1
   scatter_subset(𝑃,df,t)
end
```
"""
function test2_periodic_domain(np = 12, nq = 12)
    #domain and time parameters
    Γ = simple_periodic_domain(np, nq)
    Γ = UnitGrid(Γ.XC.grid;option="full")

    u = 0.1 ./ Γ.DXC
    v = 0.3 ./ Γ.DYC
    (u, v) = exchange(u, v, 1)

    f = (u -> MeshArrays.update_location_dpdo!(u,Γ.XC.grid))
    𝑃=FlowFields(u,u,v,v,[0.0,400.0],f)
    𝐷=(;)

    #initial conditions
    x0 = np * (0.4:0.04:0.6)
    y0 = nq * (0.4:0.04:0.6)
    x0 = vec(x0) * ones(1, length(y0))
    y0 = ones(size(x0, 1), 1) * transpose(vec(y0))
    u0 = permutedims([[x0[i];y0[i];1.0] for i in eachindex(x0)])
    du=0*u0
    
    #solve for trajectories
    prob = ODEProblem(dxdt!, u0, 𝑃.𝑇, 𝑃)
    sol = solve(prob,Euler(),dt=0.1)

    return postprocess_xy(sol, 𝑃, 𝐷),𝑃
end
