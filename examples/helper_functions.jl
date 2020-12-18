
using IndividualDisplacements, MeshArrays, OrdinaryDiffEq
p=dirname(pathof(MeshArrays))
include(joinpath(p,"../examples/Demos.jl"))

"""
    simple_flow_field(Γ::Dict,np,nz)

Set up an idealized flow field which consists of 
[rigid body rotation](https://en.wikipedia.org/wiki/Rigid_body), 
plus a convergent term, plus a sinking term.

```
Γ=simple_periodic_domain(12)
u,v,w=simple_flow_field(Γ,12,4)
```
"""
function simple_flow_field(Γ,np,nz)
    γ=Γ["XC"].grid;
    
    #Solid-body rotation around central location ...
    i=Int(np/2+1)
    u=-(Γ["YG"].-Γ["YG"][1][i,i])
    v=(Γ["XG"].-Γ["XG"][1][i,i])
    
    #... plus a convergent term to / from central location
    d=-0.01
    u=u+d*(Γ["XG"].-Γ["XG"][1][i,i])
    v=v+d*(Γ["YG"].-Γ["YG"][1][i,i])
    
    #Replicate u,v in vertical dimension
    uu=MeshArray(γ,γ.ioPrec,nz)
    [uu[k]=u[1] for k=1:nz]
    vv=MeshArray(γ,γ.ioPrec,nz)
    [vv[k]=v[1] for k=1:nz]
    
    #Vertical velocity component w    
    w=fill(-0.01,MeshArray(γ,γ.ioPrec,nz));
    
    return write(uu),write(vv),write(w)
end
    
"""
    setup_random_flow(;np=12,nq=18)

Set up a random flow field over a gridded domain of size np,nq

```
ϕ,u,v=setup_random_flow()
```
"""
function setup_random_flow(;np=12,nq=18)
    Γ=simple_periodic_domain(np,nq)
    (_,ϕ,_,_)=demo2(Γ)
    ϕ .*= 0.5

    #For the convergent / scalar potential case, ϕ is interpreted as being 
    #on center points -- hence the standard gradient function readily gives 
    #what we need
    #(u,v)=gradient(ϕ,Γ) 
    #return u[1],v[1],ϕ[1]

    #For the rotational / streamfunction case, ϕ is interpreted as being 
    #on S/W corner points -- this is ok since the grid is homegeneous, 
    #and conveniently yields an adequate substitution u,v <- -v,u; but note
    #that doing the same with gradient() would shift indices inconsistenly
    u=-(circshift(ϕ[1], (0,-1))-ϕ[1])
    v=(circshift(ϕ[1], (-1,0))-ϕ[1])
    return u,v,ϕ[1]
end

function setup_point_cloud(U::Array{T,2},V::Array{T,2};X=[],Y=[]) where T
    np,nq=size(U)
    Γ=simple_periodic_domain(np,nq)
    g=Γ["XC"].grid
    u=MeshArray(g,[U])
    v=MeshArray(g,[V])
    #vel=dxy_dt
    (u,v)=exchange(u,v,1)
    vel=dxy_dt!
    func=(u -> IndividualDisplacements.update_location_dpdo!(u,g))

    𝑃=𝑃_MeshArray2D{eltype(u)}(u,u,v,v,[0.0,10.0],func)
    pp=postprocess_xy
    isempty(X) ? X=np*rand(10) : nothing
    isempty(Y) ? Y=nq*rand(10) : nothing

    xy = permutedims([[X[i];Y[i];1.0] for i in eachindex(X)])
    tr = DataFrame(ID=Int[], x=Float64[], y=Float64[], t=Float64[])
    solv(prob) = solve(prob,Tsit5(),reltol=1e-5,abstol=1e-5)
    
    I=(position=xy,record=tr,velocity=vel,
       integration=solv,postprocessing=pp,parameters=𝑃)

    return Individuals(I)
end

"""
    setup_global_ocean(;k=10,ny=2)

Set up Global Ocean particle simulation in 2D with seasonally varying flow field.

```
𝑃=setup_global_ocean(k=1,ny=2);
```
"""
function setup_global_ocean(;k=1,ny=2)

  #k=10 #choice of vertical level
  #ny=2 #number of simulated years (20 for k>20)
  r_reset = 0.01 #fraction of the particles reset per month (0.05 for k<=10)

  #read grid and set up connections between subdomains
  p=dirname(pathof(IndividualDisplacements))
  γ=GridSpec("LatLonCap",MeshArrays.GRID_LLC90)
  Γ=GridLoad(γ)
  Γ=merge(Γ,IndividualDisplacements.NeighborTileIndices_cs(Γ))

  func=(u -> IndividualDisplacements.update_location_llc!(u,𝐷))
  Γ=merge(Γ,Dict("update_location!" => func))

  #initialize u0,u1 etc
  𝑃,𝐷=set_up_𝑃(k,0.0,Γ,ECCOclim_path);

  #add parameters for use in reset!
  tmp=(frac=r_reset, Γ=Γ)
  𝐷=merge(𝐷,tmp)

  return 𝑃,𝐷

end

##

"""
    example3_setup(;backward_in_time::Bool=false)

Define gridded variables and return result as Dictionary (`uvetc`).
"""
function example3_setup(;backward_in_time::Bool=false)

   p=dirname(pathof(IndividualDisplacements))
   dirIn=joinpath(p,"../examples/llc90_latlon/")
   γ=gcmgrid(dirIn,"PeriodicChannel",1,
                  [(360,178)], [360 178], Float32, read, write)

   Γ=Dict("XC" => read(γ.path*"XC.latlon.data",MeshArray(γ,Float32)),
   "YC" => read(γ.path*"YC.latlon.data",MeshArray(γ,Float32)),
   "XG" => read(γ.path*"XG.data",MeshArray(γ,Float32)),
   "YG" => read(γ.path*"YG.data",MeshArray(γ,Float32)),
   "DXC" => read(γ.path*"DXC.latlon.data",MeshArray(γ,Float32)),
   "DYC" => read(γ.path*"DYC.latlon.data",MeshArray(γ,Float32)) );

   file = matopen(γ.path*"uv_lonlat.mat")
   u=read(file, "u")
   v=read(file, "v")
   close(file)

   u=dropdims(mean(u,dims=3),dims=3)
   v=dropdims(mean(v,dims=3),dims=3)

   #mask out near edge values to avoid exiting domain
   u[:,1:2].=NaN
   v[:,1:2].=NaN
   u[:,end-2:end].=NaN
   v[:,end-2:end].=NaN

   u=read(u,MeshArray(γ,Float32))
   v=read(v,MeshArray(γ,Float32));

   u[findall(isnan.(u))]=0.0
   v[findall(isnan.(v))]=0.0

   backward_in_time ? s=-1.0 : s=1.0
   u0=s*u; u1=s*u;
   v0=s*v; v1=s*v;

   t0=0.0; t1=86400*366*2.0;

   u0=u0./Γ["DXC"]
   u1=u1./Γ["DXC"]
   v0=v0./Γ["DYC"]
   v1=v1./Γ["DYC"]

   #nr=50; kk=1;

   #mskW=read(γ.path*"hFacW.latlon.data",MeshArray(γ,Float32,nr))
   #mskW=1.0 .+ 0.0 * mask(mskW[:,kk],NaN,0.0)
   #mskS=read(γ.path*"hFacS.latlon.data",MeshArray(γ,Float32,nr))
   #mskS=1.0 .+ 0.0 * mask(mskS[:,kk],NaN,0.0)
   #msk=Dict("mskW" => mskW, "mskS" => mskS)

   𝑃 = (u0=u0, u1=u1, v0=v0, v1=v1, 𝑇=[t0,t1], ioSize=(360,178),
        XC=exchange(Γ["XC"]), YC=exchange(Γ["YC"]))

   return 𝑃,Γ

end

"""
    OCCA_setup(;backward_in_time::Bool=false)

Define gridded variables and return result as Dictionary (`uvetc`).
"""
function OCCA_setup(;backward_in_time::Bool=false)

   γ=GridSpec("PeriodicChannel",MeshArrays.GRID_LL360)
   Γ=GridLoad(γ)
   n=length(Γ["RC"])
   n=10

   g=Γ["XC"].grid
   func=(u -> IndividualDisplacements.update_location_dpdo!(u,g))

   delete!.(Ref(Γ), ["hFacC", "hFacW", "hFacS","DXG","DYG","RAC","RAZ","RAS"]);
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
      u[i]=u[i]./Γ["DXC"][1]
      v[i]=v[i]./Γ["DYC"][1]
   end

   for i in eachindex(u)
      u[i]=circshift(u[i],[-180 0])
      v[i]=circshift(v[i],[-180 0])
      θ[i]=circshift(θ[i],[-180 0])
#      𝑆[i]=circshift(𝑆[i],[-180 0])
   end

   for i in eachindex(w)
      w[i]=w[i]./Γ["DRC"][min(i[2]+1,n)]
      w[i]=circshift(w[i],[-180 0])
   end

   tmpx=circshift(Γ["XC"][1],[-180 0])
   tmpx[1:180,:]=tmpx[1:180,:] .- 360.0
   Γ["XC"][1]=tmpx

   tmpx=circshift(Γ["XG"][1],[-180 0])
   tmpx[1:180,:]=tmpx[1:180,:] .- 360.0
   Γ["XG"][1]=tmpx
   Γ["Depth"][1]=circshift(Γ["Depth"][1],[-180 0])

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

   𝑃=𝑃_MeshArray3D{eltype(u)}(u,u,v,v,w,w,[t0,t1],func)

   𝐷 = (θ0=θ, θ1=θ, XC=exchange(Γ["XC"]), YC=exchange(Γ["YC"]), 
   RF=Γ["RF"], RC=Γ["RC"],ioSize=(360,160,n))

   return 𝑃,𝐷,Γ

end

"""
    init_global_range(lons::Tuple = (-160.0, -150.0),lats::Tuple = (35.0, 45.0))

Randomly distribute `np` points over a lon,la region, and 
return position in grid index space (`i,j,subdomain`).
"""
function init_global_range(lons::Tuple = (-160.0, -150.0),lats::Tuple = (35.0, 45.0))
    lo0, lo1 = lons #(-160.0, -150.0)
    la0, la1 = lats #(35.0, 45.0)
    np = 100
    lon = lo0 .+ (lo1 - lo0) .* rand(np)
    lat = la0 .+ (la1 - la0) .* rand(np)
    (u0, _) = initialize_lonlat(Γ, lon, lat; msk = Γ["hFacC"][:, k])
    id=collect(1:np)
    return u0
end

"""
    init_global_randn(np ::Int , 𝑃::NamedTuple)

Randomly distribute `np` points over the Earth, within `𝑃.msk` 
region, and return position in grid index space (`i,j,subdomain`).
"""
function init_global_randn(np ::Int , 𝑃::NamedTuple)
    (lon, lat) = randn_lonlat(2*np)
    (u0, _) = initialize_lonlat(𝑃.Γ, lon, lat; msk = 𝑃.msk)
    u0[:,1:np]
end

"""
    reset_lonlat!(𝐼::Individuals)

Randomly select a fraction (𝐼.𝑃.frac) of the particles and reset their positions.
"""
function reset_lonlat!(𝐼::Individuals,𝐷::NamedTuple)
    np=length(𝐼.🆔)
    n_reset = Int(round(𝐷.frac*np))
    (lon, lat) = randn_lonlat(2*n_reset)
    (v0, _) = initialize_lonlat(𝐷.Γ, lon, lat; msk = 𝐷.msk)
    n_reset=min(n_reset,size(v0,2))
    k_reset = rand(1:np, n_reset)
    v0 = permutedims([v0[:,i] for i in 1:size(v0,2)])
    𝐼.📌[k_reset].=v0[1:n_reset]
    isempty(𝐼.🔴.ID) ? m=maximum(𝐼.🆔) : m=max(maximum(𝐼.🔴.ID),maximum(𝐼.🆔))
    𝐼.🆔[k_reset]=collect(1:n_reset) .+ m
end

##

"""
    isosurface(θ,T,z)

```
isosurface(𝐼.𝑃.θ0,15,Γ["RC"])
```    
"""
function isosurface(θ,T,z)
    d=NaN*similar(θ[:,1])
    nr=size(θ,2)
    for j=1:size(d,1)
        for k=1:nr-1
            i=findall(isnan.(d[j]).&(θ[j,k].>T).&(θ[j,k+1].<=T))
            a=(θ[j,k][i] .- T)./(θ[j,k][i] .- θ[j,k+1][i])
            d[j][i]=(1 .- a).*Γ["RC"][k] + a.*Γ["RC"][k+1]
            i=findall(isnan.(d[j]).&(θ[j,k].<=T).&(θ[j,k+1].>T))
            a=(θ[j,k+1][i] .- T)./(θ[j,k+1][i] .- θ[j,k][i])
            d[j][i]=(1 .- a).*Γ["RC"][k+1] + a.*Γ["RC"][k]
        end
    end
    return d
end
