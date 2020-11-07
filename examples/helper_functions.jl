
"""
    setup_random_flow(Γ::Dict)

Set up a random flow field over the domain specified by Γ

```
Γ=simple_periodic_domain(12)
𝑃,ϕ=setup_random_flow(Γ)
```
"""
function setup_random_flow(Γ::Dict)
  (_,ϕ,_,_)=demo2(Γ);

  (u,v)=gradient(ϕ,Γ)
  u=u./Γ["DXC"]#normalization to grid units
  v=v./Γ["DYC"]

  (u,v)=exchange(u,v,1)
  u0=-v; u1=-v;
  v0=u; v1=u;

  𝑃 = (u0=u0, u1=u1, v0=v0, v1=v1, 𝑇=[0.0,400.0], ioSize=ϕ.grid.ioSize)
  return 𝑃,ϕ

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

  #initialize u0,u1 etc
  𝑃=set_up_𝑃(k,0.0,Γ,ECCOclim_path);

  #add parameters for use in reset!
  tmp=(frac=r_reset, Γ=Γ)
  𝑃=merge(𝑃,tmp)

  return 𝑃

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

   delete!.(Ref(Γ), ["hFacC", "hFacW", "hFacS","DXG","DYG","RAC","RAZ","RAS"]);
   backward_in_time ? s=-1.0 : s=1.0

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
   θ=s*read(rd(fileIn,"theta",n),MeshArray(γ,Float32,n))

#   fileIn=OCCAclim_path*"DDsalt.0406clim.nc"
#   𝑆=s*read(rd(fileIn,"salt",n),MeshArray(γ,Float64,n))

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

#   for k=1:n
#    (tmpu,tmpv)=exchange(u[:,k],v[:,k],1)
#    u[:,k]=tmpu
#    v[:,k]=tmpv
#   end
#   for k=1:n+1
#    tmpw=exchange(w[:,k],1)
#    w[:,k]=tmpw
#   end

   𝑃 = (θ0=θ, θ1=θ, u0=u, u1=u, v0=v, v1=v, w0=w, w1=w, 𝑇=[t0,t1],
   XC=exchange(Γ["XC"]), YC=exchange(Γ["YC"]), 
   RF=Γ["RF"], RC=Γ["RC"],
   ioSize=(360,160,n))

   return 𝑃,Γ

end

##

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

function init_global_randn(np ::Int , 𝑃::NamedTuple)
    (lon, lat) = randn_lonlat(2*np)
    (u0, _) = initialize_lonlat(𝑃.Γ, lon, lat; msk = 𝑃.msk)
    u0[:,1:np]
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
