module OCCA_FlowFields

using IndividualDisplacements, Climatology, NetCDF

import IndividualDisplacements.DataFrames: DataFrame
import IndividualDisplacements.MeshArrays as MeshArrays
import IndividualDisplacements.MeshArrays: MeshArray
import IndividualDisplacements.OrdinaryDiffEq: solve, Tsit5, remake
import IndividualDisplacements.OrdinaryDiffEq: ODEProblem, EnsembleProblem

function setup(;backward_in_time::Bool=false,nmax=Inf)

   γ=MeshArrays.GridSpec("PeriodicChannel",MeshArrays.GRID_LL360)
   Γ=MeshArrays.GridLoad(γ;option="full")
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

   Climatology.get_occa_velocity_if_needed()

   fileIn=joinpath(ScratchSpaces.OCCA,"DDuvel.0406clim.nc")
   u=s*read(rd(fileIn,"u",n),MeshArray(γ,Float32,n))

   fileIn=joinpath(ScratchSpaces.OCCA,"DDvvel.0406clim.nc")
   v=s*read(rd(fileIn,"v",n),MeshArray(γ,Float32,n))

   fileIn=joinpath(ScratchSpaces.OCCA,"DDwvel.0406clim.nc")
   w=s*rd(fileIn,"w",n)
   w=-cat(w,zeros(360, 160),dims=3)
   w[:,:,1] .=0.0
   w=read(w,MeshArray(γ,Float32,n+1))

   fileIn=joinpath(ScratchSpaces.OCCA,"DDtheta.0406clim.nc")
   θ=read(rd(fileIn,"theta",n),MeshArray(γ,Float32,n))

#   fileIn=joinpath(ScratchSpaces.OCCA,"DDsalt.0406clim.nc")
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
    (tmpu,tmpv)=MeshArrays.exchange(u[:,k],v[:,k],1)
    u[:,k]=tmpu
    v[:,k]=tmpv
   end
   for k=1:n+1
    tmpw=MeshArrays.exchange(w[:,k],1)
    w[:,k]=tmpw
   end

   P=FlowFields(u,u,v,v,w,w,[t0,t1],func)

   D = (θ0=θ, θ1=θ, XC=MeshArrays.exchange(Γ.XC), YC=MeshArrays.exchange(Γ.YC), 
   RF=Γ.RF, RC=Γ.RC,ioSize=(360,160,n), Γ=Γ)

   return P,D

end

"""
    initial_positions(Γ; nf=10000, lon_rng=(-160.0,-159.0), lat_rng=(30.0,31.0))

Randomly assign initial positions in longitude,latitude ranges. Positions are 
expressed in, normalized, grid point units (x,y in the 0,nx and 0,ny range). 
To convert from longitude,latitude here we take advantage of the regularity 
of the 1 degree grid being used -- for a more general alternative, see the 
global ocean example.
"""
function initial_positions(Γ::NamedTuple, nf=10000, lon_rng=(-160.0,-159.0), lat_rng=(30.0,31.0), level=1)
   lon=lon_rng[1] .+(lon_rng[2]-lon_rng[1]).*rand(nf)
   lat=lat_rng[1] .+(lat_rng[2]-lat_rng[1]).*rand(nf)
   x=lon .+ (21. - Γ.XC[1][21,1])
   y=lat .+ (111. - Γ.YC[1][1,111])

   return DataFrame(:x => x, :y => y, :z => fill(level,nf),:f => fill(1,nf))
end

custom🔴 = DataFrame(ID=Int[], fid=Int[], x=Float64[], y=Float64[],
   k=Float64[], z=Float64[], iso=Float64[], t=Float64[],
   lon=Float64[], lat=Float64[], dlon=Float64[], dlat=Float64[], 
   year=Float64[], col=Symbol[])

function custom🔧(sol,P::uvwMeshArrays,D::NamedTuple;id=missing,T=missing)
   df=postprocess_MeshArray(sol,P,D,id=id,T=T)
   add_lonlat!(df,D.Γ.XC,D.Γ.YC)
   df.dlon=0*df.lon
   df.dlat=0*df.lat

   #add year (convenience time axis for plotting)
   df.year=df.t ./86400/365

   #add depth (i.e. the 3rd, vertical, coordinate)
   k=[[sol[i][3,1] for i in 1:size(sol,3)];[sol[i][3,end] for i in 1:size(sol,3)]]
   nz=length(D.RC)
   df.k=min.(max.(k[:],Ref(0.0)),Ref(nz)) #level
   k=Int.(floor.(df.k)); w=(df.k-k);
   df.z=D.RF[1 .+ k].*(1 .- w)+D.RF[2 .+ k].*w #depth

   #add one isotherm depth
   θ=0.5*(D.θ0+D.θ1)
   d=MeshArrays.isosurface(θ,15,D)
   d[findall(isnan.(d))].=0.
   df.iso=interp_to_xy(df,MeshArrays.exchange(d));

   #add color = f(iso-z)
   c=fill(:gold,length(df.iso))
   c[findall(df.iso.<df.z)].=:violet
   df.col=c

   #to plot e.g. Pacific Ocean transports, shift longitude convention?
   df.lon[findall(df.lon .< 0.0 )] = df.lon[findall(df.lon .< 0.0 )] .+360.0
   return df
end

function custom∫(prob)
	u0 = prob.u0
	prob_func(prob,i,repeat) = remake(prob,u0=u0[i])
	indiv_prob = ODEProblem(prob.f,u0[1],prob.tspan,prob.p)
	ensemble_prob = EnsembleProblem(indiv_prob,prob_func=prob_func,safetycopy=false)
	solve(ensemble_prob, Tsit5(), reltol=1e-8, abstol=1e-8, 
         trajectories=length(u0),saveat=365/12*86400.0)
end

end #module OCCA_FlowFields
