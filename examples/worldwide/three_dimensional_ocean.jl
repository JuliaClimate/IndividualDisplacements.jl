### A Pluto.jl notebook ###
# v0.17.3

using Markdown
using InteractiveUtils

# ╔═╡ 16eab80b-325b-43bd-8bda-6b9ed27513a8
begin
	using Pkg
	Pkg.activate()
	
	using IndividualDisplacements, OceanStateEstimation, DataFrames, MeshArrays, NetCDF
	import CairoMakie as Mkie
end

# ╔═╡ 68c92218-40d3-11ec-0397-1747ac61c311
md"""# Three Dimensional Ocean Circulation

Advect particles with climatological mean flow in three dimensions starting from a selected depth level (e.g. `k=10` for 95 m) and region using a near-global ocean state estimate ([OCCA](https://doi.org/10.1175/2009JPO4043.1) which is here repeated for two years. For additional documentation e.g. see : [1](https://JuliaClimate.github.io/MeshArrays.jl/dev/), [2](https://JuliaClimate.github.io/IndividualDisplacements.jl/dev/), [3](https://docs.juliadiffeq.org/latest/solvers/ode_solve.html), [4](https://en.wikipedia.org/wiki/Displacement_(vector))

"""
#![Three dimensional simulation](https://user-images.githubusercontent.com/20276764/94491485-595ee900-01b6-11eb-95e6-c2cacb812f46.png)

# ╔═╡ 44346351-f249-4376-b002-8147755ed489
md"""## Initialization"""

# ╔═╡ 00464caa-fab2-4fd2-b39b-177e505d6d89
md"""## Compute Displacements"""

# ╔═╡ b66f8e31-6aae-429b-9c42-bb9ea2d01eb3
md"""## Helper Functions"""

# ╔═╡ a879b36d-8536-4b9e-a22d-b3d2161e589c
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
   func=(u -> IndividualDisplacements.update_location_dpdo!(u,g))

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

# ╔═╡ 66c95828-227c-4db5-a6f1-3e3004a99785
begin
	OceanStateEstimation.get_occa_velocity_if_needed();
	𝑃,𝐷,Γ=OCCA_FlowFields(nmax=5)
	"done with flow fields"
end

# ╔═╡ 378b6547-898d-477a-a796-285a1f7e9b08
"""
    initial_positions(Γ; nf=10000, lon_rng=(-160.0,-159.0), lat_rng=(30.0,31.0))

Randomly assign initial positions in longitude,latitude ranges. Positions are
expressed in, normalized, grid point units (x,y in the 0,nx and 0,ny range).
To convert from longitude,latitude here we take advantage of the regularity
of the 1 degree grid being used -- for a more general alternative, see the
global ocean example.
"""
function initial_positions(Γ::NamedTuple, nf=10000, lon_rng=(-160.0,-159.0), lat_rng=(30.0,31.0))
   lon=lon_rng[1] .+(lon_rng[2]-lon_rng[1]).*rand(nf)
   lat=lat_rng[1] .+(lat_rng[2]-lat_rng[1]).*rand(nf)
   x=lon .+ (21. - Γ.XC[1][21,1])
   y=lat .+ (111. - Γ.YC[1][1,111])
   return x,y
end

# ╔═╡ b754f4f6-b513-4eda-b689-8e0529223417
custom🔴 = DataFrame(ID=Int[], fid=Int[], x=Float64[], y=Float64[],
   k=Float64[], z=Float64[], iso=Float64[], t=Float64[],
   lon=Float64[], lat=Float64[], year=Float64[], col=Symbol[])

# ╔═╡ c5eeb1f0-2c6f-4cd5-8a66-538da09c282d
function custom🔧(sol,𝑃::𝐹_MeshArray3D;id=missing,𝑇=missing)
   df=postprocess_MeshArray(sol,𝑃,id=id,𝑇=𝑇)
   add_lonlat!(df,𝐷.XC,𝐷.YC)

   #add year (convenience time axis for plotting)
   df.year=df.t ./86400/365

   #add depth (i.e. the 3rd, vertical, coordinate)
   k=[sol[1,i,j][3] for i in 1:size(sol,2), j in 1:size(sol,3)]
   nz=length(𝐷.RC)
   df.k=min.(max.(k[:],Ref(0.0)),Ref(nz)) #level
   k=Int.(floor.(df.k)); w=(df.k-k);
   df.z=𝐷.RF[1 .+ k].*(1 .- w)+𝐷.RF[2 .+ k].*w #depth

   #add one isotherm depth
   θ=0.5*(𝐷.θ0+𝐷.θ1)
   d=MeshArrays.isosurface(θ,15,𝐷)
   d[findall(isnan.(d))].=0.
   df.iso=interp_to_xy(df,exchange(d));

   #add color = f(iso-z)
   c=fill(:gold,length(df.iso))
   c[findall(df.iso.<df.z)].=:violet
   df.col=c

   #to plot e.g. Pacific Ocean transports, shift longitude convention?
   df.lon[findall(df.lon .< 0.0 )] = df.lon[findall(df.lon .< 0.0 )] .+360.0
   return df
end

# ╔═╡ f199f321-976a-4ccd-a003-140211aa67fe
begin
	nf=100; lo=(-160.0,-150.0); la=(30.0,40.0); kk=2.5;
	df=DataFrame(:z => fill(kk,nf),:f => fill(1,nf))
	(df.x,df.y)=initial_positions(Γ, nf, lo, la)
	
	𝐼=Individuals(𝑃,df.x,df.y,df.z,df.f,(🔴=custom🔴,🔧=custom🔧))
end

# ╔═╡ 938fdaa8-357d-477e-8fa2-e6da53806242
begin
	𝑇=(0.0,10*86400.0)
	∫!(𝐼,𝑇)
	🔴_by_t = groupby(𝐼.🔴, :t)
end

# ╔═╡ dafc7de0-d2c4-42cd-8cd9-cc152cadb33e
begin
	Mkie.set_theme!(Mkie.theme_light())
	fig=Mkie.Figure(resolution = (900, 600))
	a = Mkie.Axis(fig[1, 1],xlabel="longitude",ylabel="latitude")		
	Mkie.scatter!(a,🔴_by_t[1].lon,🔴_by_t[1].lat,color=:green2)
	Mkie.scatter!(a,🔴_by_t[end].lon,🔴_by_t[end].lat,color=:red)
	fig
end

# ╔═╡ 7f3b1f13-abfa-4fc0-8169-b69f1c3519f4
"""
    plot(𝐼::Individuals)

Plot the initial and final positions as scatter plot in x,y plane.
"""
function plot(𝐼::Individuals)
    🔴_by_t = groupby(𝐼.🔴, :t)
    if (sum(names(🔴_by_t).=="lon")==0)
        fig=scatter(🔴_by_t[1].x,🔴_by_t[1].y,c=:red,label="t0",marker = (:circle, stroke(0)))
        scatter!(🔴_by_t[end].x,🔴_by_t[end].y,c=:blue,label="t1",marker = (:circle, stroke(0)))
    else
        fig=scatter(🔴_by_t[1].lon,🔴_by_t[1].lat,c=:red,label="t0",marker = (:circle, stroke(0)))
        scatter!(🔴_by_t[end].lon,🔴_by_t[end].lat,c=:blue,label="t1",marker = (:circle, stroke(0)))
    end
    return fig
end

# ╔═╡ Cell order:
# ╟─68c92218-40d3-11ec-0397-1747ac61c311
# ╟─16eab80b-325b-43bd-8bda-6b9ed27513a8
# ╟─44346351-f249-4376-b002-8147755ed489
# ╟─66c95828-227c-4db5-a6f1-3e3004a99785
# ╟─f199f321-976a-4ccd-a003-140211aa67fe
# ╟─00464caa-fab2-4fd2-b39b-177e505d6d89
# ╟─938fdaa8-357d-477e-8fa2-e6da53806242
# ╟─dafc7de0-d2c4-42cd-8cd9-cc152cadb33e
# ╟─b66f8e31-6aae-429b-9c42-bb9ea2d01eb3
# ╟─a879b36d-8536-4b9e-a22d-b3d2161e589c
# ╟─378b6547-898d-477a-a796-285a1f7e9b08
# ╟─b754f4f6-b513-4eda-b689-8e0529223417
# ╟─c5eeb1f0-2c6f-4cd5-8a66-538da09c282d
# ╟─7f3b1f13-abfa-4fc0-8169-b69f1c3519f4
