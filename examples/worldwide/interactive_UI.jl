### A Pluto.jl notebook ###
# v0.17.3

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 747d446a-dfeb-11ea-3533-c9404fd41688
begin
	using Pkg 
	Pkg.activate()
	
	using IndividualDisplacements, MeshArrays, DataFrames, OceanStateEstimation
	using Statistics, PlutoUI, OrdinaryDiffEq, StatsPlots, NetCDF
	✓ = "😃"
	"$✓ Set up packages"
end

# ╔═╡ bf19d29c-e70e-11ea-0153-d3a49981d56c
md"""longitude,latitude of South West corner =
$(@bind lon0 NumberField(-180.0:10:180.0; default=0.0))
,
$(@bind lat0 NumberField(-90.0:5:90.0; default=-50.0))
"""

# ╔═╡ 4935fd46-e70f-11ea-386c-f9c444a20644
md"""depth level (integer), set size (nb of particles) =
$(@bind klev NumberField(1:10; default=5))
,
$(@bind npar NumberField(50:50:1000; default=100))
""" 

# ╔═╡ 9c80e722-e70f-11ea-22a6-0be2e85f3b8b
md"""time step (for plotting)= 
$(@bind tt Slider(1:26; default=13, show_value=true))
"""

# ╔═╡ e25eee9e-dfee-11ea-2a4c-3946ccb63876
begin
	lo0=lon0; lo1=lo0+5
	la0=lat0; la1=la0+5
	z_init=klev
	n_part=npar
	(particles=n_part,longitudes=(lo0,lo1),latitudes=(la0,la1),level=z_init,plot_time=tt)
end

# ╔═╡ 87c5e5d4-c343-49fb-bc7e-6b2c0e647a38
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

# ╔═╡ 1811955d-618d-49e3-9f81-63195e2632fe
begin
	𝑃,𝐷,Γ=OCCA_FlowFields(nmax=5)
	tmp=(Γ = Γ, m = "OCCA")
    𝐷=merge(𝐷,tmp)
	"$✓ Set up gridded domain"
end

# ╔═╡ f75fae30-dfee-11ea-18ef-259321acfa2f
begin	
	lon=lo0 .+(lo1-lo0).*rand(n_part)
	lat=la0 .+(la1-la0).*rand(n_part)
	(_,_,_,_,f,x,y)=InterpolationFactors(𝐷.Γ,lon,lat)
    m=findall( (f.!==0).*((!isnan).(x)) )
	df=DataFrame(x=x[m],y=y[m],f=f[m])
	
	custom🔴 = DataFrame(ID=Int[], fid=Int[], x=Float32[], y=Float32[], 
		z=Float32[], year=Float32[], t=Float32[])
	function custom🔧(sol,𝐹::𝐹_MeshArray3D;id=missing,𝑇=missing)
		df=postprocess_MeshArray(sol,𝐹,id=id,𝑇=𝑇)
		z=[sol[1,i,j][3] for i in 1:size(sol,2), j in 1:size(sol,3)]
		df.z=z[:]
		df.year=df.t ./86400/365
		return df
	end
	custom∫(prob) = solve(prob,Tsit5(),reltol=1e-8,abstol=1e-8,saveat=365/12*86400.0)
	
	𝐼=Individuals(𝑃,df.x,df.y,fill(z_init,n_part),df.f,
		(🔴=custom🔴,🔧=custom🔧,∫=custom∫))

	"$✓ Set up Individuals"
end

# ╔═╡ 9ffe84c0-dff0-11ea-2726-8924892df73a
begin
	𝑇=(0.0,𝐼.𝑃.𝑇[2])
	∫!(𝐼,𝑇)
	
	add_lonlat!(𝐼.🔴,𝐷.XC,𝐷.YC)
	
	🔴_by_ID = groupby(𝐼.🔴, :ID)
	🔴_by_t = groupby(𝐼.🔴, :t)
	
	ff(x) = x[1]
	tmp1=combine(🔴_by_ID,:lon => ff => :lo,:lat => ff => :la)
	𝐼.🔴.dlon=𝐼.🔴.lon - tmp1.lo[𝐼.🔴.ID]
	𝐼.🔴.dlat=𝐼.🔴.lat - tmp1.la[𝐼.🔴.ID]
	
	nt=length(unique(𝐼.🔴.t))

	"$✓ ∫!(𝐼,𝑇) etc"
end

# ╔═╡ f65ddffa-e63a-11ea-34a6-2fa9284e98fa
begin
	mx=20.0

	#f(x)=x[tt]-x[1]
	#f(x)=last(x).-first(x)
	#cdf = combine(🔴_by_ID,nrow,:lat => f => :dlat,:lon => f => :dlon)

	g(x)=x[tt]
	cdf = combine(🔴_by_ID,:dlat => g => :dlat,:dlon => g => :dlon)
	
	plt_hist=histogram2d(cdf.dlon,cdf.dlat,nbins = (10, 10),colorbar=false)
	#scatter(cdf.dlon,cdf.dlat,xlims=(-mx,mx),ylims=(-mx,mx))
	"$✓ dlon, dlat histogram2d"
end

# ╔═╡ 7d52252e-e006-11ea-2632-df2af831b52f
begin
	xx=vec(𝐷.Γ.XC[1][:,1])
	yy=vec(𝐷.Γ.YC[1][1,:])
	zz=transpose(log10.(𝐷.Γ.Depth[1]))
	𝐵=(x = xx, y = yy, z = zz)
	
	𝐶(g::ColorGradient) = RGB[g[z] for z=LinRange(0,1,length(🔴_by_t))]
	𝐶(t::Int) = 𝐶(cgrad(:inferno))[t]	
	"$✓ Set up plotting"
end

# ╔═╡ 0b12cf52-e6e3-11ea-1a01-dd0c49c9e641
begin
	plt_dlat = @df 🔴_by_t[1] density(:dlat, leg = :none, colour = 𝐶(1), ylims=(0,0.5))
	[@df 🔴_by_t[tt] density!(plt_dlat,:dlat, leg = :none, colour = 𝐶(tt)) for tt in 2:length(🔴_by_t)];
	density!(plt_dlat,🔴_by_t[tt].dlat, leg = :none, colour = :cyan, linewidth=4)
	"$✓ dlat density"
end

# ╔═╡ 6f70033a-e6cc-11ea-373e-6dcbaaa53d15
begin
	#clims=extrema(:t)
	#@df 𝐼.🔴 density(:lon, group = (:t), leg = :none, palette = cgrad(:ice))
	plt_dlon = @df 🔴_by_t[1] density(:dlon, leg = :none, colour = 𝐶(1), ylims=(0,0.5))
	[@df 🔴_by_t[tt] density!(plt_dlon,:dlon, leg = :none, colour = 𝐶(tt)) for tt in 2:length(🔴_by_t)];
	density!(plt_dlon,🔴_by_t[tt].dlon, leg = :none, colour = :cyan, linewidth=4)
	"$✓ dlon density"
end

# ╔═╡ a13d6ea6-dff1-11ea-0713-cb235e28cf79
begin
	plt_map=contourf(xx,yy,zz,clims=(-.5,4.),c = :ice, 
		colorbar=false, xlims=(-180.0,180.0),ylims=(-90.0,90.0))
	scatter!(plt_map,🔴_by_t[1].lon,🔴_by_t[1].lat,c=:gold,leg=:none,
		markersize=2.0, marker = (:dot, stroke(0)) )
	scatter!(plt_map,🔴_by_t[tt].lon,🔴_by_t[tt].lat,c=:red,leg=:none,
		markersize=2.0, marker = (:dot, stroke(0)) )
	
	plot(plt_map,plt_hist,plt_dlon,plt_dlat)
end

# ╔═╡ Cell order:
# ╟─bf19d29c-e70e-11ea-0153-d3a49981d56c
# ╟─4935fd46-e70f-11ea-386c-f9c444a20644
# ╟─9c80e722-e70f-11ea-22a6-0be2e85f3b8b
# ╟─e25eee9e-dfee-11ea-2a4c-3946ccb63876
# ╟─a13d6ea6-dff1-11ea-0713-cb235e28cf79
# ╟─f65ddffa-e63a-11ea-34a6-2fa9284e98fa
# ╟─0b12cf52-e6e3-11ea-1a01-dd0c49c9e641
# ╟─6f70033a-e6cc-11ea-373e-6dcbaaa53d15
# ╟─9ffe84c0-dff0-11ea-2726-8924892df73a
# ╟─7d52252e-e006-11ea-2632-df2af831b52f
# ╟─f75fae30-dfee-11ea-18ef-259321acfa2f
# ╟─747d446a-dfeb-11ea-3533-c9404fd41688
# ╟─1811955d-618d-49e3-9f81-63195e2632fe
# ╟─87c5e5d4-c343-49fb-bc7e-6b2c0e647a38
