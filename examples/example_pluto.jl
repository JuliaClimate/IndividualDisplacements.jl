### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 747d446a-dfeb-11ea-3533-c9404fd41688
begin
	using IndividualDisplacements, MeshArrays, DataFrames
	using Statistics, PlutoUI, OrdinaryDiffEq, StatsPlots

	p=dirname(pathof(IndividualDisplacements))
    include(joinpath(p,"../examples/flow_fields.jl"))
	𝑃,𝐷,Γ=OCCA_FlowFields()	
	tmp=(Γ = Γ, m = "OCCA")
    𝐷=merge(𝐷,tmp)
	✓ = "😃"
	"$✓ Set up packages, gridded domain, etc"
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
	xx=vec(𝐷.Γ["XC"][1][:,1])
	yy=vec(𝐷.Γ["YC"][1][1,:])
	zz=transpose(log10.(𝐷.Γ["Depth"][1]))
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
