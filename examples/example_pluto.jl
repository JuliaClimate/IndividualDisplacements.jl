### A Pluto.jl notebook ###
# v0.11.10

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
	using IndividualDisplacements, Plots, NetCDF, OrdinaryDiffEq, DataFrames
	using Statistics, StatsPlots, PlutoUI

	p=dirname(pathof(IndividualDisplacements))
    include(joinpath(p,"../examples/example123.jl"))
    include(joinpath(p,"../examples/helper_functions.jl"))
    include(joinpath(p,"../examples/recipes_plots.jl"))
	𝑃,Γ=OCCA_setup()	
	tmp=(Γ = Γ, m = "OCCA")
    𝑃=merge(𝑃,tmp)
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
$(@bind tt Slider(2:75; default=50, show_value=true))
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
	day=86400.0
	mon=365/12*day
	solv(prob) = solve(prob,Euler(),dt=2*day)
	
	lon=lo0 .+(lo1-lo0).*rand(n_part)
	lat=la0 .+(la1-la0).*rand(n_part)
	(xy,du)=initialize_lonlat(Γ,lon,lat)
	xy[3,:] .= z_init
	id=collect(1:size(xy,2))
	
	function solv(prob)
	  sol=solve(prob,Euler(),dt=10*86400.0)
	  nx,ny=𝑃.ioSize[1:2]
	  sol[1,:,:]=mod.(sol[1,:,:],nx)
	  sol[2,:,:]=mod.(sol[2,:,:],ny)
	  return sol
	end
	
	tr = DataFrame([fill(Int, 2) ; fill(Float64, 6)],[:ID, :fid, :x, :y, :z, :t, :lon, :lat])

	function postproc(sol,𝑃::NamedTuple;id=missing,𝑇=missing)
	  df=postprocess_lonlat(sol,𝑃,id=id,𝑇=𝑇)
	  #add third coordinate
	  z=sol[3,:,:]
	  df.z=z[:]
	  #to plot e.g. Pacific Ocean transports, shift longitude convention?
	  #df.lon[findall(df.lon .< 0.0 )] = df.lon[findall(df.lon .< 0.0 )] .+360.0
	  return df
	end
	
	"$✓ Set up Individuals"
end

# ╔═╡ 9ffe84c0-dff0-11ea-2726-8924892df73a
begin
	𝐼 = Individuals{Float64}(📌=deepcopy(xy), 🔴=deepcopy(tr), 🆔=id, 
		                     🚄 = dxyz_dt, ∫ = solv, 🔧 = postproc, 𝑃=𝑃)
	
	𝑇=(0.0,𝐼.𝑃.𝑇[2])
	∫!(𝐼,𝑇)
	
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
	mx=50.0

	#f(x)=x[tt]-x[1]
	#f(x)=last(x).-first(x)
	#cdf = combine(🔴_by_ID,nrow,:lat => f => :dlat,:lon => f => :dlon)

	g(x)=x[tt]
	cdf = combine(🔴_by_ID,:dlat => g => :dlat,:dlon => g => :dlon)
	
	plt_hist=histogram2d(cdf.dlon,cdf.dlat,nbins = (10, 10),
		xlims=(-mx,mx),ylims=(-mx,mx), colorbar=false)
	#scatter(cdf.dlon,cdf.dlat,xlims=(-mx,mx),ylims=(-mx,mx))
	"$✓ dlon, dlat histogram2d"
end

# ╔═╡ 7d52252e-e006-11ea-2632-df2af831b52f
begin
	x=vec(𝐼.𝑃.Γ["XC"][1][:,1])
	y=vec(𝐼.𝑃.Γ["YC"][1][1,:])
	z=transpose(log10.(𝐼.𝑃.Γ["Depth"][1]))
	𝐵=(x = x, y = y, z = z)
	
	𝐶(g::ColorGradient) = RGB[g[z] for z=LinRange(0,1,length(🔴_by_t))]
	𝐶(t::Int) = 𝐶(cgrad(:inferno))[t]	
	"$✓ Set up plotting"
end

# ╔═╡ 0b12cf52-e6e3-11ea-1a01-dd0c49c9e641
begin
	plt_dlat = @df 🔴_by_t[1] density(:dlat, leg = :none, colour = 𝐶(1), ylims=(0,0.5))
	[@df 🔴_by_t[tt] density!(plt_dlat,:dlat, leg = :none, colour = 𝐶(tt)) for tt in 2:length(🔴_by_t)];
	density!(plt_dlat,🔴_by_t[nt].dlat, leg = :none, colour = :cyan, linewidth=4)
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
	plt_map=contourf(x,y,z,clims=(-.5,4.),c = :ice, 
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
