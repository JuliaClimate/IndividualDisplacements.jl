### A Pluto.jl notebook ###
# v0.11.8

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
	using IndividualDisplacements, Plots, NetCDF, OrdinaryDiffEq, DataFrames, Statistics, StatsPlots

	p=dirname(pathof(IndividualDisplacements))
    include(joinpath(p,"../examples/example123.jl"))
    include(joinpath(p,"../examples/helper_functions.jl"))
    include(joinpath(p,"../examples/recipes_plots.jl"))
	𝑃,Γ=OCCA_setup()	
	tmp=(Γ = Γ, m = "OCCA")
    𝑃=merge(𝑃,tmp)
	✓ = "😃" 
end

# ╔═╡ 09f32060-dff7-11ea-1276-add3ff402f18
@bind lon0 html"<input type=range min=-180.0 max=180.0>"

# ╔═╡ 3850bd7c-dff8-11ea-3a9b-1ba4c812e614
@bind lat0 html"<input type=range min=-90.0 max=90.0>"

# ╔═╡ d79a0800-dffa-11ea-1751-9de46bd95c3c
@bind npar html"<input type=range min=0 max=1000>"

# ╔═╡ 3039d8a0-dffb-11ea-33a2-03f957d4efff
@bind klev html"<input type=range min=0 max=10>"

# ╔═╡ f1da7264-e675-11ea-0f7e-ef6ff63c9123
@bind nt html"<input type=range min=2 max=75>"

# ╔═╡ e25eee9e-dfee-11ea-2a4c-3946ccb63876
begin
	lo0=lon0; lo1=lo0+5
	la0=lat0; la1=la0+5
	z_init=klev
	n_part=npar
	(particles=n_part,longitudes=(lo0,lo1),latitudes=(la0,la1),level=z_init)
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
	
	✓
end

# ╔═╡ 9ffe84c0-dff0-11ea-2726-8924892df73a
begin
   𝐼 = Individuals{Float64}(📌=xy, 🔴=tr, 🆔=id, 🚄 = dxyz_dt, ∫ = solv, 🔧 = postproc, 𝑃=𝑃)
   𝑇=(0.0,𝐼.𝑃.𝑇[2])
   ∫!(𝐼,𝑇)
end

# ╔═╡ e59d725e-e61f-11ea-255f-0b5550795ad9
begin
	🔴_by_ID = groupby(𝐼.🔴, :ID)
	🔴_by_t = groupby(𝐼.🔴, :t)
	
	ff(x) = x[1]
	tmp1=combine(🔴_by_ID,:lon => ff => :ini_lon,:lat => ff => :ini_lat)
	𝐼.🔴.dlon=𝐼.🔴.lon - tmp1.ini_lon[𝐼.🔴.ID]
	𝐼.🔴.dlat=𝐼.🔴.lat - tmp1.ini_lat[𝐼.🔴.ID]

#	df_a = DataFrame(i=1:10, x=0.1:0.1:1.0, y='a':'j');
#	nunique(x)=length(unique(x))
#	describe(df_a,:nunique => nunique)
end

# ╔═╡ f65ddffa-e63a-11ea-34a6-2fa9284e98fa
begin
	mx=50.0
	f(x)=x[nt]-x[1]
	g(x)=x[nt]
	#f(x)=last(x).-first(x)
	#cdf = combine(🔴_by_ID,nrow,:lat => f => :dlat,:lon => f => :dlon)
	cdf = combine(🔴_by_ID,:dlat => g => :dlat,:dlon => g => :dlon)
	histogram2d(cdf.dlon,cdf.dlat,nbins = (10, 10) ,xlims=(-mx,mx),ylims=(-mx,mx))
	#scatter(cdf.dlon,cdf.dlat,xlims=(-mx,mx),ylims=(-mx,mx))
end

# ╔═╡ 6a1a6d38-e6d5-11ea-1182-f7e690710745
describe(🔴_by_ID[2])

# ╔═╡ 7d52252e-e006-11ea-2632-df2af831b52f
begin
	x=vec(𝐼.𝑃.Γ["XC"][1][:,1])
	y=vec(𝐼.𝑃.Γ["YC"][1][1,:])
	z=transpose(log10.(𝐼.𝑃.Γ["Depth"][1]))
	𝐵=(x = x, y = y, z = z)
	
	𝐶(g::ColorGradient) = RGB[g[z] for z=LinRange(0,1,length(🔴_by_t))]
	𝐶(t::Int) = 𝐶(cgrad(:inferno))[t]	
	✓
end

# ╔═╡ a13d6ea6-dff1-11ea-0713-cb235e28cf79
begin
	plt=contourf(x,y,z,clims=(-.5,4.),c = :ice, 
		colorbar=false, xlims=(-180.0,180.0),ylims=(-90.0,90.0))
	scatter!(plt,🔴_by_t[1].lon,🔴_by_t[1].lat,c=:gold,leg=:none,
		markersize=2.0, marker = (:dot, stroke(0)) )
	scatter!(plt,🔴_by_t[nt].lon,🔴_by_t[nt].lat,c=:red,leg=:none,
		markersize=2.0, marker = (:dot, stroke(0)) )
end

# ╔═╡ 0b12cf52-e6e3-11ea-1a01-dd0c49c9e641
begin
	plt2 = @df 🔴_by_t[1] density(:dlat, leg = :none, colour = 𝐶(1), ylims=(0,0.5))
	[@df 🔴_by_t[tt] density!(plt2,:dlat, leg = :none, colour = 𝐶(tt)) for tt in 2:length(🔴_by_t)];
	density!(plt2,🔴_by_t[nt].dlat, leg = :none, colour = :cyan, linewidth=4)
	plt2
end

# ╔═╡ 6f70033a-e6cc-11ea-373e-6dcbaaa53d15
begin
	#clims=extrema(:t)
	#@df 𝐼.🔴 density(:lon, group = (:t), leg = :none, palette = cgrad(:ice))
	plt_dlon = @df 🔴_by_t[1] density(:dlon, leg = :none, colour = 𝐶(1), ylims=(0,0.5))
	[@df 🔴_by_t[tt] density!(plt_dlon,:dlon, leg = :none, colour = 𝐶(tt)) for tt in 2:length(🔴_by_t)];
	density!(plt_dlon,🔴_by_t[nt].dlon, leg = :none, colour = :cyan, linewidth=4)
	plt_dlon
end

# ╔═╡ Cell order:
# ╠═09f32060-dff7-11ea-1276-add3ff402f18
# ╠═3850bd7c-dff8-11ea-3a9b-1ba4c812e614
# ╠═d79a0800-dffa-11ea-1751-9de46bd95c3c
# ╠═3039d8a0-dffb-11ea-33a2-03f957d4efff
# ╠═f1da7264-e675-11ea-0f7e-ef6ff63c9123
# ╟─e25eee9e-dfee-11ea-2a4c-3946ccb63876
# ╟─a13d6ea6-dff1-11ea-0713-cb235e28cf79
# ╟─f65ddffa-e63a-11ea-34a6-2fa9284e98fa
# ╟─0b12cf52-e6e3-11ea-1a01-dd0c49c9e641
# ╟─6f70033a-e6cc-11ea-373e-6dcbaaa53d15
# ╠═6a1a6d38-e6d5-11ea-1182-f7e690710745
# ╟─e59d725e-e61f-11ea-255f-0b5550795ad9
# ╟─9ffe84c0-dff0-11ea-2726-8924892df73a
# ╟─7d52252e-e006-11ea-2632-df2af831b52f
# ╟─f75fae30-dfee-11ea-18ef-259321acfa2f
# ╟─747d446a-dfeb-11ea-3533-c9404fd41688
