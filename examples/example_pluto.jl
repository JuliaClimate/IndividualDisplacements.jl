### A Pluto.jl notebook ###
# v0.11.7

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
	
	tr = DataFrame( ID=[], x=[], y=[], t = [], lon=[], lat=[], z=[], fid=[])
	
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
    𝐼 = Individuals{Float64}(📌=xy, 🔴=tr, 🆔=id, ⎔ = dxyz_dt, ∫ = solv, ⟁ = postproc, 𝑃=𝑃)
	start!(𝐼)
end

# ╔═╡ 7d52252e-e006-11ea-2632-df2af831b52f
begin
	x=vec(𝐼.𝑃.Γ["XC"][1][:,1])
	y=vec(𝐼.𝑃.Γ["YC"][1][1,:])
	z=transpose(log10.(𝐼.𝑃.Γ["Depth"][1]))
	#plt=contourf([x;x .+ 360],y,[z z],clims=(-.5,4.),c = :ice, 
	plt=contourf(x,y,z,clims=(-.5,4.),c = :ice, 
		colorbar=false, xlims=(-180.0,180.0),ylims=(-90.0,90.0))
	✓
end

# ╔═╡ a13d6ea6-dff1-11ea-0713-cb235e28cf79
begin	
df = 𝐼.🔴[ findall(𝐼.🔴.t .== minimum(𝐼.🔴.t)) , :]
scatter!(plt,df.lon,df.lat,markersize=2.0,c=:gold,leg=:none,
        marker = (:dot, stroke(0)))
df = 𝐼.🔴[ findall(𝐼.🔴.t .== maximum(𝐼.🔴.t)) , :]
scatter!(plt,df.lon,df.lat,markersize=2.0,c=:red,leg=:none,
        marker = (:dot, stroke(0)))
end

# ╔═╡ Cell order:
# ╟─747d446a-dfeb-11ea-3533-c9404fd41688
# ╠═09f32060-dff7-11ea-1276-add3ff402f18
# ╠═3850bd7c-dff8-11ea-3a9b-1ba4c812e614
# ╠═d79a0800-dffa-11ea-1751-9de46bd95c3c
# ╠═3039d8a0-dffb-11ea-33a2-03f957d4efff
# ╟─e25eee9e-dfee-11ea-2a4c-3946ccb63876
# ╟─a13d6ea6-dff1-11ea-0713-cb235e28cf79
# ╟─f75fae30-dfee-11ea-18ef-259321acfa2f
# ╟─9ffe84c0-dff0-11ea-2726-8924892df73a
# ╟─7d52252e-e006-11ea-2632-df2af831b52f
