### A Pluto.jl notebook ###
# v0.17.0

using Markdown
using InteractiveUtils

# ╔═╡ 3e9d08f8-3ea1-11ec-262e-bb3d43960aec
md"""# Simple Two Dimensional Flow Fields

Simulate trajectories of a particle cloud in a simple, two-dimensional, doubly-periodic, flow field. 

For additional documentation e.g. see :
[1](https://JuliaClimate.github.io/IndividualDisplacements.jl/dev/),
[2](https://JuliaClimate.github.io/MeshArrays.jl/dev/),
[3](https://docs.juliadiffeq.org/latest/solvers/ode_solve.html),
[4](https://en.wikipedia.org/wiki/Displacement_(vector))
"""

# ╔═╡ 192fc454-054c-4364-a9ed-1aa4969b612a
begin
	using IndividualDisplacements, DataFrames, MeshArrays
	using Plots, ColorSchemes, PlutoUI
	
	"done with loading packages"
end

# ╔═╡ bb710e3c-9f38-4feb-a241-f624d2fca943
TableOfContents()

# ╔═╡ 023222b6-e2ac-43cd-be10-f5d9b7ce0124
md"""## Define A Flow Field

Let's start with a simple, doubly periodic flow field defined by a streamfunction
and create the `FlowFields` data structure which will then be drive the 
individual displacement and trajectory computations.
"""

# ╔═╡ a2c57844-080e-4598-bef5-4d5dfb740a63
begin
	nx=16
	dx= π/nx
	XC = dx*(collect(1:2*nx) .- 0.5)
	YC = dx*(collect(1:nx) .- 0.5)
	
	fac=0.1
	f(x, y) = sin(x) + cos(y) #streamfunction
	ϕ = fac*[f(x, y) for x in XC,y in YC] #streamfunction
	uC = -fac*[sin(y) for x in XC,y in YC] #dphi/dy at cell center
	vC = -fac*[cos(x) for x in XC,y in YC] #-dphi/dx at cell center
	"done with defining flow field at grid cell centers"
end

# ╔═╡ 96b1d76e-1563-4f90-b821-feb75baea771
md"""
When starting with collocated velocities in m/s units (`uC,vC` at the grid cell center), one easily obtains the normalized, staggered C-grid velocities (`u,v`). The resulting  `u` (resp `v`) is staggered by `-0.5` grid point in direction `1` for `u` (`2` for `v`) relative to `uC,vC`. This staggering of variables and the normalization of velocities by corresponding grid scales (`dx` below) are  conventions that apply thrhougout `IndividualDisplacements.jl`.

These flowfields, for consecutive time steps, are then embedded in a `FlowFields` data structure, along with time range `(t0,t1)`. Between `t0` and `t1` velcocities are interpolated linearly going from `u0,v0` to `u1,v1`. By setting `u1=u0=u` and `v1=v0=v`, the flow field remains the same throughout `(t0,t1)`.
"""

# ╔═╡ aa86042f-3da6-4252-9493-9a713688b4b1
begin
	u=0.5*(circshift(uC, (1,0))+uC) /dx #staggered u converted to grid point units (m/s -> 1/s)
	v=0.5*(circshift(vC, (0,1))+vC) /dx #staggered v converted to grid point units (m/s -> 1/s)
	𝑇=(0.,10.)
	𝐹=FlowFields(u,u,v,v,𝑇)
	"done with staggered flow field definition"
end

# ╔═╡ 47c87570-ca36-468f-9b9d-5c41de490105
md"""## Initialize Individuals

Here we initialize 100 particles within a subdomainand wraps everything in the `Individuals` data structure.
"""

# ╔═╡ 53a79184-9a72-47a5-a6be-fc3c5fe7a096
begin
	np,nq=size(u)
	x=np*(0.4 .+ 0.2*rand(100))
	y=nq*(0.4 .+ 0.2*rand(100))
	𝐼=Individuals(𝐹,x,y)
end

# ╔═╡ dffe1032-a247-4008-be22-692abcbe458a
md"""## Compute Trajectories

The time period is `𝐼.𝑃.𝑇` by default, unless `∫!(𝐼,𝑇)` is called instead as done below. 

Note that the size of 🔴 is different from before -- this DataFrame is a record of the trajectories.
"""

# ╔═╡ 28a3af5d-c1b3-4d95-8d06-034e1ad4f585
begin
	∫!(𝐼,𝑇)
	𝐼
end

# ╔═╡ 6a6c82d8-b8d3-4a94-af6e-79039ceaa157
md"""## Plot Results"""

# ╔═╡ 9288b04f-0f38-4d8a-9af1-027147799079
function plot(ϕ,df)
    nx,ny=size(ϕ)
	if isa(df,GroupedDataFrame)
	    f=contour(-0.5 .+ (1:nx),-0.5 .+ (1:ny), transpose(ϕ), c = :black, colorbar=:none,
			linewidth = 0.5, xlims=(0,nx),ylims=(0,ny))	
		co=cgrad(:buda,1:length(df))
		for i in 1:length(df)
			scatter!(df[i].x,df[i].y,markersize=3.0, palette=co,
				marker = (:circle, stroke(0)), leg=:none)
		end
	else
	    f=contourf(-0.5 .+ (1:nx),-0.5 .+ (1:ny), transpose(ϕ), c = :blues, colorbar=:none,
			linewidth = 0.5, xlims=(0,nx),ylims=(0,ny))	
	    scatter!(df.x,df.y,markersize=4.0,c=:red,marker = (:circle, stroke(0)),leg=:none)
	end
	f
end

# ╔═╡ 0366f590-e132-4a7f-877b-2168566faf60
begin	
	🔴_by_t = groupby(𝐼.🔴, :t)
	plot(ϕ,🔴_by_t[end])
end

# ╔═╡ 4a5d4cf9-faa3-45e1-85fa-50132f4836ec
plot(ϕ,🔴_by_t[1:100:end])
#size(🔴_by_t[1:1000:end])

# ╔═╡ c0f0c2c6-7d74-4e2a-87f7-d96971e234de
md""" To generate a simple animation, try this for example:

```
anim = @animate for t in 1:10:length(🔴_by_t)
   plot(ϕ,🔴_by_t[t])
end

gif_file=joinpath(tempdir(),"RandomFlow.gif")
gif(anim, gif_file, fps = 20)
```
"""

# ╔═╡ eb979c4c-549e-4912-aac8-e816c93017a1
begin
	anim = @animate for t in 1:10:length(🔴_by_t)
	   plot(ϕ,🔴_by_t[t])
	end
	
	gif_file=joinpath(tempdir(),"RandomFlow.gif")
	gif(anim, gif_file, fps = 20)
end

# ╔═╡ 67c39f1e-8aa5-4fad-a0b8-a7fca1b17c36
md"""## Exercises

- change the initial distribution of particles
- increase the duration of the trajectories simulation
- treat the non-periodic domain case by padding `u,v` with zeros 
- make the flow field time variable `𝐹=FlowFields(-u,u,-v,v,(0.,10.))`
- replace `u,v` with your own two-dimensional flow fields 
"""

# ╔═╡ feb874c4-80dd-444f-beb5-90b90647d44d
md"""## Extras

Instead of using the common `Array` type it can be advantageous to use [MeshArrays.jl](https://juliaclimate.github.io/MeshArrays.jl/dev/) which provides functionalities for staggered vector fields and gridded domain decomposition. The `convert_to_FlowFields` convenience function does the conversion for you. The only other difference from the `Array` case is the need to provide a vector of subdomain indices to `Individuals`. Here this is just a vector of ones since `convert_to_FlowFields` does not decompose the gridded domain.

```
𝐹=convert_to_FlowFields(u,v,10.0)
𝐼=Individuals(𝐹,x,y,fill(1,length(x)))
```

The `random_flow_field` function, found below, provides generates random flow fields that can be used instead of the analytical formulation used above. The "Rotational Component" option is most similar to what done in the original example.

```
(U,V,Φ)=random_flow_field("Rotational Component";np=2*nx,nq=nx)
F=convert_to_FlowFields(U,V,10.0)
I=Individuals(F,x,y,fill(1,length(x)))
```

The other option, "Divergent Component", generates a purely divergent flow field instead. Try it and should you notice a qualitatively different outcome in terms of trajectories.

In general, user defined `uC, vC` fields may have both rotational and divergent components. [MeshArrays.jl](https://juliaclimate.github.io/MeshArrays.jl/dev/) provides an implementation of the [Helmholtz decomposition](https://en.wikipedia.org/wiki/Helmholtz_decomposition) to separate them out as sometimes needed.
"""

# ╔═╡ e927ba74-2c88-493e-b25b-910e23a63045
begin
	(U,V,Φ)=random_flow_field("Rotational Component";np=2*nx,nq=nx)
	F=convert_to_FlowFields(U,V,10.0)
	I=Individuals(F,x,y,fill(1,length(x)))
	∫!(I)
	I
end

# ╔═╡ ef8cbf8e-3e10-4615-9640-86cb8bc68288
function random_flow_field(option::String;np=12,nq=18)

	#define gridded domain
	Γ=simple_periodic_domain(np,nq)
	γ=Γ.XC.grid
	Γ=UnitGrid(γ;option="full")

    #initialize 2D field of random numbers
    tmp1=randn(Float64,Tuple(γ.ioSize))
    ϕ=γ.read(tmp1,MeshArray(γ,Float64))

    #apply smoother
    ϕ=smooth(ϕ,3*Γ.DXC,3*Γ.DYC,Γ);

	#derive flow field
	if option=="Divergent Component"
		#For the convergent / scalar potential case, ϕ is interpreted as being 
		#on center points -- hence the standard gradient function readily gives 
		#what we need
		(u,v)=gradient(ϕ,Γ) 
		tmp=(u[1],v[1],ϕ[1])
	elseif option=="Rotational Component"
		#For the rotational / streamfunction case, ϕ is interpreted as being 
		#on S/W corner points -- this is ok here since the grid is homegeneous, 
		#and conveniently yields an adequate substitution u,v <- -v,u; but note
		#that doing the same with gradient() would shift indices inconsistenly
		u=-(circshift(ϕ[1], (0,-1))-ϕ[1])
		v=(circshift(ϕ[1], (-1,0))-ϕ[1])
		tmp=(u,v,ϕ[1])
	else
		error("non-recognized option")
	end
	return tmp[1],tmp[2],tmp[3]
end

# ╔═╡ Cell order:
# ╟─3e9d08f8-3ea1-11ec-262e-bb3d43960aec
# ╟─192fc454-054c-4364-a9ed-1aa4969b612a
# ╟─bb710e3c-9f38-4feb-a241-f624d2fca943
# ╟─023222b6-e2ac-43cd-be10-f5d9b7ce0124
# ╠═a2c57844-080e-4598-bef5-4d5dfb740a63
# ╟─96b1d76e-1563-4f90-b821-feb75baea771
# ╠═aa86042f-3da6-4252-9493-9a713688b4b1
# ╟─47c87570-ca36-468f-9b9d-5c41de490105
# ╠═53a79184-9a72-47a5-a6be-fc3c5fe7a096
# ╟─dffe1032-a247-4008-be22-692abcbe458a
# ╠═28a3af5d-c1b3-4d95-8d06-034e1ad4f585
# ╟─6a6c82d8-b8d3-4a94-af6e-79039ceaa157
# ╟─9288b04f-0f38-4d8a-9af1-027147799079
# ╟─0366f590-e132-4a7f-877b-2168566faf60
# ╟─4a5d4cf9-faa3-45e1-85fa-50132f4836ec
# ╟─c0f0c2c6-7d74-4e2a-87f7-d96971e234de
# ╟─eb979c4c-549e-4912-aac8-e816c93017a1
# ╟─67c39f1e-8aa5-4fad-a0b8-a7fca1b17c36
# ╟─feb874c4-80dd-444f-beb5-90b90647d44d
# ╟─e927ba74-2c88-493e-b25b-910e23a63045
# ╠═ef8cbf8e-3e10-4615-9640-86cb8bc68288
