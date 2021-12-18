### A Pluto.jl notebook ###
# v0.17.3

using Markdown
using InteractiveUtils

# ╔═╡ 192fc454-054c-4364-a9ed-1aa4969b612a
begin
	using Pkg
	Pkg.activate()
	
	using IndividualDisplacements, DataFrames, MeshArrays, PlutoUI
	import CairoMakie as Mkie

	"done with loading packages"
end

# ╔═╡ 3e9d08f8-3ea1-11ec-262e-bb3d43960aec
md"""# Simple Two Dimensional Flow Fields

Simulate trajectories of a particle cloud in a simple, two-dimensional, doubly-periodic, flow field. 

For additional documentation e.g. see :
[1](https://JuliaClimate.github.io/IndividualDisplacements.jl/dev/),
[2](https://JuliaClimate.github.io/MeshArrays.jl/dev/),
[3](https://docs.juliadiffeq.org/latest/solvers/ode_solve.html),
[4](https://en.wikipedia.org/wiki/Displacement_(vector))
"""

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

# ╔═╡ 22dde07e-95e1-4a49-a3a1-fba7702dd74d
begin
	🔴_by_t = groupby(𝐼.🔴, :t)
	nt=length(🔴_by_t)

	time = Mkie.Node(nt)
	xp=Mkie.@lift( 🔴_by_t[$time].x )
	yp=Mkie.@lift( 🔴_by_t[$time].y )

	siz=size(ϕ)
	xx=-0.5 .+ collect(1:siz[1])
	yy=-0.5 .+ collect(1:siz[2])
	ll=collect(-1.0:0.2:1.0)*maximum(ϕ)
	
	Mkie.set_theme!(Mkie.theme_light())
	fig=Mkie.Figure(resolution = (900, 600))
	a = Mkie.Axis(fig[1, 1],xlabel="x",ylabel="y", title="Positions and Streamfunction")		
	Mkie.contourf!(a, xx,yy, ϕ, levels=ll, colormap=:grayC)
	Mkie.scatter!(a,🔴_by_t[1].x,🔴_by_t[1].y,color=:green2)
	Mkie.scatter!(a,xp,yp,color=:red)

	fig
end

# ╔═╡ 818b516f-c7ee-4da7-9954-8303c978bbd4
begin
	if false
		fil=joinpath(tempdir(),"random_flow_field.jl.mp4")
		Mkie.record(fig,fil, 1:10:nt; framerate = 20) do t
			time[] = t
		end
		LocalResource(fil)
	else
		"try code in this cell to generate animation"
	end
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
    ϕ=MeshArrays.smooth(ϕ,3*Γ.DXC,3*Γ.DYC,Γ);

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

# ╔═╡ e927ba74-2c88-493e-b25b-910e23a63045
begin
	(U,V,Φ)=random_flow_field("Rotational Component";np=2*nx,nq=nx)
	F=convert_to_FlowFields(U,V,10.0)
	I=Individuals(F,x,y,fill(1,length(x)))
	∫!(I)
	I
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
# ╟─22dde07e-95e1-4a49-a3a1-fba7702dd74d
# ╟─818b516f-c7ee-4da7-9954-8303c978bbd4
# ╟─67c39f1e-8aa5-4fad-a0b8-a7fca1b17c36
# ╟─feb874c4-80dd-444f-beb5-90b90647d44d
# ╟─e927ba74-2c88-493e-b25b-910e23a63045
# ╠═ef8cbf8e-3e10-4615-9640-86cb8bc68288
