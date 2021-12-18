### A Pluto.jl notebook ###
# v0.17.3

using Markdown
using InteractiveUtils

# ╔═╡ 104ce9b0-3fd1-11ec-3eff-3b029552e3d9
begin
	using Pkg
	Pkg.activate()
	
	using IndividualDisplacements, OceanStateEstimation, DataFrames, Statistics, CSV
	using MeshArrays, MITgcmTools, NetCDF, Plots
	"done with loading packages"
end

# ╔═╡ c9e9faa8-f5f0-479c-bc85-877ff7114883
md"""# Global Climatology

Advect particles with climatological monthly mean flow at selected depth level
(e.g. `k=10` for 95 m) from a global ocean state estimate ([ECCO v4 r2](https://eccov4.readthedocs.io/en/latest/) ; see also <https://ecco-group.org>)
which is here repeated for `ny` years. For additional documentation e.g. see :
[1](https://JuliaClimate.github.io/MeshArrays.jl/dev/),
[2](https://JuliaClimate.github.io/IndividualDisplacements.jl/dev/),
[3](https://docs.juliadiffeq.org/latest/solvers/ode_solve.html),
[4](https://en.wikipedia.org/wiki/Displacement_(vector))

[![simulated particle movie (5m)](https://user-images.githubusercontent.com/20276764/84766999-b801ad80-af9f-11ea-922a-610ad8a257dc.png)](https://youtu.be/W5DNqJG9jt0)
"""

# ╔═╡ 7fec71b4-849f-4369-bec2-26bfe2e00a97
md"""## 1. Grid and Velocity Files"""

# ╔═╡ 07e65622-3698-4dd8-b718-83588e116e58
begin
	#pth1=dirname(pathof(IndividualDisplacements))
	#include(joinpath(pth1,"../examples/helper_functions.jl"))
	
	OceanStateEstimation.get_ecco_velocity_if_needed();
	γ=GridSpec("LatLonCap",MeshArrays.GRID_LLC90)
	Γ=GridLoad(γ;option="full")
	Γ=merge(Γ,MeshArrays.NeighborTileIndices_cs(Γ))

	func=(u -> MeshArrays.update_location_llc!(u,Γ))
	Γ=merge(Γ,(; update_location! = func))

	"Done with Grid and Velocity Files"
end

# ╔═╡ 94ca10ae-6a8a-4038-ace0-07d7d9026712
md"""## 2. `FlowFields` Data Structure

The following parameters are used:

- select vertical level (k=1 by default; k=0 for 3D)
- select duration in years (ny=1, nm=1 by default)
- read and process grid variables
- return FlowFields (𝑃) and ancillary variables etc (𝐷) 
- read & normalize velocities (𝐷.🔄)
"""

# ╔═╡ f1215951-2eb2-490b-875a-91c1205b8f63
md"""## 3. Main Computation Loop

### 3.1 initial particle positions randomly over Global Ocean
### 3.2 initial integration from time 0 to 0.5 month
"""

# ╔═╡ 6158a5e4-89e0-4496-ab4a-044d1e3e8cc0
md""" ### 3.2 Iteration function example

In addition, `step!` is defined to provide additional flexibility around `∫!` :

- `𝐷.🔄(𝐼.𝑃,t_ϵ)` resets the velocity input streams to bracket t_ϵ=𝐼.𝑃.𝑇[2]+eps(𝐼.𝑃.𝑇[2]) 
- `reset_📌!(𝐼)` randomly selects a fraction (`𝐷.frac`) of the particles and resets their positions before each integration period. This tends to maintain homogeneous coverage of the Global Ocean by particles.
- `∫!(𝐼)` then solves for the individual trajectories over one month, with updated velocity fields (𝐼.𝑃.u0 etc), and adds diagnostics to the DataFrame used to record variables along the trajectory (𝐼.🔴).
"""

# ╔═╡ 7efadea7-4542-40cf-893a-40a75e9c52be
md"""### 3.3 Iterate For `ny*12` Months"""

# ╔═╡ 15077957-64d5-46a5-8a87-a76ad619cf38
md"""## 3.4 Compute summary statistics

See [DataFrames.jl](https://juliadata.github.io/DataFrames.jl/latest/) documentation for detail and additinal functionalities.
"""

# ╔═╡ 4b887e2f-7505-4db2-8784-400a786fba10
begin
	function CalcIntFac(Γ)
	    lon=[i for i=20.:2.0:380., j=-79.:2.0:89.]
	    lat=[j for i=20.:2.0:380., j=-79.:2.0:89.]
		(f,i,j,w,_,_,_)=InterpolationFactors(Γ,vec(lon),vec(lat))
		IntFac=(lon=lon,lat=lat,f=f,i=i,j=j,w=w)	
	end
	IntFac=CalcIntFac(Γ)
	"Done with interoikation coefficients for map"
end

# ╔═╡ af74c6c8-1859-4fdf-ae2b-5af8dccdee60
"""
    OceanDepthLog(Γ,IntFac)

Compute Ocean depth logarithm on regular grid.
"""
function OceanDepthLog(Γ,IntFac)
#    lon=[i for i=20.:2.0:380., j=-79.:2.0:89.]
#    lat=[j for i=20.:2.0:380., j=-79.:2.0:89.]
	DL=interp_to_lonlat(Γ.Depth,IntFac)
	DL[findall(DL.<0)].=0
    DL=transpose(log10.(DL))
    DL[findall((!isfinite).(DL))].=NaN
    return DL
#	return (lon=lon[:,1],lat=lat[1,:],fld=DL,rng=(1.5,5))
end

# ╔═╡ 14f7eadb-9ac4-41cd-b773-8b17d0e69a2c
"""
    read_velocities(γ::gcmgrid,t::Int,pth::String)

Read velocity components `u,v` from files in `pth`for time `t`
"""
function read_velocities(γ::gcmgrid,t::Int,pth::String)
    u=read_nctiles("$pth"*"UVELMASS/UVELMASS","UVELMASS",γ,I=(:,:,:,t))
    v=read_nctiles("$pth"*"VVELMASS/VVELMASS","VVELMASS",γ,I=(:,:,:,t))
    return u,v
end

# ╔═╡ 11ea0fe5-b713-453f-ab66-77c75fd74ea4
begin
	"""
	    update_FlowFields!(𝑃::𝐹_MeshArray2D,𝐷::NamedTuple,t::Float64)
	
	Update flow field arrays (in 𝑃), 𝑃.𝑇, and ancillary variables (in 𝐷) 
	according to the chosen time `t` (in `seconds`). 
	
	_Note: for now, it is assumed that (1) the time interval `dt` between 
	consecutive records is diff(𝑃.𝑇), (2) monthly climatologies are used 
	with a periodicity of 12 months, (3) vertical 𝑃.k is selected_
	"""
	function update_FlowFields!(𝑃::𝐹_MeshArray2D,𝐷::NamedTuple,t::Float64)
	    dt=𝑃.𝑇[2]-𝑃.𝑇[1]
	
	    m0=Int(floor((t+dt/2.0)/dt))
	    m1=m0+1
	    t0=m0*dt-dt/2.0
	    t1=m1*dt-dt/2.0
	
	    m0=mod(m0,12)
	    m0==0 ? m0=12 : nothing
	    m1=mod(m1,12)
	    m1==0 ? m1=12 : nothing
	
	    (U,V)=read_velocities(𝑃.u0.grid,m0,𝐷.pth)
	    u0=U[:,𝐷.k]; v0=V[:,𝐷.k]
	    u0[findall(isnan.(u0))]=0.0; v0[findall(isnan.(v0))]=0.0 #mask with 0s rather than NaNs
	    u0=u0.*𝐷.iDXC; v0=v0.*𝐷.iDYC; #normalize to grid units
	    (u0,v0)=exchange(u0,v0,1) #add 1 point at each edge for u and v
	
	    (U,V)=read_velocities(𝑃.u0.grid,m1,𝐷.pth)
	    u1=U[:,𝐷.k]; v1=V[:,𝐷.k]
	    u1[findall(isnan.(u1))]=0.0; v1[findall(isnan.(v1))]=0.0 #mask with 0s rather than NaNs
	    u1=u1.*𝐷.iDXC; v1=v1.*𝐷.iDYC; #normalize to grid units
	    (u1,v1)=exchange(u1,v1,1) #add 1 point at each edge for u and v
	
	    𝑃.u0[:]=u0[:]
	    𝑃.u1[:]=u1[:]
	    𝑃.v0[:]=v0[:]
	    𝑃.v1[:]=v1[:]
	    𝑃.𝑇[:]=[t0,t1]
	
	end
end

# ╔═╡ b9b561f8-da40-423a-a7e0-2bf9eafc6e57

"""
    update_FlowFields!(𝑃::𝐹_MeshArray3D,𝐷::NamedTuple,t::Float64)

Update flow field arrays (in 𝑃), 𝑃.𝑇, and ancillary variables (in 𝐷) 
according to the chosen time `t` (in `seconds`). 

_Note: for now, it is assumed that (1) the time interval `dt` between 
consecutive records is diff(𝑃.𝑇), (2) monthly climatologies are used 
with a periodicity of 12 months, (3) vertical 𝑃.k is selected_
"""
function update_FlowFields!(𝑃::𝐹_MeshArray3D,𝐷::NamedTuple,t::Float64)
    dt=𝑃.𝑇[2]-𝑃.𝑇[1]

    m0=Int(floor((t+dt/2.0)/dt))
    m1=m0+1
    t0=m0*dt-dt/2.0
    t1=m1*dt-dt/2.0

    m0=mod(m0,12)
    m0==0 ? m0=12 : nothing
    m1=mod(m1,12)
    m1==0 ? m1=12 : nothing

    (_,nr)=size(𝐷.Γ.hFacC)

    (U,V)=read_velocities(𝑃.u0.grid,m0,𝐷.pth)
    u0=U; v0=V
    u0[findall(isnan.(u0))]=0.0; v0[findall(isnan.(v0))]=0.0 #mask with 0s rather than NaNs
    for k=1:nr
        u0[:,k]=u0[:,k].*𝐷.iDXC; v0[:,k]=v0[:,k].*𝐷.iDYC; #normalize to grid units
        (tmpu,tmpv)=exchange(u0[:,k],v0[:,k],1) #add 1 point at each edge for u and v
        u0[:,k]=tmpu
        v0[:,k]=tmpv
    end
    w0=read_nctiles(𝐷.pth*"WVELMASS/WVELMASS","WVELMASS",𝑃.u0.grid,I=(:,:,:,m0))
    w0[findall(isnan.(w0))]=0.0 #mask with 0s rather than NaNs

    (U,V)=read_velocities(𝑃.u0.grid,m1,𝐷.pth)
    u1=U; v1=V
    u1[findall(isnan.(u1))]=0.0; v1[findall(isnan.(v1))]=0.0 #mask with 0s rather than NaNs
    for k=1:nr
        u1[:,k]=u1[:,k].*𝐷.iDXC; v1[:,k]=v1[:,k].*𝐷.iDYC; #normalize to grid units
        (tmpu,tmpv)=exchange(u1[:,k],v1[:,k],1) #add 1 point at each edge for u and v
        u1[:,k]=tmpu
        v1[:,k]=tmpv
    end
    w1=read_nctiles(𝐷.pth*"WVELMASS/WVELMASS","WVELMASS",𝑃.u0.grid,I=(:,:,:,m1))
    w1[findall(isnan.(w1))]=0.0 #mask with 0s rather than NaNs

    𝑃.u0[:,:]=u0[:,:]
    𝑃.u1[:,:]=u1[:,:]
    𝑃.v0[:,:]=v0[:,:]
    𝑃.v1[:,:]=v1[:,:]
    for k=1:nr
        tmpw=exchange(-w0[:,k],1)
        𝑃.w0[:,k]=tmpw./𝐷.Γ.DRC[k]
        tmpw=exchange(-w1[:,k],1)
        𝑃.w1[:,k]=tmpw./𝐷.Γ.DRC[k]
    end
    𝑃.w0[:,1]=0*exchange(-w0[:,1],1)
    𝑃.w1[:,1]=0*exchange(-w1[:,1],1)
    𝑃.w0[:,nr+1]=0*exchange(-w0[:,1],1)
    𝑃.w1[:,nr+1]=0*exchange(-w1[:,1],1)

    #θ0=read_nctiles(𝐷.pth*"THETA/THETA","THETA",𝑃.u0.grid,I=(:,:,:,m0))
    #θ0[findall(isnan.(θ0))]=0.0 #mask with 0s rather than NaNs
    #𝐷.θ0[:,:]=θ0[:,:]

    #θ1=read_nctiles(𝐷.pth*"THETA/THETA","THETA",𝑃.u0.grid,I=(:,:,:,m1))
    #θ1[findall(isnan.(θ1))]=0.0 #mask with 0s rather than NaNs
    #𝐷.θ1[:,:]=θ1[:,:]

    𝑃.𝑇[:]=[t0,t1]
end

# ╔═╡ d466146a-f5b2-41c7-9415-da4a24a61209
"""
    set_up_FlowFields(k::Int,Γ::NamedTuple,pth::String)

Define `FlowFields` data structure (𝑃) for the specified grid (`Γ` dictionary), 
vertical level (`k`), and  file location (`pth`).

_Note: the initial implementation approximates month durations to 
365 days / 12 months for simplicity and sets 𝑃.𝑇 to [-mon/2,mon/2]_
"""
function set_up_FlowFields(k::Int,Γ::NamedTuple,pth::String)
    XC=exchange(Γ.XC) #add 1 lon point at each edge
    YC=exchange(Γ.YC) #add 1 lat point at each edge
    iDXC=1. ./Γ.DXC
    iDYC=1. ./Γ.DYC
    γ=Γ.XC.grid
    mon=86400.0*365.0/12.0
    func=Γ.update_location!

    if k==0
        msk=Γ.hFacC
        (_,nr)=size(msk)
        𝑃=FlowFields(MeshArray(γ,Float64,nr),MeshArray(γ,Float64,nr),
        MeshArray(γ,Float64,nr),MeshArray(γ,Float64,nr),
        MeshArray(γ,Float64,nr+1),MeshArray(γ,Float64,nr+1),
        [-mon/2,mon/2],func)
    else
        msk=Γ.hFacC[:, k]
        𝑃=FlowFields(MeshArray(γ,Float64),MeshArray(γ,Float64),
        MeshArray(γ,Float64),MeshArray(γ,Float64),[-mon/2,mon/2],func)    
    end

	𝐷 = (🔄 = update_FlowFields!, pth=pth,
	 XC=XC, YC=YC, iDXC=iDXC, iDYC=iDYC, 
	 k=k, msk=msk, θ0=similar(msk), θ1=similar(msk))

    𝐷 = merge(𝐷 , MeshArrays.NeighborTileIndices_cs(Γ))
    
    return 𝑃,𝐷
end

# ╔═╡ 218b9beb-68f2-4498-a96d-08e0719b4cff
begin
	#func=(u -> update_location_llc!(u,𝐷))
	#Γ=merge(Γ,(; update_location! = func))

	ny=1
	nm=1
	k=1

	𝑃,𝐷=set_up_FlowFields(k,Γ,ECCOclim_path)

	#add parameters for use in reset! and grid variables
    frac=0.01 #fraction of the particles reset per month (0.05 for k<=10)
	tmp=(frac=frac, Γ=Γ)
	𝐷=merge(𝐷,tmp)
	
	𝐷.🔄(𝑃,𝐷,0.0)
end

# ╔═╡ f727992f-b72a-45bc-93f1-cc8daf89af0f
begin
	np=500
	
	#xy = init_global_randn(np,𝐷)
	#df=DataFrame(x=xy[1,:],y=xy[2,:],f=xy[3,:])
	
	p=dirname(pathof(IndividualDisplacements))
	fil=joinpath(p,"../examples/worldwide/global_ocean_circulation.csv")
	df=DataFrame(CSV.File(fil))

	if !(k==0)
		𝐼=Individuals(𝑃,df.x[1:np],df.y[1:np],df.f[1:np])
	else
		kk=2.5
		𝐼=Individuals(𝑃,df.x[1:np],df.y[1:np],fill(kk,np),df.f[1:np])
	end
	fieldnames(typeof(𝐼))
end

# ╔═╡ 1495fda9-e46b-424e-922a-3b823f3fe200
𝐼

# ╔═╡ cc7cb4a8-86ea-42b0-bbb9-ca78469ad4ad
df

# ╔═╡ a3e45927-5d53-42be-b7b7-489d6e7a6fe5
begin
	📌ini=deepcopy(𝐼.📌)
	𝑇=(0.0,𝐼.𝑃.𝑇[2])
	∫!(𝐼,𝑇)
	✔1="done"
end

# ╔═╡ c57f60b8-cec6-4ef0-bb63-0201c18c9ece
"""
    reset_📌!(𝐼::Individuals,frac::Number,📌::Array)

Randomly select a fraction (frac) of the particles and reset 
their positions (𝐼.📌) to a random subset of the specificed 📌.
"""
function reset_📌!(𝐼::Individuals,frac::Number,📌::Array)
    np=length(𝐼.🆔)
    n_reset = Int(round(𝐷.frac*np))
    k_reset = rand(1:np, n_reset)
    l_reset = rand(1:np, n_reset)
    𝐼.📌[k_reset]=deepcopy(📌[l_reset])
    isempty(𝐼.🔴.ID) ? m=maximum(𝐼.🆔) : m=max(maximum(𝐼.🔴.ID),maximum(𝐼.🆔))
    𝐼.🆔[k_reset]=collect(1:n_reset) .+ m
end

# ╔═╡ a2375720-f599-43b9-a7fb-af17956309b6
function step!(𝐼::Individuals)
    t_ϵ=𝐼.𝑃.𝑇[2]+eps(𝐼.𝑃.𝑇[2])
    𝐷.🔄(𝐼.𝑃,𝐷,t_ϵ)
    reset_📌!(𝐼,𝐷.frac,📌ini)
    ∫!(𝐼)
end

# ╔═╡ 1044c5aa-1a56-45b6-a4c6-63d24eea878d
begin
	✔1
	[step!(𝐼) for y=1:ny, m=1:nm]
	add_lonlat!(𝐼.🔴,𝐷.XC,𝐷.YC)
	✔2="done"
end

# ╔═╡ 6e43a2af-bf01-4f42-a4ba-1874a8cf4885
begin
	✔2
	gdf = groupby(𝐼.🔴, :ID)
	sgdf= combine(gdf,nrow,:lat => mean)
end

# ╔═╡ e1cdcac9-c3cc-4ce4-a477-452ca460a3d5
begin
	fig=plot(;xlims=(-180,180),ylims=(-90,90),legend=:none)
	p!(x,y)=scatter!(fig,x,y,markersize=1.1,markerstrokewidth=0)
	[p!(gdf[i].lon,gdf[i].lat) for i in rand(collect(1:length(gdf)),10)]
	fig
end

# ╔═╡ 4a7ba3ff-449a-44e1-ad10-1de15a6d31cc
"""
    map(𝐼::Individuals,background::NamedTuple)

Plot initial and final positions, superimposed on a map of ocean depth log.
"""
function map(𝐼::Individuals,𝐵::NamedTuple)
    xlims=extrema(𝐵.lon)
    ylims=extrema(𝐵.lat)
    plt=contourf(𝐵.lon,𝐵.lat,𝐵.fld,clims=𝐵.rng,c = :ice, 
    colorbar=false, xlims=xlims,ylims=ylims)

    🔴_by_t = groupby(𝐼.🔴, :t)
    lo=deepcopy(🔴_by_t[1].lon); lo[findall(lo.<xlims[1])]=lo[findall(lo.<xlims[1])].+360
    scatter!(lo,🔴_by_t[1].lat,markersize=2.5,c=:red,leg=:none,marker = (:circle, stroke(0)))
    lo=deepcopy(🔴_by_t[end].lon); lo[findall(lo.<xlims[1])]=lo[findall(lo.<xlims[1])].+360
    scatter!(lo,🔴_by_t[end].lat,markersize=2.5,c=:yellow,leg=:none,marker = (:dot, stroke(0)))

    return plt
end

# ╔═╡ b4841dc0-c257-45e0-8657-79121f2c9ce8
begin
	DL=(lon=IntFac.lon[:,1],lat=IntFac.lat[1,:],fld=OceanDepthLog(Γ,IntFac),rng=(1.5,5))
	map(𝐼,DL)
end

# ╔═╡ Cell order:
# ╟─c9e9faa8-f5f0-479c-bc85-877ff7114883
# ╟─104ce9b0-3fd1-11ec-3eff-3b029552e3d9
# ╟─7fec71b4-849f-4369-bec2-26bfe2e00a97
# ╟─07e65622-3698-4dd8-b718-83588e116e58
# ╟─94ca10ae-6a8a-4038-ace0-07d7d9026712
# ╟─218b9beb-68f2-4498-a96d-08e0719b4cff
# ╟─f1215951-2eb2-490b-875a-91c1205b8f63
# ╟─f727992f-b72a-45bc-93f1-cc8daf89af0f
# ╟─1495fda9-e46b-424e-922a-3b823f3fe200
# ╟─cc7cb4a8-86ea-42b0-bbb9-ca78469ad4ad
# ╟─a3e45927-5d53-42be-b7b7-489d6e7a6fe5
# ╟─6158a5e4-89e0-4496-ab4a-044d1e3e8cc0
# ╟─a2375720-f599-43b9-a7fb-af17956309b6
# ╟─7efadea7-4542-40cf-893a-40a75e9c52be
# ╟─1044c5aa-1a56-45b6-a4c6-63d24eea878d
# ╟─15077957-64d5-46a5-8a87-a76ad619cf38
# ╟─6e43a2af-bf01-4f42-a4ba-1874a8cf4885
# ╟─e1cdcac9-c3cc-4ce4-a477-452ca460a3d5
# ╟─b4841dc0-c257-45e0-8657-79121f2c9ce8
# ╟─4b887e2f-7505-4db2-8784-400a786fba10
# ╟─af74c6c8-1859-4fdf-ae2b-5af8dccdee60
# ╟─d466146a-f5b2-41c7-9415-da4a24a61209
# ╟─11ea0fe5-b713-453f-ab66-77c75fd74ea4
# ╟─b9b561f8-da40-423a-a7e0-2bf9eafc6e57
# ╟─14f7eadb-9ac4-41cd-b773-8b17d0e69a2c
# ╟─c57f60b8-cec6-4ef0-bb63-0201c18c9ece
# ╟─4a7ba3ff-449a-44e1-ad10-1de15a6d31cc
