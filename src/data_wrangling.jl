"""
    convert_to_FlowFields(U::Array{T,2},V::Array{T,2},t1::T) where T

Convert a pair of U,V arrays (staggered C-grid velocity field in 2D) to
a `𝐹_MeshArray2D` struct ready for integration of individual displacements
from time `t0=0` to time `t1`.

```
_,u,v=random_flow_field()
𝐹=convert_to_FlowFields(u,v,10.0)
```
"""
function convert_to_FlowFields(U::Array{T,2},V::Array{T,2},t1::T) where T
    np,nq=size(U)
    Γ=simple_periodic_domain(np,nq)

    g=Γ["XC"].grid
    u=MeshArray(g,[U])
    v=MeshArray(g,[V])
    (u,v)=exchange(u,v,1)
    func=(u -> IndividualDisplacements.update_location_dpdo!(u,g))

    𝐹_MeshArray2D{eltype(u)}(u,u,v,v,[0,t1],func)
end

"""
    postprocess_MeshArray(sol,𝑃::FlowFields; id=missing, 𝑇=missing)

Copy `sol` to a `DataFrame` & map position to lon,lat coordinates
using "exchanged" 𝐷.XC, 𝐷.YC via `add_lonlat!`
"""
function postprocess_MeshArray(sol::ODESolution,𝑃::FlowFields; id=missing, 𝑇=missing)
    ismissing(id) ? id=collect(1:size(sol,2)) : nothing
    ismissing(𝑇) ? 𝑇=𝑃.𝑇 : nothing

    nd=length(size(sol))
    nt=size(sol,nd)
    nf=size(sol,nd-1)
    id=id*ones(1,size(sol,nd))
    if (size(sol,1)>1)&&(nd>2)
        x=sol[1,:,:]
        y=sol[2,:,:]
        fIndex=sol[end,:,:]
    elseif (nd>2)
        x=[sol[1,i,j][1] for i in 1:nf, j in 1:nt]
        y=[sol[1,i,j][2] for i in 1:nf, j in 1:nt]
        fIndex=[sol[1,i,j][end] for i in 1:nf, j in 1:nt]
    else
        x=sol[1,:]
        y=sol[2,:]
        fIndex=sol[end,:]
        nf=1
    end

    𝑃.u0.grid.nFaces==1 ? fIndex=ones(size(x)) : nothing

    t=[ceil(i/nf)-1 for i in 1:nt*nf]
    t=𝑇[1] .+ (𝑇[2]-𝑇[1])/t[end].*t
    
    df = DataFrame(ID=Int.(id[:]), x=x[:], y=y[:], fid=Int.(fIndex[:]), t=t[:])
    return df
end

"""
    add_lonlat!(df::DataFrame,XC,YC)

Add lon & lat to dataframe using "exchanged" XC, YC
"""
function add_lonlat!(df::DataFrame,XC,YC)
    x=df[!,:x];
    y=df[!,:y];
    f=Int.(df[!,:fid]);
    dx,dy=(x - floor.(x) .+ 0.5,y - floor.(y) .+ 0.5);
    i_c = Int32.(floor.(x)) .+ 1;
    j_c = Int32.(floor.(y)) .+ 1;

    lon=zeros(length(x),4)
    [lon[k,:]=XC[f[k]][i_c[k]:i_c[k]+1,j_c[k]:j_c[k]+1][:] for k in 1:length(i_c)]
    lat=zeros(length(x),4)
    [lat[k,:]=YC[f[k]][i_c[k]:i_c[k]+1,j_c[k]:j_c[k]+1][:] for k in 1:length(i_c)]

    k=findall(vec(maximum(lon,dims=2)-minimum(lon,dims=2)) .> 180.0)
    tmp=view(lon,k,:)
    tmp[findall(tmp.<0.0)]=tmp[findall(tmp.<0.0)] .+ 360.0

    df.lon=(1.0 .-dx).*(1.0 .-dy).*lon[:,1]+dx.*(1.0 .-dy).*lon[:,2] +
         (1.0 .-dx).*dy.*lon[:,3]+dx.*dy.*lon[:,4]
    df.lat=(1.0 .-dx).*(1.0 .-dy).*lat[:,1]+dx.*(1.0 .-dy).*lat[:,2] +
         (1.0 .-dx).*dy.*lat[:,3]+dx.*dy.*lat[:,4]

    return df
end

"""
    postprocess_xy()

Copy `sol` to a `DataFrame` & map position to x,y coordinates,
and define time axis for a simple doubly periodic domain
"""
function postprocess_xy(sol,𝑃::FlowFields; id=missing, 𝑇=missing)
    ismissing(id) ? id=collect(1:size(sol,2)) : nothing
    ismissing(𝑇) ? 𝑇=𝑃.𝑇 : nothing

    nf=size(sol,2)
    nt=size(sol,3)

    isa(𝑃.u0,MeshArray) ? (nx,ny)=𝑃.u0.grid.ioSize[1:2] : (nx,ny)=size(𝑃.u0)[1:2]
    nd=length(size(sol))

    id=id*ones(1,size(sol,nd))
    if (size(sol,1)>1)&&(nd>2)
        x=mod.(sol[1,:,:],Ref(nx))
        y=mod.(sol[2,:,:],Ref(ny))
    elseif (nd>2)
        x=[mod(sol[1,i,j][1],nx) for i in 1:nf, j in 1:nt]
        y=[mod(sol[1,i,j][2],ny) for i in 1:nf, j in 1:nt]
    else
        x=mod.(sol[1,:],Ref(nx))
        y=mod.(sol[2,:],Ref(ny))
    end
    t=[ceil(i/nf)-1 for i in 1:nt*nf]
    #size(𝐷.XC,1)>1 ? fIndex=sol[3,:,:] : fIndex=fill(1.0,size(x))
    t=𝑇[1] .+ (𝑇[2]-𝑇[1])/t[end].*t

    return DataFrame(ID=Int.(id[:]), t=t[:], x=x[:], y=y[:])
end

"""
    randn_lonlat(nn=1,seed=missing)

Randomly distributed longitude, latitude positions on the sphere.
"""
function randn_lonlat(nn=1;seed=missing)
    !ismissing(seed) ? rng = MersenneTwister(1234) : rng = MersenneTwister()
    tmp = randn(rng, Float64, (nn, 3))
    tmpn = tmp ./ sqrt.(sum(tmp.^2, dims=2))
    lon = rad2deg.(atan.(tmpn[:,2], tmpn[:,1]))
    lat = 90.0 .- rad2deg.(acos.(tmpn[:,3]))
    return lon, lat
end

"""
    nearest_to_xy(α::MeshArray,x,y,f)

Value of α at eachindex of the grid cell center nearest to `x,y` on subdomain array / facet `f`
"""
nearest_to_xy(α::MeshArray,x,y,f) = [α[Int(f[i]),1][ Int(round(x[i] .+ 0.5)), Int(round(y[i] .+ 0.5)) ] for i in eachindex(x)]

"""
    nearest_to_xy(α::Array,x,y)

Value of α at eachindex of the grid cell center nearest to `x,y`
"""
nearest_to_xy(α::Array,x,y) = [α[ Int(round(x[i] .+ 0.5)), Int(round(y[i] .+ 0.5)) ] for i in eachindex(x)]

"""
    interp_to_lonlat

Use MeshArrays.Interpolate() to interpolate to e.g. a regular grid (e.g. maps for plotting purposes).

```jldoctest
using IndividualDisplacements
p=dirname(pathof(IndividualDisplacements))
include(joinpath(p,"../examples/helper_functions.jl"))
𝑃,𝐷=global_ocean_circulation(k=1,ny=2);

lon=[i for i=20.:20.0:380., j=-70.:10.0:70.]
lat=[j for i=20.:20.0:380., j=-70.:10.0:70.]
(f,i,j,w,_,_,_)=InterpolationFactors(𝐷.Γ,vec(lon),vec(lat))
IntFac=(lon=lon,lat=lat,f=f,i=i,j=j,w=w)

tmp1=interp_to_lonlat(𝐷.Γ["Depth"],𝐷.Γ,lon,lat)
tmp2=interp_to_lonlat(𝐷.Γ["Depth"],IntFac)

ref=[5896. 5896.]
prod(isapprox.([maximum(tmp1) maximum(tmp2)],ref,atol=1.0))

# output

true
```
"""
function interp_to_lonlat(X::MeshArray,Γ::Dict,lon,lat)
    (f,i,j,w,_,_,_)=InterpolationFactors(Γ,vec(lon),vec(lat))
    return reshape(Interpolate(X,f,i,j,w),size(lon))
end

function interp_to_lonlat(X::MeshArray,IntFac::NamedTuple)
    @unpack f,i,j,w,lon,lat = IntFac
    return reshape(Interpolate(X,f,i,j,w),size(lon))
end

"""
    interp_to_xy(df::DataFrame,Zin::MeshArray)

Interpolate "exchanged" / "hallo-included" Zin to df[!,:x], df[!,:y] on df[!,:fid]
"""
function interp_to_xy(df::DataFrame,Zin::MeshArray)
    x=df[!,:x];
    y=df[!,:y];
    f=Int.(df[!,:fid]);
    dx,dy=(x - floor.(x) .+ 0.5,y - floor.(y) .+ 0.5);
    i_c = Int32.(floor.(x)) .+ 1;
    j_c = Int32.(floor.(y)) .+ 1;

    Z=zeros(length(x),4)
    [Z[k,:]=Zin[f[k]][i_c[k]:i_c[k]+1,j_c[k]:j_c[k]+1][:] for k in 1:length(i_c)]

    return (1.0 .-dx).*(1.0 .-dy).*Z[:,1]+dx.*(1.0 .-dy).*Z[:,2] +
           (1.0 .-dx).*dy.*Z[:,3]+dx.*dy.*Z[:,4]
end
