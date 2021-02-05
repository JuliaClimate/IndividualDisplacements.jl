#needed within OrdinaryDiffEq somehow:
import Base: zero
zero(tmp::Array) = zero.(tmp)

#needed to avoid allocations:
Flo=Union{Float32,Float64}
mydt(tim::Flo,𝑇::Array{Float32,1})=(tim-𝑇[1])/(𝑇[2]-𝑇[1])
mydt(tim::Flo,𝑇::Array{Float64,1})=(tim-𝑇[1])/(𝑇[2]-𝑇[1])

#needed to avoid allocations:
fSize(f::Array{Tuple{Int,Int}},i::Int) = (f[i][1],f[i][2])
fview(f::Array{Array{Float32,2},1},i::Int) = view(f[i],:,:)
fview(f::Array{Array{Float64,2},1},i::Int) = view(f[i],:,:)
fview(f::Array{Array{Float32,2},2},i::Int,j::Int) = view(f[i,j],:,:)
fview(f::Array{Array{Float64,2},2},i::Int,j::Int) = view(f[i,j],:,:)

"""
    dxyz_dt!(du,u,p::𝐹_MeshArray3D,tim)

Interpolate velocity from gridded fields (3D; with halos) to position `u`
(`x,y,z,fIndex`) to compute the derivative of position v time  `du_dt`.

```jldoctest
using IndividualDisplacements, Statistics
p=dirname(pathof(IndividualDisplacements))
include(joinpath(p,"../examples/worldwide/three_dimensional_ocean.jl"))
ref=[211. 34. -70.]
prod(isapprox.([mean(𝐼.🔴.lon) mean(𝐼.🔴.lat) mean(𝐼.🔴.z)],ref,atol=50.0))

# output

true
```
"""
function dxyz_dt!(du::Array{T,1},u::Array{T,1},𝑃::𝐹_MeshArray3D,tim) where T
    dt=mydt(tim,𝑃.𝑇)
    g=𝑃.u0.grid
    #
    while location_is_out(u,g)
        𝑃.update_location!(u)
    end

    x=u[1]
    y=u[2]
    z=u[3]
    fIndex = Int(u[4])
    (nx,ny)=fSize(g.fSize,fIndex)
    nz=size(𝑃.u0,2)
    #
    (dx,dy,dz)=(x - floor(x),y - floor(y),z - floor(z))
    i_c = Int32(floor(x)) + 2
    j_c = Int32(floor(y)) + 2
    k_c = Int32(floor(z)) + 1
    #
    (i_w,i_e)=(i_c,i_c+1)
    (j_s,j_n)=(j_c,j_c+1)
    (k_l,k_r)=(k_c,k_c+1)
    #
    k_c=min(k_c,nz)
    k_l=min(k_l,nz)
    k_r=min(k_r,nz)
    k_c=max(k_c,1)
    k_l=max(k_l,1)
    k_r=max(k_r,1)
    #
    onemdt=(1.0-dt)
    onemdx=(1.0-dx)
    onemdy=(1.0-dy)
    onemdz=(1.0-dz)
    #interpolate u to position and time
    u0=fview(𝑃.u0.f,fIndex,k_c)
    u1=fview(𝑃.u1.f,fIndex,k_c)
    du[1]=onemdx*onemdt*u0[i_w,j_c]+
    dx*onemdt*u0[i_e,j_c]+
    onemdx*dt*u1[i_w,j_c]+
    dx*dt*u1[i_e,j_c]
    #interpolate v to position and time
    v0=fview(𝑃.v0.f,fIndex,k_c)
    v1=fview(𝑃.v1.f,fIndex,k_c)
    du[2]=onemdy*onemdt*v0[i_c,j_s]+
    dy*onemdt*v0[i_c,j_n]+
    onemdy*dt*v1[i_c,j_s]+
    dy*dt*v1[i_c,j_n]
    #interpolate w to position and time
    du[3]=onemdz*onemdt*𝑃.w0.f[fIndex,k_l][i_c,j_c]+
    dz*onemdt*𝑃.w0.f[fIndex,k_r][i_c,j_c]+
    onemdz*dt*𝑃.w1.f[fIndex,k_l][i_c,j_c]+
    dz*dt*𝑃.w1.f[fIndex,k_r][i_c,j_c]
    #leave face index unchanged
    du[4]=0.0
    #
    return du
end

function dxyz_dt!(du::Array{T,2},u::Array{T,2},𝑃::𝐹_MeshArray3D,tim) where T
    [dxyz_dt!(du[i],u[i],𝑃,tim) for i=1:size(u,2)]
end

"""
    dxy_dt!(du,u,p::𝐹_MeshArray2D,tim)

Interpolate velocity from gridded fields (2D; with halos) to position `u`
(`x,y,fIndex`) to compute the derivative of position v time  `du_dt`.

```jldoctest
using IndividualDisplacements, Statistics
p=dirname(pathof(IndividualDisplacements))
include(joinpath(p,"../examples/flow_fields.jl"));

u,v,ϕ=random_flow_field()
#𝐹=𝐹_Array2D(u,u,v,v,[0.,10.])
𝐹=convert_to_FlowFields(u,v,10.0)

np,nq=size(u)
x=np*(0.4 .+ 0.2*rand(100))
y=nq*(0.4 .+ 0.2*rand(100))

a=ones(size(x))
𝐼=Individuals(𝐹,x,y,a)
∫!(𝐼)

ref=size(u) ./2
prod(isapprox.([mean(𝐼.🔴.x) mean(𝐼.🔴.y)],ref,atol=10.0))

# output

true
```
"""
function dxy_dt!(du::Array{T,1},u::Array{T,1},𝑃::𝐹_MeshArray2D,tim) where T
    #compute positions in index units
    dt=mydt(tim,𝑃.𝑇)
    g=𝑃.u0.grid
    #
    while location_is_out(u,g)
        𝑃.update_location!(u)
    end

    x=u[1]
    y=u[2]
    fIndex = Int(u[3])
    (nx,ny)=fSize(g.fSize,fIndex)
    #
    (dx,dy)=(x - floor(x),y - floor(y))
    i_c = Int32(floor(x)) + 2
    j_c = Int32(floor(y)) + 2
    #
    (i_w,i_e)=(i_c,i_c+1)
    (j_s,j_n)=(j_c,j_c+1)
    #interpolate u to position and time
    u0=fview(𝑃.u0.f,fIndex)
    u1=fview(𝑃.u1.f,fIndex)
    du[1]=(1.0-dx)*(1.0-dt)*u0[i_w,j_c]+
    dx*(1.0-dt)*u0[i_e,j_c]+
    (1.0-dx)*dt*u1[i_w,j_c]+
    dx*dt*u1[i_e,j_c]
    #interpolate v to position and time
    v0=fview(𝑃.v0.f,fIndex)
    v1=fview(𝑃.v1.f,fIndex)
    du[2]=(1.0-dy)*(1.0-dt)*v0[i_c,j_s]+
    dy*(1.0-dt)*v0[i_c,j_n]+
    (1.0-dy)*dt*v1[i_c,j_s]+
    dy*dt*v1[i_c,j_n]
    #leave face index unchanged
    du[3]=0.0
    #
    return du
end

function dxy_dt!(du::Array{T,2},u::Array{T,2},𝑃::𝐹_MeshArray2D,tim) where T
    [dxy_dt!(du[i],u[i],𝑃,tim) for i=1:size(u,2)]
end

"""
    dxyz_dt(du,u,𝑃::𝐹_Array3D,tim)

Interpolate velocity from gridded fields (3D; NO halos) to position `u`
(`x,y,z`) to compute the derivative of position v time  `du_dt`.

```jldoctest
using IndividualDisplacements
p=dirname(pathof(IndividualDisplacements))
include(joinpath(p,"../examples/basics/solid_body_rotation.jl"))
ref=[7.767441577479032 9.513402495574852 0.7065855989421701]
prod(isapprox.(𝐼.📌',ref,atol=1.0))

# output

true
```
"""
function dxyz_dt(du::Array{T,1},u::Array{T,1},𝑃::𝐹_Array3D,tim) where T
    #compute positions in index units
    dt=mydt(tim,𝑃.𝑇)
    #
    (nx,ny,nz) = size(𝑃.u0)
    x=mod(u[1],nx)
    y=mod(u[2],ny)
    z=mod(u[3],nz)
    #
    dx=x - floor(x)
    dy=y - floor(y)
    dz=z - floor(z)
    #
    i_c = Int32(floor(x)) + 1
    j_c = Int32(floor(y)) + 1
    k_c = Int32(floor(z)) + 1
    #
    i_c==0 ? i_c=nx : nothing
    j_c==0 ? j_c=ny : nothing
    #k_c==0 ? k_c=nz : nothing
    #
    (i_w,i_e)=(i_c,i_c+1)
    x>=nx-1 ? (i_w,i_e)=(nx,1) : nothing
    (j_s,j_n)=(j_c,j_c+1)
    y>=ny-1 ? (j_s,j_n)=(ny,1) : nothing
    (k_l,k_r)=(k_c,k_c+1)
    #z>=nz-1 ? (k_l,k_r)=(nz,1) : nothing

    #interpolate u to position and time
    du[1]=(1.0-dx)*(1.0-dt)*𝑃.u0[i_w,j_c,k_c]+
    dx*(1.0-dt)*𝑃.u0[i_e,j_c,k_c]+
    (1.0-dx)*dt*𝑃.u1[i_w,j_c,k_c]+
    dx*dt*𝑃.u1[i_e,j_c,k_c]
    #interpolate v to position and time
    du[2]=(1.0-dy)*(1.0-dt)*𝑃.v0[i_c,j_s,k_c]+
    dy*(1.0-dt)*𝑃.v0[i_c,j_n,k_c]+
    (1.0-dy)*dt*𝑃.v1[i_c,j_s,k_c]+
    dy*dt*𝑃.v1[i_c,j_n,k_c]
    #interpolate w to position and time
    du[3]=(1.0-dz)*(1.0-dt)*𝑃.w0[i_c,j_c,k_l]+
    dz*(1.0-dt)*𝑃.w0[i_c,j_c,k_r]+
    (1.0-dz)*dt*𝑃.w1[i_c,j_c,k_l]+
    dz*dt*𝑃.w1[i_c,j_c,k_r]
    #
    return du
end

function dxyz_dt(du::Array{T,2},u::Array{T,2},𝑃::𝐹_Array3D,tim) where T
    [dxyz_dt(du[i],u[i],𝑃,tim) for i=1:size(u,2)]
end

"""
    dxy_dt(du,u,𝑃::𝐹_Array2D,tim)

Interpolate velocity from gridded fields (2D; NO halos) to position `u`
(`x,y`) to compute the derivative of position v time  `du_dt`.

```jldoctest
using IndividualDisplacements, Statistics
p=dirname(pathof(IndividualDisplacements))
include(joinpath(p,"../examples/basics/particle_cloud.jl"))
ref=[29.381183342468674  19.890831699436823]
prod(isapprox.([mean(𝐼.🔴.x) mean(𝐼.🔴.y)],ref,atol=1.0))

# output

true
```
"""
function dxy_dt(du::Array{T,1},u::Array{T,1},𝑃::𝐹_Array2D,tim) where T
    dt=mydt(tim,𝑃.𝑇)
    #
    (nx,ny) = size(𝑃.u0)
    x=mod(u[1],nx)
    y=mod(u[2],ny)
    #
    dx=x - floor(x)
    dy=y - floor(y)
    #
    i_c = Int32(floor(x)) + 1
    j_c = Int32(floor(y)) + 1
    #
    i_c==0 ? i_c=nx : nothing
    j_c==0 ? j_c=ny : nothing
    #
    (i_w,i_e)=(i_c,i_c+1)
    x>=nx-1 ? (i_w,i_e)=(nx,1) : nothing
    (j_s,j_n)=(j_c,j_c+1)
    y>=ny-1 ? (j_s,j_n)=(ny,1) : nothing
    #interpolate u to position and time
    du[1]=(1.0-dx)*(1.0-dt)*𝑃.u0[i_w,j_c]+
    dx*(1.0-dt)*𝑃.u0[i_e,j_c]+
    (1.0-dx)*dt*𝑃.u1[i_w,j_c]+
    dx*dt*𝑃.u1[i_e,j_c]
    #interpolate v to position and time
    du[2]=(1.0-dy)*(1.0-dt)*𝑃.v0[i_c,j_s]+
    dy*(1.0-dt)*𝑃.v0[i_c,j_n]+
    (1.0-dy)*dt*𝑃.v1[i_c,j_s]+
    dy*dt*𝑃.v1[i_c,j_n]
    #
    return du
end

function dxy_dt(du::Array{T,2},u::Array{T,2},𝑃::𝐹_Array2D,tim) where T
    [dxy_dt(du[i],u[i],𝑃,tim) for i=1:size(u,2)]
end

"""
    dxy_dt_CyclicArray(du,u,𝑃::NamedTuple,tim)

Nearest neighbor (?) velocity from gridded fields (2D; NO halos but
not needed when CyclicArrays is used to extend valid indice ranges).

_notes:_ spatial interpolation & temporal interpolation are lacking
"""
function dxy_dt_CyclicArray(du::Array{T,2},u::Array{T,2},𝑃::NamedTuple,tim) where T
    np=size(du,2)
    xi,yi=(u[1,:],u[2,:])

    #not needed:
    #xi=floor.(𝑃.xg[1,Int.(sign.(xi).*floor.(abs.(xi)))])+rem.(xi,1)
    #yi=floor.(𝑃.yg[Int.(sign.(yi).*floor.(abs.(yi))),1])+rem.(yi,1)

    i=Int.(floor.(𝑃.xg[1,Int.(floor.(xi))]))
    j=Int.(floor.(𝑃.yg[Int.(floor.(yi)),1]))
    du[1,:]=[𝑃.u[i[ii],j[ii]] for ii in 1:np]
    du[2,:]=[𝑃.v[i[ii],j[ii]] for ii in 1:np]
    return du
end

"""
    dict_to_nt(tmp::Dict)

Attempt to convert from Dict to NamedTuple
"""
dict_to_nt(tmp::Dict) = tmp=(; zip(Symbol.(keys(tmp)), values(tmp))...)

dxy_dt!(du,u,𝑃::Dict,tim) = dxy_dt!(du,u,dict_to_nt(𝑃),tim)
dxy_dt(du,u,𝑃::Dict,tim) = dxy_dt(du,u,dict_to_nt(𝑃),tim)
dxyz_dt(du,u,𝑃::Dict,tim) = dxyz_dt(du,u,dict_to_nt(𝑃),tim)

"""
    dxy_dt_replay(du,u,p::DataFrame,t)

Interpolate velocity from MITgcm float_trajectories output and return
position increment `du`.
"""
function dxy_dt_replay(du,u,p::DataFrame,t)
    tt=t/3600.0
    tt0=Int32(floor(tt))
    w=tt-tt0
    du[1]=(1.0-w)*p[tt0+1,:uVel]+w*p[tt0+2,:uVel]
    du[2]=(1.0-w)*p[tt0+1,:vVel]+w*p[tt0+2,:vVel]
end
