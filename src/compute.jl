"""
    dxyz_dt!(du,u,p::NamedTuple,tim)

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
function dxyz_dt!(du::Array{Float64,1},u::Array{Float64,1},𝑃::NamedTuple,tim)
    #compute positions in index units
    dt=(tim-𝑃.𝑇[1])/(𝑃.𝑇[2]-𝑃.𝑇[1])
    dt>1.0 ? error("dt>1.0") : nothing
    dt<0.0 ? error("dt>0.0") : nothing
    g=𝑃.u0.grid
    #
    g.class=="PeriodicDomain" ? update_location_dpdo!(u,g) : nothing
    g.class=="CubeSphere" ? update_location_cs!(u,𝑃) : nothing
    g.class=="LatLonCap" ? update_location_cs!(u,𝑃) : nothing

    x,y,z = u[1:3]
    fIndex = Int(u[4])
    nx,ny=g.fSize[fIndex]
    nz=size(𝑃.u0,2)
    #
    dx,dy,dz=[x - floor(x),y - floor(y),z - floor(z)]
    i_c,j_c = Int32.(floor.([x y])) .+ 2
    k_c = Int32.(floor.(z)) .+ 1
    #
    i_w,i_e=[i_c i_c+1]
    j_s,j_n=[j_c j_c+1]
    k_l,k_r=[k_c k_c+1]
    #interpolate u to position and time
    du[1]=(1.0-dx)*(1.0-dt)*𝑃.u0.f[fIndex,k_c][i_w,j_c]+
    dx*(1.0-dt)*𝑃.u0.f[fIndex,k_c][i_e,j_c]+
    (1.0-dx)*dt*𝑃.u1.f[fIndex,k_c][i_w,j_c]+
    dx*dt*𝑃.u1.f[fIndex,k_c][i_e,j_c]
    #interpolate v to position and time
    du[2]=(1.0-dy)*(1.0-dt)*𝑃.v0.f[fIndex,k_c][i_c,j_s]+
    dy*(1.0-dt)*𝑃.v0.f[fIndex,k_c][i_c,j_n]+
    (1.0-dy)*dt*𝑃.v1.f[fIndex,k_c][i_c,j_s]+
    dy*dt*𝑃.v1.f[fIndex,k_c][i_c,j_n]
    #interpolate w to position and time
    du[3]=(1.0-dz)*(1.0-dt)*𝑃.w0.f[fIndex,k_l][i_c,j_c]+
    dz*(1.0-dt)*𝑃.w0.f[fIndex,k_r][i_c,j_c]+
    (1.0-dz)*dt*𝑃.w1.f[fIndex,k_l][i_c,j_c]+
    dz*dt*𝑃.w1.f[fIndex,k_r][i_c,j_c]
    
    #leave face index unchanged
    du[4]=0.0
    #
    return du
end

function dxyz_dt!(du::Array{Float64,2},u::Array{Float64,2},𝑃::NamedTuple,tim)
    for i=1:size(u,2)
        tmpdu=du[1:4,i]
        tmpu=u[1:4,i]
        dxyz_dt!(tmpdu,tmpu,𝑃,tim)
        du[1:4,i]=tmpdu
        u[1:4,i]=tmpu
    end
    return du
end

"""
    dxy_dt!(du,u,p::NamedTuple,tim)

Interpolate velocity from gridded fields (2D; with halos) to position `u`
(`x,y,fIndex`) to compute the derivative of position v time  `du_dt`.

```jldoctest
using IndividualDisplacements, Statistics
p=dirname(pathof(IndividualDisplacements))
include(joinpath(p,"../examples/basics/random_flow_field.jl"))
ref=[4. 6.]
prod(isapprox.([mean(𝐼.🔴.x) mean(𝐼.🔴.y)],ref,atol=10.0))

# output

true
```

```jldoctest
using IndividualDisplacements, Statistics
p=dirname(pathof(IndividualDisplacements))
include(joinpath(p,"../examples/worldwide/global_ocean_circulation.jl"))
ref=[78. 88.]
prod(isapprox.([mean(𝐼.🔴.x) mean(𝐼.🔴.y)],ref,atol=10.0))

# output

true
```
"""
function dxy_dt!(du::Array{Float64,1},u::Array{Float64,1},𝑃::NamedTuple,tim)
    #compute positions in index units
    dt=(tim-𝑃.𝑇[1])/(𝑃.𝑇[2]-𝑃.𝑇[1])
    dt>1.0 ? error("dt>1.0") : nothing
    dt<0.0 ? error("dt>0.0") : nothing
    g=𝑃.u0.grid
    #
    g.class=="PeriodicDomain" ? update_location_dpdo!(u,g) : nothing
    g.class=="CubeSphere" ? update_location_cs!(u,𝑃) : nothing
    g.class=="LatLonCap" ? update_location_cs!(u,𝑃) : nothing

    x,y = u[1:2]
    fIndex = Int(u[3])
    nx,ny=g.fSize[fIndex]
    #
    dx,dy=[x - floor(x),y - floor(y)]
    i_c,j_c = Int32.(floor.([x y])) .+ 2
    #
    i_w,i_e=[i_c i_c+1]
    j_s,j_n=[j_c j_c+1]
    #interpolate u to position and time
    du[1]=(1.0-dx)*(1.0-dt)*𝑃.u0.f[fIndex][i_w,j_c]+
    dx*(1.0-dt)*𝑃.u0.f[fIndex][i_e,j_c]+
    (1.0-dx)*dt*𝑃.u1.f[fIndex][i_w,j_c]+
    dx*dt*𝑃.u1.f[fIndex][i_e,j_c]
    #interpolate v to position and time
    du[2]=(1.0-dy)*(1.0-dt)*𝑃.v0.f[fIndex][i_c,j_s]+
    dy*(1.0-dt)*𝑃.v0.f[fIndex][i_c,j_n]+
    (1.0-dy)*dt*𝑃.v1.f[fIndex][i_c,j_s]+
    dy*dt*𝑃.v1.f[fIndex][i_c,j_n]
    #leave face index unchanged
    du[3]=0.0
    #
    return du
end

function dxy_dt!(du::Array{Float64,2},u::Array{Float64,2},𝑃::NamedTuple,tim)
    for i=1:size(u,2)
        tmpdu=du[1:3,i]
        tmpu=u[1:3,i]
        dxy_dt!(tmpdu,tmpu,𝑃,tim)
        du[1:3,i]=tmpdu
        u[1:3,i]=tmpu
    end
    return du
end

"""
    dxyz_dt(du,u,𝑃::NamedTuple,tim)

Interpolate velocity from gridded fields (3D; NO halos) to position `u`
(`x,y,z`) to compute the derivative of position v time  `du_dt`.

```jldoctest
using IndividualDisplacements
p=dirname(pathof(IndividualDisplacements))
include(joinpath(p,"../examples/basics/solid_body_rotation.jl"))
ref=[7.767441577479032 9.513402495574852 0.7065855989421701]
prod(isapprox.(𝐼.📌',ref))

# output

true
```
"""
function dxyz_dt(du::Array{Float64,1},u::Array{Float64,1},𝑃::NamedTuple,tim)
    #compute positions in index units
    dt=(tim-𝑃.𝑇[1])/(𝑃.𝑇[2]-𝑃.𝑇[1])
    #
    x,y,z = u[1:3]
    nx,ny,nz = 𝑃.ioSize
    x,y=[mod(x,nx),mod(y,ny)]
    z=mod(z,nz)
    #
    dx,dy,dz=[x - floor(x),y - floor(y),z - floor(z)]
    i_c,j_c,k_c = Int32.(floor.([x y z])) .+ 1
    i_c==0 ? i_c=nx : nothing
    j_c==0 ? j_c=ny : nothing
    #k_c==0 ? k_c=nz : nothing
    #
    i_w,i_e=[i_c i_c+1]
    x>=nx-1 ? (i_w,i_e)=(nx,1) : nothing
    j_s,j_n=[j_c j_c+1]
    y>=ny-1 ? (j_s,j_n)=(ny,1) : nothing
    k_l,k_r=[k_c k_c+1]
    #z>=nz-1 ? (k_l,k_r)=(nz,1) : nothing

    #interpolate u to position and time
    du[1]=(1.0-dx)*(1.0-dt)*𝑃.u0.f[1,k_c][i_w,j_c]+
    dx*(1.0-dt)*𝑃.u0.f[1,k_c][i_e,j_c]+
    (1.0-dx)*dt*𝑃.u1.f[1,k_c][i_w,j_c]+
    dx*dt*𝑃.u1.f[1,k_c][i_e,j_c]
    #interpolate v to position and time
    du[2]=(1.0-dy)*(1.0-dt)*𝑃.v0.f[1,k_c][i_c,j_s]+
    dy*(1.0-dt)*𝑃.v0.f[1,k_c][i_c,j_n]+
    (1.0-dy)*dt*𝑃.v1.f[1,k_c][i_c,j_s]+
    dy*dt*𝑃.v1.f[1,k_c][i_c,j_n]
    #interpolate w to position and time
    du[3]=(1.0-dz)*(1.0-dt)*𝑃.w0.f[1,k_l][i_c,j_c]+
    dz*(1.0-dt)*𝑃.w0.f[1,k_r][i_c,j_c]+
    (1.0-dz)*dt*𝑃.w1.f[1,k_l][i_c,j_c]+
    dz*dt*𝑃.w1.f[1,k_r][i_c,j_c]
    #
    return du
end

function dxyz_dt(du::Array{Float64,2},u::Array{Float64,2},𝑃::NamedTuple,tim)
    for i=1:size(u,2)
        tmpdu=du[1:3,i]
        dxyz_dt(tmpdu,u[1:3,i],𝑃,tim)
        du[1:3,i]=tmpdu
    end
    return du
end

"""
    dxy_dt(du,u,𝑃::NamedTuple,tim)

Interpolate velocity from gridded fields (2D; NO halos) to position `u`
(`x,y`) to compute the derivative of position v time  `du_dt`.

```jldoctest
using IndividualDisplacements, Statistics
p=dirname(pathof(IndividualDisplacements))
include(joinpath(p,"../examples/basics/particle_cloud.jl"))
ref=[29.381183342468674  19.890831699436823]
prod(isapprox.([mean(𝐼.🔴.x) mean(𝐼.🔴.y)],ref))

# output

true
```
"""
function dxy_dt(du::Array{Float64,1},u::Array{Float64,1},𝑃::NamedTuple,tim)
    #compute positions in index units
    dt=(tim-𝑃.𝑇[1])/(𝑃.𝑇[2]-𝑃.𝑇[1])
    #
    x,y = u[1:2]
    nx,ny=𝑃.ioSize
    x,y=[mod(x,nx),mod(y,ny)]
    #
    dx,dy=[x - floor(x),y - floor(y)]
    i_c,j_c = Int32.(floor.([x y])) .+ 1
    i_c==0 ? i_c=nx : nothing
    j_c==0 ? j_c=ny : nothing
    #
    i_w,i_e=[i_c i_c+1]
    x>=nx-1 ? (i_w,i_e)=(nx,1) : nothing
    j_s,j_n=[j_c j_c+1]
    y>=ny-1 ? (j_s,j_n)=(ny,1) : nothing
    #interpolate u to position and time
    du[1]=(1.0-dx)*(1.0-dt)*𝑃.u0.f[1][i_w,j_c]+
    dx*(1.0-dt)*𝑃.u0.f[1][i_e,j_c]+
    (1.0-dx)*dt*𝑃.u1.f[1][i_w,j_c]+
    dx*dt*𝑃.u1.f[1][i_e,j_c]
    #interpolate v to position and time
    du[2]=(1.0-dy)*(1.0-dt)*𝑃.v0.f[1][i_c,j_s]+
    dy*(1.0-dt)*𝑃.v0.f[1][i_c,j_n]+
    (1.0-dy)*dt*𝑃.v1.f[1][i_c,j_s]+
    dy*dt*𝑃.v1.f[1][i_c,j_n]
    #
    return du
end

function dxy_dt(du::Array{Float64,2},u::Array{Float64,2},𝑃::NamedTuple,tim)
    for i=1:size(u,2)
        tmpdu=du[1:2,i]
        dxy_dt(tmpdu,u[1:2,i],𝑃,tim)
        du[1:2,i]=tmpdu
    end
    return du
end

"""
    dxy_dt_CyclicArray(du,u,𝑃::NamedTuple,tim)

Nearest neighbor (?) velocity from gridded fields (2D; NO halos but
not needed when CyclicArrays is used to extend valid indice ranges).

_notes:_ spatial interpolation & temporal interpolation are lacking
"""
function dxy_dt_CyclicArray(du::Array{Float64,2},u::Array{Float64,2},𝑃::NamedTuple,tim)
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
