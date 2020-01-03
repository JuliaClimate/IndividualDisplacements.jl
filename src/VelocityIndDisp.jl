
using DataFrames

"""
    NeighborTileIndices_dpdo(ni::Int,nj::Int)

List of W, E, S, N neighbor tile IDs in the case of a doubly
periodic domain with ni x nj tiles.
"""
function NeighborTileIndices_dpdo(ni::Int,nj::Int)
    tmp=fill(0,ni*nj,4)
    for i=1:ni
        for j=1:nj
            k=i+ni*(j-1)
            kS=j-1; kS==0 ? kS=nj : nothing; kS=i+ni*(kS-1)
            kN=j+1; kN==nj+1 ? kN=1 : nothing; kN=i+ni*(kN-1)
            kW=i-1; kW==0 ? kW=ni : nothing; kW=kW+ni*(j-1)
            kE=i+1; kE==ni+1 ? kE=1 : nothing; kE=kE+ni*(j-1)
            tmp[k,1]=kW
            tmp[k,2]=kE
            tmp[k,3]=kS
            tmp[k,4]=kN
        end
    end
    return tmp
end

"""
    UpdateLocation!

Update location (x,y,fIndex) when out of domain. Note: initially, this
only works for the `dpdo` grid type provided by `MeshArrays.jl`.
"""
function UpdateLocation!(u::Array{Float64,1},grid::gcmgrid)
    x,y = u[1:2]
    fIndex = Int(u[3])
    #
    nx,ny=grid.fSize[fIndex]
    ni,nj=Int.(transpose(grid.ioSize)./grid.fSize[1])
    WESN=NeighborTileIndices_dpdo(ni,nj)
    #
    if x<0
        x=x+nx
        u[1]=x
        fIndex=WESN[fIndex,1]
        u[3]=fIndex
    elseif x>=nx
        x=x-nx
        u[1]=x
        fIndex=WESN[fIndex,2]
        u[3]=fIndex
    end
    #
    if y<0
        y=y+ny
        u[2]=y
        fIndex=WESN[fIndex,3]
        u[3]=fIndex
    elseif y>=ny
        y=y-ny
        u[2]=y
        fIndex=WESN[fIndex,4]
        u[3]=fIndex
    end
    #
    return u
end


"""
    VelComp!(du,u,p::Dict,tim)

Interpolate velocity from gridded fields and return position increment `du`
"""
function VelComp!(du::Array{Float64,1},u::Array{Float64,1},p::Dict,tim)
    #compute positions in index units
    dt=(tim-p["t0"])/(p["t1"]-p["t0"])
    #
    UpdateLocation!(u,p["u0"].grid)
    x,y = u[1:2]
    fIndex = Int(u[3])
    nx,ny=p["u0"].grid.fSize[fIndex]
    #debugging stuff
    if (false & (mod(x,nx)!=x)|(mod(y,ny)!=y))
        println("crossing domain edge"*"$x and $y")
    end
    #
    dx,dy=[x - floor(x),y - floor(y)]
    i_c,j_c = Int32.(floor.([x y])) .+ 2
    #
    i_w,i_e=[i_c i_c+1]
    j_s,j_n=[j_c j_c+1]
    #interpolate u to position and time
    du[1]=(1.0-dx)*(1.0-dt)*p["u0"].f[fIndex][i_w,j_c]+
    dx*(1.0-dt)*p["u0"].f[fIndex][i_e,j_c]+
    (1.0-dx)*dt*p["u1"].f[fIndex][i_w,j_c]+
    dx*dt*p["u1"].f[fIndex][i_e,j_c]
    #interpolate v to position and time
    du[2]=(1.0-dy)*(1.0-dt)*p["v0"].f[fIndex][i_c,j_s]+
    dy*(1.0-dt)*p["v0"].f[fIndex][i_c,j_n]+
    (1.0-dy)*dt*p["v1"].f[fIndex][i_c,j_s]+
    dy*dt*p["v1"].f[fIndex][i_c,j_n]
    #leave face index unchanged
    du[3]=0.0
    #
    return du
end

function VelComp!(du::Array{Float64,2},u::Array{Float64,2},p::Dict,tim)
    for i=1:size(u,2)
        tmpdu=du[1:3,i]
        tmpu=u[1:3,i]
        VelComp!(tmpdu,tmpu,p,tim)
        du[1:3,i]=tmpdu
        u[1:3,i]=tmpu
    end
    return du
end

"""
    VelComp(du,u,p::Dict,tim)

Interpolate velocity from gridded fields and return position increment `du`
"""
function VelComp(du::Array{Float64,1},u::Array{Float64,1},p::Dict,tim)
    #compute positions in index units
    dt=(tim-p["t0"])/(p["t1"]-p["t0"])
    #
    x,y = u[1:2]
    nx,ny=p["u0"].grid.ioSize
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
    #debugging stuff
    if false
        println((x,y,i_c,j_c))
        println((i_w,i_e,j_s,j_n))
        println((dx,dy,dt))
        #println(du)
    end
    #interpolate u to position and time
    du[1]=(1.0-dx)*(1.0-dt)*p["u0"].f[1][i_w,j_c]+
    dx*(1.0-dt)*p["u0"].f[1][i_e,j_c]+
    (1.0-dx)*dt*p["u1"].f[1][i_w,j_c]+
    dx*dt*p["u1"].f[1][i_e,j_c]
    #interpolate v to position and time
    du[2]=(1.0-dy)*(1.0-dt)*p["v0"].f[1][i_c,j_s]+
    dy*(1.0-dt)*p["v0"].f[1][i_c,j_n]+
    (1.0-dy)*dt*p["v1"].f[1][i_c,j_s]+
    dy*dt*p["v1"].f[1][i_c,j_n]
    #
    return du
end

function VelComp(du::Array{Float64,2},u::Array{Float64,2},p::Dict,tim)
    for i=1:size(u,2)
        tmpdu=du[1:2,i]
        VelComp(tmpdu,u[1:2,i],p,tim)
        du[1:2,i]=tmpdu
    end
    return du
end

# VelCopy.jl

"""
    VelCopy(du,u,p::DataFrame,t)

Interpolate velocity from MITgcm float_trajectories output and return
position increment `du`.
"""
function VelCopy(du,u,p::DataFrame,t)
    tt=t/3600.0
    tt0=Int32(floor(tt))
    w=tt-tt0
    du[1]=(1.0-w)*p[tt0+1,:uVel]+w*p[tt0+2,:uVel]
    du[2]=(1.0-w)*p[tt0+1,:vVel]+w*p[tt0+2,:vVel]
end
