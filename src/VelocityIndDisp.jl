
using DataFrames

#VelCopy(du,uInit,tmp,3600.0)
#VelComp(du,uInit,uvetc,3600.0)


"""
    getNeighborList(ni::Int,nj::Int)

Get list of E, W, S, N neighbors in the doubly perdiodic domain case
"""
function getNeighborList(ni::Int,nj::Int)

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
    VelComp!(du,u,p::Dict,tim)

Interpolate velocity from gridded fields and return position increment `du`
"""
function VelComp!(du::Array{Float64,1},u::Array{Float64,1},p::Dict,tim)
    #compute positions in index units
    dt=(tim-p["t0"])/(p["t1"]-p["t0"])
    #
    x,y = u[1:2]
    fIndex = Int(u[3])
    nx,ny=p["u0"].grid.fSize[fIndex]
    #
    ni,nj=Int.(transpose(p["u0"].grid.ioSize)./p["u0"].grid.fSize[1])
    NeighborList=getNeighborList(ni,nj)
    #debugging stuff
    if false
        println(("a",x,y,fIndex))
    end
    #
    if x<0
        x=x+nx
        u[1]=x
        fIndex=NeighborList[fIndex,1]
        u[3]=fIndex
    elseif x>=nx
        x=x-nx
        u[1]=x
        fIndex=NeighborList[fIndex,2]
        u[3]=fIndex
    end
    #
    if y<0
        y=y+ny
        u[2]=y
        fIndex=NeighborList[fIndex,3]
        u[3]=fIndex
    elseif y>=ny
        y=y-ny
        u[2]=y
        fIndex=NeighborList[fIndex,4]
        u[3]=fIndex
    end
    #debugging stuff
    if false
        println(("b",x,y,fIndex))
    end
    #
    if (mod(x,nx)!=x) | (mod(y,ny)!=y)
        println("crossing domain edge"*"$x and $y")
    end
    #
    dx,dy=[x - floor(x),y - floor(y)]
    i_c,j_c = Int32.(floor.([x y])) .+ 2
    #
    i_w,i_e=[i_c i_c+1]
    j_s,j_n=[j_c j_c+1]
    #debugging stuff
    if i_w>22|i_c>22|i_e>22|j_s>22|j_c>22|j_n>22
        println((x,y,fIndex,i_c,j_c))
        println((i_w,i_e,j_s,j_n))
        println((dx,dy,dt))
        #println(du)
    end
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