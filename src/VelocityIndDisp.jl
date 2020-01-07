
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
    UpdateLocation_cs!

Update location (x,y,fIndex) when out of domain. Note: initially, this
only works for the `dpdo` grid type provided by `MeshArrays.jl`.
"""
function UpdateLocation_cs!(u::Array{Float64,1},grid::Dict)
    x,y = u[1:2]
    fIndex = Int(u[3])
    nx,ny=grid["XC"].fSize[fIndex]
    if x<0||x>nx||y<0||y>ny
        j = 0
        x<0 ? j=grid["aW"][fIndex] : nothing
        x>nx ? j=grid["aE"][fIndex] : nothing
        y<0 ? j=grid["aS"][fIndex] : nothing
        y>ny ? j=grid["aN"][fIndex] : nothing
        (x,y)=grid["RelocFunctions"][j,fIndex](x,y)
        u[1]=x
        u[2]=y
        u[3]=j
    end
    #
    return u
end

"""
    UpdateLocation_dpdo!

Update location (x,y,fIndex) when out of domain. Note: initially, this
only works for the `dpdo` grid type provided by `MeshArrays.jl`.
"""
function UpdateLocation_dpdo!(u::Array{Float64,1},grid::gcmgrid)
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
    NeighborTileIndices_cs(grid::Dict)

Derive list of neighboring tile indices for a cs or llc grid + functions that
convert indices from one tile to another. Returns a Dict to merge later.
"""
function NeighborTileIndices_cs(grid::Dict)
    s = grid["XC"].fSize
    nFaces = length(s)
    nFaces == 5 ? s = vcat(s, s[3]) : nothing
    aW=Array{Int,1}(undef,nFaces)
    aE=similar(aW); aS=similar(aW); aN=similar(aW);
    for i = 1:nFaces
        (aW[i], aE[i], aS[i], aN[i], _, _, _, _) = MeshArrays.exch_cs_sources(i, s, 1)
    end
    RelocFunctions=RelocationFunctions_cs(grid["XC"])
    return Dict("aW" => aW, "aE" => aE, "aS" => aS, "aN" => aN, "RelocFunctions" => RelocFunctions)
end

"""
    RelocationFunctions_cs(xmpl)

Define matrix of functions to convert indices across neighboring tiles
"""
function RelocationFunctions_cs(xmpl::MeshArray)

# f1 : 0-n,0-n => -n-0,0-n     for 1->2, 3->4, 5->6
# f2 : 0-n,0-n => n-0,-n-0     for 2->4, 4->6, 6->2
# f3 : 0-n,0-n => 0-n,-n-0     for 2->3, 4->5, 6->1
# f4 : 0-n,0-n => -n-0,n-0     for 1->3, 3->5, 5->1
# g1, g2, g3, g4 : the reverse connections

    f1(x, y, nx, ny) = (x .- Float64(nx), y)
    f2(x, y, nx, ny) = (Float64(ny) .- y .- 1.0, x .- Float64(nx))
    f3(x, y, nx, ny) = (x, y .- Float64(ny))
    f4(x, y, nx, ny) = (y .- Float64(ny), Float64(nx) .- x .- 1.0)

    g1(x, y, nx, ny) = (x .+ Float64(nx), y)
    g2(x, y, nx, ny) = (y .+ Float64(ny), Float64(nx) .- x .- 1.0)
    g3(x, y, nx, ny) = (x, y .+ Float64(ny))
    g4(x, y, nx, ny) = (Float64(ny) .- y .- 1.0, x .+ Float64(nx))

#

    s = size.(xmpl.f)
    nFaces = length(s)
    tmp = Array{Function,2}(undef, 6, 6)

# f1, f2, f3, f4 : always get nx & ny from the source tile

    tmp[2, 1] = (x, y) -> f1(x, y, s[1][1], s[1][2])
    tmp[4, 3] = (x, y) -> f1(x, y, s[3][1], s[3][2])
    tmp[6, 5] = (x, y) -> f1(x, y, s[5][1], s[5][2])

    tmp[4, 2] = (x, y) -> f2(x, y, s[2][1], s[2][2])
    tmp[6, 4] = (x, y) -> f2(x, y, s[4][1], s[4][2])
    tmp[2, 6] = (x, y) -> f2(x, y, s[6][1], s[6][2])

    tmp[3, 2] = (x, y) -> f3(x, y, s[2][1], s[2][2])
    tmp[5, 4] = (x, y) -> f3(x, y, s[4][1], s[4][2])
    tmp[1, 6] = (x, y) -> f3(x, y, s[6][1], s[6][2])

    tmp[3, 1] = (x, y) -> f4(x, y, s[1][1], s[1][2])
    tmp[5, 3] = (x, y) -> f4(x, y, s[3][1], s[3][2])
    tmp[1, 5] = (x, y) -> f4(x, y, s[5][1], s[5][2])

# g1, g2, g3, g4 : nx or ny can come from source or target + notice nx/ny flips

    tmp[1, 2] = (x, y) -> g1(x, y, s[1][1], s[2][2])
    tmp[3, 4] = (x, y) -> g1(x, y, s[3][1], s[4][2])
    tmp[5, 6] = (x, y) -> g1(x, y, s[5][1], s[6][2])

    tmp[2, 4] = (x, y) -> g2(x, y, s[4][1], s[2][1])
    tmp[4, 6] = (x, y) -> g2(x, y, s[6][1], s[4][1])
    tmp[6, 2] = (x, y) -> g2(x, y, s[2][1], s[6][1])

    tmp[2, 3] = (x, y) -> g3(x, y, s[3][1], s[2][2])
    tmp[4, 5] = (x, y) -> g3(x, y, s[5][1], s[4][2])
    tmp[6, 1] = (x, y) -> g3(x, y, s[1][1], s[6][2])

    tmp[1, 3] = (x, y) -> g4(x, y, s[1][2], s[3][2])
    tmp[3, 5] = (x, y) -> g4(x, y, s[3][2], s[5][2])
    tmp[5, 1] = (x, y) -> g4(x, y, s[5][2], s[1][2])

    return tmp

end

"""
    RelocationFunctions_cs_check(xmpl,RF,trgt)

Visualize that RelocationFunctions_cs behaves as expected
"""
function RelocationFunctions_cs_check(
    xmpl::MeshArray,
    RF::Array{Function,2},
    trgt::Int,
)

    s = size.(xmpl.f)
    nFaces = length(s)
    nFaces == 5 ? s = vcat(s, s[3]) : nothing

    (aW, aE, aS, aN, iW, iE, iS, iN) = MeshArrays.exch_cs_sources(trgt, s, 1)
    nx, ny = s[trgt]
    p = plot([0.0, nx], [0.0, ny], color = :black)
    plot!([0.0, nx], [ny, 0.0], color = :black)
    for i = 1:nFaces
        (nx, ny) = s[i]
        x = [i - 0.5 for i = 1:nx, j = 1:ny]
        y = [j - 0.5 for i = 1:nx, j = 1:ny]
        c = missing
        if aW == i
            println("source West=$i")
            c = :red
        end
        if aE == i
            println("source East=$i")
            c = :orange
        end
        if aS == i
            println("source South=$i")
            c = :blue
        end
        if aN == i
            println("source North=$i")
            c = :cyan
        end
        if !ismissing(c)
            (x, y) = RF[trgt, i](x, y)
            p = scatter!(
                x,
                y,
                color = c,
                legend = false,
                marker = :rect,
                markerstrokewidth = 0.0,
                markersize = 1.0,
            )
        end
    end
    return p
end

"""
    VelComp!(du,u,p::Dict,tim)

Interpolate velocity from gridded fields (after exchange on u0,v0)
and return position increment `du` (i.e. `x,y,fIndex`).
"""
function VelComp!(du::Array{Float64,1},u::Array{Float64,1},p::Dict,tim)
    #compute positions in index units
    dt=(tim-p["t0"])/(p["t1"]-p["t0"])
    g=p["u0"].grid
    #
    g.class=="dpdo" ? UpdateLocation_dpdo!(u,g) : nothing
    g.class=="cs" ? UpdateLocation_cs!(u,p) : nothing
    g.class=="llc" ? UpdateLocation_cs!(u,p) : nothing

    x,y = u[1:2]
    fIndex = Int(u[3])
    nx,ny=g.fSize[fIndex]
    #
    dx,dy=[x - floor(x),y - floor(y)]
    i_c,j_c = Int32.(floor.([x y])) .+ 2
    #
    i_w,i_e=[i_c i_c+1]
    j_s,j_n=[j_c j_c+1]
    #debugging stuff
    if false && (i_e>nx+2||j_n>ny+2)
        println("nx=$nx")
        println("ny=$ny")
        println("u=$u")
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
