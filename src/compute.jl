"""
    dxy_dt!(du,u,p::Dict,tim)

Interpolate velocity from gridded fields (2D; with halos) to position `u`
(`x,y,fIndex`) to compute the derivative of position v time  `du_dt`.
"""
function dxy_dt!(du::Array{Float64,1},u::Array{Float64,1},p::Dict,tim)
    #compute positions in index units
    dt=(tim-p["t0"])/(p["t1"]-p["t0"])
    dt>1.0 ? error("dt>1.0") : nothing
    dt<0.0 ? error("dt>0.0") : nothing
    g=p["u0"].grid
    #
    g.class=="PeriodicDomain" ? update_location_dpdo!(u,g) : nothing
    g.class=="CubeSphere" ? update_location_cs!(u,p) : nothing
    g.class=="LatLonCap" ? update_location_cs!(u,p) : nothing

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

function dxy_dt!(du::Array{Float64,2},u::Array{Float64,2},p::Dict,tim)
    for i=1:size(u,2)
        tmpdu=du[1:3,i]
        tmpu=u[1:3,i]
        dxy_dt!(tmpdu,tmpu,p,tim)
        du[1:3,i]=tmpdu
        u[1:3,i]=tmpu
    end
    return du
end

"""
    dxyz_dt(du,u,p::Dict,tim)

Interpolate velocity from gridded fields (3D; NO halos) to position `u`
(`x,y,z`) to compute the derivative of position v time  `du_dt`.
"""
function dxyz_dt(du::Array{Float64,1},u::Array{Float64,1},p::Dict,tim)
    #compute positions in index units
    dt=(tim-p["t0"])/(p["t1"]-p["t0"])
    #
    x,y,z = u[1:3]
    nx,ny=p["u0"].grid.ioSize
    x,y=[mod(x,nx),mod(y,ny)]
    nz=size(p["u0"],2)
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
    du[1]=(1.0-dx)*(1.0-dt)*p["u0"].f[1,k_c][i_w,j_c]+
    dx*(1.0-dt)*p["u0"].f[1,k_c][i_e,j_c]+
    (1.0-dx)*dt*p["u1"].f[1,k_c][i_w,j_c]+
    dx*dt*p["u1"].f[1,k_c][i_e,j_c]
    #interpolate v to position and time
    du[2]=(1.0-dy)*(1.0-dt)*p["v0"].f[1,k_c][i_c,j_s]+
    dy*(1.0-dt)*p["v0"].f[1,k_c][i_c,j_n]+
    (1.0-dy)*dt*p["v1"].f[1,k_c][i_c,j_s]+
    dy*dt*p["v1"].f[1,k_c][i_c,j_n]
    #interpolate w to position and time
    du[3]=(1.0-dz)*(1.0-dt)*p["w0"].f[1,k_l][i_c,j_c]+
    dz*(1.0-dt)*p["w0"].f[1,k_r][i_c,j_c]+
    (1.0-dz)*dt*p["w1"].f[1,k_l][i_c,j_c]+
    dz*dt*p["w1"].f[1,k_r][i_c,j_c]
    #
    return du
end

function dxyz_dt(du::Array{Float64,2},u::Array{Float64,2},p::Dict,tim)
    for i=1:size(u,2)
        tmpdu=du[1:3,i]
        dxyz_dt(tmpdu,u[1:3,i],p,tim)
        du[1:3,i]=tmpdu
    end
    return du
end

"""
    dxy_dt(du,u,p::Dict,tim)

Interpolate velocity from gridded fields (2D; NO halos) to position `u`
(`x,y`) to compute the derivative of position v time  `du_dt`.
"""
function dxy_dt(du::Array{Float64,1},u::Array{Float64,1},p::Dict,tim)
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

function dxy_dt(du::Array{Float64,2},u::Array{Float64,2},p::Dict,tim)
    for i=1:size(u,2)
        tmpdu=du[1:2,i]
        dxy_dt(tmpdu,u[1:2,i],p,tim)
        du[1:2,i]=tmpdu
    end
    return du
end

"""
    dxy_dt_CyclicArray(du,u,p::NamedTuple,tim)

Nearest neighbor (?) velocity from gridded fields (2D; NO halos but
not needed when CyclicArrays is used to extend valid indice ranges).

_notes:_ spatial interpolation & temporal interpolation are lacking
"""
function dxy_dt_CyclicArray(du::Array{Float64,2},u::Array{Float64,2},p::NamedTuple,tim)
    np=size(du,2)
    xi,yi=(u[1,:],u[2,:])
    @unpack xg,yg,u,v = p

    #not needed:
    #xi=floor.(xg[1,Int.(sign.(xi).*floor.(abs.(xi)))])+rem.(xi,1)
    #yi=floor.(yg[Int.(sign.(yi).*floor.(abs.(yi))),1])+rem.(yi,1)

    i=Int.(floor.(xg[1,Int.(floor.(xi))]))
    j=Int.(floor.(yg[Int.(floor.(yi)),1]))
    du[1,:]=[p.u[i[ii],j[ii]] for ii in 1:np]
    du[2,:]=[p.v[i[ii],j[ii]] for ii in 1:np]
    return du
end
