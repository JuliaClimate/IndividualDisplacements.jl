
using DataFrames

#VelCopy(du,uInit,tmp,3600.0)
#VelComp(du,uInit,uvetc,3600.0)

"""
    VelComp(du,u,p::Dict,tim)

Interpolate velocity from gridded fields and return position increment `du`
"""
function VelComp(du,u,p::Dict,tim)
    #compute positions in index units
    dt=(tim-p["t0"])/(p["t1"]-p["t0"])
    #
    x,y = u ./ p["dx"]
    x,y=[mod(x,80.0),mod(y,42.0)]
    #
    dx,dy=[x - floor(x),y - floor(y)]
    i_c,j_c = Int32.(floor.([x y])) .+ 1
    i_c==0 ? i_c=80 : nothing
    j_c==0 ? j_c=42 : nothing
    #
    i_w,i_e=[i_c i_c+1]
    x>=79.0 ? (i_w,i_e)=(80,1) : nothing
    j_s,j_n=[j_c j_c+1]
    y>=41.0 ? (j_s,j_n)=(42,1) : nothing
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
    #show(p["u0"].f[1][x0,y0,1])
    ##du[2]=0.4 * du[2]
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
