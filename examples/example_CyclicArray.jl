
using CyclicArrays
using UnPack

function dxy_dt_CyclicArray(du::Array{Float64,2},u::Array{Float64,2},p::NamedTuple,tim)
    np=size(du,2)
    xi,yi=(u[1,:],u[2,:])
    @unpack xg,yg,u,v = p

    #xi=floor.(xg[1,Int.(sign.(xi).*floor.(abs.(xi)))])+rem.(xi,1)
    #yi=floor.(yg[Int.(sign.(yi).*floor.(abs.(yi))),1])+rem.(yi,1)

    i=Int.(floor.(xg[1,Int.(floor.(xi))]))
    j=Int.(floor.(yg[Int.(floor.(yi)),1]))
    du[1,:]=[p.u[i[ii],j[ii]] for ii in 1:np]
    du[2,:]=[p.v[i[ii],j[ii]] for ii in 1:np]
    return du
end

function cyclicarray_example()
    𝑃=cyclicarray_setup()

    x0=collect(0:20:490).%100 .+ 170
    y0=reshape(repeat(collect(0:20:80),1,5)',25) .+ 160
    uInit=Float64.([x0';y0'])

    N=1500
    tspan = (0.0,N*𝑃.dt)

    prob = ODEProblem(dxy_dt_CyclicArray,uInit,tspan,𝑃)
    #sol = solve(prob,Tsit5(),reltol=1e-5,abstol=1e-5)
    sol = solve(prob,Euler(),dt=𝑃.dt)
end

function cyclicarray_setup()
    nx=ny=100
    nfaces=1;

    faces=zeros(nfaces,2,2,4);
    faces[1,1,1,:]=[1,1,2,0];
    faces[1,1,2,:]=[1,1,1,0];
    faces[1,2,1,:]=[1,2,2,0];
    faces[1,2,2,:]=[1,2,1,0];

    g=CyclicArray(faces);

    ind=range(1, 360, length = 360)
    x = y = ind/180*pi
    f(y,x) = sin(x) + sin(y)

    X = repeat(reshape(x, 1, :), length(y), 1)
    Y = repeat(y, 1, length(x))

    F=f.(Y,X)

    #f = sin.(2*x .+ 10 .* sin.(x)) .+ sin.(3*y)

    phi=CyclicArray(F,g)
    contourf(phi.data, aspect_ratio=1,xlim=(0,360),ylim=(0,360))
    xg=CyclicArray(repeat(ind, 1, length(x))',g)
    yg=CyclicArray(repeat(ind, 1, length(x)),g)
    u=diff(phi,dims=2)
    v=-diff(phi,dims=1)

    (u = u, v = v, xg = xg, yg = yg, dt = 100)
end
