"""
    test1_setup()

Call `gcmgrid`, initialize a single point,
rely on `dxy_dt`, and just output `sol` at the end.

```
using IndividualDisplacements, MeshArrays, OrdinaryDiffEq
𝑃,sol=test1_setup()
```
"""
function test1_setup()

    mygrid=gcmgrid("flt_example/","ll",1,[(80,42)], [80 42], Float32, read, write)
    XC=MeshArray(mygrid,Float32); XC[1]=vec(2500.:5000.:397500.0)*ones(1,42);
    XG=MeshArray(mygrid,Float32); XG[1]=vec(0.:5000.:395000.0)*ones(1,42);
    YC=MeshArray(mygrid,Float32); YC[1]=ones(80,1)*transpose(vec(2500.:5000.:207500.0));
    YG=MeshArray(mygrid,Float32); YG[1]=ones(80,1)*transpose(vec(0.:5000.:205000.0));

    dx=5000.0
    t0=0.0; t1=18001.0*3600.0
    u=-(YG.-YC[1][40,21])/2000000.
    v=(XG.-XC[1][40,21])/2000000.
    u0=u./dx; u1=u./dx
    v0=v./dx; v1=v./dx

    𝑃=𝑃_Array2D{eltype(u)}(u0[1], u1[1], v0[1], v1[1], [t0,t1])
    
    u0=[200000.0;0.0]./dx
    du=fill(0.0,2);
    prob = ODEProblem(dxy_dt,u0,[0.0,2998*3600.0],𝑃)
    sol = solve(prob,Tsit5(),reltol=1e-8,abstol=1e-8)

    return 𝑃,sol
end

"""
    test2_periodic_domain(np = 12, nq = 12)

Call `simple_periodic_domain`, initialize 6x6 point cloud,
rely on `dxy_dt!`, and call `postprocess_xy` at the end.

```
using IndividualDisplacements, MeshArrays, OrdinaryDiffEq
df,𝑃=test2_periodic_domain()

using Plots
@gif for t in 𝑃.t0:1.0:𝑃.t1
   scatter_subset(𝑃,df,t)
end
```
"""
function test2_periodic_domain(np = 12, nq = 12)
    #domain and time parameters
    Γ = simple_periodic_domain(np, nq)
    Γ = IndividualDisplacements.dict_to_nt(Γ)

    u = 0.1 ./ Γ.DXC
    v = 0.3 ./ Γ.DYC
    (u, v) = exchange(u, v, 1)

    f = (u -> IndividualDisplacements.update_location_dpdo!(u,Γ.XC.grid))
    𝑃=𝑃_MeshArray2D{eltype(u)}(u,u,v,v,[0.0,400.0],f)

    #initial conditions
    x0 = np * (0.4:0.04:0.6)
    y0 = nq * (0.4:0.04:0.6)
    x0 = vec(x0) * ones(1, length(y0))
    y0 = ones(size(x0, 1), 1) * transpose(vec(y0))
    u0 = permutedims([[x0[i];y0[i];1.0] for i in eachindex(x0)])
    du=0*u0
    
    #solve for trajectories
    prob = ODEProblem(dxy_dt!, u0, 𝑃.𝑇, 𝑃)
    sol = solve(prob,Euler(),dt=0.1)

    return postprocess_xy(sol, 𝑃),𝑃
end

"""
    scatter_subset(𝑃,df,t=missing,dt=1.0)

```
@gif for t in 𝑃.t0:1.0:𝑃.t1
   scatter_subset(𝑃,df,t)
end
```
"""
function scatter_subset(𝑃,df,t=missing,dt=1.0)
    ismissing(t) ? t=maximum(df[!,:t]) : nothing
    df_t = df[ (df.t.>t-dt).&(df.t.<=t) , :]
    nx,ny=size(𝑃.XC[1])
    scatter(df_t.x,df_t.y,markersize=2.0,c=:red,
    xlims=(0,nx),ylims=(0,ny),leg=:none,marker = (:circle, stroke(0)))
end
