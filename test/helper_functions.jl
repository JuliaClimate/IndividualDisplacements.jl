function test1_setup()

    mygrid=gcmgrid("flt_example/","ll",1,[(80,42)], [80 42], Float32, read, write)
    XC=MeshArray(mygrid,Float32); XC[1]=vec(2500.:5000.:397500.0)*ones(1,42);
    XG=MeshArray(mygrid,Float32); XG[1]=vec(0.:5000.:395000.0)*ones(1,42);
    YC=MeshArray(mygrid,Float32); YC[1]=ones(80,1)*transpose(vec(2500.:5000.:207500.0));
    YG=MeshArray(mygrid,Float32); YG[1]=ones(80,1)*transpose(vec(0.:5000.:205000.0));
    GridVariables=Dict("XC" => XC,"YC" => YC,"XG" => XG,"YG" => YG,"dx" => 5000.0);

    t0=0.0; t1=18001.0*3600.0
    u0=-(YG.-YC[1][40,21])/2000000.; u1=u0
    v0=(XG.-XC[1][40,21])/2000000.; v1=v0

    u0=u0./GridVariables["dx"]
    u1=u1./GridVariables["dx"]
    v0=v0./GridVariables["dx"]
    v1=v1./GridVariables["dx"]

    uvt = Dict("u0" => u0, "u1" => u1, "v0" => v0, "v1" => v1, "t0" => t0, "t1" => t1)
    return merge(uvt,GridVariables)
end


function test2_periodic_domain(np = 12, nq = 12)
    #domain and time parameters
    Γ = simple_periodic_domain(np, nq)

    u = 0.1 * Γ["DXC"]
    v = 0.3 * Γ["DYC"]
    (u, v) = exchange(u, v, 1)
    𝑃 = Dict(
        "u0" => u,
        "u1" => u,
        "v0" => v,
        "v1" => v,
        "t0" => 0.0,
        "t1" => 400.0,
        "dt" => 0.1,
    )
    𝑃 = merge(𝑃, Γ)#add grid variables

    #initial conditions
    x0 = np * (0.4:0.04:0.6)
    y0 = nq * (0.4:0.04:0.6)
    x0 = vec(x0) * ones(1, length(y0))
    y0 = ones(size(x0, 1), 1) * transpose(vec(y0))
    u0 = transpose([x0[:] y0[:] ones(size(x0[:]))])

    #solve for trajectories
    𝑇 = (𝑃["t0"], 𝑃["t1"])
    prob = ODEProblem(⬡!, u0, 𝑇, 𝑃)
    sol = solve(prob,Euler(),dt=𝑃["dt"])

    return postprocess_xy(sol, 𝑃)
end

"""
    scatter_subset(Γ,df,t=missing,dt=1.0)

```
@gif for t in 𝑃["t0"]:1.0:𝑃["t1"]
   scatter_subset(Γ,df,t)
end
```
"""
function scatter_subset(Γ,df,t=missing,dt=1.0)
    ismissing(t) ? t=maximum(df[!,:t]) : nothing
    df_t = df[ (df.t.>t-dt).&(df.t.<=t) , :]
    nx,ny=size(Γ["XC"][1])
    scatter(df_t.x,df_t.y,markersize=2.0,c=:red,
    xlims=(0,nx),ylims=(0,ny),leg=:none,marker = (:circle, stroke(0)))
end
