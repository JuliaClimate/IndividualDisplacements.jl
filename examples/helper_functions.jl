
"""
    get_grid_if_needed()

Download global `MITgcm` grid and transport output to `examples/GRID_LLC90`
"""
function get_grid_if_needed()
  p=dirname(pathof(IndividualDisplacements))
  p=joinpath(p,"../examples/GRID_LLC90")
  r="https://github.com/gaelforget/GRID_LLC90"
  !isdir(p) ? run(`git clone $r $p`) : nothing
end

"""
    get_velocity_if_needed()

Download `MITgcm` transport output to `examples/nctiles_climatology` if needed
"""
function get_velocity_if_needed()
    p=dirname(pathof(IndividualDisplacements))
    pth="$p/../examples/nctiles_climatology/"
    !isdir("$pth") ? mkdir("$pth") : nothing
    !isdir("$pth"*"UVELMASS") ? get_from_dataverse("UVELMASS",pth) : nothing
    !isdir("$pth"*"VVELMASS") ? get_from_dataverse("VVELMASS",pth) : nothing
end

"""
    get_flt_ex_if_needed()

Download simple grid, velocity, & trajectory output from `MITgcm/pkg/flt`
to `examples/flt_example`
"""
function get_flt_ex_if_needed()
    p=dirname(pathof(IndividualDisplacements))
    p=joinpath(p,"../examples/flt_example")
    r="https://github.com/gaelforget/flt_example"
    !isdir(p) ? run(`git clone $r $p`) : nothing
end

"""
    setup_random_flow(Γ::Dict)

Set up a random flow field over the domain specified by Γ

```
Γ=simple_periodic_domain(12)
𝑃,ϕ=setup_random_flow(Γ)
```
"""
function setup_random_flow(Γ::Dict)
  (_,ϕ,_,_)=demo2(Γ);

  (u,v)=gradient(ϕ,Γ)
  u=u./Γ["DXC"]#normalization to grid units
  v=v./Γ["DYC"]

  (u,v)=exchange(u,v,1)
  u0=-v; u1=-v;
  v0=u; v1=u;

  𝑃 = (u0=u0, u1=u1, v0=v0, v1=v1, 𝑇=[0.0,400.0], ioSize=ϕ.grid.ioSize)
  return 𝑃,ϕ

end

"""
    setup_global_ocean(k::Int)

Set up Global Ocean particle simulation in 2D with seasonally varying flow field.

```
𝑃=setup_global_ocean(10);
```
"""
function setup_global_ocean(k::Int)

  #k=10 #choice of vertical level
  ny=2 #number of simulated years (20 for k>20)
  r_reset = 0.01 #fraction of the particles reset per month (0.05 for k<=10)

  #read grid and set up connections between subdomains
  γ=GridSpec("LatLonCap",joinpath(p,"../examples/GRID_LLC90/"))
  Γ=GridLoad(γ)
  Γ=merge(Γ,IndividualDisplacements.NeighborTileIndices_cs(Γ))

  #initialize u0,u1 etc
  𝑃=set_up_𝑃(k,0.0,Γ,joinpath(p,"../examples/nctiles_climatology/"));

  #add parameters for use in reset!
  tmp=(frac=r_reset, Γ=Γ)
  𝑃=merge(𝑃,tmp)

  return 𝑃

end

function init_global_range(lons::Tuple = (-160.0, -150.0),lats::Tuple = (35.0, 45.0))
    lo0, lo1 = lons #(-160.0, -150.0)
    la0, la1 = lats #(35.0, 45.0)
    np = 100
    lon = lo0 .+ (lo1 - lo0) .* rand(np)
    lat = la0 .+ (la1 - la0) .* rand(np)
    (u0, _) = initialize_lonlat(Γ, lon, lat; msk = Γ["hFacC"][:, k])
    id=collect(1:np)
    return u0
end

function init_global_randn(np ::Int , 𝑃::NamedTuple)
    (lon, lat) = randn_lonlat(2*np)
    (u0, _) = initialize_lonlat(𝑃.Γ, lon, lat; msk = 𝑃.msk)
    u0[:,1:np]
end
