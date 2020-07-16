
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

  𝑃 = Dict( "u0" => u0, "u1" => u1, "v0" => v0, "v1" => v1,
            "t0" => 0.0, "t1" => 400.0, "dt" => 0.1)
  𝑃=merge(𝑃,Γ)#add grid variables

  return 𝑃,ϕ
end
