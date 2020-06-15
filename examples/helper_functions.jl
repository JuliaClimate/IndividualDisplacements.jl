
"""
    SetupPeriodicDomain(np::Integer=16)

Set up a periodic domain of size np x np

```
np=16 #domain size is np x np
Γ=SetPeriodicDomain(np)
```
"""
function SetupPeriodicDomain(np::Integer=16)
    γ,Γ=GridOfOnes("PeriodicDomain",1,np)
    Γ["XC"][1]=vec(0.5:1.0:np-0.5)*ones(1,np)
    Γ["XG"][1]=vec(0.0:1.0:np-1.0)*ones(1,np)
    Γ["YC"][1]=ones(np,1)*transpose(vec(0.5:1.0:np-0.5))
    Γ["YG"][1]=ones(np,1)*transpose(vec(0.0:1.0:np-1.0))
    return Γ
end


"""
    SetupRandomFlow(Γ::Dict)

Set up a random flow field over the domain specified by Γ

```
Γ=SetPeriodicDomain(16)
𝑃,ϕ=SetupRandomFlow(Γ)
```
"""
function SetupRandomFlow(Γ::Dict)
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

