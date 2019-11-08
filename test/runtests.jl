using Test
using IndividualDisplacements, MeshArrays, DifferentialEquations

@testset "MeshArrays tests:" begin

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
uvetc=merge(uvt,GridVariables);

uInit=[200000.0;0.0]./uvetc["dx"]
nSteps=3000-2
du=fill(0.0,2);
tspan = (0.0,nSteps*3600.0)
prob = ODEProblem(IndividualDisplacements.VelComp,uInit,tspan,uvetc)
sol = solve(prob,Tsit5(),reltol=1e-8,abstol=1e-8)

@test isapprox(sol[1,end],117237.0./uvetc["dx"]; atol=100.)
@test isapprox(sol[2,end],40448.0./uvetc["dx"]; atol=100.)

end
