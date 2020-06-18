using Test
using IndividualDisplacements, MeshArrays, OrdinaryDiffEq
include("helper_functions.jl")

@testset "MeshArrays tests:" begin

    uvetc=test1_setup()

    uInit=[200000.0;0.0]./uvetc["dx"]
    nSteps=3000-2
    du=fill(0.0,2);
    tspan = (0.0,nSteps*3600.0)
    prob = ODEProblem(⬡,uInit,tspan,uvetc)
    sol = solve(prob,Tsit5(),reltol=1e-8,abstol=1e-8)

    @test isapprox(sol[1,end],117237.0./uvetc["dx"]; atol=100.)
    @test isapprox(sol[2,end],40448.0./uvetc["dx"]; atol=100.)

end
