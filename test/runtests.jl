using Test, Documenter
using IndividualDisplacements, MeshArrays, OrdinaryDiffEq, OceanStateEstimation
include(joinpath(dirname(pathof(IndividualDisplacements)),"../examples/helper_functions.jl"))

IndividualDisplacements.flt_example_download()
OceanStateEstimation.get_ecco_velocity_if_needed()
OceanStateEstimation.get_occa_velocity_if_needed()
GridLoad(GridSpec("LatLonCap",MeshArrays.GRID_LLC90))
GridLoad(GridSpec("PeriodicChannel",MeshArrays.GRID_LL360))

@testset "test1" begin
    uvetc,sol=test1_setup()
    @test isapprox(sol[1,end],23.4474; atol=0.01)
    @test isapprox(sol[2,end],8.0896; atol=0.01)
end

@testset "test2" begin
    df,𝑃=test2_periodic_domain()
    @test prod(isapprox.(df[end-35:6:end,:y],12*(0.4:0.04:0.6)))
    @test prod(isapprox.(df[end-5:end,:x],12*(0.4:0.04:0.6).+4.0))
end

@testset "test3" begin
    p=dirname(pathof(IndividualDisplacements))
    include(joinpath(p,"../examples/basics/random_flow_field.jl"))

    show(𝐼)
    diff(𝐼)
    size(𝐼)
    𝐽=similar(𝐼)
    @test isa(𝐽,Individuals)

    tmp1=randn_lonlat(10)
    𝐺=convert_to_FlowFields(u,v,10.0)
    tmp2=nearest_to_xy(𝐺.u0,3.,3.,1.)
    @test isa(tmp2,Array)
    tmp3=nearest_to_xy(𝐹.u0,3.,3.)
    @test isa(tmp3,Array)
end

@testset "doctests" begin
    doctest(IndividualDisplacements; manual = false)
end
