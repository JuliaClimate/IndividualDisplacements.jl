using Test, Documenter
using IndividualDisplacements, OceanStateEstimation
import IndividualDisplacements.MeshArrays as MeshArrays

OceanStateEstimation.get_ecco_velocity_if_needed()
OceanStateEstimation.get_occa_velocity_if_needed()
MeshArrays.GridLoad(MeshArrays.GridSpec("LatLonCap",MeshArrays.GRID_LLC90))
MeshArrays.GridLoad(MeshArrays.GridSpec("PeriodicChannel",MeshArrays.GRID_LL360))
IndividualDisplacements.flt_example_download()

@testset "test3" begin
    p=dirname(pathof(IndividualDisplacements))
    include(joinpath(p,"../examples/jupyter/random_flow_field.jl"))

    tmp1=randn_lonlat(10)

    show(𝐼)
    diff(𝐼)
    size(𝐼)
    𝐽=similar(𝐼)
    @test isa(𝐽,Individuals)

    (U,V,Φ)=IndividualDisplacements.random_flow_field("Rotational Component")
    𝐺=convert_to_FlowFields(u,v,10.0)
    tmp2=nearest_to_xy(𝐺.u0,3.,3.,1.)
    @test isa(tmp2,Array)
    tmp3=nearest_to_xy(𝐹.u0,3.,3.)
    @test isa(tmp3,Array)
end

@testset "doctests" begin
    doctest(IndividualDisplacements; manual = false)
end
