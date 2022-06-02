using DRRs
using Test

@testset "DRRs.jl" begin
    @info "Testing I/O submodule"
    volume, ΔX, ΔY, ΔZ = read_dicom("../data/cxr")
    @test size(volume) == (512, 512, 133)
    @show ΔX, ΔY, ΔZ

    @info "Testing camera submodule"
    detector = Detector([2, 1, 0.2], [-1, 0, 1], 101, 101, 0.4, 0.4)
    @test size(make_plane(detector)) == (101, 101)
end
