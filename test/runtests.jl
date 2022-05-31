using DRRs
using Test

@testset "DRRs.jl" begin
    @info "Testing I/O submodule"
    volume, ΔX, ΔY, ΔZ = read_dicom("../data/cxr")
    @test size(volume) == (512, 512, 133)
    @show ΔX, ΔY, ΔZ
end
