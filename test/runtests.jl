using DRRs
using Test

@testset "DRRs.jl" begin
    @info "Testing I/O submodule..."
    ct = read_dicom("../data/cxr")
    @test size(ct.volume) == (512, 512, 133)
    @show ct.ΔX, ct.ΔY, ct.ΔZ

    @info "Testing camera submodule..."
    camera = Camera([0, 0, 0])
    detector = Detector([2, 1, 0.2], 101, 101, 0.4, 0.4)
    @test size(make_xrays(camera, detector)) == (101, 101)

    @info "Testing trilinear interpolation..."
    x0, y0, z0, x1, y1, z1 = randn(6)
    M = make_coordinate_matrix(x0, y0, z0, x1, y1, z1)
    Minv = make_inverse_coordinate_matrix(x0, y0, z0, x1, y1, z1)
    @test size(M) == size(Minv)
    @test inv(M) ≈ Minv

    @info "Testing DRR generator..."
    drr = sampler(ct, camera, detector; ray2pix=:trilinear, spacing=0.1)
    @test size(drr) == (101, 101)
    drr = sampler(ct, camera, detector; ray2pix=:sample, spacing=0.1)
    @test size(drr) == (101, 101)
    drr = siddon(ct, camera, detector)

    @info "Testing DRR visualization..."
    # @test plot(ct, camera, detector)
    # @test plot(ct, camera, detector, drr)

end
