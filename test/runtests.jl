using DRRs
using Test

@testset "DRRs.jl" begin
    @info "Testing I/O submodule..."
    volume, ΔX, ΔY, ΔZ = read_dicom("../data/cxr")
    @test size(volume) == (512, 512, 133)
    @show ΔX, ΔY, ΔZ

    @info "Testing camera submodule..."
    camera = Camera([0, 0, 0])
    detector = Detector([2, 1, 0.2], [-1, 0, 1], 101, 101, 0.4, 0.4)
    @test size(make_plane(detector)) == (101, 101)
    @test size(get_rays(camera, detector)) == (101, 101)

    @info "Testing trilinear interpolation..."
    x0, y0, z0, x1, y1, z1 = randn(6)
    M = make_coordinate_matrix(x0, y0, z0, x1, y1, z1)
    Minv = make_inverse_coordinate_matrix(x0, y0, z0, x1, y1, z1)
    @test size(M) == size(Minv)
    @test inv(M) ≈ Minv

    @info "Testing DRR generator..."
    drr = make_drr(volume, ΔX, ΔY, ΔZ, detector, camera, 0.1, trilinear)
    @test size(drr) == (101, 101)
    drr = make_drr(volume, ΔX, ΔY, ΔZ, detector, camera, 0.1, sample)
    @test size(drr) == (101, 101)

end
