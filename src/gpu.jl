using Adapt
using CUDA


# Redefine DRR structs for the GPU
struct CTGPU
    volume::CuArray{T,3} where {T<:Real}
    ΔX::Float64
    ΔY::Float64
    ΔZ::Float64
    X₀::Float64
    Y₀::Float64
    Z₀::Float64
end


struct CameraGPU
    center::CuArray{Float32,1}
end


struct DetectorGPU
    center::CuArray{Float32,1}
    height::Int64
    width::Int64
    Δx::Float32
    Δy::Float32
end


struct RayGPU
    origin::CuArray{Float32,1}
    target::CuArray{Float32,1}
end


# Define CPU -> GPU struct conversion functions
gpu(ct::CT) = CTGPU(cu(ct.volume), ct.ΔX, ct.ΔY, ct.ΔZ, ct.X₀, ct.Y₀, ct.Z₀)
gpu(camera::Camera) = CameraGPU(cu(camera.center))
gpu(detector::Detector) = DetectorGPU(cu(detector.center), detector.height, detector.width, detector.Δx, detector.Δy)
gpu(ray::Ray) = RayGPU(cu(ray.origin), cu(ray.target))

