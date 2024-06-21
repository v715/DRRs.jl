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



struct RayGPU
    origin::CuArray{Float32,1}
    target::CuArray{Float32,1}
end
trace(t::Float32; ray::RayGPU) = ray.origin + (ray.target - ray.origin) * t


# Define CPU -> GPU struct conversion functions
gpu(ct::CT) = CTGPU(cu(ct.volume), ct.ΔX, ct.ΔY, ct.ΔZ, ct.X₀, ct.Y₀, ct.Z₀)
gpu(ray::Ray) = RayGPU(cu(ray.origin), cu(ray.target))

# Add converter for the rays to trace
gpu(projector::Matrix{Ray}) = [gpu(ray) for ray in projector]
