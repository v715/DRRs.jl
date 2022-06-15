using LinearAlgebra
using StaticArrays

import Base: product

"""
    Camera
"""
mutable struct Camera
    center::MVector{3,Float64}
end

"""
    Detector Plane
"""
mutable struct Detector
    center::MVector{3,Float64}
    height::Int64
    width::Int64
    Δx::Float64
    Δy::Float64
end

"""
    Ray
"""
mutable struct Ray
    origin::SVector{3,Float64}
    target::SVector{3,Float64}
end

trace(t::Float64; ray::Ray) = ray.origin + (ray.target - ray.origin) * t

"""
    make_xrays

Given a camera and detector, construct the rays eminating from the camera center
and the detector plane which they hit.
"""
function make_xrays(camera::Camera, detector::Detector)
    # Get the detector plane normal vector
    n̂ = camera.center - detector.center |> normalize
    d = detector.center' * n̂

    # Construct the detector plane
    xs = ((-detector.height÷2:1:detector.height÷2) .* detector.Δx) .+ detector.center[1]
    ys = ((-detector.width÷2:1:detector.width÷2) * detector.Δy) .+ detector.center[2]
    xys = product(xs, ys)
    zs = [findz(xy..., n̂..., d) for xy in xys]
    plane = [SVector(xy..., z) for (xy, z) in zip(xys, zs)]

    # Get the rays
    rays = [Ray(camera.center, target) for target in plane]

    return rays
end

findz(x, y, a, b, c, d) = (d - a * x - b * y) / c