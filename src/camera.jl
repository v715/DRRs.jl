using LinearAlgebra
using StaticArrays

import Base: product

"""
    Camera
"""
mutable struct Camera
    center::SVector{3,Float64}
end

"""
    Detector Plane
"""
mutable struct Detector
    center::SVector{3,Float64}
    normal::SVector{3,Float64}
    height::Int64
    width::Int64
    Δx::Float64
    Δy::Float64
end

function make_plane(detector::Detector)
    d = detector.center' * normalize(detector.normal)
    xs = ((-detector.height÷2:1:detector.height÷2) .* detector.Δx) .+ detector.center[1]
    ys = ((-detector.width÷2:1:detector.width÷2) * detector.Δy) .+ detector.center[2]
    xys = product(xs, ys)
    zs = [findz(xy..., normalize(detector.normal)..., d) for xy in xys]
    return [SVector(xy..., z) for (xy, z) in zip(xys, zs)]
end

findz(x, y, a, b, c, d) = (d - a * x - b * y) / c

"""
    Ray
"""
mutable struct Ray
    origin::SVector{3,Float64}
    direction::SVector{3,Float64}
end

function get_rays(camera::Camera, detector::Detector)
    directions = make_plane(detector) .|> pixel -> pixel - camera.center
    return [Ray(camera.center, direction) for direction in directions]
end

trace(t::Float64; ray::Ray) = ray.origin + ray.direction * t