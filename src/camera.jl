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
    origin::SVector{3,Float64}
    normal::SVector{3,Float64}
    height::Int64
    width::Int64
    Δx::Float64
    Δy::Float64
end

function make_plane(detector::Detector)
    d = d = detector.origin' * detector.normal
    xs = (-detector.height÷2:1:detector.height÷2) * detector.Δx
    ys = (-detector.width÷2:1:detector.width÷2) * detector.Δy
    xys = product(xs, ys)
    zs = [findz(xy..., detector.normal..., d) for xy in xys]
    return [SVector(xy..., z) for (xy, z) in zip(xys, zs)]
end

findz(x, y, a, b, c, d) = (d - a * x - b * y) / c