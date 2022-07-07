using LinearAlgebra
using StaticArrays

import Base: product, length

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
length(ray::Ray) = norm(ray.target - ray.origin)

"""
    CT
"""

"""
    make_xrays

Given a camera and detector, construct the rays eminating from the camera center
and the detector plane which they hit.
"""
function make_xrays(camera::Camera, detector::Detector)

    # Get the detector plane normal vector
    n̂ = camera.center - detector.center |> normalize
    u, v, w = get_basis(n̂)
    @assert w ≈ n̂

    # Construct the detector plane
    t = (-detector.height÷2:1:detector.height÷2) .* detector.Δx
    s = (-detector.width÷2:1:detector.width÷2) * detector.Δy
    get_coord(u, v, coef) = (coef .* (u, v)) |> sum |> x -> x + detector.center
    plane = [get_coord(u, v, coef) for coef in product(t, s)]

    # Get the rays
    rays = [Ray(camera.center, target) for target in plane]

    return rays

end