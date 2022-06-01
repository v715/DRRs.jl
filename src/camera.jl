using StaticArrays

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