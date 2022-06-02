using LinearAlgebra
using StaticArrays


function make_coordinate_matrix(x0, y0, z0, x1, y1, z1)
    M = SMatrix{8,8}([
        1 x0 y0 z0 x0*y0 x0*z0 y0*z0 x0*y0*z0
        1 x1 y0 z0 x1*y0 x1*z0 y0*z0 x1*y0*z0
        1 x0 y1 z0 x0*y1 x0*z0 y1*z0 x0*y1*z0
        1 x0 y1 z1 x0*y1 x0*z1 y1*z1 x0*y1*z1
        1 x0 y0 z1 x0*y0 x0*z1 y0*z1 x0*y0*z1
        1 x1 y0 z1 x1*y0 x1*z1 y0*z1 x1*y0*z1
        1 x1 y1 z0 x1*y1 x1*z0 y1*z0 x1*y1*z0
        1 x1 y1 z1 x1*y1 x1*z1 y1*z1 x1*y1*z1
    ])
    return M
end


function make_inverse_coordinate_matrix(x0, y0, z0, x1, y1, z1)
    volume = (x1 - x0) * (y1 - y0) * (z1 - z0)
    Minv = SMatrix{8,8}([
        x1*y1*z1 -x0*y1*z1 -x1*y0*z1 x1*y0*z0 -x1*y1*z0 x0*y1*z0 x0*y0*z1 -x0*y0*z0
        -y1*z1 y1*z1 y0*z1 -y0*z0 y1*z0 -y1*z0 -y0*z1 y0*z0
        -x1*z1 x0*z1 x1*z1 -x1*z0 x1*z0 -x0*z0 -x0*z1 x0*z0
        -x1*y1 x0*y1 x1*y0 -x1*y0 x1*y1 -x0*y1 -x0*y0 x0*y0
        z1 -z1 -z1 z0 -z0 z0 z1 -z0
        y1 -y1 -y0 y0 -y1 y1 y0 -y0
        x1 -x0 -x1 x1 -x1 x0 x0 -x0
        -1 1 1 -1 1 -1 -1 1
    ])
    return Minv / volume
end


function get_colors(volume, xidx, yidx, zidx)
    c000 = volume[xidx, yidx, zidx]
    c100 = volume[xidx+1, yidx, zidx]
    c010 = volume[xidx, yidx+1, zidx]
    c110 = volume[xidx+1, yidx+1, zidx]
    c001 = volume[xidx, yidx, zidx+1]
    c101 = volume[xidx+1, yidx, zidx+1]
    c110 = volume[xidx+1, yidx+1, zidx]
    c111 = volume[xidx+1, yidx+1, zidx+1]
    return SVector{8}([c000 c100 c010 c110 c001 c101 c110 c111])
end


function trilinear(x::Float64, y::Float64, z::Float64; volume::AbstractArray, xs::StepRangeLen, ys::StepRangeLen, zs::StepRangeLen)

    # Find the indices of the lower left corner of the cube we're inside of
    xidx = findlast(xs .<= x)
    yidx = findlast(ys .<= y)
    zidx = findlast(zs .<= z)
    if any(map(x -> isnothing(x), (xidx, yidx, zidx)))
        return -1024.0
    end
    if xidx == length(xs) || yidx == length(ys) || zidx == length(zs)
        return -1024.0
    end

    # Get the coordinate values of the lower left and upper right corners
    x0, y0, z0 = xs[xidx], ys[yidx], zs[zidx]
    x1, y1, z1 = xs[xidx+1], ys[yidx+1], zs[zidx+1]

    # Get the coordinate matrices
    Minv = make_inverse_coordinate_matrix(x0, y0, z0, x1, y1, z1)

    # Get the colors of the corners
    c = get_colors(volume, xidx, yidx, zidx)

    # Get the component vector
    p = SVector{8}([1 x y z x*y x*z y*z x*y*z])

    return p' * Minv * c

end
