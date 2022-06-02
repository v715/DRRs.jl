function sample(x::Float64, y::Float64, z::Float64; volume::AbstractArray, xs::StepRangeLen, ys::StepRangeLen, zs::StepRangeLen)
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

    return volume[xidx, yidx, zidx]
end