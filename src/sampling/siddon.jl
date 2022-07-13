struct Siddon{ArrFloat,ArrInt}
    origin::ArrFloat
    target::ArrFloat
    spacing::ArrFloat
    isocenter::ArrFloat
    dims::ArrInt
end


function get_α(i::Int, j::Int, k::Int; sid::Siddon)
    planes = [i, j, k]
    return @. (sid.isocenter + planes * sid.spacing - sid.origin) / (sid.target - sid.origin)
end


function get_α_minmax(sid::Siddon)
    αx0, αy0, αz0 = get_α(0, 0, 0; sid)
    αx1, αy1, αz1 = (sid.dims .- 1) |> nxyz -> get_α(nxyz...; sid)
    αxmin, αxmax = minmax(αx0, αx1)
    αymin, αymax = minmax(αy0, αy1)
    αzmin, αzmax = minmax(αz0, αz1)
    αmin = max(αxmin, αymin, αzmin)
    αmax = min(αxmax, αymax, αzmax)
    return αxmin, αxmax, αymin, αymax, αzmin, αzmax, αmin, αmax
end


function get_φ(α::Float64; sid::Siddon)
    pxyz = @. sid.origin + α * (sid.target - sid.origin)  # Trace the ray
    return @. (pxyz - sid.isocenter) / sid.spacing
end


function get_idx_minmax(
    αmin::Float64, αmax::Float64, αxmin::Float64, αxmax::Float64,
    ixmin::Float64, ixmax::Float64, p1::Float64, p2::Float64, nx::Int64
)
    if p1 ≤ p2
        imin = αmin == αxmin ? 1 : trunc(Int, ixmin + 1)
        imax = αmax == αxmax ? nx - 1 : trunc(Int, ixmax)
    else
        imin = αmax == αxmax ? 1 : trunc(Int, ixmax + 1)
        imax = αmin == αxmin ? nx - 2 : trunc(Int, ixmin)
    end
    return imin, imax
end


function initialize(sid::Siddon)
    αxmin, αxmax, αymin, αymax, αzmin, αzmax, αmin, αmax = get_α_minmax(sid)
    ixmin, jxmin, kxmin = get_φ(αmin; sid)
    ixmax, jxmax, kxmax = get_φ(αmax; sid)
    imin, imax = get_idx_minmax(αmin, αmax, αxmin, αxmax, ixmin, ixmax, sid.origin[1], sid.target[1], sid.dims[1])
    jmin, jmax = get_idx_minmax(αmin, αmax, αymin, αymax, jxmin, jxmax, sid.origin[2], sid.target[2], sid.dims[2])
    kmin, kmax = get_idx_minmax(αmin, αmax, αzmin, αzmax, kxmin, kxmax, sid.origin[3], sid.target[3], sid.dims[3])
    return αmin, αmax, imin, imax, jmin, jmax, kmin, kmax
end


# TODO: When to convert to CartesianIndex?
function get_voxel_idx(α::Float64; sid::Siddon)
    xidxs = get_φ(α; sid)
    idxs = trunc.(Int, xidxs) .+ 1
    return idxs
    # return CartesianIndex(idxs...)
end


function (sid::Siddon)(volume)

    # Get update conditions
    iu = sid.origin[1] ≤ sid.target[1] ? 1 : -1
    ju = sid.origin[2] ≤ sid.target[2] ? 1 : -1
    ku = sid.origin[3] ≤ sid.target[3] ? 1 : -1
    update_ijk = [iu, ju, ku]
    update_α = @. sid.spacing / abs(sid.target - sid.origin)

    # Initialize the loop
    αmin, αmax, imin, imax, jmin, jmax, kmin, kmax = initialize(sid)
    αcurr = αmin

    steps = get_α(imin, jmin, kmin; sid)  # Get potential next steps in xyz planes
    αnext, idx = findmin(steps)  # Find the smallest step (i.e., the next plane)

    αmid = (αcurr + αnext) / 2
    voxel = get_voxel_idx(αmid; sid)  # Get the voxel at the midpoint between αcurr and αnext

    step_len = αnext - αcurr
    d12 = @views step_len * volume[voxel...]
    αcurr = αnext

    # Loop over all voxels that the ray passes through
    ⪅(a, b) = (a < b) && !(a ≈ b)  # Test if a is less than b, but not approximately equal to b
    while αcurr ⪅ αmax             # Necessary because floating point errors means that αcurr ≈ αmax at some point 
        voxel[idx] += update_ijk[idx]
        steps[idx] += update_α[idx]
        αnext, idx = findmin(steps)
        step_len = αnext - αcurr
        d12 += @views step_len * volume[voxel...]
        αcurr = αnext
    end
    return d12
end


function siddon(ray::Ray, ct::CT)
    spacing = [ct.ΔX, ct.ΔY, ct.ΔZ]
    isocenter = [ct.X₀, ct.Y₀, ct.Z₀]
    dims = size(ct.volume) .+ 1
    sid = Siddon(
        Array(ray.origin),
        Array(ray.target),
        spacing,
        isocenter,
        dims
    )
    volume = ct.volume[end:-1:1, :, :] |> vol -> permutedims(vol, (2, 1, 3))  # Reverse the rows and swap rows/cols
    radiologic_path_length = sid(volume)
    return length(ray) * radiologic_path_length
end
