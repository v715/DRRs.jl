using DRRs

# Get XYZ coordinates of ray at distance α from source
# use trace(α; ray)

# Get the spatial coordinate of a plane, assuming the image isocenter is the origin
# Planes are indexed from 0:n where n is the number of slices in a given direction
X_plane(i::Int64, ct::CT; X₀::Float64=0.0) = X₀ + i * ct.ΔX
Y_plane(j::Int64, ct::CT; Y₀::Float64=0.0) = Y₀ + j * ct.ΔY
Z_plane(k::Int64, ct::CT; Z₀::Float64=0.0) = Z₀ + k * ct.ΔZ


# Get the α values for each plane
function get_α(i::Int64, j::Int64, k::Int64; ray::Ray, ct::CT)
    x1, y1, z1 = ray.origin
    x2, y2, z2 = ray.target
    αx = (X_plane(i, ct) - x1) / (x2 - x1)
    αy = (Y_plane(j, ct) - y1) / (y2 - y1)
    αz = (Z_plane(k, ct) - z1) / (z2 - z1)
    return αx, αy, αz
end

get_αx(i::Int64; ray::Ray, ct::CT) = (X_plane(i, ct) - ray.origin[1]) / (ray.target[1] - ray.origin[1])
get_αy(j::Int64; ray::Ray, ct::CT) = (Y_plane(j, ct) - ray.origin[2]) / (ray.target[2] - ray.origin[2])
get_αz(k::Int64; ray::Ray, ct::CT) = (Z_plane(k, ct) - ray.origin[3]) / (ray.target[3] - ray.origin[3])


function get_α_minmax(ray::Ray, ct::CT)
    nx, ny, nz = size(ct.volume)
    αx0, αy0, αz0 = get_α(0, 0, 0; ray, ct)
    αx1, αy1, αz1 = get_α(nx, ny, nz; ray, ct)
    αmin = max(0, min(αx0, αx1), min(αy0, αy1), min(αz0, αz1))
    αmax = min(1, max(αx0, αx1), max(αy0, αy1), max(αz0, αz1))
    return αmin, αmax
end


# Get the first and last intersections of a ray with orthogonal planes
function get_xidx_minmax(x1::Float64, x2::Float64, ΔX::Float64, nx::Int64, αmin::Float64, αmax::Float64)
    if x2 - x1 >= 0
        imin = nx - (X_plane(nx, ct) - αmin * (x2 - x1) - x1) / ΔX
        imax = 1 + (x1 + αmax * (x2 - x1) - X_plane(0, ct)) / ΔX
    else
        imin = nx - (X_plane(nx, ct) - αmax * (x2 - x1) - x2) / ΔX
        imax = 1 + (x1 + αmin * (x2 - x1) - X_plane(0, ct)) / ΔX
    end
    return Int(floor(imin)), Int(floor(imax))
end


function get_idx_minmax(αmin::Float64, αmax::Float64, ray::Ray, ct::CT)
    x1, y1, z1 = ray.origin
    x2, y2, z2 = ray.target
    nx, ny, nz = size(ct.volume)
    imin, imax = get_xidx_minmax(x1, x2, ct.ΔX, nx, αmin, αmax)
    jmin, jmax = get_xidx_minmax(y1, y2, ct.ΔY, ny, αmin, αmax)
    kmin, kmax = get_xidx_minmax(z1, z2, ct.ΔZ, nz, αmin, αmax)
    return imin, imax, jmin, jmax, kmin, kmax
end


# Get the list of α's
function get_merged_α(ray::Ray, ct::CT)
    αmin, αmax = get_α_minmax(ray, ct)
    imin, imax, jmin, jmax, kmin, kmax = get_idx_minmax(αmin, αmax, ray, ct)
    α̲x = get_αx.(imin:imax; ray, ct)
    α̲y = get_αy.(jmin:jmax; ray, ct)
    α̲z = get_αz.(kmin:kmax; ray, ct)
    α̲ = vcat(αmin, α̲x, α̲y, α̲z, αmax) |> sort
    return α̲
end


# Get the voxel characterized by two adjacent α values
function get_weighted_voxel(m::Int64, α̲::Vector{Float64}, ct::CT; X₀::Float64=0.0, Y₀::Float64=0.0, Z₀::Float64=0.0)
    αmid = (α̲[m] + α̲[m-1]) / 2
    x1, y1, z1 = ray.origin
    x2, y2, z2 = ray.target
    i = floor((x1 + αmid * (x2 - x1) - X₀) / (ct.ΔX))
    j = floor((y1 + αmid * (y2 - y1) - Y₀) / (ct.ΔY))
    k = floor((z1 + αmid * (z2 - z1) - Z₀) / (ct.ΔZ))
    return (α̲[m+1] - α̲[m]) * ct.volume[i, j, k]
end


function siddon(ray::Ray, ct::CT)
    α̲ = get_merged_α(ray, ct)
    n = length(α̲)
    radiologic_path_length = 0.0
    for m in 1:n
        radiologic_path_length += get_weighted_voxel(m, α̲, ct)
    end
    d12 = norm(ray.target - ray.origin)
    return d12 * radiologic_path_length
end