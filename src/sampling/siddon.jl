# Get the spatial coordinate of a plane, assuming the image isocenter is the origin
# Planes are indexed from 1:(n+1) where n is the number of voxels in a given direction
X_plane(i::Int64, ct::CT) = ct.X₀ + (i - 1) * ct.ΔX
Y_plane(j::Int64, ct::CT) = ct.Y₀ + (j - 1) * ct.ΔY
Z_plane(k::Int64, ct::CT) = ct.Z₀ + (k - 1) * ct.ΔZ


# Get the α values for a plane
get_αx(i::Int64; ray::Ray, ct::CT) = (X_plane(i, ct) - ray.origin[1]) / (ray.target[1] - ray.origin[1])
get_αy(j::Int64; ray::Ray, ct::CT) = (Y_plane(j, ct) - ray.origin[2]) / (ray.target[2] - ray.origin[2])
get_αz(k::Int64; ray::Ray, ct::CT) = (Z_plane(k, ct) - ray.origin[3]) / (ray.target[3] - ray.origin[3])

function get_α(i::Int64, j::Int64, k::Int64; ray::Ray, ct::CT)
    αx = get_αx(i; ray, ct)
    αy = get_αy(j; ray, ct)
    αz = get_αz(k; ray, ct)
    return αx, αy, αz
end


# Get the min/max α values for a ray
function get_α_minmax(ray::Ray, ct::CT)
    nx, ny, nz = size(ct.volume) .+ 1  # Number of x,y,z planes
    αx0, αy0, αz0 = get_α(1, 1, 1; ray, ct)
    αx1, αy1, αz1 = get_α(nx, ny, nz; ray, ct)
    αmin = max(0, min(αx0, αx1), min(αy0, αy1), min(αz0, αz1))
    αmax = min(1, max(αx0, αx1), max(αy0, αy1), max(αz0, αz1))
    return αmin, αmax
end


# Get the first and last intersections of a ray with orthogonal planes
function get_idx_minmax(lookup::Function, p1::Float64, p2::Float64, Δ::Float64, n::Int64; αmin::Float64, αmax::Float64, ct::CT)
    if p2 - p1 ≥ 0
        pmin = n - (lookup(n, ct) - αmin * (p2 - p1) - p1) / Δ |> ceil |> Int
        pmax = 1 + (p1 + αmax * (p2 - p1) - lookup(1, ct)) / Δ |> floor |> Int
    else
        pmin = n - (lookup(n, ct) - αmax * (p2 - p1) - p1) / Δ |> ceil |> Int
        pmax = 1 + (p1 + αmin * (p2 - p1) - lookup(1, ct)) / Δ |> floor |> Int
    end
    return pmin, pmax
end

function get_idx_minmax(αmin::Float64, αmax::Float64, ray::Ray, ct::CT)
    x1, y1, z1 = ray.origin
    x2, y2, z2 = ray.target
    nx, ny, nz = size(ct.volume)
    imin, imax = get_idx_minmax(X_plane, x1, x2, ct.ΔX, nx; αmin, αmax, ct)
    jmin, jmax = get_idx_minmax(Y_plane, y1, y2, ct.ΔY, ny; αmin, αmax, ct)
    kmin, kmax = get_idx_minmax(Z_plane, z1, z2, ct.ΔZ, nz; αmin, αmax, ct)
    return imin, imax, jmin, jmax, kmin, kmax
end


# Get the merged list of α's
function get_merged_αs(ray::Ray, ct::CT)
    αmin, αmax = get_α_minmax(ray, ct)
    imin, imax, jmin, jmax, kmin, kmax = get_idx_minmax(αmin, αmax, ray, ct)
    α̲x = get_αx.(imin:imax; ray, ct)
    α̲y = get_αy.(jmin:jmax; ray, ct)
    α̲z = get_αz.(kmin:kmax; ray, ct)
    α̲ = vcat(α̲x, α̲y, α̲z) |> sort
    return α̲
end


# Get the voxel characterized by two adjacent α values
function get_weighted_voxel(ray::Ray, m::Int64, α̲::Vector{Float64}, ct::CT)
    αmid = (α̲[m] + α̲[m-1]) / 2
    x1, y1, z1 = ray.origin
    x2, y2, z2 = ray.target
    i = 1 + (x1 + αmid * (x2 - x1) - ct.X₀) / ct.ΔX |> floor |> Int
    j = 1 + (y1 + αmid * (y2 - y1) - ct.Y₀) / ct.ΔY |> floor |> Int
    k = 1 + (z1 + αmid * (z2 - z1) - ct.Z₀) / ct.ΔZ |> floor |> Int
    return (α̲[m] - α̲[m-1]) * ct.volume[i, j, k]
end


# Main function to get pixel value from ray
function siddon(ray::Ray, ct::CT)
    α̲ = get_merged_αs(ray, ct)
    n = length(α̲)
    radiologic_path_length = 0.0
    for m in 2:n
        radiologic_path_length += get_weighted_voxel(ray, m, α̲, ct)
    end
    return length(ray) * radiologic_path_length
end