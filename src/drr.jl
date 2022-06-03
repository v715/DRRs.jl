function DRR(volume, ΔX, ΔY, ΔZ, detector, camera, spacing, sampling_function)

    # Make the rays in the projector
    projector = get_rays(camera, detector)

    # Get the spatial dimensions of the volume
    nx, ny, nz = size(volume)
    xs = 0:ΔX:(nx-1)*ΔX
    ys = 0:ΔY:(ny-1)*ΔY
    zs = 0:ΔZ:(nz-1)*ΔZ

    # Trace rays through the voxel grid
    t = 0:spacing:1
    drr = [sampler(ray, t, sampling_function; volume, xs, ys, zs) for ray in projector]
    return drr

end


function sampler(ray, t, sampling_function; volume, xs, ys, zs)
    out = trace.(t; ray) .|> xyz -> sampling_function(xyz...; volume, xs, ys, zs)
    return sum(out) / length(out)
end