function DRR(ct::CT, detector::Detector, camera::Camera, spacing::Float64, sampling_function)

    # Make the rays in the projector
    projector = get_rays(camera, detector)

    # Get the spatial dimensions of the volume
    nx, ny, nz = size(CT.volume)
    xs = 0:CT.ΔX:(nx-1)*CT.ΔX
    ys = 0:CT.ΔY:(ny-1)*CT.ΔY
    zs = 0:CT.ΔZ:(nz-1)*CT.ΔZ

    # Trace rays through the voxel grid
    t = 0:spacing:1
    drr = [sampler(ray, t, sampling_function; ct.volume, xs, ys, zs) for ray in projector]
    return drr

end


function sampler(ray, t, sampling_function; volume, xs, ys, zs)
    out = trace.(t; ray) .|> xyz -> sampling_function(xyz...; volume, xs, ys, zs)
    return sum(out) / length(out)
end