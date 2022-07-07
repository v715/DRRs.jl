function DRR(ct::CT, detector::Detector, camera::Camera; sampler::String, spacing::Union{Float64,Nothing}=nothing)

    # Make the rays in the projector
    projector = make_xrays(camera, detector)

    # Get the spatial dimensions of the volume
    nx, ny, nz = size(ct.volume)
    xs = 0:ct.ΔX:(nx-1)*ct.ΔX
    ys = 0:ct.ΔY:(ny-1)*ct.ΔY
    zs = 0:ct.ΔZ:(nz-1)*ct.ΔZ

    # Trace rays through the voxel grid
    if sampler == "sample" || sampler == "trilinear"
        t = 0:spacing:1
        sampling_function = sampler == "sample" ? sample : trilinear
        drr = [sampler(ray, t, sampling_function; ct.volume, xs, ys, zs) for ray in projector]
    else if sampler == "siddon"
        drr = [siddon(ray, ct) for ray in projector]
    else
        error("Unknown sampler")
    end

    return drr

end


function sampler(ray, t, sampling_function; volume, xs, ys, zs)
    out = trace.(t; ray) .|> xyz -> sampling_function(xyz...; volume, xs, ys, zs)
    return sum(out) / length(out)
end