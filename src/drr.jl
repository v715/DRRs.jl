function DRR(ct::CT, detector::Detector, camera::Camera; ray2pix::Symbol, spacing::Union{Float64,Nothing}=nothing)

    # Make the rays in the projector
    projector = make_xrays(camera, detector)

    # Trace rays through the voxel grid
    if ray2pix == :sample || ray2pix == :trilinear
        t = 0:spacing:1
        nx, ny, nz = size(ct.volume)
        xs = 0:ct.ΔX:(nx-1)*ct.ΔX
        ys = 0:ct.ΔY:(ny-1)*ct.ΔY
        zs = 0:ct.ΔZ:(nz-1)*ct.ΔZ
        sampling_function = ray2pix == :sample ? sample : trilinear
        drr = [sampler(ray, t, sampling_function; ct.volume, xs, ys, zs) for ray in projector]
    elseif ray2pix == :siddon
        drr = [siddon(ray, ct) for ray in projector]
    else
        error("Unknown sampling function: ", ray2pix)
    end

    return drr

end


function sampler(ray, t, sampling_function; volume, xs, ys, zs)
    out = trace.(t; ray) .|> xyz -> sampling_function(xyz...; volume, xs, ys, zs)
    return sum(out) / length(out)
end