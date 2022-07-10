using Base.Threads


function sampler(ct::CT, camera::Camera, detector::Detector; ray2pix::Symbol, spacing::Float64)

    # Parse the sampling function specification
    if ray2pix == :sample
        sampling_function = sample
    elseif ray2pix == :trilinear
        sampling_function = trilinear
    else
        error("Unkown ray2pix method: :$ray2pix. Must be in {:sample, :trilinear}.")
    end

    # Get the model parameters
    t = 0.0:spacing:1.0
    nx, ny, nz = size(ct.volume)
    xs = 0:ct.ΔX:(nx-1)*ct.ΔX
    ys = 0:ct.ΔY:(ny-1)*ct.ΔY
    zs = 0:ct.ΔZ:(nz-1)*ct.ΔZ

    # Make the DRR
    projector = make_xrays(camera, detector)
    drr = zeros(Float64, detector.height, detector.width)

    @threads for i in eachindex(projector)
        @inbounds ray = projector[i]
        out = trace.(t; ray) .|> xyz -> sampling_function(xyz...; ct.volume, xs, ys, zs)
        intensity = sum(out) / length(out)  # Normalize by the number of samples
        @inbounds drr[i] = intensity
    end

    return drr
end


function siddon(ct::CT, camera::Camera, detector::Detector)
    projector = make_xrays(camera, detector)
    drr = zeros(Float64, detector.height, detector.width)
    @threads for i in eachindex(projector)
        @inbounds drr[i] = siddon(projector[i], ct)
    end
    return drr
end
