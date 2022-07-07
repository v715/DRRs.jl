using PlotlyJS

import PlotlyJS: plot


"""
    get_slice_{x,y,z}

Get a slice of a CT volume along the x (coronal), y (sagittal) or z (transverse) axis. 
The slice is plotted as an isosurface in Plotly.
"""
function get_slice_x(i::Int64; vol, xs, ys, zs, ctkwargs...)
    X, Y, Z = mgrid([xs[i]], ys, zs)
    value = vol[i, :, :]
    return isosurface(x=X[:], y=Y[:], z=Z[:], value=value[:], name="Coronal"; ctkwargs...)
end


function get_slice_y(i::Int64; vol, xs, ys, zs, ctkwargs...)
    X, Y, Z = mgrid(xs, [ys[i]], zs)
    value = vol[:, i, :]
    return isosurface(x=X[:], y=Y[:], z=Z[:], value=value[:], name="Sagittal"; ctkwargs...)
end


function get_slice_z(i::Int64; vol, xs, ys, zs, ctkwargs...)
    X, Y, Z = mgrid(xs, ys, [zs[i]])
    value = vol[:, :, i]
    return isosurface(x=X[:], y=Y[:], z=Z[:], value=value[:], name="Transverse"; ctkwargs...)
end

"""
    plot_ct

Make a Plotly trace for the CT volume.
`ctkwargs...` is a dictionary of keyword arguments to pass to the isosurface function.
"""
function plot_ct(ct::CT; x::Int64=256, y::Int64=256, z::Int64=66, ctkwargs...)
    # Get the coordinate spacings
    nx, ny, nz = size(ct.volume)
    xs = 0:ct.ΔX:(nx-1)*ct.ΔX
    ys = 0:ct.ΔY:(ny-1)*ct.ΔY
    zs = 0:ct.ΔZ:(nz-1)*ct.ΔZ
    return [
        get_slice_x(x; vol=ct.volume, xs, ys, zs, ctkwargs..., showscale=false),
        get_slice_y(y; vol=ct.volume, xs, ys, zs, ctkwargs..., showscale=false),
        get_slice_z(z; vol=ct.volume, xs, ys, zs, ctkwargs..., showscale=true),
    ]
end


"""
    plot_camera

Make a Plotly trace for the camera center a single point in a 3D scatterplot.
"""
function plot_camera(camera::Camera)
    return scatter(
        name="Camera",
        x=[camera.center[1]],
        y=[camera.center[2]],
        z=[camera.center[3]],
        mode="markers",
        type="scatter3d",
        showlegend=false,
    )
end


"""
    plot_detector

Make two traces: one for the detector plane with a 3D mesh, and another
for the rays emanating from the camera. The rays are created using multiple
3D lineplots. The rays are subsampled (every 10th ray) to reduce overhead.
"""
function plot_detector(camera::Camera, detector::Detector; skips::Int64=20)

    # Make the detector plane
    rays = make_xrays(camera, detector)
    pts = rays[[1, detector.height, end - detector.height + 1, end]]
    pts = [pt.target for pt in pts]
    plane = mesh3d(
        name="Detector",
        x=[pt[1] for pt in pts],
        y=[pt[2] for pt in pts],
        z=[pt[3] for pt in pts],
        i=[0, 0],
        j=[1, 3],
        k=[3, 2],
    )

    # Make the rays
    function plot_ray(ray::DRRs.Ray)
        origin, target = ray.origin, ray.target
        return scatter(
            x=[origin[1], target[1]],
            y=[origin[2], target[2]],
            z=[origin[3], target[3]],
            type="scatter3d",
            mode="lines",
            line=attr(color="darkblue", width=2),
            opacity=0.25,
            showlegend=false,
        )
    end
    rays = [plot_ray(ray) for ray in rays[1:skips:detector.height, 1:skips:detector.width]][:]

    return [plane, rays...]

end


"""
    plot_image(drr::Matrix{T}; name::String="DRR") where {T<:Real}

Plot a DRR as a 2D heatmap.
"""
function plot_image(drr::Matrix{T}; name::String="DRR") where {T<:Real}
    trace = heatmap(z=drr, colorscale="Greys", name=name)
    height, width = size(drr)
    layout = Layout(
        showlegend=false,
        width=width,
        height=height,
        autosize=false
    )
    plot(trace, layout)
end


"""
    plot(ct::CT, camera::Camera, detector::Detector)

Overload the plot function to render the projector geometry.
"""
function plot(ct::CT, camera::Camera, detector::Detector; axlim=[-500, 1100], ctkwargs...)
    traces = [
        plot_ct(ct; ctkwargs...)...,
        plot_camera(camera),
        plot_detector(camera, detector)...,
    ]
    layout = Layout(scene=attr(
        xaxis=attr(range=axlim),
        yaxis=attr(range=axlim),
        zaxis=attr(range=axlim),
        aspectratio=(x=1, y=1, z=1)
    ))
    plot(traces, layout)
end


"""
    plot(ct::CT, camera::Camera, detector::Detector, drr::Matrix{Float64})

Plot the geometry and the resulting DRR. 
NOTE: Very slow.
"""
function plot(ct::CT, camera::Camera, detector::Detector, drr::Matrix{T}; axlim=[-500, 1100], ctkwargs...) where {T<:Real}

    fig = make_subplots(
        rows=1, cols=2,
        specs=[Spec(kind="scene") Spec(kind="xy")]
    )

    # Plot the 3D geometric setup
    traces = [
        plot_ct(ct; ctkwargs...)...,
        plot_camera(camera),
        plot_detector(camera, detector)...,
    ]
    layout = Layout(scene=attr(
        xaxis=attr(range=axlim),
        yaxis=attr(range=axlim),
        zaxis=attr(range=axlim),
        aspectratio=(x=1, y=1, z=1)
    ))
    for trace in traces
        add_trace!(fig, trace, row=1, col=1)
    end

    # Plot the DRR
    p = heatmap(z=drr, colorscale="Greys")
    add_trace!(fig, p, row=1, col=2)

    relayout!(fig)
    relayout!(layout)
    fig
end
