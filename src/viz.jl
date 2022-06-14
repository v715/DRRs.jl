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
    nx, ny, nz = size(CT.volume)
    xs = 0:CT.ΔX:(nx-1)*CT.ΔX
    ys = 0:CT.ΔY:(ny-1)*CT.ΔY
    zs = 0:CT.ΔZ:(nz-1)*CT.ΔZ
    return [
        get_slice_x(x; CT.volume, xs, ys, zs, ctkwargs...),
        get_slice_y(y; CT.volume, xs, ys, zs, ctkwargs...),
        get_slice_z(z; CT.volume, xs, ys, zs, ctkwargs...),
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

    Make a Plotly trace for the detector plane with a 3D mesh.
"""
function plot_detector(detector::Detector)
    plane = make_plane(detector)
    pts = plane[[1, 201, end - 200, end]]
    return mesh3d(
        name="Detector",
        x=[pt[1] for pt in pts],
        y=[pt[2] for pt in pts],
        z=[pt[3] for pt in pts],
        i=[0, 1],
        j=[1, 2],
        k=[2, 3],
    )
end


"""
plot_rays

    Make a trace for the rays using multiple 3D lineplots.
    The rays are subsampled (every 10th ray) to reduce overhead.
"""
function plot_rays(camera::Camera, detector::Detector)
    function plot_ray(ray::DRRs.Ray)
        origin, dest = ray.origin, ray.direction
        return scatter(
            x=[origin[1], dest[1]],
            y=[origin[2], dest[2]],
            z=[origin[3], dest[3]],
            type="scatter3d",
            mode="lines",
            line=attr(color="darkblue", width=2),
            opacity=0.25,
            showlegend=false,
        )
    end
    return [plot_ray(ray) for ray in get_rays(camera, detector)[1:20:201, 1:20:201]][:]
end


"""
plot(ct::CT, camera::Camera, detector::Detector)

    Overload the plot function to render the projector geometry.
"""
function plot(ct::CT, camera::Camera, detector::Detector)
    traces = [plot_rays(camera, detector)..., plot_ct(ct; ctkwargs...)..., plot_camera(camera), plot_detector(detector)]
    layout = Layout(scene=attr(
        xaxis=attr(range=[-500, 1100]),
        yaxis=attr(range=[-500, 1100]),
        zaxis=attr(range=[-500, 1100]),
        aspectratio=(x=1, y=1, z=1)
    ))
    plot(traces, layout)
end


"""
plot_drr

    Plot the geometry and the resulting DRR.
    NOTE: Very slow.
"""
function plot_drr(vol, ΔX, ΔY, ΔZ, camera::Camera, detector::Detector)

    fig = make_subplots(
        rows=1, cols=2,
        specs=[Spec(kind="scene") Spec(kind="xy")]
    )

    traces = [plot_rays(camera, detector)..., plot_ct(vol, ΔX, ΔY, ΔZ; ctkwargs...)..., plot_camera(camera), plot_detector(detector)]
    layout = Layout(scene=attr(
        xaxis=attr(range=[-500, 1100]),
        yaxis=attr(range=[-500, 1100]),
        zaxis=attr(range=[-500, 1100]),
        aspectratio=(x=1, y=1, z=1)
    ))
    for trace in traces
        add_trace!(fig, trace, row=1, col=1)
    end

    p = heatmap(z=DRR(vol, ΔX, ΔY, ΔZ, detector, camera, 0.1, trilinear), colorscale="Greys")
    add_trace!(fig, p, row=1, col=2)

    relayout!(fig)
    relayout!(layout)
    fig
end
