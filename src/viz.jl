using PlotlyJS

"""
get_slice_{x,y,z}

    Get a slice of a CT volume along the x (coronal), y (sagittal) or z (transverse) axis.
    The slice is plotted as an isosurface in Plotly.
"""
function get_slice_x(i::Int64; ctkwargs...)
    X, Y, Z = mgrid([xs[i]], ys, zs)
    value = vol[i, :, :]
    return isosurface(x=X[:], y=Y[:], z=Z[:], value=value[:], name="Coronal"; ctkwargs...)
end


function get_slice_y(i::Int64; ctkwargs...)
    X, Y, Z = mgrid(xs, [ys[i]], zs)
    value = vol[:, i, :]
    return isosurface(x=X[:], y=Y[:], z=Z[:], value=value[:], name="Sagittal"; ctkwargs...)
end


function get_slice_z(i::Int64; ctkwargs...)
    X, Y, Z = mgrid(xs, ys, [zs[i]])
    value = vol[:, :, i]
    return isosurface(x=X[:], y=Y[:], z=Z[:], value=value[:], name="Transverse"; ctkwargs...)
end