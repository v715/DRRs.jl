using DICOM


struct CT
    volume::Array{T,3} where {T<:Real}  # 3D image volume
    ΔX::Float64                         # X-direction spacing
    ΔY::Float64                         # Y-direction spacing
    ΔZ::Float64                         # Z-direction spacing
    X₀::Float64                         # X-coordinate isocenter
    Y₀::Float64                         # Y-coordinate isocenter
    Z₀::Float64                         # Z-coordinate isocenter
end


function read_dicom(path::String)

    dcm_data_array = dcmdir_parse(path)

    # Get the dimensions of the volume (image coordinates)
    n_slices = length(dcm_data_array)
    nx, ny = size(dcm_data_array[1][tag"PixelData"])
    volume = Array{Int16}(undef, (nx, ny, n_slices))

    # Read the volume and get the x-, y-, and z-directional voxel spacing
    ΔX, ΔY = dcm_data_array[1][tag"PixelSpacing"]
    Zs = Array{Float64}(undef, n_slices)

    for i in 1:n_slices
        slice = dcm_data_array[end-i+1]  # Read the slices in reverse order
        volume[:, :, i] = slice[tag"PixelData"][end:-1:1, :]  # Reverse the slice
        Zs[i] = slice[tag"ImagePositionPatient"][3]
    end

    ΔZ = diff(Zs)
    ΔZ = all(z ≈ ΔZ[1] for z in ΔZ) ? ΔZ[1] : ΔZ

    # Get the origin of the volume (image coordinates)
    # For now, hard-code the origin to be the center of the volume
    # TODO: Make this settable
    X₀, Y₀, Z₀ = 0.0, 0.0, 0.0

    return CT(volume, ΔX, ΔY, ΔZ, X₀, Y₀, Z₀)

end