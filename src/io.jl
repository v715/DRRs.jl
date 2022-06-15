using DICOM


struct CT
    volume::Array{T,3} where {T<:Real}
    ΔX::Float64
    ΔY::Float64
    ΔZ::Float64
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

    return CT(volume, ΔX, ΔY, ΔZ)

end