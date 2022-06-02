module DRRs

export read_dicom
include("io.jl")

export Camera, Detector, make_plane, get_rays, trace
include("camera.jl")

end
