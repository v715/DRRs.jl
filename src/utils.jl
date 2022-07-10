using LinearAlgebra
using StaticArrays


"""
    get_basis

Given a new vector, construct an orthonormal basis.
"""
function get_basis(a::MVector{3,Float64})
    w = normalize(a)
    t = get_noncollinear_vector(w) |> normalize
    u = cross(t, w)
    v = cross(u, w)
    return u, v, w
end


function get_noncollinear_vector(w)
    t = deepcopy(w)
    i = argmin(abs.(w))
    t[i] = 1
    return t
end


"""
    adjsum

Add adjacent pairs of elements in an array. Summation analog of `diff`.
"""
function adjsum(arr)
    arr1 = @view arr[begin:end-1]
    arr2 = @view arr[begin+1:end]
    arr1 + arr2
end