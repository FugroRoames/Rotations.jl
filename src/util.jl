"""
    perpendicular_vector(vec)

Compute a vector perpendicular to `vec` by switching the two elements with
largest absolute value, flipping the sign of the second largest, and setting the
remaining element to zero.
"""
function perpendicular_vector(vec::SVector{3})
    T = eltype(vec)

    # find indices of the two elements of vec with the largest absolute values:
    absvec = abs.(vec)
    ind1 = argmax(absvec) # index of largest element
    tmin = typemin(T)
    @inbounds absvec2 = @SVector [ifelse(i == ind1, tmin, absvec[i]) for i = 1 : 3] # set largest element to typemin(T)
    ind2 = argmax(absvec2) # index of second-largest element

    # perp[ind1] = -vec[ind2], perp[ind2] = vec[ind1], set remaining element to zero:
    @inbounds perpind1 = -vec[ind2]
    @inbounds perpind2 = vec[ind1]
    tzero = zero(T)
    perp = @SVector [ifelse(i == ind1, perpind1, ifelse(i == ind2, perpind2, tzero)) for i = 1 : 3]
end

@noinline length_error(v, len) =
    throw(DimensionMismatch("Expected length $len, got length $(length(v))"))

@inline function check_length(v, len)
    if length(v) != len
        length_error(v, len)
    end
end

function skew(v::AbstractVector)
    check_length(v, 3)
    @SMatrix [0   -v[3]  v[2];
              v[3] 0    -v[1];
             -v[2] v[1]  0]
end

function vee(S::AbstractMatrix)
    return @SVector [S[3,2], S[1,3], S[2,1]]
end

"""
The element type for a rotation matrix with a given angle type is composed of
trigonometric functions of that type.
"""
Base.@pure rot_eltype(angle_type) = typeof(sin(zero(angle_type)))