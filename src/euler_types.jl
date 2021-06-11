# We have rotations along one, two or three axes (e.g. RotX, RotXY and RotXYZ).
# They compose together nicely, so the user can make Euler angles by:
#
#     RotX(θx) * RotY(θy) * RotZ(θz) -> RotXYZ(θx, θy, θz)
#
# and never get confused by the order of application, etc.

#########################
# Single axis rotations #
#########################

# The element type for a rotation matrix with a given angle type is composed of
# trigonometric functions of that type.
Base.@pure rot_eltype(angle_type) = typeof(sin(zero(angle_type)))

for axis in [:X, :Y, :Z]
    RotType = Symbol("Rot" * string(axis))
    @eval begin
        struct $RotType{T} <: Rotation{3,T}
            theta::T
            $RotType{T}(theta) where {T} = new{T}(theta)
            $RotType{T}(r::$RotType) where {T} = new{T}(r.theta)
        end

        @inline function $RotType(theta)
            $RotType{rot_eltype(typeof(theta))}(theta)
        end
        @inline $RotType(r::$RotType{T}) where {T} = $RotType{T}(r)

        @inline $RotType(::NTuple{9}) = error("Cannot construct a cardinal axis rotation from a matrix")
        @inline $RotType{T}(::NTuple{9}) where T = error("Cannot construct a cardinal axis rotation from a matrix")

        @inline Base.:*(r1::$RotType, r2::$RotType) = $RotType(r1.theta + r2.theta)

        @inline Base.inv(r::$RotType) = $RotType(-r.theta)

        # define null rotations for convenience
        @inline Base.one(::Type{$RotType}) = $RotType(0.0)
        @inline Base.one(::Type{$RotType{T}}) where {T} = $RotType{T}(zero(T))

        params(r::$RotType) = SVector{1}(r.theta)
    end
end

function Base.rand(::Type{R}) where R <: Union{RotX,RotY,RotZ}
    T = eltype(R)
    if T == Any
        T = Float64
    end

    return R(2*pi*rand(T))
end


"""
    struct RotX{T} <: Rotation{3,T}
    RotX(theta)

A 3×3 rotation matrix which represents a rotation by `theta` about the X axis.
"""
RotX

@inline function Base.getindex(r::RotX{T}, i::Int) where T
    T2 = Base.promote_op(sin, T)
    if i == 1
        one(T2)
    elseif i < 5
        zero(T2)
    elseif i == 5
        cos(r.theta)
    elseif i == 6
        sin(r.theta)
    elseif i == 7
        zero(T2)
    elseif i == 8
        -sin(r.theta)
    elseif i == 9
        cos(r.theta)
    else
        throw(BoundsError(r,i))
    end
end

@inline function Base.Tuple(r::RotX{T}) where T
    s, c = sincos(r.theta)
    o = one(s)
    z = zero(s)
    (o,  z,  z,   # transposed representation
     z,  c,  s,
     z, -s,  c)
end

@inline function Base.:*(r::RotX, v::StaticVector)
    if length(v) != 3
        throw("Dimension mismatch: cannot rotate a vector of length $(length(v))")
    end

    st, ct = sincos(r.theta)
    T = Base.promote_op(*, typeof(st), eltype(v))
    return similar_type(v,T)(v[1],
                             v[2] * ct - v[3] * st,
                             v[3] * ct + v[2] * st)
end


"""
    struct RotY{T} <: Rotation{3,T}
    RotY(theta)

A 3×3 rotation matrix which represents a rotation by `theta` about the Y axis.
"""
RotY

@inline function Base.getindex(r::RotY{T}, i::Int) where T
    T2 = Base.promote_op(sin, T)
    if i == 1
        cos(r.theta)
    elseif i == 2
        zero(T2)
    elseif i == 3
        -sin(r.theta)
    elseif i == 4
        zero(T2)
    elseif i == 5
        one(T2)
    elseif i == 6
        zero(T2)
    elseif i == 7
        sin(r.theta)
    elseif i == 8
        zero(T2)
    elseif i == 9
        cos(r.theta)
    else
        throw(BoundsError(r,i))
    end
end

@inline function Base.Tuple(r::RotY{T}) where T
    s, c = sincos(r.theta)
    o = one(s)
    z = zero(s)
    (c,  z, -s,   # transposed representation
     z,  o,  z,
     s,  z,  c)
end

@inline function Base.:*(r::RotY, v::StaticVector)
    if length(v) != 3
        throw("Dimension mismatch: cannot rotate a vector of length $(length(v))")
    end

    st, ct = sincos(r.theta)
    T = Base.promote_op(*, typeof(st), eltype(v))
    return similar_type(v,T)(v[1] * ct + v[3] * st,
                             v[2],
                             v[3] * ct - v[1] * st)
end


"""
    struct RotZ{T} <: Rotation{3,T}
    RotZ(theta)

A 3×3 rotation matrix which represents a rotation by `theta` about the Z axis.
"""
RotZ

@inline function Base.getindex(r::RotZ{T}, i::Int) where T
    T2 = Base.promote_op(sin, T)
    if i == 1
        cos(r.theta)
    elseif i == 2
        sin(r.theta)
    elseif i == 3
        zero(T2)
    elseif i == 4
        -sin(r.theta)
    elseif i == 5
        cos(r.theta)
    elseif i < 9
        zero(T2)
    elseif i == 9
        one(T2)
    else
        throw(BoundsError(r,i))
    end
end

@inline function Base.Tuple(r::RotZ{T}) where T
    s, c = sincos(r.theta)
    o = one(s)
    z = zero(s)
    ( c, s, z,   # transposed representation
     -s, c, z,
      z, z, o)
end

@inline function Base.:*(r::RotZ, v::StaticVector)
    if length(v) != 3
        throw("Dimension mismatch: cannot rotate a vector of length $(length(v))")
    end

    st, ct = sincos(r.theta)
    T = Base.promote_op(*, typeof(st), eltype(v))
    return similar_type(v,T)(v[1] * ct - v[2] * st,
                             v[2] * ct + v[1] * st,
                             v[3])
end


################################################################################
################################################################################

######################
# Two axis rotations #
######################

for axis1 in [:X, :Y, :Z]
    Rot1Type = Symbol("Rot" * string(axis1))
    for axis2 in filter(axis -> axis != axis1, [:X, :Y, :Z])
        Rot2Type = Symbol("Rot" * string(axis2))
        RotType = Symbol("Rot" * string(axis1) * string(axis2))
        InvRotType = Symbol("Rot" * string(axis2) * string(axis1))

        @eval begin
            struct $RotType{T} <: Rotation{3,T}
                theta1::T
                theta2::T
                $RotType{T}(theta1, theta2) where {T} = new{T}(theta1, theta2)
                $RotType{T}(r::$RotType) where {T} = new{T}(r.theta1, r.theta2)
            end

            @inline function $RotType(theta1, theta2)
                ts = promote(theta1, theta2)
                $RotType{rot_eltype(eltype(ts))}(ts...)
            end
            @inline $RotType(r::$RotType{T}) where {T} = $RotType{T}(r)

            @inline function Base.getindex(r::$RotType{T}, i::Int) where T
                Tuple(r)[i] # Slow...
            end

            @inline (::Type{R})(t::NTuple{9}) where {R<:$RotType} = error("Cannot construct a two-axis rotation from a matrix")

            # Composing single-axis rotations to obtain a two-axis rotation:
            @inline Base.:*(r1::$Rot1Type, r2::$Rot2Type) = $RotType(r1.theta, r2.theta)

            # Composing single-axis rotations with two-axis rotations:
            @inline Base.:*(r1::$RotType, r2::$Rot2Type) = $RotType(r1.theta1, r1.theta2 + r2.theta)
            @inline Base.:*(r1::$Rot1Type, r2::$RotType) = $RotType(r1.theta + r2.theta1, r2.theta2)

            @inline Base.inv(r::$RotType) = $InvRotType(-r.theta2, -r.theta1)

            # define null rotations for convenience
            @inline Base.one(::Type{$RotType}) = $RotType(0.0, 0.0)
            @inline Base.one(::Type{$RotType{T}}) where {T} = $RotType{T}(zero(T), zero(T))

            params(r::$RotType) = SVector{2}(r.theta1, r.theta2)
        end
    end
end

function Base.rand(::Type{R}) where R <: Union{RotXY,RotYZ,RotZX, RotXZ, RotYX, RotZY}
    T = eltype(R)
    if T == Any
        T = Float64
    end

    # Not really sure what this distribution is, but it's also not clear what
    # it should be! rand(RotXY) *is* invariant to pre-rotations by a RotX and
    # post-rotations by a RotY...
    return R(2*pi*rand(T), 2*pi*rand(T))
end


"""
    struct RotXY{T} <: Rotation{3,T}
    RotXY(theta_x, theta_y)

A 3×3 rotation matrix which represents a rotation by `theta_y` about the Y axis,
followed by a rotation by `theta_x` about the X axis.
"""
RotXY

@inline function Base.Tuple(r::RotXY{T}) where T
    sinθ₁, cosθ₁ = sincos(r.theta1)
    sinθ₂, cosθ₂ = sincos(r.theta2)
    z = zero(sinθ₁)

    # transposed representation
    (cosθ₂,  sinθ₁*sinθ₂,    cosθ₁*-sinθ₂,
     z,      cosθ₁,          sinθ₁,
     sinθ₂,  -sinθ₁*cosθ₂,   cosθ₁*cosθ₂)
end

@inline function Base.:*(r::RotXY, v::StaticVector)
    if length(v) != 3
        throw("Dimension mismatch: cannot rotate a vector of length $(length(v))")
    end

    sinθ₁, cosθ₁ = sincos(r.theta1)
    sinθ₂, cosθ₂ = sincos(r.theta2)

    T = Base.promote_op(*, typeof(sinθ₁), eltype(v))

    return similar_type(v,T)(cosθ₂*v[1] + sinθ₂*v[3],
                             sinθ₁*sinθ₂*v[1] + cosθ₁*v[2] + -sinθ₁*cosθ₂*v[3],
                             cosθ₁*-sinθ₂*v[1] + sinθ₁*v[2] + cosθ₁*cosθ₂*v[3])
end


"""
    struct RotYX{T} <: Rotation{3,T}
    RotYX(theta_y, theta_x)

A 3×3 rotation matrix which represents a rotation by `theta_x` about the X axis,
followed by a rotation by `theta_y` about the Y axis.
"""
RotYX

@inline function Base.Tuple(r::RotYX{T}) where T
    sinθ₁, cosθ₁ = sincos(r.theta1)
    sinθ₂, cosθ₂ = sincos(r.theta2)
    z = zero(sinθ₁)

    # transposed representation
    (cosθ₁,        z,       -sinθ₁,
     sinθ₁*sinθ₂,  cosθ₂,   cosθ₁*sinθ₂,
     sinθ₁*cosθ₂,  -sinθ₂,  cosθ₁*cosθ₂)
end

@inline function Base.:*(r::RotYX, v::StaticVector)
    if length(v) != 3
        throw("Dimension mismatch: cannot rotate a vector of length $(length(v))")
    end

    sinθ₁, cosθ₁ = sincos(r.theta1)
    sinθ₂, cosθ₂ = sincos(r.theta2)

    T = Base.promote_op(*, typeof(sinθ₁), eltype(v))

    return similar_type(v,T)(cosθ₁*v[1] + sinθ₁*sinθ₂*v[2] + sinθ₁*cosθ₂*v[3],
                             cosθ₂*v[2] + -sinθ₂*v[3],
                             -sinθ₁*v[1] + cosθ₁*sinθ₂*v[2] + cosθ₁*cosθ₂*v[3])
end


"""
    struct RotXZ{T} <: Rotation{3,T}
    RotXZ(theta_x, theta_z)

A 3×3 rotation matrix which represents a rotation by `theta_z` about the Z axis,
followed by a rotation by `theta_x` about the X axis.
"""
RotXZ

@inline function Base.Tuple(r::RotXZ{T}) where T
    sinθ₁, cosθ₁ = sincos(r.theta1)
    sinθ₂, cosθ₂ = sincos(r.theta2)
    z = zero(sinθ₁)

    # transposed representation
    (cosθ₂,   cosθ₁*sinθ₂,  sinθ₁*sinθ₂,
     -sinθ₂,  cosθ₁*cosθ₂,  sinθ₁*cosθ₂,
     z,       -sinθ₁,       cosθ₁)
end

@inline function Base.:*(r::RotXZ, v::StaticVector)
    if length(v) != 3
        throw("Dimension mismatch: cannot rotate a vector of length $(length(v))")
    end

    sinθ₁, cosθ₁ = sincos(r.theta1)
    sinθ₂, cosθ₂ = sincos(r.theta2)

    T = Base.promote_op(*, typeof(sinθ₁), eltype(v))

    return similar_type(v,T)(cosθ₂*v[1] + -sinθ₂*v[2],
                             cosθ₁*sinθ₂*v[1] + cosθ₁*cosθ₂*v[2] + -sinθ₁*v[3],
                             sinθ₁*sinθ₂*v[1] + sinθ₁*cosθ₂*v[2] + cosθ₁*v[3])
end


"""
    struct RotZX{T} <: Rotation{3,T}
    RotZX(theta_z, theta_x)

A 3×3 rotation matrix which represents a rotation by `theta_x` about the X axis,
followed by a rotation by `theta_z` about the Z axis.
"""
RotZX

@inline function Base.Tuple(r::RotZX{T}) where T
    sinθ₁, cosθ₁ = sincos(r.theta1)
    sinθ₂, cosθ₂ = sincos(r.theta2)
    z = zero(sinθ₁)

    # transposed representation
    ( cosθ₁,         sinθ₁,         z,
     -sinθ₁*cosθ₂,   cosθ₁*cosθ₂,   sinθ₂,
      sinθ₁*sinθ₂,   cosθ₁*-sinθ₂,  cosθ₂)
end

@inline function Base.:*(r::RotZX, v::StaticVector)
    if length(v) != 3
        throw("Dimension mismatch: cannot rotate a vector of length $(length(v))")
    end

    sinθ₁, cosθ₁ = sincos(r.theta1)
    sinθ₂, cosθ₂ = sincos(r.theta2)

    T = Base.promote_op(*, typeof(sinθ₁), eltype(v))

    return similar_type(v,T)(cosθ₁*v[1] + -sinθ₁*cosθ₂*v[2] + sinθ₁*sinθ₂*v[3],
                             sinθ₁*v[1] + cosθ₁*cosθ₂*v[2] + cosθ₁*-sinθ₂*v[3],
                             sinθ₂*v[2] + cosθ₂*v[3])
end


"""
    struct RotZY{T} <: Rotation{3,T}
    RotZY(theta_z, theta_y)

A 3×3 rotation matrix which represents a rotation by `theta_y` about the Y axis,
followed by a rotation by `theta_z` about the Z axis.
"""
RotZY

@inline function Base.Tuple(r::RotZY{T}) where T
    sinθ₁, cosθ₁ = sincos(r.theta1)
    sinθ₂, cosθ₂ = sincos(r.theta2)
    z = zero(sinθ₁)

    # transposed representation
    ( cosθ₁*cosθ₂,  sinθ₁*cosθ₂, -sinθ₂,
     -sinθ₁,        cosθ₁,        z,
      cosθ₁*sinθ₂,  sinθ₁*sinθ₂,  cosθ₂)
end

@inline function Base.:*(r::RotZY, v::StaticVector)
    if length(v) != 3
        throw("Dimension mismatch: cannot rotate a vector of length $(length(v))")
    end

    sinθ₁, cosθ₁ = sincos(r.theta1)
    sinθ₂, cosθ₂ = sincos(r.theta2)

    T = Base.promote_op(*, typeof(sinθ₁), eltype(v))

    return similar_type(v,T)(cosθ₁*cosθ₂*v[1] + -sinθ₁*v[2] + cosθ₁*sinθ₂*v[3],
                             sinθ₁*cosθ₂*v[1] + cosθ₁*v[2] + sinθ₁*sinθ₂*v[3],
                             -sinθ₂*v[1] + cosθ₂*v[3])
end


"""
    struct RotYZ{T} <: Rotation{3,T}
    RotYZ(theta_y, theta_z)

A 3×3 rotation matrix which represents a rotation by `theta_z` about the Z axis,
followed by a rotation by `theta_y` about the Y axis.
"""
RotYZ

@inline function Base.Tuple(r::RotYZ{T}) where T
    sinθ₁, cosθ₁ = sincos(r.theta1)
    sinθ₂, cosθ₂ = sincos(r.theta2)
    z = zero(sinθ₁)

    # transposed representation
    (cosθ₁*cosθ₂,   sinθ₂,    -sinθ₁*cosθ₂,
     cosθ₁*-sinθ₂,  cosθ₂,     sinθ₁*sinθ₂,
     sinθ₁,         z,         cosθ₁)
end

@inline function Base.:*(r::RotYZ, v::StaticVector)
    if length(v) != 3
        throw("Dimension mismatch: cannot rotate a vector of length $(length(v))")
    end

    sinθ₁, cosθ₁ = sincos(r.theta1)
    sinθ₂, cosθ₂ = sincos(r.theta2)

    T = Base.promote_op(*, typeof(sinθ₁), eltype(v))

    return similar_type(v,T)(cosθ₁*cosθ₂*v[1] + cosθ₁*-sinθ₂*v[2] + sinθ₁*v[3],
                             sinθ₂*v[1] + cosθ₂*v[2],
                             -sinθ₁*cosθ₂*v[1] + sinθ₁*sinθ₂*v[2] + cosθ₁*v[3])
end

################################################################################
################################################################################

########################
# Three axis Rotations #
########################

for axis1 in [:X, :Y, :Z]
    Rot1Type = Symbol("Rot" * string(axis1))
    for axis2 in filter(axis -> axis != axis1, [:X, :Y, :Z])
        Rot2Type = Symbol("Rot" * string(axis2))
        Rot12Type = Symbol("Rot" * string(axis1) * string(axis2))
        for axis3 in filter(axis -> axis != axis2, [:X, :Y, :Z])
            Rot3Type = Symbol("Rot" * string(axis3))
            Rot23Type = Symbol("Rot" * string(axis2) * string(axis3))
            RotType = Symbol("Rot" * string(axis1) * string(axis2) * string(axis3))
            InvRotType = Symbol("Rot" * string(axis3) * string(axis2) * string(axis1))

            @eval begin
                struct $RotType{T} <: Rotation{3,T}
                    theta1::T
                    theta2::T
                    theta3::T
                    $RotType{T}(theta1, theta2, theta3) where {T} = new{T}(theta1, theta2, theta3)
                    $RotType{T}(r::$RotType) where {T} = new{T}(r.theta1, r.theta2, r.theta3)
                end

                @inline function $RotType(theta1, theta2, theta3)
                    ts = promote(theta1, theta2, theta3)
                    $RotType{rot_eltype(eltype(ts))}(ts...)
                end
                @inline $RotType(r::$RotType{T}) where {T} = $RotType{T}(r)

                @inline function Base.getindex(r::$RotType{T}, i::Int) where T
                    Tuple(r)[i] # Slow...
                end

                # Composing single-axis rotations with two-axis rotations:
                @inline Base.:*(r1::$Rot1Type, r2::$Rot23Type) = $RotType(r1.theta, r2.theta1, r2.theta2)
                @inline Base.:*(r1::$Rot12Type, r2::$Rot3Type) = $RotType(r1.theta1, r1.theta2, r2.theta)

                # Composing with single-axis rotations:
                @inline Base.:*(r1::$RotType, r2::$Rot3Type) = $RotType(r1.theta1, r1.theta2, r1.theta3 + r2.theta)
                @inline Base.:*(r1::$Rot1Type, r2::$RotType) = $RotType(r1.theta + r2.theta1, r2.theta2, r2.theta3)

                @inline Base.inv(r::$RotType) = $InvRotType(-r.theta3, -r.theta2, -r.theta1)

                # define null rotations for convenience
                @inline Base.one(::Type{$RotType}) = $RotType(0.0, 0.0, 0.0)
                @inline Base.one(::Type{$RotType{T}}) where {T} = $RotType{T}(zero(T), zero(T), zero(T))

                params(r::$RotType) = SVector{3}(r.theta1, r.theta2, r.theta3)
            end
        end
    end
end


################################################################################
################################################################################

##########################
# Proper Euler Rotations #
##########################

"""
    struct RotXYX{T} <: Rotation{3,T}
    RotXYX(theta1, theta2, theta3)

A 3×3 rotation matrix parameterized by the "proper" XYX Euler angle convention,
consisting of first a rotation about the X axis by `theta3`, followed by a
rotation about the Y axis by `theta2`, and finally a rotation about the X axis
by `theta1`.
"""
RotXYX

@inline function (::Type{Rot})(t::NTuple{9}) where Rot <: RotXYX
    R = SMatrix{3,3}(t)

    t1 = atan(R[2, 1], (-R[3, 1] + eps(t[1])) - eps(t[1]))  # TODO: handle denormal numbers better, as atan(0,0) != atan(0,-0)
    st1, ct1 = sincos(t1)

    Rot(t1,
        atan((R[1, 2] * R[1, 2] + R[1, 3] * R[1, 3])^(1/2), R[1, 1]),
        atan(- R[2, 3]*ct1 - R[3, 3]*st1, R[2, 2]*ct1 + R[3, 2]*st1))
end

@inline function Base.Tuple(r::RotXYX{T}) where T
    sinθ₁, cosθ₁ = sincos(r.theta1)
    sinθ₂, cosθ₂ = sincos(r.theta2)
    sinθ₃, cosθ₃ = sincos(r.theta3)

    # transposed representation
    (cosθ₂,        sinθ₁*sinθ₂,                        cosθ₁*-sinθ₂,
     sinθ₂*sinθ₃,  cosθ₁*cosθ₃ + -sinθ₁*cosθ₂*sinθ₃,   sinθ₁*cosθ₃ + cosθ₁*cosθ₂*sinθ₃,
     sinθ₂*cosθ₃,  cosθ₁*-sinθ₃ + -sinθ₁*cosθ₂*cosθ₃,  sinθ₁*-sinθ₃ + cosθ₁*cosθ₂*cosθ₃)
end

@inline function Base.:*(r::RotXYX, v::StaticVector)
    if length(v) != 3
        throw("Dimension mismatch: cannot rotate a vector of length $(length(v))")
    end

    sinθ₁, cosθ₁ = sincos(r.theta1)
    sinθ₂, cosθ₂ = sincos(r.theta2)
    sinθ₃, cosθ₃ = sincos(r.theta3)

    T = Base.promote_op(*, typeof(sinθ₁), eltype(v))

    return similar_type(v,T)(
        cosθ₂*v[1] + sinθ₂*sinθ₃*v[2] + sinθ₂*cosθ₃*v[3],
        sinθ₁*sinθ₂*v[1] + (cosθ₁*cosθ₃ + -sinθ₁*cosθ₂*sinθ₃)*v[2] + (cosθ₁*-sinθ₃ + -sinθ₁*cosθ₂*cosθ₃)*v[3],
        cosθ₁*-sinθ₂*v[1] + (sinθ₁*cosθ₃ + cosθ₁*cosθ₂*sinθ₃)*v[2] + (sinθ₁*-sinθ₃ + cosθ₁*cosθ₂*cosθ₃)*v[3])
end


"""
    struct RotXZX{T} <: Rotation{3,T}
    RotXZX(theta1, theta2, theta3)

A 3×3 rotation matrix parameterized by the "proper" XZX Euler angle convention,
consisting of first a rotation about the X axis by `theta3`, followed by a
rotation about the Z axis by `theta2`, and finally a rotation about the X axis
by `theta1`.
"""
RotXZX

@inline function (::Type{Rot})(t::NTuple{9}) where Rot <: RotXZX
    R = SMatrix{3,3}(t)

    t1 = atan(R[3, 1], R[2, 1])
    st1, ct1 = sincos(t1)

    Rot(t1,
        atan((R[1, 2] * R[1, 2] + R[1, 3] * R[1, 3])^(1/2), R[1, 1]),
        atan(R[3, 2]*ct1 - R[2, 2]*st1, R[3, 3]*ct1 - R[2, 3]*st1))
end

@inline function Base.Tuple(r::RotXZX{T}) where T
    sinθ₁, cosθ₁ = sincos(r.theta1)
    sinθ₂, cosθ₂ = sincos(r.theta2)
    sinθ₃, cosθ₃ = sincos(r.theta3)

    # transposed representation
    (cosθ₂,         cosθ₁*sinθ₂,                        sinθ₁*sinθ₂,
     -sinθ₂*cosθ₃,  cosθ₁*cosθ₂*cosθ₃ + -sinθ₁*sinθ₃,   sinθ₁*cosθ₂*cosθ₃ + cosθ₁*sinθ₃,
     sinθ₂*sinθ₃,   cosθ₁*cosθ₂*-sinθ₃ + -sinθ₁*cosθ₃,  sinθ₁*cosθ₂*-sinθ₃ + cosθ₁*cosθ₃)
end

@inline function Base.:*(r::RotXZX, v::StaticVector)
    if length(v) != 3
        throw("Dimension mismatch: cannot rotate a vector of length $(length(v))")
    end

    sinθ₁, cosθ₁ = sincos(r.theta1)
    sinθ₂, cosθ₂ = sincos(r.theta2)
    sinθ₃, cosθ₃ = sincos(r.theta3)

    T = Base.promote_op(*, typeof(sinθ₁), eltype(v))

    return similar_type(v,T)(
        cosθ₂*v[1] + -sinθ₂*cosθ₃*v[2] + sinθ₂*sinθ₃*v[3],
        cosθ₁*sinθ₂*v[1] + (cosθ₁*cosθ₂*cosθ₃ + -sinθ₁*sinθ₃)*v[2] + (cosθ₁*cosθ₂*-sinθ₃ + -sinθ₁*cosθ₃)*v[3],
        sinθ₁*sinθ₂*v[1] + (sinθ₁*cosθ₂*cosθ₃ + cosθ₁*sinθ₃)*v[2] + (sinθ₁*cosθ₂*-sinθ₃ + cosθ₁*cosθ₃)*v[3])
end


"""
    struct RotYXY{T} <: Rotation{3,T}
    RotYXY(theta1, theta2, theta3)

A 3×3 rotation matrix parameterized by the "proper" YXY Euler angle convention,
consisting of first a rotation about the Y axis by `theta3`, followed by a
rotation about the X axis by `theta2`, and finally a rotation about the Y axis
by `theta1`.
"""
RotYXY

@inline function (::Type{Rot})(t::NTuple{9}) where Rot <: RotYXY
    R = SMatrix{3,3}(t)

    t1 = atan(R[1, 2], R[3, 2])
    st1, ct1 = sincos(t1)

    Rot(t1,
        atan((R[2, 1] * R[2, 1] + R[2, 3] * R[2, 3])^(1/2), R[2, 2]),
        atan(R[1, 3]*ct1 - R[3, 3]*st1, R[1, 1]*ct1 - R[3, 1]*st1))
end

@inline function Base.Tuple(r::RotYXY{T}) where T
    sinθ₁, cosθ₁ = sincos(r.theta1)
    sinθ₂, cosθ₂ = sincos(r.theta2)
    sinθ₃, cosθ₃ = sincos(r.theta3)

    # transposed representation
    (cosθ₁*cosθ₃ + sinθ₁*cosθ₂*-sinθ₃,  sinθ₂*sinθ₃,   -sinθ₁*cosθ₃ + cosθ₁*cosθ₂*-sinθ₃,
     sinθ₁*sinθ₂,                       cosθ₂,         cosθ₁*sinθ₂,
     cosθ₁*sinθ₃ + sinθ₁*cosθ₂*cosθ₃,   -sinθ₂*cosθ₃,  -sinθ₁*sinθ₃ + cosθ₁*cosθ₂*cosθ₃)
end

@inline function Base.:*(r::RotYXY, v::StaticVector)
    if length(v) != 3
        throw("Dimension mismatch: cannot rotate a vector of length $(length(v))")
    end

    sinθ₁, cosθ₁ = sincos(r.theta1)
    sinθ₂, cosθ₂ = sincos(r.theta2)
    sinθ₃, cosθ₃ = sincos(r.theta3)

    T = Base.promote_op(*, typeof(sinθ₁), eltype(v))

    return similar_type(v,T)(
        (cosθ₁*cosθ₃ + sinθ₁*cosθ₂*-sinθ₃)*v[1] + sinθ₁*sinθ₂*v[2] + (cosθ₁*sinθ₃ + sinθ₁*cosθ₂*cosθ₃)*v[3],
        sinθ₂*sinθ₃*v[1] + cosθ₂*v[2] + -sinθ₂*cosθ₃*v[3],
        (-sinθ₁*cosθ₃ + cosθ₁*cosθ₂*-sinθ₃)*v[1] + cosθ₁*sinθ₂*v[2] + (-sinθ₁*sinθ₃ + cosθ₁*cosθ₂*cosθ₃)*v[3])
end


"""
    struct RotYZY{T} <: Rotation{3,T}
    RotYZY(theta1, theta2, theta3)

A 3×3 rotation matrix parameterized by the "proper" YXY Euler angle convention,
consisting of first a rotation about the Y axis by `theta3`, followed by a
rotation about the Z axis by `theta2`, and finally a rotation about the Y axis
by `theta1`.
"""
RotYZY

@inline function (::Type{Rot})(t::NTuple{9}) where Rot <: RotYZY
    R = SMatrix{3,3}(t)

    t1 = atan(R[3, 2], -R[1, 2])  # TODO: handle denormal numbers better, as atan(0,0) != atan(0,-0)
    st1, ct1 = sincos(t1)

    Rot(t1,
        atan((R[2, 1] * R[2, 1] + R[2, 3] * R[2, 3])^(1/2), R[2, 2]),
        atan(- R[3, 1]*ct1 - R[1, 1]*st1, R[3, 3]*ct1 + R[1, 3]*st1))
end

@inline function Base.Tuple(r::RotYZY{T}) where T
    sinθ₁, cosθ₁ = sincos(r.theta1)
    sinθ₂, cosθ₂ = sincos(r.theta2)
    sinθ₃, cosθ₃ = sincos(r.theta3)

    # transposed representation
    (cosθ₁*cosθ₂*cosθ₃ + sinθ₁*-sinθ₃,  sinθ₂*cosθ₃,  -sinθ₁*cosθ₂*cosθ₃ + cosθ₁*-sinθ₃,
     cosθ₁*-sinθ₂,                      cosθ₂,        sinθ₁*sinθ₂,
     cosθ₁*cosθ₂*sinθ₃ + sinθ₁*cosθ₃,   sinθ₂*sinθ₃,  -sinθ₁*cosθ₂*sinθ₃ + cosθ₁*cosθ₃)
end

@inline function Base.:*(r::RotYZY, v::StaticVector)
    if length(v) != 3
        throw("Dimension mismatch: cannot rotate a vector of length $(length(v))")
    end

    sinθ₁, cosθ₁ = sincos(r.theta1)
    sinθ₂, cosθ₂ = sincos(r.theta2)
    sinθ₃, cosθ₃ = sincos(r.theta3)

    T = Base.promote_op(*, typeof(sinθ₁), eltype(v))

    return similar_type(v,T)(
        (cosθ₁*cosθ₂*cosθ₃ + sinθ₁*-sinθ₃)*v[1] + cosθ₁*-sinθ₂*v[2] + (cosθ₁*cosθ₂*sinθ₃ + sinθ₁*cosθ₃)*v[3],
        sinθ₂*cosθ₃*v[1] + cosθ₂*v[2] + sinθ₂*sinθ₃*v[3],
        (-sinθ₁*cosθ₂*cosθ₃ + cosθ₁*-sinθ₃)*v[1] + sinθ₁*sinθ₂*v[2] + (-sinθ₁*cosθ₂*sinθ₃ + cosθ₁*cosθ₃)*v[3])
end


"""
    struct RotZXZ{T} <: Rotation{3,T}
    RotZXZ(theta1, theta2, theta3)

A 3×3 rotation matrix parameterized by the "proper" ZXZ Euler angle convention,
consisting of first a rotation about the Z axis by `theta3`, followed by a
rotation about the X axis by `theta2`, and finally a rotation about the Z axis
by `theta1`.
"""
RotZXZ

@inline function (::Type{Rot})(t::NTuple{9}) where Rot <: RotZXZ
    R = SMatrix{3,3}(t)

    t1 = atan(R[1, 3], (-R[2, 3] + eps()) - eps())  # TODO: handle denormal numbers better, as atan(0,0) != atan(0,-0)
    st1, ct1 = sincos(t1)

    Rot(t1,
        atan((R[3, 1] * R[3, 1] + R[3, 2] * R[3, 2])^(1/2), R[3, 3]),
        atan(- R[1, 2]*ct1 - R[2, 2]*st1, R[1, 1]*ct1 + R[2, 1]*st1))
end

@inline function Base.Tuple(r::RotZXZ{T}) where T
    sinθ₁, cosθ₁ = sincos(r.theta1)
    sinθ₂, cosθ₂ = sincos(r.theta2)
    sinθ₃, cosθ₃ = sincos(r.theta3)

    # transposed representation
    (cosθ₁*cosθ₃ + -sinθ₁*cosθ₂*sinθ₃,   sinθ₁*cosθ₃ + cosθ₁*cosθ₂*sinθ₃,   sinθ₂*sinθ₃,
     cosθ₁*-sinθ₃ + -sinθ₁*cosθ₂*cosθ₃,  sinθ₁*-sinθ₃ + cosθ₁*cosθ₂*cosθ₃,  sinθ₂*cosθ₃,
     sinθ₁*sinθ₂,                        cosθ₁*-sinθ₂,                      cosθ₂)
end

@inline function Base.:*(r::RotZXZ, v::StaticVector)
    if length(v) != 3
        throw("Dimension mismatch: cannot rotate a vector of length $(length(v))")
    end

    sinθ₁, cosθ₁ = sincos(r.theta1)
    sinθ₂, cosθ₂ = sincos(r.theta2)
    sinθ₃, cosθ₃ = sincos(r.theta3)

    T = Base.promote_op(*, typeof(sinθ₁), eltype(v))

    return similar_type(v,T)(
         (cosθ₁*cosθ₃ + -sinθ₁*cosθ₂*sinθ₃)*v[1] + (cosθ₁*-sinθ₃ + -sinθ₁*cosθ₂*cosθ₃)*v[2] + sinθ₁*sinθ₂*v[3],
         (sinθ₁*cosθ₃ + cosθ₁*cosθ₂*sinθ₃)*v[1] + (sinθ₁*-sinθ₃ + cosθ₁*cosθ₂*cosθ₃)*v[2] + cosθ₁*-sinθ₂*v[3],
         sinθ₂*sinθ₃*v[1] + sinθ₂*cosθ₃*v[2] + cosθ₂*v[3])
end


"""
    struct RotZYZ{T} <: Rotation{3,T}
    RotZYZ(theta1, theta2, theta3)

A 3×3 rotation matrix parameterized by the "proper" ZXZ Euler angle convention,
consisting of first a rotation about the Z axis by `theta3`, followed by a
rotation about the Y axis by `theta2`, and finally a rotation about the Z axis
by `theta1`.
"""
RotZYZ

@inline function (::Type{Rot})(t::NTuple{9}) where Rot <: RotZYZ
    R = SMatrix{3,3}(t)

    t1 = atan(R[2, 3], R[1, 3])
    st1, ct1 = sincos(t1)

    Rot(t1,
        atan((R[3, 1] * R[3, 1] + R[3, 2] * R[3, 2])^(1/2), R[3, 3]),
        atan(R[2, 1]*ct1 - R[1, 1]*st1, R[2, 2]*ct1 - R[1, 2]*st1))
end

@inline function Base.Tuple(r::RotZYZ{T}) where T
    sinθ₁, cosθ₁ = sincos(r.theta1)
    sinθ₂, cosθ₂ = sincos(r.theta2)
    sinθ₃, cosθ₃ = sincos(r.theta3)

    # transposed representation
    (cosθ₁*cosθ₂*cosθ₃ + -sinθ₁*sinθ₃,   sinθ₁*cosθ₂*cosθ₃ + cosθ₁*sinθ₃,   -sinθ₂*cosθ₃,
     cosθ₁*cosθ₂*-sinθ₃ + -sinθ₁*cosθ₃,  sinθ₁*cosθ₂*-sinθ₃ + cosθ₁*cosθ₃,  sinθ₂*sinθ₃,
     cosθ₁*sinθ₂,                        sinθ₁*sinθ₂,                       cosθ₂)
end

@inline function Base.:*(r::RotZYZ, v::StaticVector)
    if length(v) != 3
        throw("Dimension mismatch: cannot rotate a vector of length $(length(v))")
    end

    sinθ₁, cosθ₁ = sincos(r.theta1)
    sinθ₂, cosθ₂ = sincos(r.theta2)
    sinθ₃, cosθ₃ = sincos(r.theta3)

    T = Base.promote_op(*, typeof(sinθ₁), eltype(v))

    return similar_type(v,T)(
        (cosθ₁*cosθ₂*cosθ₃ + -sinθ₁*sinθ₃)*v[1] + (cosθ₁*cosθ₂*-sinθ₃ + -sinθ₁*cosθ₃)*v[2] + cosθ₁*sinθ₂*v[3],
        (sinθ₁*cosθ₂*cosθ₃ + cosθ₁*sinθ₃)*v[1] + (sinθ₁*cosθ₂*-sinθ₃ + cosθ₁*cosθ₃)*v[2] + sinθ₁*sinθ₂*v[3],
        -sinθ₂*cosθ₃*v[1] + sinθ₂*sinθ₃*v[2] + cosθ₂*v[3])
end

###############################
# Tait-Bryant Euler Rotations #
###############################

"""
    struct RotXYZ{T} <: Rotation{3,T}
    RotXYZ(theta1, theta2, theta3)
    RotXYZ(roll=r, pitch=p, yaw=y)

A 3×3 rotation matrix parameterized by the "Tait-Bryant" XYZ Euler angle
convention, consisting of first a rotation about the Z axis by `theta3`,
followed by a rotation about the Y axis by `theta2`, and finally a rotation
about the X axis by `theta1`.

The keyword argument form applies roll, pitch and yaw to the X, Y and Z axes
respectively, in XYZ order. (Because it is a right-handed coordinate system,
note that positive pitch is heading in the negative Z axis).
"""
RotXYZ

@inline (::Type{Rot})(; roll=0, pitch=0, yaw=0) where {Rot<:RotXYZ} = Rot(roll, pitch, yaw)

@inline function (::Type{Rot})(t::NTuple{9}) where Rot <: RotXYZ
    R = SMatrix{3,3}(t)

    t1 = atan(-R[2, 3], R[3, 3])
    st1, ct1 = sincos(t1)

    Rot(t1,
        atan(R[1, 3], (R[1, 1] * R[1, 1] + R[1, 2] * R[1, 2])^(1/2)),
        atan(R[2, 1]*ct1 + R[3, 1]*st1, R[2, 2]*ct1 + R[3, 2]*st1))
end

@inline function Base.Tuple(r::RotXYZ{T}) where T
    sinθ₁, cosθ₁ = sincos(r.theta1)
    sinθ₂, cosθ₂ = sincos(r.theta2)
    sinθ₃, cosθ₃ = sincos(r.theta3)

    # transposed representation
    (cosθ₂*cosθ₃,   sinθ₁*sinθ₂*cosθ₃ + cosθ₁*sinθ₃,   cosθ₁*-sinθ₂*cosθ₃ + sinθ₁*sinθ₃,
     cosθ₂*-sinθ₃,  sinθ₁*sinθ₂*-sinθ₃ + cosθ₁*cosθ₃,  cosθ₁*-sinθ₂*-sinθ₃ + sinθ₁*cosθ₃,
     sinθ₂,         -sinθ₁*cosθ₂,                      cosθ₁*cosθ₂)
end

@inline function Base.:*(r::RotXYZ, v::StaticVector)
    if length(v) != 3
        throw("Dimension mismatch: cannot rotate a vector of length $(length(v))")
    end

    sinθ₁, cosθ₁ = sincos(r.theta1)
    sinθ₂, cosθ₂ = sincos(r.theta2)
    sinθ₃, cosθ₃ = sincos(r.theta3)

    T = Base.promote_op(*, typeof(sinθ₁), eltype(v))

    return similar_type(v,T)(
        cosθ₂*cosθ₃*v[1] + cosθ₂*-sinθ₃*v[2] + sinθ₂*v[3],
        (sinθ₁*sinθ₂*cosθ₃ + cosθ₁*sinθ₃)*v[1] + (sinθ₁*sinθ₂*-sinθ₃ + cosθ₁*cosθ₃)*v[2] + -sinθ₁*cosθ₂*v[3],
        (cosθ₁*-sinθ₂*cosθ₃ + sinθ₁*sinθ₃)*v[1] + (cosθ₁*sinθ₂*sinθ₃ + sinθ₁*cosθ₃)*v[2] + cosθ₁*cosθ₂*v[3])
end


"""
    struct RotZYX{T} <: Rotation{3,T}
    RotZYX(theta1, theta2, theta3)
    RotZYX(roll=r, pitch=p, yaw=y)

A 3×3 rotation matrix parameterized by the "Tait-Bryant" ZYX Euler angle
convention, consisting of first a rotation about the X axis by `theta3`,
followed by a rotation about the Y axis by `theta2`, and finally a rotation
about the Z axis by `theta1`.

The keyword argument form applies roll, pitch and yaw to the X, Y and Z axes
respectively, in ZYX order. (Because it is a right-handed coordinate system,
note that positive pitch is heading in the negative Z axis).
"""
RotZYX

@inline (::Type{Rot})(; roll=0, pitch=0, yaw=0) where {Rot<:RotZYX} = Rot(yaw, pitch, roll)

@inline function (::Type{Rot})(t::NTuple{9}) where Rot <: RotZYX
    R = SMatrix{3,3}(t)

    t1 = atan(R[2, 1], R[1, 1])
    st1, ct1 = sincos(t1)

    Rot(t1,
        atan(-R[3, 1], (R[3, 2] * R[3, 2] + R[3, 3] * R[3, 3])^(1/2)),
        atan(R[1, 3]*st1 - R[2, 3]*ct1, R[2, 2]*ct1 - R[1, 2]*st1))
end

@inline function Base.Tuple(r::RotZYX{T}) where T
    sinθ₁, cosθ₁ = sincos(r.theta1)
    sinθ₂, cosθ₂ = sincos(r.theta2)
    sinθ₃, cosθ₃ = sincos(r.theta3)

    # transposed representation
    ( cosθ₁*cosθ₂,                       sinθ₁*cosθ₂,                      -sinθ₂,
     -sinθ₁*cosθ₃ + cosθ₁*sinθ₂*sinθ₃,   cosθ₁*cosθ₃ + sinθ₁*sinθ₂*sinθ₃,   cosθ₂*sinθ₃,
      sinθ₁*sinθ₃ + cosθ₁*sinθ₂*cosθ₃,   cosθ₁*-sinθ₃ + sinθ₁*sinθ₂*cosθ₃,  cosθ₂*cosθ₃)
end

@inline function Base.:*(r::RotZYX, v::StaticVector)
    if length(v) != 3
        throw("Dimension mismatch: cannot rotate a vector of length $(length(v))")
    end

    sinθ₁, cosθ₁ = sincos(r.theta1)
    sinθ₂, cosθ₂ = sincos(r.theta2)
    sinθ₃, cosθ₃ = sincos(r.theta3)

    T = Base.promote_op(*, typeof(sinθ₁), eltype(v))

    return similar_type(v,T)(
        cosθ₁*cosθ₂*v[1] + (-sinθ₁*cosθ₃ + cosθ₁*sinθ₂*sinθ₃)*v[2] + (sinθ₁*sinθ₃ + cosθ₁*sinθ₂*cosθ₃)*v[3],
        sinθ₁*cosθ₂*v[1] + (cosθ₁*cosθ₃ + sinθ₁*sinθ₂*sinθ₃)*v[2] + (cosθ₁*-sinθ₃ + sinθ₁*sinθ₂*cosθ₃)*v[3],
        -sinθ₂*v[1] + cosθ₂*sinθ₃*v[2] + cosθ₂*cosθ₃*v[3])
end


"""
    struct RotXZY{T} <: Rotation{3,T}
    RotXZY(theta1, theta2, theta3)
    RotXZY(roll=r, pitch=p, yaw=y)

A 3×3 rotation matrix parameterized by the "Tait-Bryant" XZY Euler angle
convention, consisting of first a rotation about the Y axis by `theta3`,
followed by a rotation about the Z axis by `theta2`, and finally a rotation
about the X axis by `theta1`.

The keyword argument form applies roll, pitch and yaw to the X, Y and Z axes
respectively, in XZY order. (Because it is a right-handed coordinate system,
note that positive pitch is heading in the negative Z axis).
"""
RotXZY

@inline (::Type{Rot})(; roll=0, pitch=0, yaw=0) where {Rot<:RotXZY} = Rot(roll, yaw, pitch)

@inline function (::Type{Rot})(t::NTuple{9}) where Rot <: RotXZY
    R = SMatrix{3,3}(t)

    t1 = atan(R[3, 2], R[2, 2])
    st1, ct1 = sincos(t1)

    Rot(t1,
        atan(-R[1, 2], (R[1, 1] * R[1, 1] + R[1, 3] * R[1, 3])^(1/2)),
        atan(R[2, 1]*st1 - R[3, 1]*ct1, R[3, 3]*ct1 - R[2, 3]*st1))
end

@inline function Base.Tuple(r::RotXZY{T}) where T
    sinθ₁, cosθ₁ = sincos(r.theta1)
    sinθ₂, cosθ₂ = sincos(r.theta2)
    sinθ₃, cosθ₃ = sincos(r.theta3)

    # transposed representation
    ( cosθ₂*cosθ₃,  cosθ₁*sinθ₂*cosθ₃ + sinθ₁*sinθ₃,   sinθ₁*sinθ₂*cosθ₃ + cosθ₁*-sinθ₃,
     -sinθ₂,        cosθ₁*cosθ₂,                       sinθ₁*cosθ₂,
      cosθ₂*sinθ₃,  cosθ₁*sinθ₂*sinθ₃ + -sinθ₁*cosθ₃,  sinθ₁*sinθ₂*sinθ₃ + cosθ₁*cosθ₃)
end

@inline function Base.:*(r::RotXZY, v::StaticVector)
    if length(v) != 3
        throw("Dimension mismatch: cannot rotate a vector of length $(length(v))")
    end

    sinθ₁, cosθ₁ = sincos(r.theta1)
    sinθ₂, cosθ₂ = sincos(r.theta2)
    sinθ₃, cosθ₃ = sincos(r.theta3)

    T = Base.promote_op(*, typeof(sinθ₁), eltype(v))

    return similar_type(v,T)(
        cosθ₂*cosθ₃*v[1] + -sinθ₂*v[2] + cosθ₂*sinθ₃*v[3],
        (cosθ₁*sinθ₂*cosθ₃ + sinθ₁*sinθ₃)*v[1] + cosθ₁*cosθ₂*v[2] + (cosθ₁*sinθ₂*sinθ₃ + -sinθ₁*cosθ₃)*v[3],
        (sinθ₁*sinθ₂*cosθ₃ + cosθ₁*-sinθ₃)*v[1] + sinθ₁*cosθ₂*v[2] + (sinθ₁*sinθ₂*sinθ₃ + cosθ₁*cosθ₃)*v[3])
end


"""
    struct RotYZX{T} <: Rotation{3,T}
    RotYZX(theta1, theta2, theta3)
    RotYZX(roll=r, pitch=p, yaw=y)

A 3×3 rotation matrix parameterized by the "Tait-Bryant" YZX Euler angle
convention, consisting of first a rotation about the X axis by `theta3`,
followed by a rotation about the Z axis by `theta2`, and finally a rotation
about the Y axis by `theta1`.

The keyword argument form applies roll, pitch and yaw to the X, Y and Z axes
respectively, in YZX order. (Because it is a right-handed coordinate system,
note that positive pitch is heading in the negative Z axis).
"""
RotYZX

@inline (::Type{Rot})(; roll=0, pitch=0, yaw=0) where {Rot<:RotYZX} = Rot(pitch, yaw, roll)

@inline function (::Type{Rot})(t::NTuple{9}) where Rot <: RotYZX
    R = SMatrix{3,3}(t)

    t1 = atan(-R[3, 1], R[1, 1])
    st1, ct1 = sincos(t1)

    Rot(t1,
        atan(R[2, 1], (R[2, 2] * R[2, 2] + R[2, 3] * R[2, 3])^(1/2)),
        atan(R[3, 2]*ct1 + R[1, 2]*st1, R[3, 3]*ct1 + R[1, 3]*st1))
end

@inline function Base.Tuple(r::RotYZX{T}) where T
    sinθ₁, cosθ₁ = sincos(r.theta1)
    sinθ₂, cosθ₂ = sincos(r.theta2)
    sinθ₃, cosθ₃ = sincos(r.theta3)

    # transposed representation
    (cosθ₁*cosθ₂,                        sinθ₂,         -sinθ₁*cosθ₂,
     cosθ₁*-sinθ₂*cosθ₃ + sinθ₁*sinθ₃,   cosθ₂*cosθ₃,   sinθ₁*sinθ₂*cosθ₃ + cosθ₁*sinθ₃,
     cosθ₁*sinθ₂*sinθ₃ + sinθ₁*cosθ₃,    cosθ₂*-sinθ₃,  sinθ₁*sinθ₂*-sinθ₃ + cosθ₁*cosθ₃)
end

@inline function Base.:*(r::RotYZX, v::StaticVector)
    if length(v) != 3
        throw("Dimension mismatch: cannot rotate a vector of length $(length(v))")
    end

    sinθ₁, cosθ₁ = sincos(r.theta1)
    sinθ₂, cosθ₂ = sincos(r.theta2)
    sinθ₃, cosθ₃ = sincos(r.theta3)

    T = Base.promote_op(*, typeof(sinθ₁), eltype(v))

    return similar_type(v,T)(
        cosθ₁*cosθ₂*v[1] + (cosθ₁*-sinθ₂*cosθ₃ + sinθ₁*sinθ₃)*v[2] + (cosθ₁*-sinθ₂*-sinθ₃ + sinθ₁*cosθ₃)*v[3],
        sinθ₂*v[1] + cosθ₂*cosθ₃*v[2] + cosθ₂*-sinθ₃*v[3],
        -sinθ₁*cosθ₂*v[1] + (sinθ₁*sinθ₂*cosθ₃ + cosθ₁*sinθ₃)*v[2] + (sinθ₁*sinθ₂*-sinθ₃ + cosθ₁*cosθ₃)*v[3])
end


"""
    struct RotYXZ{T} <: Rotation{3,T}
    RotYXZ(theta1, theta2, theta3)
    RotYXZ(roll=r, pitch=p, yaw=y)

A 3×3 rotation matrix parameterized by the "Tait-Bryant" YXZ Euler angle
convention, consisting of first a rotation about the Z axis by `theta3`,
followed by a rotation about the X axis by `theta2`, and finally a rotation
about the Y axis by `theta1`.

The keyword argument form applies roll, pitch and yaw to the X, Y and Z axes
respectively, in YXZ order. (Because it is a right-handed coordinate system,
note that positive pitch is heading in the negative Z axis).
"""
RotYXZ

@inline (::Type{Rot})(; roll=0, pitch=0, yaw=0) where {Rot<:RotYXZ} = Rot(pitch, roll, yaw)

@inline function (::Type{Rot})(t::NTuple{9}) where Rot <: RotYXZ
    R = SMatrix{3,3}(t)

    t1 = atan(R[1, 3], R[3, 3])
    st1, ct1 = sincos(t1)

    Rot(t1,
        atan(-R[2, 3], (R[2, 1] * R[2, 1] + R[2, 2] * R[2, 2])^(1/2)),
        atan(R[3, 2]*st1 - R[1, 2]*ct1, R[1, 1]*ct1 - R[3, 1]*st1))
end

@inline function Base.Tuple(r::RotYXZ{T}) where T
    sinθ₁, cosθ₁ = sincos(r.theta1)
    sinθ₂, cosθ₂ = sincos(r.theta2)
    sinθ₃, cosθ₃ = sincos(r.theta3)

    # transposed representation
    (cosθ₁*cosθ₃ + sinθ₁*sinθ₂*sinθ₃,   cosθ₂*sinθ₃,  -sinθ₁*cosθ₃ + cosθ₁*sinθ₂*sinθ₃,
     cosθ₁*-sinθ₃ + sinθ₁*sinθ₂*cosθ₃,  cosθ₂*cosθ₃,  sinθ₁*sinθ₃ + cosθ₁*sinθ₂*cosθ₃,
     sinθ₁*cosθ₂,                       -sinθ₂,       cosθ₁*cosθ₂)
end

@inline function Base.:*(r::RotYXZ, v::StaticVector)
    if length(v) != 3
        throw("Dimension mismatch: cannot rotate a vector of length $(length(v))")
    end

    sinθ₁, cosθ₁ = sincos(r.theta1)
    sinθ₂, cosθ₂ = sincos(r.theta2)
    sinθ₃, cosθ₃ = sincos(r.theta3)

    T = Base.promote_op(*, typeof(sinθ₁), eltype(v))

    return similar_type(v,T)(
        (cosθ₁*cosθ₃ + sinθ₁*sinθ₂*sinθ₃)*v[1] + (cosθ₁*-sinθ₃ + sinθ₁*sinθ₂*cosθ₃)*v[2] + sinθ₁*cosθ₂*v[3],
        cosθ₂*sinθ₃*v[1] + cosθ₂*cosθ₃*v[2] + -sinθ₂*v[3],
        (-sinθ₁*cosθ₃ + cosθ₁*sinθ₂*sinθ₃)*v[1] + (sinθ₁*sinθ₃ + cosθ₁*sinθ₂*cosθ₃)*v[2] + cosθ₁*cosθ₂*v[3])
end


"""
    struct RotZXY{T} <: Rotation{3,T}
    RotZXY(theta1, theta2, theta3)
    RotZXY(roll=r, pitch=p, yaw=y)

A 3×3 rotation matrix parameterized by the "Tait-Bryant" ZXY Euler angle
convention, consisting of first a rotation about the Y axis by `theta3`,
followed by a rotation about the X axis by `theta2`, and finally a rotation
about the Z axis by `theta1`.

The keyword argument form applies roll, pitch and yaw to the X, Y and Z axes
respectively, in ZXY order. (Because it is a right-handed coordinate system,
note that positive pitch is heading in the negative Z axis).
"""
RotZXY

@inline (::Type{Rot})(; roll=0, pitch=0, yaw=0) where {Rot<:RotZXY} = Rot(yaw, roll, pitch)

@inline function (::Type{Rot})(t::NTuple{9}) where Rot <: RotZXY
    R = SMatrix{3,3}(t)

    t1 = atan(-R[1, 2], R[2, 2])
    st1, ct1 = sincos(t1)

    Rot(t1,
        atan(R[3, 2], (R[3, 1] * R[3, 1] + R[3, 3] * R[3, 3])^(1/2)),
        atan(R[1, 3]*ct1 + R[2, 3]*st1, R[1, 1]*ct1 + R[2, 1]*st1))
end

@inline function Base.Tuple(r::RotZXY{T}) where T
    sinθ₁, cosθ₁ = sincos(r.theta1)
    sinθ₂, cosθ₂ = sincos(r.theta2)
    sinθ₃, cosθ₃ = sincos(r.theta3)

    # transposed representation
    ( cosθ₁*cosθ₃ + sinθ₁*sinθ₂*-sinθ₃,  sinθ₁*cosθ₃ + cosθ₁*-sinθ₂*-sinθ₃,  cosθ₂*-sinθ₃,
     -sinθ₁*cosθ₂,                       cosθ₁*cosθ₂,                        sinθ₂,
      cosθ₁*sinθ₃ + sinθ₁*sinθ₂*cosθ₃,   sinθ₁*sinθ₃ + cosθ₁*-sinθ₂*cosθ₃,   cosθ₂*cosθ₃)
end

@inline function Base.:*(r::RotZXY, v::StaticVector)
    if length(v) != 3
        throw("Dimension mismatch: cannot rotate a vector of length $(length(v))")
    end

    sinθ₁, cosθ₁ = sincos(r.theta1)
    sinθ₂, cosθ₂ = sincos(r.theta2)
    sinθ₃, cosθ₃ = sincos(r.theta3)

    T = Base.promote_op(*, typeof(sinθ₁), eltype(v))

    return similar_type(v,T)(
        (cosθ₁*cosθ₃ + sinθ₁*sinθ₂*-sinθ₃)*v[1] + -sinθ₁*cosθ₂*v[2] + (cosθ₁*sinθ₃ + sinθ₁*sinθ₂*cosθ₃)*v[3],
        (sinθ₁*cosθ₃ + cosθ₁*-sinθ₂*-sinθ₃)*v[1] + cosθ₁*cosθ₂*v[2] + (sinθ₁*sinθ₃ + cosθ₁*-sinθ₂*cosθ₃)*v[3],
         cosθ₂*-sinθ₃*v[1] + sinθ₂*v[2] + cosθ₂*cosθ₃*v[3])
end
