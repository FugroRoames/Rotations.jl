################################################################################
################################################################################
"""
    struct AngleAxis{T} <: Rotation{3,T}
    AngleAxis(Θ, x, y, z)

A 3×3 rotation matrix parameterized by a 3D rotation by angle θ about an
arbitrary axis `[x, y, z]`.

Note that the axis is not unique for θ = 0, and that this parameterization does
not continuously map the neighbourhood of the null rotation (and therefore
might not be suitable for autodifferentation and optimization purposes).

Note: by default, the constructor will renormalize the input so that the axis
has length 1 (x² + y² + z² = 1).

Renormalization can be skipped by passing `false` as an additional constructor
argument, in which case the user provides the guarantee that the input arguments
represent a normalized rotation axis. Operations on an `AngleAxis` with a rotation
axis that does not have unit norm, created by skipping renormalization in this fashion,
are not guaranteed to do anything sensible.
"""
struct AngleAxis{T} <: Rotation{3,T}
    theta::T
    axis_x::T
    axis_y::T
    axis_z::T

    @inline function AngleAxis{T}(θ, x, y, z, normalize::Bool = true) where {T}
        if normalize
            # Not sure what to do with theta?? Should it become theta * norm ?
            norm = sqrt(x*x + y*y + z*z)
            new(θ, x/norm, y/norm, z/norm)
        else
            new(θ, x, y, z)
        end
    end
end

params(aa::AngleAxis) = SVector{4}(aa.theta, aa.axis_x, aa.axis_y, aa.axis_z)

# StaticArrays will take over *all* the constructors and put everything in a tuple...
# but this isn't quite what we mean when we have 4 inputs (not 9).
@inline function AngleAxis(θ::Θ, x::X, y::Y, z::Z, normalize::Bool = true) where {Θ,X,Y,Z}
    AngleAxis{promote_type(promote_type(promote_type(Θ, X), Y), Z)}(θ, x, y, z, normalize)
end

# These functions are enough to satisfy the entire StaticArrays interface:
@inline (::Type{AA})(t::NTuple{9}) where {AA <: AngleAxis} = AA(UnitQuaternion(t)) # TODO: consider going directly from tuple (RotMatrix) to AngleAxis
@inline Base.getindex(aa::AngleAxis, i::Int) = UnitQuaternion(aa)[i]

@inline function Base.Tuple(aa::AngleAxis{T}) where T
    # Rodrigues' rotation formula.
    s, c = sincos(aa.theta)
    c1 = one(T) - c

    c1x2 = c1 * aa.axis_x^2
    c1y2 = c1 * aa.axis_y^2
    c1z2 = c1 * aa.axis_z^2

    c1xy = c1 * aa.axis_x * aa.axis_y
    c1xz = c1 * aa.axis_x * aa.axis_z
    c1yz = c1 * aa.axis_y * aa.axis_z

    sx = s * aa.axis_x
    sy = s * aa.axis_y
    sz = s * aa.axis_z

    # Note that the RotMatrix constructor argument order makes this look transposed:
    (one(T) - c1y2 - c1z2, c1xy + sz, c1xz - sy,
        c1xy - sz, one(T) - c1x2 - c1z2, c1yz + sx,
        c1xz + sy, c1yz - sx, one(T) - c1x2 - c1y2)
end

@inline function (::Type{Q})(aa::AngleAxis) where Q <: UnitQuaternion
    s, c = sincos(aa.theta / 2)
    return Q(c, s * aa.axis_x, s * aa.axis_y, s * aa.axis_z, false)
end

@inline function (::Type{AA})(q::UnitQuaternion) where AA <: AngleAxis
    s2 = q.x * q.x + q.y * q.y + q.z * q.z
    sin_t2 = sqrt(s2)
    theta = 2 * atan(sin_t2, q.w)
    num_pert = eps(typeof(theta))^4
    inv_sin_t2 = 1 / (sin_t2 + num_pert)
    return principal_value(AA(theta, inv_sin_t2 * (q.x + num_pert), inv_sin_t2 * q.y, inv_sin_t2 * q.z, false))
end

# Using Rodrigues formula on an AngleAxis parameterization (assume unit axis length) to do the rotation
# (implementation from: https://ceres-solver.googlesource.com/ceres-solver/+/1.10.0/include/ceres/rotation.h)
function Base.:*(aa::AngleAxis, v::StaticVector)
    if length(v) != 3
        throw("Dimension mismatch: cannot rotate a vector of length $(length(v))")
    end

    w = rotation_axis(aa)
    st, ct = sincos(aa.theta)
    w_cross_pt = cross(w, v)
    m = dot(v, w) * (one(w_cross_pt[1]) - ct)
    T = promote_type(eltype(aa), eltype(v))
    return similar_type(v,T)(v[1] * ct + w_cross_pt[1] * st + w[1] * m,
                             v[2] * ct + w_cross_pt[2] * st + w[2] * m,
                             v[3] * ct + w_cross_pt[3] * st + w[3] * m)
end

@inline Base.:*(aa::AngleAxis, r::Rotation) = UnitQuaternion(aa) * r
@inline Base.:*(aa::AngleAxis, r::RotMatrix) = UnitQuaternion(aa) * r
@inline Base.:*(aa::AngleAxis, r::MRP) = UnitQuaternion(aa) * r
@inline Base.:*(r::Rotation, aa::AngleAxis) = r * UnitQuaternion(aa)
@inline Base.:*(r::RotMatrix, aa::AngleAxis) = r * UnitQuaternion(aa)
@inline Base.:*(r::MRP, aa::AngleAxis) = r * UnitQuaternion(aa)
@inline Base.:*(aa1::AngleAxis, aa2::AngleAxis) = UnitQuaternion(aa1) * UnitQuaternion(aa2)

@inline Base.inv(aa::AngleAxis) = AngleAxis(-aa.theta, aa.axis_x, aa.axis_y, aa.axis_z)
@inline Base.:^(aa::AngleAxis, t::Real) = AngleAxis(aa.theta*t, aa.axis_x, aa.axis_y, aa.axis_z)
@inline Base.:^(aa::AngleAxis, t::Integer) = AngleAxis(aa.theta*t, aa.axis_x, aa.axis_y, aa.axis_z) # to avoid ambiguity


# define null rotations for convenience
@inline Base.one(::Type{AngleAxis}) = AngleAxis(0.0, 1.0, 0.0, 0.0)
@inline Base.one(::Type{AngleAxis{T}}) where {T} = AngleAxis{T}(zero(T), one(T), zero(T), zero(T))

# accessors
@inline rotation_angle(aa::AngleAxis) = aa.theta #  - floor((aa.theta+pi) / (2*pi)) * 2*pi
@inline rotation_axis(aa::AngleAxis) = SVector(aa.axis_x, aa.axis_y, aa.axis_z)


################################################################################
################################################################################
"""
    struct RotationVec{T} <: Rotation{3,T}
    RotationVec(sx, sy, sz)

Rodrigues vector parameterization of a 3×3 rotation matrix. The direction of the
vector [sx, sy, sz] defines the axis of rotation, and the rotation angle is
given by its norm.
"""
struct RotationVec{T} <: Rotation{3,T}
    sx::T
    sy::T
    sz::T
end

params(aa::RotationVec) = SVector{3}(aa.sx, aa.sy, aa.sz)

# StaticArrays will take over *all* the constructors and put everything in a tuple...
# but this isn't quite what we mean when we have 4 inputs (not 9).
@inline RotationVec(x::X, y::Y, z::Z) where {X,Y,Z} = RotationVec{promote_type(promote_type(X, Y), Z)}(x, y, z)

# These functions are enough to satisfy the entire StaticArrays interface:
@inline (::Type{RV})(t::NTuple{9}) where {RV <: RotationVec} = RV(UnitQuaternion(t)) # TODO: go through AngleAxis once it's faster
@inline Base.getindex(aa::RotationVec, i::Int) = UnitQuaternion(aa)[i]
@inline Base.Tuple(rv::RotationVec) = Tuple(UnitQuaternion(rv))

function (::Type{AA})(rv::RotationVec) where AA <: AngleAxis
    # TODO: consider how to deal with derivative near theta = 0. There should be a first-order expansion here.
    theta = rotation_angle(rv)
    return theta > 0 ? AA(theta, rv.sx / theta, rv.sy / theta, rv.sz / theta, false) : AA(zero(theta), one(theta), zero(theta), zero(theta), false)
end

function (::Type{RV})(aa::AngleAxis) where RV <: RotationVec
    return RV(aa.theta * aa.axis_x, aa.theta * aa.axis_y, aa.theta * aa.axis_z)
end

function (::Type{Q})(rv::RotationVec) where Q <: UnitQuaternion
    return UnitQuaternion(AngleAxis(rv))
end

(::Type{RV})(q::UnitQuaternion) where {RV <: RotationVec} = RV(AngleAxis(q))

function Base.:*(rv::RotationVec{T1}, v::StaticVector{3, T2}) where {T1,T2}
    theta = rotation_angle(rv)
    if (theta > eps(T1)) # use eps here because we have the 1st order series expansion defined
        return AngleAxis(rv) * v
    else
        return similar_type(typeof(v), promote_type(T1,T2))(
                    v[1] + rv.sy * v[3] - rv.sz * v[2],
                    v[2] + rv.sz * v[1] - rv.sx * v[3],
                    v[3] + rv.sx * v[2] - rv.sy * v[1])
    end
end

@inline Base.:*(rv::RotationVec, r::Rotation) = UnitQuaternion(rv) * r
@inline Base.:*(rv::RotationVec, r::RotMatrix) = UnitQuaternion(rv) * r
@inline Base.:*(rv::RotationVec, r::MRP) = UnitQuaternion(rv) * r
@inline Base.:*(rv::RotationVec, r::AngleAxis) = UnitQuaternion(rv) * r
@inline Base.:*(r::Rotation, rv::RotationVec) = r * UnitQuaternion(rv)
@inline Base.:*(r::RotMatrix, rv::RotationVec) = r * UnitQuaternion(rv)
@inline Base.:*(r::MRP, rv::RotationVec) = r * UnitQuaternion(rv)
@inline Base.:*(r::AngleAxis, rv::RotationVec) = r * UnitQuaternion(rv)
@inline Base.:*(rv1::RotationVec, rv2::RotationVec) = UnitQuaternion(rv1) * UnitQuaternion(rv2)

@inline Base.inv(rv::RotationVec) = RotationVec(-rv.sx, -rv.sy, -rv.sz)
@inline Base.:^(rv::RotationVec, t::Real) = RotationVec(rv.sx*t, rv.sy*t, rv.sz*t)
@inline Base.:^(rv::RotationVec, t::Integer) = RotationVec(rv.sx*t, rv.sy*t, rv.sz*t) # to avoid ambiguity

# rotation properties
@inline rotation_angle(rv::RotationVec) = sqrt(rv.sx * rv.sx + rv.sy * rv.sy + rv.sz * rv.sz)
function rotation_axis(rv::RotationVec)     # what should this return for theta = 0?
    theta = rotation_angle(rv)
    return (theta > 0 ? SVector(rv.sx / theta, rv.sy / theta, rv.sz / theta) : SVector(one(theta), zero(theta), zero(theta)))
end

# define null rotations for convenience
@inline Base.one(::Type{RotationVec}) = RotationVec(0.0, 0.0, 0.0)
@inline Base.one(::Type{RotationVec{T}}) where {T} = RotationVec{T}(zero(T), zero(T), zero(T))
