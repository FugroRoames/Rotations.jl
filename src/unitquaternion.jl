import Base: +, -, *, /, \, exp, log, ≈, ==, inv, conj

"""
    UnitQuaternion{T} <: Rotation

4-parameter attitute representation that is singularity-free. Quaternions with unit norm
represent a double-cover of SO(3). The `UnitQuaternion` does NOT strictly enforce the unit
norm constraint, but certain methods will assume you have a unit quaternion.
Follows the Hamilton convention for quaternions.

# Constructors
```julia
UnitQuaternion(w,x,y,z)
UnitQuaternion(q::AbstractVector)
```
where `w` is the scalar (real) part, `x`,`y`, and `z` are the vector (imaginary) part,
and `q = [w,x,y,z]`.
"""
struct UnitQuaternion{T} <: Rotation{3,T}
    w::T
    x::T
    y::T
    z::T

    @inline function UnitQuaternion{T}(w, x, y, z, normalize::Bool = true) where T
        if normalize
            inorm = inv(sqrt(w*w + x*x + y*y + z*z))
            new{T}(w*inorm, x*inorm, y*inorm, z*inorm)
        else
            new{T}(w, x, y, z)
        end
    end

    UnitQuaternion{T}(q::UnitQuaternion) where T = new{T}(q.w, q.x, q.y, q.z)
end

# ~~~~~~~~~~~~~~~ Constructors ~~~~~~~~~~~~~~~ #
# Use default map
function UnitQuaternion(w,x,y,z, normalize::Bool = true)
    types = promote(w,x,y,z)
    UnitQuaternion{eltype(types)}(w,x,y,z, normalize)
end

# Pass in Vectors
@inline function (::Type{Q})(q::AbstractVector, normalize::Bool = true) where Q <: UnitQuaternion
    check_length(q, 4)
    Q(q[1], q[2], q[3], q[4], normalize)
end
@inline (::Type{Q})(q::StaticVector{4}, normalize::Bool = true) where Q <: UnitQuaternion =
    Q(q[1], q[2], q[3], q[4], normalize)

# Copy constructors
UnitQuaternion(q::UnitQuaternion) = q

# UnitQuaternion <=> Quat
# (::Type{Q})(q::Quat) where Q <: UnitQuaternion = Q(q.w, q.x, q.y, q.z, false)
# (::Type{Q})(q::UnitQuaternion) where Q <: Quat = Q(q.w, q.x, q.y, q.z, false)
# const AllQuats{T} = Union{<:Quat{T}, <:UnitQuaternion{T}}


# ~~~~~~~~~~~~~~~ StaticArrays Interface ~~~~~~~~~~~~~~~ #
function (::Type{Q})(t::NTuple{9}) where Q<:UnitQuaternion
    #=
    This function solves the system of equations in Section 3.1
    of https://arxiv.org/pdf/math/0701759.pdf. This cheap method
    only works for matrices that are already orthonormal (orthogonal
    and unit length columns). The nearest orthonormal matrix can
    be found by solving Wahba's problem:
    https://en.wikipedia.org/wiki/Wahba%27s_problem as shown below.

    not_orthogonal = randn(3,3)
    u,s,v = svd(not_orthogonal)
    is_orthogoral = u * diagm([1, 1, sign(det(u * transpose(v)))]) * transpose(v)
    =#

    a = 1 + t[1] + t[5] + t[9]
    b = 1 + t[1] - t[5] - t[9]
    c = 1 - t[1] + t[5] - t[9]
    d = 1 - t[1] - t[5] + t[9]
    max_abcd = max(a, b, c, d)
    if a == max_abcd
        b = t[6] - t[8]
        c = t[7] - t[3]
        d = t[2] - t[4]
    elseif b == max_abcd
        a = t[6] - t[8]
        c = t[2] + t[4]
        d = t[7] + t[3]
    elseif c == max_abcd
        a = t[7] - t[3]
        b = t[2] + t[4]
        d = t[6] + t[8]
    else
        a = t[2] - t[4]
        b = t[7] + t[3]
        c = t[6] + t[8]
    end
    return Q(a, b, c, d)
end


function Base.getindex(q::UnitQuaternion, i::Int)
    if i == 1
        ww = (q.w * q.w)
        xx = (q.x * q.x)
        yy = (q.y * q.y)
        zz = (q.z * q.z)

        ww + xx - yy - zz
    elseif i == 2
        xy = (q.x * q.y)
        zw = (q.w * q.z)

        2 * (xy + zw)
    elseif i == 3
        xz = (q.x * q.z)
        yw = (q.y * q.w)

        2 * (xz - yw)
    elseif i == 4
        xy = (q.x * q.y)
        zw = (q.w * q.z)

        2 * (xy - zw)
    elseif i == 5
        ww = (q.w * q.w)
        xx = (q.x * q.x)
        yy = (q.y * q.y)
        zz = (q.z * q.z)

        ww - xx + yy - zz
    elseif i == 6
        yz = (q.y * q.z)
        xw = (q.w * q.x)

        2 * (yz + xw)
    elseif i == 7
        xz = (q.x * q.z)
        yw = (q.y * q.w)

        2 * (xz + yw)
    elseif i == 8
        yz = (q.y * q.z)
        xw = (q.w * q.x)

        2 * (yz - xw)
    elseif i == 9
        ww = (q.w * q.w)
        xx = (q.x * q.x)
        yy = (q.y * q.y)
        zz = (q.z * q.z)

        ww - xx - yy + zz
    else
        throw(BoundsError(r,i))
    end
end

function Base.Tuple(q::UnitQuaternion)
    ww = (q.w * q.w)
    xx = (q.x * q.x)
    yy = (q.y * q.y)
    zz = (q.z * q.z)
    xy = (q.x * q.y)
    zw = (q.w * q.z)
    xz = (q.x * q.z)
    yw = (q.y * q.w)
    yz = (q.y * q.z)
    xw = (q.w * q.x)

    # initialize rotation part
    return (ww + xx - yy - zz,
            2 * (xy + zw),
            2 * (xz - yw),
            2 * (xy - zw),
            ww - xx + yy - zz,
            2 * (yz + xw),
            2 * (xz + yw),
            2 * (yz - xw),
            ww - xx - yy + zz)
end

# ~~~~~~~~~~~~~~~ Getters ~~~~~~~~~~~~~~~ #
@inline scalar(q::UnitQuaternion) = q.w
@inline vector(q::UnitQuaternion) = SVector{3}(q.x, q.y, q.z)

"""
    params(R::Rotation)

Return an `SVector` of the underlying parameters used by the rotation representation.

# Example
```julia
p = MRP(1.0, 2.0, 3.0)
Rotations.params(p) == @SVector [1.0, 2.0, 3.0]  # true
```
"""
@inline params(q::UnitQuaternion) = SVector{4}(q.w, q.x, q.y, q.z)

# ~~~~~~~~~~~~~~~ Initializers ~~~~~~~~~~~~~~~ #
Base.rand(::Type{<:UnitQuaternion{T}}) where T =
    normalize(UnitQuaternion{T}(randn(T), randn(T), randn(T), randn(T)))
Base.rand(::Type{UnitQuaternion}) = Base.rand(UnitQuaternion{Float64})
@inline Base.one(::Type{Q}) where Q <: UnitQuaternion = Q(1.0, 0.0, 0.0, 0.0)


# ~~~~~~~~~~~~~~~ Math Operations ~~~~~~~~~~~~~~~ #

# Inverses
conj(q::Q) where Q <: UnitQuaternion = Q(q.w, -q.x, -q.y, -q.z, false)
inv(q::UnitQuaternion) = conj(q)
(-)(q::Q) where Q <: UnitQuaternion = Q(-q.w, -q.x, -q.y, -q.z, false)

# Norms
LinearAlgebra.norm(q::UnitQuaternion) = sqrt(q.w^2 + q.x^2 + q.y^2 + q.z^2)
vecnorm(q::UnitQuaternion) = sqrt(q.x^2 + q.y^2 + q.z^2)

function LinearAlgebra.normalize(q::Q) where Q <: UnitQuaternion
    n = inv(norm(q))
    Q(q.w*n, q.x*n, q.y*n, q.z*n)
end

# Identity
(::Type{Q})(I::UniformScaling) where Q <: UnitQuaternion = one(Q)

# Exponentials and Logarithms
"""
    pure_quaternion(v::AbstractVector)
    pure_quaternion(x, y, z)

Create a `UnitQuaternion` with zero scalar part (i.e. `q.w == 0`).
"""
function pure_quaternion(v::AbstractVector)
    check_length(v, 3)
    UnitQuaternion(zero(eltype(v)), v[1], v[2], v[3], false)
end

@inline pure_quaternion(x::Real, y::Real, z::Real) =
    UnitQuaternion(zero(x), x, y, z, false)

function exp(q::Q) where Q <: UnitQuaternion
    θ = vecnorm(q)
    sθ,cθ = sincos(θ)
    es = exp(q.w)
    M = es*sθ/θ
    Q(es*cθ, q.x*M, q.y*M, q.z*M, false)
end

function expm(ϕ::AbstractVector)
    check_length(ϕ, 3)
    θ = norm(ϕ)
    sθ,cθ = sincos(θ/2)
    M = 1//2 *sinc(θ/π/2)
    UnitQuaternion(cθ, ϕ[1]*M, ϕ[2]*M, ϕ[3]*M, false)
end

function log(q::Q, eps=1e-6) where Q <: UnitQuaternion
    # Assumes unit quaternion
    θ = vecnorm(q)
    if θ > eps
        M = atan(θ, q.w)/θ
    else
        M = (1-(θ^2/(3q.w^2)))/q.w
    end
    pure_quaternion(M*vector(q))
end

function logm(q::UnitQuaternion)
    # Assumes unit quaternion
    2*vector(log(q))
end

# Composition
"""
    (*)(q::UnitQuaternion, w::UnitQuaternion)

Quternion Composition

Equivalent to
```julia
lmult(q) * SVector(w)
rmult(w) * SVector(q)
```

Sets the output mapping equal to the mapping of `w`
"""
function (*)(q::UnitQuaternion, w::UnitQuaternion)
    UnitQuaternion(q.w * w.w - q.x * w.x - q.y * w.y - q.z * w.z,
                   q.w * w.x + q.x * w.w + q.y * w.z - q.z * w.y,
                   q.w * w.y - q.x * w.z + q.y * w.w + q.z * w.x,
                   q.w * w.z + q.x * w.y - q.y * w.x + q.z * w.w, false)
end

"""
    (*)(q::UnitQuaternion, r::StaticVector)

Rotate a vector

Equivalent to `hmat()' lmult(q) * rmult(q)' hmat() * r`
"""
function Base.:*(q::UnitQuaternion, r::StaticVector)  # must be StaticVector to avoid ambiguity
    check_length(r, 3)
    w = q.w
    v = vector(q)
    (w^2 - v'v)*r + 2*v*(v'r) + 2*w*cross(v,r)
end

"""
    (*)(q::UnitQuaternion, w::Real)

Scalar multiplication of a quaternion. Breaks unit norm.
"""
function (*)(q::Q, w::Real) where Q<:UnitQuaternion
    return Q(q.w*w, q.x*w, q.y*w, q.z*w, false)
end
(*)(w::Real, q::UnitQuaternion) = q*w



(\)(q1::UnitQuaternion, q2::UnitQuaternion) = conj(q1)*q2  # Equivalent to inv(q1)*q2
(/)(q1::UnitQuaternion, q2::UnitQuaternion) = q1*conj(q2)  # Equivalent to q1*inv(q2)

(\)(q::UnitQuaternion, r::SVector{3}) = conj(q)*r          # Equivalent to inv(q)*r

"""
    rotation_between(from, to)

Compute the quaternion that rotates vector `from` so that it aligns with vector
`to`, along the geodesic (shortest path).
"""
rotation_between(from::AbstractVector, to::AbstractVector) = rotation_between(SVector{3}(from), SVector{3}(to))
function rotation_between(from::SVector{3}, to::SVector{3})
    # Robustified version of implementation from https://www.gamedev.net/topic/429507-finding-the-quaternion-betwee-two-vectors/#entry3856228
    normprod = sqrt(dot(from, from) * dot(to, to))
    T = typeof(normprod)
    normprod < eps(T) && throw(ArgumentError("Input vectors must be nonzero."))
    w = normprod + dot(from, to)
    v = abs(w) < 100 * eps(T) ? perpendicular_vector(from) : cross(from, to)
    @inbounds return UnitQuaternion(w, v[1], v[2], v[3]) # relies on normalization in constructor
end

# ~~~~~~~~~~~~~~~ Kinematics ~~~~~~~~~~~~~~~ $
"""
    kinematics(R::Rotation{3}, ω::AbstractVector)

The time derivative of the rotation R, according to the definition

``Ṙ = \\lim_{Δt → 0} \\frac{R(t + Δt) - R(t)}{Δt}``

where `ω` is the angular velocity. This is equivalent to

``Ṙ = \\lim_{Δt → 0} \\frac{R δR - R}{Δt}``

where ``δR`` is some small rotation, parameterized by a small rotation ``δθ`` about
an axis ``r``, such that ``lim_{Δt → 0} \\frac{δθ r}{Δt} = ω``

The kinematics are extremely useful when computing the dynamics of rigid bodies, since
`Ṙ = kinematics(R,ω)` is the first-order ODE for the evolution of the attitude dynamics.

See "Fundamentals of Spacecraft Attitude Determination and Control" by Markley and Crassidis
Sections 3.1-3.2 for more details.
"""
function kinematics(q::Q, ω::AbstractVector) where Q <: UnitQuaternion
    1//2 * params(q*Q(0.0, ω[1], ω[2], ω[3], false))
end

# ~~~~~~~~~~~~~~~ Linear Algebraic Conversions ~~~~~~~~~~~~~~~ #
"""
    lmult(q::UnitQuaternion)
    lmult(q::StaticVector{4})

`lmult(q2)*params(q1)` returns a vector equivalent to `q2*q1` (quaternion composition)
"""
function lmult(q::UnitQuaternion)
    SA[
        q.w -q.x -q.y -q.z;
        q.x  q.w -q.z  q.y;
        q.y  q.z  q.w -q.x;
        q.z -q.y  q.x  q.w;
    ]
end
lmult(q::StaticVector{4}) = lmult(UnitQuaternion(q, false))

"""
    rmult(q::UnitQuaternion)
    rmult(q::StaticVector{4})

`rmult(q1)*params(q2)` return a vector equivalent to `q2*q1` (quaternion composition)
"""
function rmult(q::UnitQuaternion)
    SA[
        q.w -q.x -q.y -q.z;
        q.x  q.w  q.z -q.y;
        q.y -q.z  q.w  q.x;
        q.z  q.y -q.x  q.w;
    ]
end
rmult(q::SVector{4}) = rmult(UnitQuaternion(q, false))

"""
    tmat()

`tmat()*params(q)`return a vector equivalent to `inv(q)`, where `q` is a `UnitQuaternion`
"""
function tmat(::Type{T}=Float64) where T
    SA{T}[
        1  0  0  0;
        0 -1  0  0;
        0  0 -1  0;
        0  0  0 -1;
    ]
end

"""
    vmat()

`vmat()*params(q)`` returns the imaginary
    (vector) part of the quaternion `q` (equivalent to `vector(q)``)
"""
function vmat(::Type{T}=Float64) where T
    SA{T}[
        0 1 0 0;
        0 0 1 0;
        0 0 0 1
    ]
end

"""
    hmat()
    hmat(r::AbstractVector)

`hmat()*r` or `hmat(r)` converts `r` into a pure quaternion, where `r` is 3-dimensional.

`hmat() == vmat()'`
"""
function hmat(::Type{T}=Float64) where T
    SA{T}[
        0 0 0;
        1 0 0;
        0 1 0;
        0 0 1.;
    ]
end

function hmat(r)
    @assert length(r) == 3
    SA[0, r[1], r[2], r[3]]
end


# ~~~~~~~~~~~~~~~ Useful Jacobians ~~~~~~~~~~~~~~~ #
"""
    ∇differential(q::UnitQuaternion)

Jacobian of `lmult(q) QuatMap(ϕ)`, when ϕ is near zero.

Useful for converting Jacobians from R⁴ to R³ and
    correctly account for unit norm constraint. Jacobians for different
    differential quaternion parameterization are the same up to a constant.
"""
function ∇differential(q::UnitQuaternion)
    SA[
        -q.x -q.y -q.z;
         q.w -q.z  q.y;
         q.z  q.w -q.x;
        -q.y  q.x  q.w;
    ]
end

"""
    ∇²differential(q::UnitQuaternion, b::AbstractVector)

Jacobian of `(∂/∂ϕ lmult(q) QuatMap(ϕ))`b, evaluated at ϕ=0, and `b` has length 4.
"""
function ∇²differential(q::UnitQuaternion, b::AbstractVector)
    check_length(b, 4)
    b1 = -params(q)'b
    Diagonal(@SVector fill(b1,3))
end

"""
    ∇rotate(R::Rotation{3}, r::AbstractVector)

Jacobian of `R*r` with respect to the rotation
"""
function ∇rotate(q::UnitQuaternion, r::AbstractVector)
    check_length(r, 3)
    rhat = UnitQuaternion(zero(eltype(r)), r[1], r[2], r[3], false)
    2vmat()*rmult(q)'rmult(rhat)
end

"""
    ∇composition1(R2::Rotation{3}, R1::Rotation{3})

Jacobian of `R2*R1` with respect to `R1`
"""
function ∇composition1(q2::UnitQuaternion, q1::UnitQuaternion)
    lmult(q2)
end

"""
    ∇composition2(R2::Rotation{3}, R1::Rotation{3})

Jacobian of `R2*R1` with respect to `R2`
"""
function ∇composition2(q2::UnitQuaternion, q1::UnitQuaternion)
    rmult(q1)
end
