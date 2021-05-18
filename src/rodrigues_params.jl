
"""
    RodriguesParam{T}

Rodrigues parameters are a three-dimensional parameterization of rotations.
They have a singularity at 180° but do not inherit the sign ambiguities of quaternions
or MRPs
"""
struct RodriguesParam{T} <: Rotation{3,T}
    x::T
    y::T
    z::T
    @inline function RodriguesParam{T}(x, y, z) where T
        new{T}(x, y, z)
    end
end

# ~~~~~~~~~~~~~~~ Constructors ~~~~~~~~~~~~~~~ #
RodriguesParam(x::X, y::Y, z::Z) where {X,Y,Z} = RodriguesParam{promote_type(X,Y,Z)}(x, y, z)
(::Type{RP})(g::StaticVector) where RP<:RodriguesParam = RP(g[1], g[2], g[3])

# ~~~~~~~~~~~~~~~ Conversions ~~~~~~~~~~~~~~~ #
params(g::RodriguesParam) = SVector{3}(g.x, g.y, g.z)

# ~~~~~~~~~~~~~~~ Initializers ~~~~~~~~~~~~~~~ #
@inline function Random.rand(rng::AbstractRNG, ::Random.SamplerType{RP}) where RP <: RodriguesParam
    RP(rand(rng, UnitQuaternion))
end
@inline Base.one(::Type{RP}) where RP <: RodriguesParam = RP(0.0, 0.0, 0.0)

# ~~~~~~~~~~~~~~~ Quaternion <=> RP ~~~~~~~~~~~~~~~~~~ #
@inline function (::Type{Q})(g::RodriguesParam) where Q<:UnitQuaternion
    M = 1/sqrt(1 + g.x*g.x + g.y*g.y + g.z*g.z)
    q = Q(M, M*g.x, M*g.y, M*g.z, false)
end

@inline function (::Type{G})(q::UnitQuaternion) where G<:RodriguesParam
    M = 1/q.w
    G(q.x*M, q.y*M, q.z*M)
end

# ~~~~~~~~~~~~~~~ StaticArrays Interface ~~~~~~~~~~~~~~~ #
@inline (::Type{RP})(t::NTuple{9}) where RP<:RodriguesParam = convert(RP, UnitQuaternion(t))
@inline Base.getindex(rp::RodriguesParam, i::Int) = convert(UnitQuaternion, rp)[i]
@inline Base.Tuple(rp::RodriguesParam) = Tuple(convert(UnitQuaternion, rp))

# ~~~~~~~~~~~~~~~ Math Operations ~~~~~~~~~~~~~~~ #
LinearAlgebra.norm(g::RodriguesParam) = sqrt(g.x^2 + g.y^2 + g.z^2)


# ~~~~~~~~~~~~~~~ Composition ~~~~~~~~~~~~~~~ #
function (*)(g2::RodriguesParam, g1::RodriguesParam)
    g2 = params(g2)
    g1 = params(g1)
    RodriguesParam((g2+g1 + g2 × g1)/(1-g2'g1))
end

function (\)(g1::RodriguesParam, g2::RodriguesParam)
    g2 = params(g2)
    g1 = params(g1)
    RodriguesParam((g2-g1 + g2 × g1)/(1+g1'g2))
end

function (/)(g1::RodriguesParam, g2::RodriguesParam)
    g2 = params(g2)
    g1 = params(g1)
    RodriguesParam((g1-g2 + g2 × g1)/(1+g1'g2))
end

# ~~~~~~~~~~~~~~~ Rotation ~~~~~~~~~~~~~~~ #
(*)(g::RodriguesParam, r::StaticVector) = UnitQuaternion(g)*r
(\)(g::RodriguesParam, r::StaticVector) = inv(UnitQuaternion(g))*r


# ~~~~~~~~~~~~~~~ Kinematics ~~~~~~~~~~~~~~~ #
function kinematics(g::RodriguesParam, ω)
    check_length(ω, 3)
    g = params(g)
    0.5*(I + skew(g) + g*g')*ω
end


function ∇rotate(g0::RodriguesParam, r)
    check_length(r, 3)
    g = params(g0)
    ghat = skew(g)
    n1 = 1/(1 + g'g)
    gxr = cross(g,r) + r
    d1 = ghat*gxr * -2*n1^2 * g'
    d2 = -(ghat*skew(r) + skew(gxr))*n1
    return 2d1 + 2d2
end

function ∇composition1(g2::RodriguesParam, g1::RodriguesParam)
    g2 = params(g2)
    g1 = params(g1)

    N = g2 + g1 + g2 × g1
    D = 1/(1 - g2'g1)
    (I + skew(g2) + D*N*g2')*D
end

function ∇²composition1(g2::RodriguesParam, g1::RodriguesParam, b::AbstractVector)
    check_length(b, 3)
    g2 = params(g2)
    g1 = params(g1)

    N = g2 + g1 + g2 × g1  # 3x1
    D = 1/(1 - g2'g1)  # scalar
    dN = I + skew(g2)
    dD = D^2*g2'
    return g2*b'*(N*(2*D*dD) + D^2*dN) + (I - skew(g2))*b*dD
end

function ∇composition2(g2::RodriguesParam, g1::RodriguesParam)
    g2 = params(g2)
    g1 = params(g1)

    N = g2 + g1 + g2 × g1
    D = 1/(1 - g2'g1)
    (I - skew(g1) + D*N*g1')*D
end

function ∇differential(g::RodriguesParam)
    g = params(g)
    (I + skew(g) + g*g')
end

function ∇²differential(g::RodriguesParam, b::AbstractVector)
    check_length(b, 3)
    g = params(g)
    return g*b'*(2g*g' + I + skew(g)) + (I - skew(g))*b*g'
end
