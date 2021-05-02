"""
    MRP{T} <: Rotation

Modified Rodrigues Parameter. Is a 3D parameterization of attitude, and is a sterographic
projection of the 4D unit sphere onto the plane tangent to the negative real pole. They
have a singularity at θ = ±360°.

# Constructors
MRP(x, y, z)
MRP(r::AbstractVector)
"""
struct MRP{T} <: Rotation{3,T}
    x::T
    y::T
    z::T
    @inline MRP{T}(x,y,z) where T = new{T}(x,y,z)
    @inline MRP{T}(p::MRP) where T = new{T}(p.x, p.y, p.z)
end

# ~~~~~~~~~~~~~~~ Constructors ~~~~~~~~~~~~~~~ #
MRP(x::X, y::Y, z::Z) where {X,Y,Z} = MRP{promote_type(X,Y,Z)}(x, y, z)
(::Type{M})(g::StaticVector) where M<:MRP= M(g[1], g[2], g[3])

# ~~~~~~~~~~~~~~~ Conversions ~~~~~~~~~~~~~~~ #
params(g::MRP) = SVector{3}(g.x, g.y, g.z)

# ~~~~~~~~~~~~~~~ Initializers ~~~~~~~~~~~~~~~ #
Base.rand(::Type{RP}) where RP <: MRP = RP(rand(UnitQuaternion))
Base.one(::Type{RP}) where RP <: MRP = RP(0.0, 0.0, 0.0)

# ~~~~~~~~~~~~~~~~ Quaternion <=> MRP ~~~~~~~~~~~~~~~~~~ #
@inline function (::Type{Q})(p::MRP) where Q <: UnitQuaternion
    n2 = p.x * p.x + p.y * p.y + p.z * p.z
    M = 2/(1+n2)
    q = Q((1-n2)/(1+n2), M*p.x, M*p.y, M*p.z, false)
end

@inline function (::Type{RP})(q::UnitQuaternion) where RP<:MRP
    M = 1/(1+q.w)
    RP(q.x*M, q.y*M, q.z*M)
end

# ~~~~~~~~~~~~~~~ StaticArrays Interface ~~~~~~~~~~~~~~~ #
@inline (::Type{RP})(t::NTuple{9}) where RP<:MRP = convert(RP, UnitQuaternion(t))
@inline Base.getindex(p::MRP, i::Int) = convert(UnitQuaternion, p)[i]
@inline Base.Tuple(p::MRP) = Tuple(convert(UnitQuaternion, p))

# ~~~~~~~~~~~~~~~ Math Operations ~~~~~~~~~~~~~~~ #
LinearAlgebra.norm(p::MRP) = sqrt(p.x^2 + p.y^2 + p.z^2)
Base.inv(p::MRP) = MRP(-p.x, -p.y, -p.z)

# ~~~~~~~~~~~~~~~ Composition ~~~~~~~~~~~~~~~ #
function (*)(p2::MRP, p1::MRP)
    p2, p1 = params(p2), params(p1)
    MRP(((1-p2'p2)*p1 + (1-p1'p1)*p2 - cross(2p1, p2) ) / (1+p1'p1*p2'p2 - 2p1'p2))
end

function (\)(p1::MRP, p2::MRP)
    q1,q2 = params(p1), params(p2)
    n1,n2 = q1'q1, q2'q2
    θ = 1/((1+n1)*(1+n2))
    s1,s2 = (1-n1), (1-n2)
    v1 = -2q1
    v2 =  2q2
    s = s1*s2 - v1'v2
    v = s1*v2 + s2*v1 + v1 × v2

    M = θ/(1+θ*s)
    return MRP(v*M)
end

function (/)(p1::MRP, p2::MRP)
    q1,q2 = params(p1), params(p2)
    n1,n2 = q1'q1, q2'q2
    θ = 1/((1+n1)*(1+n2))
    s1,s2 = (1-n1), (1-n2)
    v1 =  2q1
    v2 = -2q2
    s = s1*s2 - v1'v2
    v = s1*v2 + s2*v1 + v1 × v2

    M = θ/(1+θ*s)
    return MRP(v*M)
end


# ~~~~~~~~~~~~~~~ Rotation ~~~~~~~~~~~~~~~ #
@inline (*)(p::MRP, r::StaticVector) = UnitQuaternion(p)*r
@inline (\)(p::MRP, r::StaticVector) = UnitQuaternion(p)\r


# ~~~~~~~~~~~~~~~ Kinematics ~~~~~~~~~~~~~~~ #
function kinematics(p::MRP, ω)
    # From Fundamentals of Spacecraft Attitude Determination and Control, eq (3.24)
    #     here A = (I(3) + 2*(skew(p)^2 + skew(p))/(1+p'p)) * (1+p'p)
    #            = √R * (1 + p'p)
    #     where √R*√R = (√R')*(√R') = RotMatrix(p)
    p = params(p)
    A = SA[  #
        1 + p[1]^2 - p[2]^2 - p[3]^2  2(p[1]*p[2] - p[3])      2(p[1]*p[3] + p[2]);
        2(p[2]*p[1] + p[3])            1-p[1]^2+p[2]^2-p[3]^2   2(p[2]*p[3] - p[1]);
        2(p[3]*p[1] - p[2])            2(p[3]*p[2] + p[1])      1-p[1]^2-p[2]^2+p[3]^2]
    0.25*A*ω
end


# ~~~~~~~~~~~~~~~ Useful Jacobians ~~~~~~~~~~~~~~~ #
function ∇rotate(p::MRP, r)
    p = params(p)
    4( (1-p'p)*skew(r)*(4p*p'/(1+p'p) - I) - (4/(1+p'p)*skew(p) + I)*2*skew(p)*r*p'
      - 2*(skew(p)*skew(r) + skew(skew(p)*r)))/(1+p'p)^2
end

function ∇composition2(p2::MRP, p1::MRP)
    p2,p1 = params(p2), params(p1)
    n1 = p1'p1
    n2 = p2'p2
    D = 1 / (1+n1*n2 - 2p1'p2)
    d1 = (-2p1*p2' + (1-n1)*I - skew(2p1) ) * D
    d2 = -((1-n2)*p1 + (1-n1)*p2 - cross(2p1, p2) ) * D^2 *
        (2p2*n1 - 2p1)'
    return d1 + d2
end

function ∇composition1(p2::MRP, p1::MRP)
    p2,p1 = params(p2), params(p1)
    n1 = p1'p1
    n2 = p2'p2
    D = 1 / (1+n1*n2 - 2p1'p2)
    d1 = ((1-n2)*I + -2p2*p1' + 2skew(p2) ) * D
    d2 = -((1-n2)*p1 + (1-n1)*p2 - cross(2p1, p2) ) * D^2 *
        (2p1*n2 - 2p2)'
    d1 + d2
end

function ∇²composition1(p2::MRP, p1::MRP, b::AbstractVector)
    check_length(b, 3)
    p2,p1 = params(p2), params(p1)
    n1 = p1'p1
    n2 = p2'p2
    D = 1 / (1+n1*n2 - 2p1'p2)  # scalar
    dD = -D^2 * (n2*2p1 - 2p2)  # 3x1 (correct)
    A = -((1-n2)*p1 + (1-n1)*p2 - cross(2p1, p2) )  # 3x1
    dA = -I*(1-n2) + 2p2*p1' - 2skew(p2)  # 3x3 (correct)
    B = 2(p1*n2 -  p2)  # 3x1
    dB = 2n2*b  # 3x1

    # d1b = ((1-n2)*I + -2p1*p2' - 2skew(p2) ) * D * b
    dd1b = -2I*(p2'b)*D + ((1-n2)*I + -2p1*p2' - 2skew(p2) )*b*dD' # 3x3 (correct)

    dd2 = D^2 * B * b' * dA +
         D^2 * A'b*2*(n2*I) +
         B*A'b* 2D * dD' # 3x3 (correct)
    dd1b + dd2
end

function ∇differential(p::MRP)
    p = params(p)
    n = p'p
    # p = params(p)
    # (1-n)I + 2(skew(p) + p*p')
    # @SMatrix [n + 2p.x^2      2(p.x*p.y-p.z)  2(p.x*p.z+p.y);
    #           2(p.y*p.x+p.z)  n + 2p.y^2      2(p.y*p.z-p.x);
    #           2(p.z*p.x-p.y)  2(p.z*p.y+p.x)  n + 2p.z^2]
    #
    # p2 = params(p)
    # n2 = p2'p2
    return (1-n)*I + 2(skew(p) + p*p')
end


function ∇²differential(p2::MRP, b::AbstractVector)
    check_length(b, 3)
    p2 = params(p2)
    n2 = p2'p2
    A = -p2  # 3x1
    B = -2p2  # 3x1
    D = 1
    dD = 2p2

    dA = -I*(1-n2) - 2skew(p2)  # 3x3
    dB = 2n2*I  # 3x3

    d1 = (-2p2'b*I*D) - (dA'b * dD')
    d2 = dB * A'b * D^2 +
        (-(1-n2)*b + (skew(2p2))*b)*B' * D^2 +
        B*A'b * 2D * dD'
    d1 + d2'
end

"""
    ∇err(p1::MRP, p2::MRP)

Jacobian of `p1\\p2` wrt `p2`
"""
function ∇err(p1::MRP, p2::MRP)
    q1,q2 = params(p1), params(p2)
    n1,n2 = q1'q1,   q2'q2
    θ = 1/((1+n1)*(1+n2))
    s1,s2 = (1-n1), (1-n2)
    v1 = -2q1
    v2 =  2q2
    s = s1*s2 - v1'v2
    v = s1*v2 + s2*v1 + v1 × v2
    dsdp = -2s1*q2 - 2v1
    dvdp = 2s1*I + -2v1*q2' + 2skew(v1)
    dθdp = -θ^2*(1+n1)*2q2

    M = θ/(1+θ*s)
    dMdp = 1/(1+θ*s)*dθdp - θ/(1+θ*s)^2*(dθdp*s + θ*dsdp)
    return dvdp*M + v*dMdp'
end

"""
    ∇²err(p1::MRP, p2::MRP, b::StaticVector{3})

Jacobian of `(∂/∂p p1\\p2)'b` wrt `p2`
"""
function ∇²err(p1::MRP, p2::MRP, b::AbstractVector)
    check_length(b, 3)
    q1,q2 = params(p1), params(p2)
    n1,n2 = q1'q1, q2'q2
    θ = 1/((1+n1)*(1+n2))
    s1,s2 = (1-n1), (1-n2)
    v1 = -2q1
    v2 =  2q2
    s = s1*s2 - v1'v2
    v = s1*v2 + s2*v1 + v1 × v2

    dsdp = -2s1*q2 - 2v1  # 3x1
    dsdp2 = -2s1*I  # 3x3

    dvdp = 2s1*b + -2q2*v1'b - 2skew(v1)*b
    dvdp2 = -I*2v1'b

    dθdp = -θ^2*(1+n1)*2q2  # 3x1
    dθdp2 = -2θ*(1+n1)*2q2*dθdp' - θ^2*(1+n1)*2I # 3x3

    M = θ/(1+θ*s)  # scalar
    dMdp = 1/(1+θ*s)*dθdp - θ/(1+θ*s)^2*(dθdp*s + θ*dsdp) # 3x1
    dM2 = θ/(1+θ*s)^2  # scalar
    dM3 = dθdp*s + θ*dsdp  # 3x1
    dM2dp = dθdp'/(1+θ*s)^2 - 2θ/(1+θ*s)^3 * (dθdp*s + θ*dsdp)'
    dM3dp = dθdp2*s + dθdp*dsdp' + dsdp*dθdp' + θ*dsdp2

    dMdp2 = -1/(1+θ*s)^2*dθdp*(dθdp*s + dsdp*θ)' + 1/(1+θ*s)*dθdp2  # good
    dMdp2 -= dM3*dM2dp + dM2*dM3dp

    vb = s1*v2'b + s2*v1'b + b'skew(v1)*v2  # scalar
    vpdp = s1*2b' - 2q2' * (v1'b)  + b'skew(v1)*2 # good
    # vpdp = s1*2b' - 2q2 * (v1'b)

    d1 = M*dvdp2 + dvdp*dMdp'
    d2 = dMdp2*vb + dMdp*vpdp
    return d1 + d2
end
