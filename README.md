# Rotations.jl

*3D rotations made easy in Julia*

[![Build Status](https://github.com/JuliaGeometry/Rotations.jl/workflows/CI/badge.svg)](https://github.com/JuliaGeometry/Rotations.jl/actions?query=workflow%3ACI)
[![Coverage](https://codecov.io/gh/JuliaGeometry/Rotations.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaGeometry/Rotations.jl)

This package implements various 3D rotation parameterizations and defines
conversions between them. At their heart, each rotation parameterization *is*
a 3×3 unitary (orthogonal) matrix (based on the [StaticArrays.jl package](https://github.com/JuliaArrays/StaticArrays.jl)),
and acts to rotate a 3-vector about the origin through matrix-vector multiplication.

While the `RotMatrix` type is a dense representation of a `3×3` matrix, we also
have sparse (or computed, rather) representations such as quaternions,
angle-axis parameterizations, and Euler angles.

For composing rotations about the origin with other transformations, see
the [CoordinateTransformations.jl](https://github.com/JuliaGeometry/CoordinateTransformations.jl)
package.

### Interface
The following operations are supported by all of the implemented rotation parameterizations.

#### Composition
Any two rotations of the same type can be composed with simple multiplication:
```julia
q3 = q2*q1
```
Rotations can be composed with the opposite (or inverse) rotation with the appropriate
division operation
```julia
q1 = q2\q3
q2 = q3/q1
```

#### Rotation
Any rotation can operate on a 3D vector (represented as a `SVector{3}`), again through
simple multiplication:
```julia
r2 = q*r
```
which also supports multiplication by the opposite rotation
```julia
r = q\r2
```

#### Rotation Angle / Axis
The rotation angle and axis can be obtained for any rotation parameterization using
```julia
rotation_axis(r::Rotation)
rotation_angle(r::Rotation)
```

#### Initialization
All rotation types support `one(R)` to construct the identity rotation for the desired parameterization. A random rotation, uniformly sampled over the space of rotations, can be sampled using `rand(R)`. For example:
```julia
r = one(UnitQuaternion)  # equivalent to Quat(1.0, 0.0, 0.0, 0.0)
q = rand(UnitQuaternion)
p = rand(MRP{Float32})
```

#### Conversion
All rotatations can be converted to another parameterization by simply calling the constructor for the desired parameterization. For example:
```julia
q = rand(UnitQuaternion)
aa = AngleAxis(q)
r = RotMatrix(aa)
```

### Example Usage

```julia
using Rotations, StaticArrays

# create the null rotation (identity matrix)
id = one(RotMatrix{3, Float64})

# create a random rotation matrix (uniformly distributed over all 3D rotations)
r = rand(RotMatrix{3}) # uses Float64 by default

# create a point
p = SVector(1.0, 2.0, 3.0) # from StaticArrays.jl, but could use any AbstractVector...

# convert to a quaternion (Quat) and rotate the point
q = UnitQuaternion(r)
p_rotated = q * p

# Compose rotations
q2 = rand(UnitQuaternion)
q3 = q * q2

# Take the inverse (equivalent to transpose)
q_inv = transpose(q)
q_inv == inv(q)
p ≈ q_inv * (q * p)
q4 = q3 / q2  # q4 = q3 * inv(q2)
q5 = q3 \ q2  # q5 = inv(q3) * q2

# convert to a Modified Rodrigues Parameter (aka Stereographic quaternion projection, recommended for applications with differentiation)
spq = MRP(r)

# convert to the Angle-axis parameterization, or related Rotation vector
aa = AngleAxis(r)
rv = RotationVec(r)
ϕ = rotation_angle(r)
v = rotation_axis(r)

# convert to Euler angles, composed of X/Y/Z axis rotations (Z applied first)
# (all combinations of "RotABC" are defined)
r_xyz = RotXYZ(r)

# Rotation about the X axis by 0.1 radians
r_x = RotX(0.1)

# Composing axis rotations together automatically results in Euler parameterization
RotX(0.1) * RotY(0.2) * RotZ(0.3) === RotXYZ(0.1, 0.2, 0.3)

# Can calculate Jacobian - derivatives of rotations with respect to parameters
j1 = Rotations.jacobian(RotMatrix, q) # How does the matrix change w.r.t the 4 Quat parameters?
j2 = Rotations.jacobian(q, p) # How does the rotated point q*p change w.r.t. the 4 Quat parameters?
# ... all Jacobian's involving RotMatrix, MRP and Quat are implemented
# (MRP is ideal for optimization purposes - no constaints/singularities)
```

### Rotation Parameterizations

1. **Rotation Matrix** `RotMatrix{N, T}`

    An N x N rotation matrix storing the rotation.  This is a simple wrapper for
    a [StaticArrays](https://github.com/JuliaArrays/StaticArrays.jl) `SMatrix{N,N,T}`.
    A rotation matrix `R` should have the property `I = R * R'`, but this isn't
    enforced by the constructor. On the other hand, all the types below are
    guaranteed to be "proper" rotations for all input parameters (equivalently:
    parity conserving, in *SO(3)*, `det(r) = 1`, or a rotation without
    reflection).

2. **Arbitrary Axis Rotation** `AngleAxis{T}`

    A 3D rotation with fields `theta`, `axis_x`, `axis_y`, and
    `axis_z` to store the rotation angle and axis of the rotation.
    Like all other types in this package, once it is constructed it acts and
    behaves as a 3×3 `AbstractMatrix`. The axis will be automatically
    renormalized by the constructor to be a unit vector, so that `theta` always
    represents the rotation angle in radians.

3. **Quaternions** `UnitQuaternion{T}`

    A 3D rotation parameterized by a unit quaternion. Note that the constructor
    will renormalize the quaternion to be a unit quaternion, and that although
    they follow the same multiplicative *algebra* as quaternions, it is better
    to think of `UnitQuaternion` as a 3×3 matrix rather than as a quaternion *number*.

    Previously `Quat`.

4. **Rotation Vector** `RotationVec{T}`

    A 3D rotation encoded by an angle-axis representation as `angle * axis`.
    This type is used in packages such as [OpenCV](http://docs.opencv.org/2.4/modules/calib3d/doc/camera_calibration_and_3d_reconstruction.html#void%20Rodrigues%28InputArray%20src,%20OutputArray%20dst,%20OutputArray%20jacobian%29).

    Note: If you're differentiating a Rodrigues Vector check the result is what
    you expect at theta = 0.  The first derivative of the rotation *should*
    behave, but higher-order derivatives of it (as well as parameterization
    conversions) should be tested.  The Stereographic Quaternion Projection (`MRP`) is
    the recommended three parameter format for differentiation.

    Previously `RodriguesVec`.

6. **Rodrigues Parameters** `RodriguesParam{T}`
    A 3-parameter representation of 3D rotations that has a singularity at 180 degrees. They can be interpreted as a projection of the unit quaternion onto the plane tangent to the quaternion identity. They are computationally efficient and do not have a sign ambiguity.

7. **Modified Rodrigues Parameter** `MRP{T}`
    A 3D rotation encoded by the stereographic projection of a unit quaternion.  This projection can be visualized as a pin hole camera, with the pin hole matching the quaternion `[-1,0,0,0]` and the image plane containing the origin and having normal direction `[1,0,0,0]`.  The "null rotation" `Quaternion(1.0,0,0,0)` then maps to the `MRP(0,0,0)`

    These are similar to the Rodrigues vector in that the axis direction is stored in an unnormalized form, and the rotation angle is encoded in the length of the axis.  This type has the nice property that the derivatives of the rotation matrix w.r.t. the `MRP` parameters are rational functions, making the `MRP` type a good choice for differentiation / optimization.

    They are frequently used in aerospace applications since they are a 3-parameter representation whose singularity happens at 360 degrees. In practice, the singularity can be avoided with some switching logic between one of two equivalent MRPs (obtained by projecting the negated quaternion).

    Previously `SPQuat`.

8. **Cardinal axis rotations** `RotX{T}`, `RotY{T}`, `RotZ{T}`

    Sparse representations of 3D rotations about the X, Y, or Z axis, respectively.

9. **Two-axis rotations** `RotXY{T}`, etc

    Conceptually, these are compositions of two of the cardinal axis rotations above,
    so that `RotXY(x, y) == RotX(x) * RotY(y)` (note that the order of application to
    a vector is right-to-left, as-in matrix-matrix-vector multiplication: `RotXY(x, y) * v == RotX(x) * (RotY(y) * v)`).

10. **Euler Angles - Three-axis rotations** `RotXYZ{T}`, `RotXYX{T}`, etc

    A composition of 3 cardinal axis rotations is typically known as a Euler
    angle parameterization of a 3D rotation. The rotations with 3 unique axes,
    such as `RotXYZ`, are said to follow the [**Tait Bryan**](https://en.wikipedia.org/wiki/Euler_angles#Tait.E2.80.93Bryan_angles) angle ordering,
    while those which repeat (e.g. `EulerXYX`) are said to use [**Proper Euler**](https://en.wikipedia.org/wiki/Euler_angles#Conventions) angle ordering.

    Like the two-angle versions, the order of application to a vector is right-to-left, so that `RotXYZ(x, y, z) * v == RotX(x) * (RotY(y) * (RotZ(z) * v))`.  This may be interpreted as an "extrinsic" rotation about the Z, Y, and X axes or as an "intrinsic" rotation about the X, Y, and Z axes.  Similarly, `RotZYX(z, y, x)` may be interpreted as an "extrinsic" rotation about the X, Y, and Z axes or an "intrinsic" rotation about the Z, Y, and X axes.

### The Rotation Error state and Linearization
It is often convenient to consider perturbations or errors about a particular 3D rotation, such as applications in state estimation or optimization-based control. Intuitively, we expect these errors to live in three-dimensional space. For globally non-singular parameterizations such as unit quaternions, we need a way to map between the four parameters of the quaternion to this three-dimensional plane tangent to the four-dimensional hypersphere on which quaternions live.

There are several of these maps provided by Rotations.jl:
* `ExponentialMap`: A very common mapping that uses the quaternion
exponential and the quaternion logarithm. The quaternion logarithm
converts a 3D rotation vector (i.e. axis-angle vector) to a unit quaternion.
It tends to be the most computationally expensive mapping.

* `CayleyMap`: Represents the differential quaternion using Rodrigues
parameters. This parameterization goes singular at 180° but does not
inherit the sign ambiguities of the unit quaternion. It offers an
excellent combination of cheap computation and good behavior.

* `MRPMap`: Uses Modified Rodrigues Parameters (MRPs) to represent the
differential unit quaternion. This mapping goes singular at 360°.

* `QuatVecMap`: Uses the vector part of the unit quaternion as the
differential unit quaternion. This mapping also goes singular at 180° but is
the computationally cheapest map and often performs well.

Rotations.jl provides the `RotationError` type for representing rotation errors, that act just like a `SVector{3}` but carry the nonlinear map used to compute the error, which can also be used to convert the error back to a `UnitQuaternion` (and, by extention, any other 3D rotation parameterization). The following methods are useful for computing `RotationError`s and adding them back to the nominal rotation:
```julia
rotation_error(R1::Rotation, R2::Rotation, error_map::ErrorMap)  # compute the error between `R1` and `R2` using `error_map`
add_error(R::Roation, err::RotationError)  # "adds" the error to `R` by converting back a `UnitQuaterion` and composing with `R`
```
or their aliases
```julia
R1 ⊖ R2   # caclulate the error using the default error map
R1 ⊕ err  # alias for `add_error(R1, err)`
```

For a practical application of these ideas, see the quatenrion-multiplicative Extended Kalman Filter (MEKF). [This article](https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20040037784.pdf) provides a good description.

When taking derivatives with respect to quaternions we need to account both for these mappings and the fact that local perturbations to a rotation act through composition instead of addition, as they do in vector space (e.g. `q * dq` vs `x + dx`). The following methods are useful for computing these Jacobians for `UnitQuaternion`, `RodriguesParam` or `MRP`
* `∇rotate(q,r)`: Jacobian of the `q*r` with respect to the rotation
* `∇composition1(q2,q1)`: Jacobian of `q2*q1` with respect to q1
* `∇composition2(q2,q1,b)`: Jacobian of `q2*q1` with respect to q2
* `∇²composition1(q2,q1)`: Jacobian of `∇composition1(q2,q2)'b` where b is an arbitrary vector
* `∇differential(q)`: Jacobian of composing the rotation with an infinitesimal rotation, with
respect to the infinitesimal rotation. For unit quaternions, this is a 4x3 matrix.
* `∇²differential(q,b)`: Jacobian of `∇differential(q)'b` for some vector b.


### Import / Export

All parameterizations can be converted to and from (mutable or immutable)
3×3 matrices, e.g.

```julia
using StaticArrays, Rotations

# export
q = UnitQuaternion(1.0,0,0,0)
matrix_mutable = Array(q)
matrix_immutable = SMatrix{3,3}(q)

# import
q2 = Quaternion(matrix_mutable)
q3 = Quaternion(matrix_immutable)
```

### Notes

This package assumes [active (right handed) rotations](https://en.wikipedia.org/wiki/Active_and_passive_transformation) where applicable.


### Why use immutables / StaticArrays?

They're faster (Julia's `Array` and BLAS aren't great for 3x3 matrices) and
don't need preallocating or garbage collection. For example, see this benchmark
case where we get a 20× speedup:

```julia
julia> cd(Pkg.dir("Rotations") * "/test")

julia> include("benchmark.jl")

julia > BenchMarkRotations.benchmark_mutable()
Rotating using mutables (Base.Matrix and Base.Vector):
  0.124035 seconds (2 allocations: 224 bytes)
Rotating using immutables (Rotations.RotMatrix and StaticArrays.SVector):
  0.006006 seconds
```

## Acknowledgements

[![FugroRoames](https://avatars.githubusercontent.com/FugroRoames?s=150)](https://github.com/FugroRoames)
