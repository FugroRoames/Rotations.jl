#=
The following tests check the distribution of `rand(RotType)`.
To visualize the distribution, try the following script.

    using Rotations, Plots
    RotType = RotMatrix3{Float64}

    # Sampling
    n = 10000
    rs = rand(RotType,n)
    q = rand(RotType)
    angles = rotation_angle.([r/q for r in rs])/π

    # Cumulative distribution function
    cdf_sampled(t) = count(<(t), angles)/n
    cdf_true(t) = t-sinpi(t)/π

    # Probability density function and histogram of sampling distribution
    histogram(angles, normed = true, label=false)
    plot!(t->1-cospi(t), 0, 1, label=false)

    # Cumulative distribution function
    plot(cdf_sampled, 0, 1, label=false)
    plot!(cdf_true, 0, 1, label=false)
=#

@testset "Distribution" begin

    # TODO: add RodriguesParam (#154)
    # TODO: consider to remove Euler rotations (#155)
    Type_SO3 =  (RotMatrix{3}, AngleAxis, RotationVec,
                UnitQuaternion, MRP,
                RotXYZ, RotYZX, RotZXY, RotXZY, RotYXZ, RotZYX,
                RotXYX, RotYZY, RotZXZ, RotXZX, RotYXY, RotZYZ)

    # TODO: add Angle2d and RotMatrix2
    Type_SO2 = (RotX, RotY, RotZ)

    # Number of sampling
    N = 100000

    @testset "Distribution $(RotType)" for RotType in Type_SO3
        # Sampling
        rs = rand(RotType, N)
        q = rand(RotType)
        angles = rotation_angle.([r*inv(q) for r in rs])/π

        # Check sampled rotations are rotations
        @test all(isrotation.(rs))
        @test isrotation(q)

        # Cumulative distribution function
        cdf_sampled(t) = count(<(t), angles)/N
        cdf_true(t) = t-sinpi(t)/π

        # Check CDF functions are approximately equal
        ts = 0:0.01:1
        @test norm(cdf_true.(ts) - cdf_sampled.(ts), Inf) < 0.01
    end

    @testset "distribution $(RotType)" for RotType in Type_SO2
        # Sampling
        rs = rand(RotType, N)
        q = rand(RotType)
        angles = rotation_angle.([r*inv(q) for r in rs])/π

        # Check sampled rotations are rotations
        @test all(isrotation.(rs))
        @test isrotation(q)

        # Cumulative distribution function
        cdf_sampled(t) = count(<(t), angles)/N
        cdf_true(t) = t

        # Check CDF functions are approximately equal
        ts = 0:0.01:1
        @test norm(cdf_true.(ts) - cdf_sampled.(ts), Inf) < 0.01
    end

end
