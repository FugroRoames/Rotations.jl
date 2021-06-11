using Rotations, StaticArrays, Test

@testset "2d Rotations" begin

    #################################
    # Traits, types, and construction
    #################################

    @testset "Core" begin
        r = one(RotMatrix{2,Float32})
        @test RotMatrix((1,0,0,1)) == RotMatrix(@SMatrix [1 0; 0 1])
        @test_throws DimensionMismatch RotMatrix((1,0,0,1,0))
    end

    ###############################
    # Check fixed relationships
    ###############################

    @testset "Identity rotation checks" begin
        I = one(SMatrix{2,2,Float64})
        I32 = one(SMatrix{2,2,Float32})
        for R in [RotMatrix{2}, Angle2d]
            @test @inferred(size(R)) == (2,2)
            @test @inferred(size(R{Float32})) == (2,2)
            @test one(R)::R == I
            @test one(R{Float32})::R{Float32} == I32
        end
    end

    ################################
    # check on the inverse function
    ################################

    @testset "Testing inverse()" begin
        repeats = 100
        for R in [RotMatrix{2,Float64}, Angle2d{Float64}]
            I = one(R)
            Random.seed!(0)
            for i = 1:repeats
                r = rand(R)
                @test isrotation(r)
                @test inv(r) == adjoint(r)
                @test inv(r) == transpose(r)
                @test inv(r)*r ≈ I
                @test r*inv(r) ≈ I
            end
        end
    end

    #########################################################################
    # Rotate some stuff
    #########################################################################

    # a random rotation of a random point
    @testset "Rotate Points" begin
        repeats = 100
        for R in [RotMatrix{2}, Angle2d]
            Random.seed!(0)
            for i = 1:repeats
                r = rand(R)
                m = SMatrix(r)
                v = randn(SVector{2})

                @test r*v ≈ m*v
            end

            # Test Base.Vector also
            r = rand(R)
            m = SMatrix(r)
            v = randn(2)

            @test r*v ≈ m*v
        end
    end

    # compose two random rotations
    @testset "Compose rotations" begin
        repeats = 100
        for R in [RotMatrix{2}, Angle2d]
            Random.seed!(0)
            for i = 1:repeats
                r1 = rand(R)
                m1 = SMatrix(r1)

                r2 = rand(R)
                m2 = SMatrix(r2)

                @test r1*r2 ≈ m1*m2
                θ1, θ2 = atan(r1[2,1],r1[1,1]), atan(r2[2,1],r2[1,1])
                @test r1*r2 ≈ RotMatrix(θ1+θ2)
                @test r1/r2 ≈ RotMatrix(θ1-θ2)
                @test r1\r2 ≈ RotMatrix(θ2-θ1)
            end
        end
    end


    #########################################################################
    # Test conversions between rotation types
    #########################################################################
    @testset "Convert rotations" begin
        repeats = 100
        for R in [RotMatrix{2}, Angle2d]
            Random.seed!(0)
            for _ in 1:repeats
                r1 = rand(R)
                @test R(r1) == r1
            end
        end
    end

    @testset "Types and products" begin
        for (R,T) in ((RotMatrix(pi/4), Float64),
                      (RotMatrix(Float32(pi/4)), Float32),
                      (RotMatrix{2,Float32}(pi/4), Float32))
            @test eltype(R) == T
            @test size(R) == (2,2)
            @test R == T[cos(pi/4) -sin(pi/4); sin(pi/4) cos(pi/4)]
            @test R * R ≈ T[0 -1; 1 0]
        end
    end

    @testset "RotMatrix{2} vs Angle2d" begin
        repeats = 100
        Random.seed!(0)
        for _ in 1:repeats
            theta = randn()
            r1 = RotMatrix{2}(theta)
            r2 = Angle2d(theta)
            v = randn(2)
            @test r1 ≈ r2
            @test r1 * v ≈ r2 * v
            @test RotMatrix{2}(r2) ≈ r1
            @test Angle2d(r1) ≈ r2
        end
    end

    @testset "Angle2d" begin
        repeats = 100
        Random.seed!(0)
        for _ in 1:repeats
            r = rand(Angle2d)
            @test r^(-1) ≈ inv(r)
            @test r^(-1.0) ≈ inv(r)
            @test r^1 ≈ r
            @test r^1.0 ≈ r
            @test r^2 ≈ r*r
            @test r^2.0 ≈ r*r

            x = randn()
            y = randn()
            @test r^(x+y) ≈ r^x * r^y
        end
    end

    @testset "angle" begin
        repeats = 100
        Random.seed!(0)
        for _ in 1:repeats
            theta = 2pi*rand() - pi
            a2d = Angle2d(theta)
            r2 = RotMatrix{2}(theta)
            @test rotation_angle(r2) ≈ theta
            @test rotation_angle(a2d) == theta
            @test rotation_angle(a2d) == Rotations.params(a2d)[1]
        end
    end
end

nothing
