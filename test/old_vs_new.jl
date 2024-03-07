include("rhs_old.jl")
@testset "test rewrite" begin
    NSamples = 2
    ptest = collect(LinRange(1,3,5))
    t = 0.1
    for LL in [10,20,30]
        ptest[1] = LL


        Q,P,LC,LK,LIm = MFDecoupling.gen_helpers(LL)
        for i in 1:NSamples
            X0test = randn(ComplexF64, (LC+LK+10))
            X0test[1] = real(X0test[1])
            X0test[5] = real(X0test[5])
            X0test[8] = real(X0test[8])
            X0test[10] = real(X0test[10])

            for i in 1:LL
                Qii=floor(Int,(LL+1)*i - (i+1)*i/2 + i + 10)
                X0test[Qii] = real(X0test[Qii])
            end
            @test !any(isnan.(X0test))

            dX0test_1 = similar(X0test)
            dX0test_2 = similar(X0test)
            fill!(dX0test_1, NaN)
            fill!(dX0test_2, NaN)

            MFDecoupling.rhs!(dX0test_1, X0test, ptest, t, LC, LK, LIm, Q, P)
            test!(dX0test_2, X0test, ptest, t)

            @test !any(isnan.(dX0test_1))
            @test !any(isnan.(dX0test_2))
            @test all(dX0test_1 .≈ dX0test_2)
        end
    end
end





@testset "test real" begin
    NSamples = 2
    ptest = collect(LinRange(1,3,5))
    t = 0.1
    for LL in [10,20,30]
        ptest[1] = LL
        Q,P,LC,LK,LIm = MFDecoupling.gen_helpers(LL)


        for i in 1:NSamples
            X0test = randn(Float64, 2*(LC+LK+10))


            dX0test_1 = similar(X0test)
            dX0test_2 = similar(X0test)
            fill!(dX0test_1, NaN)
            fill!(dX0test_2, NaN)

            rhs_real_test1!(dX0test_1, X0test, ptest, t)
            MFDecoupling.rhs_real!(dX0test_2, X0test, ptest, t, LC, LK, LIm, Q, P)


           @test !any(isnan.(dX0test_1))
           @test !any(isnan.(dX0test_2))
           @test all(dX0test_1 .≈ dX0test_2)
        end
    end
end
