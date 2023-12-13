@testset "test rewrite" begin
    NSamples = 2
    ptest = collect(LinRange(1,3,5))
    t = 0.1
    LL = 10
    ptest[1] = 10


    for i in 1:NSamples
        X0test = randn(ComplexF64, 120)
        X0test[1] = real(X0test[1])
        X0test[5] = real(X0test[5])
        X0test[8] = real(X0test[8])
        X0test[10] = real(X0test[10])

        for ii in 1:LL
            Qii=floor(Int,(10+1)*i - (i+1)*i/2 + i +10)
            X0test[Qii] = real(X0test[Qii])
        end

        dX0test_1 = similar(X0test)
        dX0test_2 = similar(X0test)
        fill!(dX0test_1, NaN)
        fill!(dX0test_2, NaN)

        MFDecoupling.rhs!(dX0test_1, X0test, ptest, t)
        MFDecoupling.test!(dX0test_2, X0test, ptest, t)

        @test !any(isnan.(dX0test_1))
        @test !any(isnan.(dX0test_2))
        @test all(dX0test_1 .≈ dX0test_2)
    end
end





@testset "test real" begin
    NSamples = 2
    ptest = collect(LinRange(1,3,5))
    t = 0.1
    LL = 10
    ptest[1] = 10


    for i in 1:NSamples
        X0test = randn(Float64, 240)


        dX0test_1 = similar(X0test)
        dX0test_2 = similar(X0test)
        fill!(dX0test_1, NaN)
        fill!(dX0test_2, NaN)

        MFDecoupling.rhs_real_test1!(dX0test_1, X0test, ptest, t)
        MFDecoupling.rhs_real!(dX0test_2, X0test, ptest, t)


        """
        println("rho")
        for index in 1:10
            println(index,"   ",dX0test_1[index] - dX0test_2[index],"  ",dX0test_1[index+120] - dX0test_2[index+120])
        end
        println("C")
        for index in 11:75
            println(index,"   ",dX0test_1[index] - dX0test_2[index],"  ",dX0test_1[index+120] - dX0test_2[index+120])
        end

        println("K")
        for index in 76:120
            println(index,"   ",dX0test_1[index] - dX0test_2[index],"  ",dX0test_1[index+120] - dX0test_2[index+120])
        end
        """


       @test !any(isnan.(dX0test_1))
       @test !any(isnan.(dX0test_2))
       @test all(dX0test_1 .≈ dX0test_2)
    end
end
