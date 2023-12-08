@testset "test rewrite" begin
    NSamples = 2
    ptest = collect(LinRange(1,3,5))
    t = 0.1
    LL = 10
    ptest[1] = 10


    for i in 1:NSamples
        X0test = randn(ComplexF64, 120)
        X0test[1] = real(X0test[1]) .+ 1im * 0
        X0test[5] = real(X0test[5]) .+ 1im * 0
        X0test[8] = real(X0test[8]) .+ 1im * 0
        X0test[10] = real(X0test[10]) .+ 1im * 0

        for ii in 1:LL
            Qii=floor(Int,(10+1)*i - (i+1)*i/2 + i +10)
            X0test[Qii] = real(X0test[Qii]) .+ 1im * 0
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
