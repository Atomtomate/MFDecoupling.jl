@testset "test rewrite" begin
    NSamples = 2
    ptest = collect(LinRange(1,3,5))
    t = 0.1
    LL = 10
    ptest[1] = 10

    for i in 1:NSamples
        X0test = randn(ComplexF64, 120)
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
