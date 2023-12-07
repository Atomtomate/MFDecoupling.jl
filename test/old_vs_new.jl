@testset "test rewrite" begin
    NSamples = 5
    ptest = LinRange(1,3,5)
    t = 0.1
    LL = 10
    ptest[1] = 10

    for i in 1:NSamples
        X0test = randn(185)
        dX0test_1 = similar(X0test)
        dX0test_2 = similar(X0test)

        MFDecoupling.rhs!(dX0test_1, X0test, ptest, t)
        MFDecoupling.test!(dX0test_2, X0test, ptest, t)
        @test all(dX0test_1 .≈ dX0test_2)
    end
end
