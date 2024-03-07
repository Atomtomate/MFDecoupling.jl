using Pkg
Pkg.activate(joinpath(@__DIR__,".."))
using  MFDecoupling

fp1 = abspath("/home/julisn/Codes/MFDecoupling_TestData/t1/CK_U0.1V0.5.dat") #ARGS[1]
fp2 = abspath("/home/julisn/Codes/MFDecoupling_TestData/t1/rho_U0.1V0.5.dat")

LL = 1000
Uin = 0.1
Vin = 0.5
UU = 4.0
VV = 0.5
tmin = 0.0
tmax = 3.0


tspan = (tmin,tmax)

X0, rhsf! = setup_calculation(fp1, fp2, LL; mode=:real)
X0_c, rhsf_c! = setup_calculation(fp1, fp2, LL; mode=:complex)
# solve

linsolve = MFDecoupling.KrylovJL_GMRES()
    #SSPSDIRK2(autodiff=false)
#TODO: test various preconditioners (e.g. algebraicmultigrid), see advanced_ode/#stiff
#TODO: test CVODE_BDF, KenCarp47
#TODO: stiff-switch threshold
#TODO: steady state
#TODO: test if Jocobian is sparse

alg_impl1 = MFDecoupling.AutoTsit5(MFDecoupling.KenCarp47(linsolve = linsolve))
alg_impl2 = MFDecoupling.AutoTsit5(MFDecoupling.ImplicitEuler(linsolve = linsolve))
alg_impl1_c = MFDecoupling.AutoTsit5(MFDecoupling.KenCarp47(linsolve = linsolve, autodiff=false))
alg_impl2_c = MFDecoupling.AutoTsit5(MFDecoupling.ImplicitEuler(linsolve = linsolve, autodiff=false))
    #KenCarp47(linsolve = KrylovJL_GMRES(), autodiff=false)
    #

p_0  = [LL, UU, VV, 0.0, 0.0];

prob_c = MFDecoupling.ODEProblem(rhsf_c!,X0_c,tspan,p_0,
    progress = false,
    progress_steps = 0)
prob = MFDecoupling.ODEProblem(rhsf!,X0,tspan,p_0,
    progress = false,
    progress_steps = 0)
p_0_2  = [LL, UU - 1.0, VV, 0.0, 0.0];
prob2 = MFDecoupling.ODEProblem(rhsf!,X0,tspan,p_0_2,
    progress = false,
    progress_steps = 0)

@time sol1 = MFDecoupling.solve(prob, alg_impl1; save_everystep = false, abstol=1e-6, reltol=1e-6);
@time sol2 = MFDecoupling.solve(prob2, alg_impl2; save_everystep = false, abstol=1e-6, reltol=1e-6);
@time sol1_c = MFDecoupling.solve(prob_c, alg_impl1_c; save_everystep = false, abstol=1e-6, reltol=1e-6);
@time sol2_c = MFDecoupling.solve(prob_c, alg_impl2_c; save_everystep = false, abstol=1e-6, reltol=1e-6);

# @time sol2 = MFDecoupling.solve(prob1, alg_impl2; save_everystep = true, abstol=1e-8, reltol=1e-8);

println("Errors:\n",sol1.errors)
println(" =========================== ")
println("Algorithm Details:\n",sol1.alg)
println(" =========================== ")
println("Solution stats:\n", sol1.stats)


println("Errors:\n",sol2.errors)
println(" =========================== ")
println("Algorithm Details:\n",sol2.alg)
println(" =========================== ")
println("Solution stats:\n", sol2.stats)


println("Errors:\n",sol1_c.errors)
println(" =========================== ")
println("Algorithm Details:\n",sol1_c.alg)
println(" =========================== ")
println("Solution stats:\n", sol1_c.stats)


println("Errors:\n",sol2_c.errors)
println(" =========================== ")
println("Algorithm Details:\n",sol2_c.alg)
println(" =========================== ")
println("Solution stats:\n", sol2_c.stats)


