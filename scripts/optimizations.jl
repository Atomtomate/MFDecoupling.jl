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
tmax = 10.0


tspan = (tmin,tmax)

X0, Q, P, S, LC, LK = read_inputs(fp1, fp2, LL)
# solve

linsolve = MFDecoupling.KrylovJL_GMRES()
    #SSPSDIRK2(autodiff=false)
#TODO: autodiff
#TODO: test various preconditioners (e.g. algebraicmultigrid), see advanced_ode/#stiff
#TODO: test CVODE_BDF, KenCarp47
#TODO: stiff-switch threshold
#TODO: steady state
#TODO: test if Jocobian is sparse
#TDO: save_everystep=false

alg_impl1 = MFDecoupling.AutoTsit5(MFDecoupling.KenCarp47(autodiff=false, linsolve = linsolve))
alg_impl2 = MFDecoupling.AutoTsit5(MFDecoupling.ImplicitEuler(autodiff=false, linsolve = linsolve))
    #KenCarp47(linsolve = KrylovJL_GMRES(), autodiff=false)
    #

p_0  = [LL, UU, VV, 0.0, 0.0];

prob1 = MFDecoupling.ODEProblem(MFDecoupling.test!,X0,tspan,p_0,
    progress = false,
    progress_steps = 0)
prob2 = MFDecoupling.ODEProblem(MFDecoupling.test2!,X0,tspan,p_0,
    progress = false,
    progress_steps = 0)

@time sol1_2 = MFDecoupling.solve(prob2, alg_impl1; save_everystep = false, abstol=1e-7, reltol=1e-7);
@time sol1_1 = MFDecoupling.solve(prob1, alg_impl1; save_everystep = false, abstol=1e-7, reltol=1e-7);

# @time sol2 = MFDecoupling.solve(prob1, alg_impl2; save_everystep = true, abstol=1e-8, reltol=1e-8);

println("Errors:\n",sol1_1.errors)
println(" =========================== ")
println("Algorithm Details:\n",sol1_1.alg)
println(" =========================== ")
println("Solution stats:\n", sol1_1.stats)


println("Errors:\n",sol1_2.errors)
println(" =========================== ")
println("Algorithm Details:\n",sol1_2.alg)
println(" =========================== ")
println("Solution stats:\n", sol1_2.stats)

