using Pkg
Pkg.activate(joinpath(@__DIR__,".."))
using  MFDecoupling
using JLD2

fp1 = abspath("/home/localadmin/Nextcloud/Fuer_Julian/Initial_For_Time_Evolution/U0.1V0.5/CK_U0.1V0.5.dat") #ARGS[1]
fp2 = abspath("/home/localadmin/Nextcloud/Fuer_Julian/Initial_For_Time_Evolution/U0.1V0.5/rho_U0.1V0.5.dat")

outputf = abspath("/home/localadmin/Nextcloud/TeX/SIAM_MF_decoupling/Results/Time_evolution_Julia/")


LL = 1000
Uin = 0.1
Vin = 0.5
UU = 4.0
VV = 0.5
tmin = 0.0
tmax = 500.0


tspan = (tmin,tmax)

X0_real, rhsf! = setup_calculation(fp1, fp2, LL)
X0, rhsf_c! = setup_calculation(fp1, fp2, LL; mode=:complex)
# solve

linsolve = MFDecoupling.KrylovJL_GMRES()
    #SSPSDIRK2(autodiff=false)
#TODO: autodiff
#TODO: test various preconditioners (e.g. algebraicmultigrid), see advanced_ode/#stiff
#TODO: test CVODE_BDF, KenCarp47
#TODO: stiff-switch threshold
#TODO: steady state
#TODO: test if Jocobian is sparse
#TODO: save_everystep=false

#   maxstiffstep = 10, maxnonstiffstep = 3,
#    nonstifftol::T = 9 // 10, stifftol::T = 9 // 10,
#    dtfac = 2.0, stiffalgfirst = false
use_real = true
alg_impl1 = MFDecoupling.AutoTsit5(MFDecoupling.KenCarp47(autodiff=use_real, linsolve = linsolve); maxnonstiffstep=3, stiffalgfirst=true)
alg_impl2 = MFDecoupling.AutoTsit5(MFDecoupling.ImplicitEuler(autodiff=use_real, linsolve = linsolve))
alg_impl3 = MFDecoupling.KenCarp47(autodiff=use_real, linsolve = linsolve);
#TODO: concrete, but sparse jacobian : alg_impl4 = MFDecoupling.RadauIIA3(autodiff=false, linsolve = linsolve)
alg_impl4 = MFDecoupling.TRBDF2(autodiff=use_real, linsolve = linsolve);

    #KenCarp47(linsolve = KrylovJL_GMRES(), autodiff=false)
    #

p_0  = [LL, UU, VV, 0.0, 0.0];

prob = if use_real
    MFDecoupling.ODEProblem(rhsf!,X0_real,tspan,p_0,
        progress = false,
        progress_steps = 0)
else
    MFDecoupling.ODEProblem(rhsf_c!,X0,tspan,p_0,
        progress = false,
        progress_steps = 0)
end

#idxs_list = union(1:10, [10 + floor(Int,(LL+1)*i - (i+1)*i/2 + i+1) for i in 0:LL-1])
idxs_list = collect(1:11)
tlist = LinRange(tmin,tmax,10000)
#ODER: save_everystep=true
#@time sol1_1 = MFDecoupling.solve(prob1, alg_impl1; save_everystep = true, ab1e-7);
# @time sol2 = MFDecoupling.solve(prob1, alg_impl2; save_everystep = true, abstol=1e-8, reltol=1e-8);

#println("Errors:\n",sol1_1.errors)
#println(" =========================== ")
#println("Algorithm Details:\n",sol1_1.alg)
#println(" =========================== ")
#println("Solution stats:\n", sol1_1.stats)
#tol_list = [1e-5,1e-6,1e-7,1e-8,1e-9,1e-10,1e-11]
tol_list = [1e-10]
for tol in tol_list
    #@time sol2_2 = MFDecoupling.solve(prob, alg_impl4; abstol=tol, reltol=tol, save_idxs=idxs_list, save_everystep=true, dt=1e-4) #, saveat=tlist);
    @time sol1_2 = MFDecoupling.solve(prob, alg_impl5; abstol=tol, reltol=tol, save_idxs=idxs_list, save_everystep=true) #, saveat=tlist);
    println("Errors:\n",sol1_2.errors)
    println(" =========================== ")
    println("Algorithm Details:\n",sol1_2.alg)
    println(" =========================== ")
    println("Solution stats:\n", sol1_2.stats)
    #jldopen(outputf * "Q_Uin($Uin)Vin($Vin)_Ufi($UU)Vfi($VV)_tol$tol.jld2","w") do f
    jldopen(outputf * "TEST_NEW.jld2","w") do f
      f["solution"] = sol1_2.u
      f["t"] = sol1_2.t
      f["Uin"] = Uin
      f["Vin"] = Vin
      f["UU"] = UU
      f["VV"]=VV
      f["L"]=LL
    end
end
