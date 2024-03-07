using Pkg
Pkg.activate(joinpath(@__DIR__,".."))
using MFDecoupling

fp1 = ARGS[1]
fp2 = ARGS[2]

LL = parse(Int, ARGS[3])
Uin = parse(Float64, ARGS[4]);
Vin = parse(Float64, ARGS[5]);
UU = parse(Float64, ARGS[6])
VV = parse(Float64, ARGS[7]);
tmin = parse(Float64, ARGS[8]);
tmax = parse(Float64, ARGS[9]);
outputf = ARGS[10]


tspan = (tmin,tmax)

X0, rhsf! = setup_calculation(fp1, fp2, LL)
# solve

linsolve = KrylovJL_GMRES()
alg_expl = AutoTsit5(Rosenbrock23(autodiff=false))
    #SSPSDIRK2(autodiff=false)
alg_impl = AutoTsit5(ImplicitEuler(autodiff=false, linsolve = KrylovJL_GMRES()))
alg_impl_2 = MFDecoupling.AutoTsit5(MFDecoupling.KenCarp47(autodiff=use_real, linsolve = linsolve); maxnonstiffstep=3, stiffalgfirst=true)
    #KenCarp47(linsolve = KrylovJL_GMRES(), autodiff=false)
    #

p_0  = [LL, UU, VV, 0.0, 0.0];

prob = ODEProblem{true}(rhsf!,X0,tspan,p_0,
    progress = true,
    progress_steps = 1)

sol = solve(prob, alg_impl; save_everystep = true, abstol=1e-8, reltol=1e-8);

println("Errors:\n",sol.errors)
println(" =========================== ")
println("Algorithm Details:\n",sol.alg)
println(" =========================== ")
println("Solution stats:\n", sol.stats)

jldopen(outputf,"w") do f
  f["solution"] = sol.u
  f["t"] = sol.t
  f["Uin"] = Uin
  f["Vin"] = Vin
  f["UU"] = UU
  f["VV"]=VV
  f["L"]=LL
end

