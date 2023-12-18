using Pkg
Pkg.activate(joinpath(@__DIR__,".."))
using MFDecoupling
using JLD2


fp = ARGS[1]
fpout = ARGS[2]
fp1 = joinpath(fp,"CK_U0.1V0.5.dat") #ARGS[1]
fp2 = joinpath(fp,"rho_U0.1V0.5.dat")

#TODO: read from file
const LL::Int = 1000
const Uin::Float64 = 0.1
const Vin::Float64 = 0.5
const VV::Float64 = 0.5
const tmin::Float64 = 0.0
const tmax::Float64 = 1.0
tspan = (tmin,tmax)



UList = LinRange(3.3,4.3,5)



X0, Q, P, LC, LK = read_inputs(fp1, fp2, LL)

function solve_time_evolution(U::Float64, X0, Q, P, LC, LK, L, V, tspan)
    LIm= trunc(Int, 10 + LC +LK)
    X0, rhsf! = setup_calculation(X0, Q, P, LC, LK, L; mode=:real)
    p_0  = [L, U, V, 0.0, 0.0];
    linsolve = MFDecoupling.KrylovJL_GMRES()
    alg= MFDecoupling.AutoTsit5(MFDecoupling.KenCarp47(linsolve = linsolve))
    prob = MFDecoupling.ODEProblem(rhsf!,X0,tspan,p_0,
        progress = false,
        progress_steps = 0)
    tlist = LinRange(350.0,500.0,5000)
    idxs_list = union(collect(1:11), LIm.+ collect(1:11)) 
    @time sol = MFDecoupling.solve(prob, alg; save_idxs=idxs_list, saveat=tlist,
                                    save_everystep = false, abstol=1e-2, reltol=1e-2);
    return U, sol
end


results = []
for Ui in UList
    push!(results, solve_time_evolution(Ui, X0, Q, P, LC, LK, LL, VV, tspan))
end

jldopen(fpout, "w") do f
    for (i,res_i) in enumerate(results)
        Ui, res = res_i
        f["res_$i/U"] = Ui
        f["res_$i/sol"] = res
    end
end

