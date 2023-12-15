using Distributed
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

linsolve = MFDecoupling.KrylovJL_GMRES()
alg_impl1 = MFDecoupling.AutoTsit5(MFDecoupling.KenCarp47(linsolve = linsolve))


UList = LinRange(3.3,4.3,5)


@everywhere function solve_time_evolution(U::Float64)
    LC = (LL+3)*LL/2
    LK = (LL-1)*LL/2
    LIm= 10 + LC +LK
    X0, rhsf! = setup_calculation(fp1, fp2, LL; mode=:real)
    p_0  = [LL, U, VV, 0.0, 0.0];
    prob2 = MFDecoupling.ODEProblem(rhsf!,X0,tspan,p_0,
        progress = false,
        progress_steps = 0)
    idxs_list = union(collect(1:11),Lim .+ collect(1:11)) 
    @time sol = MFDecoupling.solve(prob2, alg_impl1; save_idxs=idxs_list, save_everystep = true, abstol=1e-12, reltol=1e-12);
    return sol
end

futures = []
for (i,wi) in enumerate(workers())
    push!(futures, remotecall(, wi, fullParamList[batch_indices[i]]))
end

jldopen(fout) do f
    for (i,wi) in enumerate(workers())
        res = fetch(futures[i])
        Ui = UList[i]
        f["res_$i/U"] = Ui
        f["res_$i/sol"] = res
    end
end
