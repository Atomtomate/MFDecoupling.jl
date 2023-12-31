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
tsave = LinRange(300,500,2000)



UList = LinRange(3.3,4.3,5)



X0 = read_inputs(fp1, fp2, LL)

function solve_time_evolution(U::Float64, V::Float64, X0::Vector, L::Int, tspan, tsave, index, fpout)

    Q,P,LC,LK,LIm = MFDecoupling.gen_helpers(L)
    X0, rhsf! = setup_calculation(X0, L; mode=:real)
    p_0  = [L, U, V, 0.0, 0.0];

    linsolve = MFDecoupling.KrylovJL_GMRES()
    alg = MFDecoupling.AutoTsit5(MFDecoupling.KenCarp47(linsolve = linsolve))
    prob = MFDecoupling.ODEProblem(rhsf!,X0,tspan,p_0,
        progress = false,
        progress_steps = 0)
    idxs_list = union(collect(1:11),LIm .+ collect(1:11)) 
    @time sol = MFDecoupling.solve(prob, alg; save_idxs=idxs_list, saveat=tsave, abstol=1e-10, reltol=1e-10);
    println("DONE with $U")
    flush(stdout)
    # jldopen(joinpath(fpout,"res_$index.jld2"), "w") do f
    #     f["res_$index/U"] = U
    #     f["res_$index/sol"] = sol
    # end
    return nothing #U, sol
end
solve_time_evolution(UList[1], VV, X0, LL, tspan, tsave, 1, fpout)
