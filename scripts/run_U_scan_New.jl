using Distributed

#addprocs(96, topology=:master_worker, restrict=true, exeflags=["-J/scratch/projects/hhp00048/MFDecoupling/UScan/MFDecouplingSysimage_03.so", "--check-bounds=no"])
addprocs(90, topology=:master_worker, restrict=true, exeflags=["-J/scratch/projects/hhp00048/MFDecoupling/UScan/MFDecoupling.so", "--check-bounds=no"])
println("addprocs done")
flush(stdout)

@everywhere using Pkg
@everywhere Pkg.activate(joinpath(@__DIR__,".."))
@everywhere using MFDecoupling
@everywhere using JLD2

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
const tmax::Float64 = 400.0
tspan = (tmin,tmax)
tsave = LinRange(100.0,400.0,1000)



UList = LinRange(0.0,4.0,90)
X0 = read_inputs(fp1, fp2, LL)

@everywhere function solve_time_evolution(U::Float64, V::Float64, X0::Vector, L::Int, tspan, tsave, index, fpout)
    Q,P,LC,LK,LIm = MFDecoupling.gen_helpers(L)
    X0, rhsf! = setup_calculation(X0, L; mode=:real)
    p_0  = [L, U, V, 0.0, 0.0];

    linsolve = MFDecoupling.KrylovJL_GMRES()
    alg = MFDecoupling.AutoTsit5(MFDecoupling.ImplicitEuler(autodiff=false, linsolve = KrylovJL_GMRES()))
    prob = MFDecoupling.ODEProblem(rhsf!,X0,tspan,p_0,
        progress = false,
        progress_steps = 0)
    idxs_list = union(collect(1:11),LIm .+ collect(1:11)) 
    @time sol = MFDecoupling.solve(prob, alg; save_idxs=idxs_list, saveat=tsave, abstol=1e-6, reltol=1e-6);
    println("DONE with $U")
    flush(stdout)
    jldopen(joinpath(fpout,"res_$index.jld2"), "w") do f
        f["res_$index/U"] = U
        f["res_$index/sol"] = sol[:,:]
        f["res_$index/t"] = sol.t
    end
    return nothing #U, sol
end

wp = default_worker_pool()

futures = []
for (i,Ui) in enumerate(UList)
    push!(futures, remotecall(solve_time_evolution, wp, Ui, VV, X0, LL, tspan, tsave, i, fpout))
end

for fi in futures
    wait(fi)
end

# jldopen(fpout, "w") do f
#     for (i,fi) in enumerate(futures)
#         Ui, res = fetch(fi)
#         f["res_$i/U"] = Ui
#         f["res_$i/sol"] = res
#     end
# end
