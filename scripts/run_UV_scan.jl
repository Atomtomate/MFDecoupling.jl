using Distributed

#addprocs(96, topology=:master_worker, restrict=true, exeflags=["-J/scratch/projects/hhp00048/MFDecoupling/UScan/MFDecouplingSysimage_03.so", "--check-bounds=no"])
addprocs(90, topology=:master_worker, restrict=true)#, exeflags=["-J/scratch/projects/hhp00048/MFDecoupling/UScan/MFDecoupling.so", "--check-bounds=no"])
println("addprocs done")
flush(stdout)

@everywhere using Pkg
@everywhere Pkg.activate(joinpath(@__DIR__,".."))
@everywhere using MFDecoupling
@everywhere using JLD2

fp = ARGS[1]
fpout = ARGS[2]
use_real = parse(Bool, ARGS[3])

fp1 = joinpath(fp,"CK_U0.1V0.5.dat") #ARGS[1]
fp2 = joinpath(fp,"rho_U0.1V0.5.dat")

#TODO: read from file
const LL::Int = 1000
const Vin::Float64 = 0.5
const tmin::Float64 = 0.0
const tmax::Float64 = 500.0
tspan = (tmin,tmax)
tsave = LinRange(0.0,500.0,2000)



UList = union(LinRange(0.0,12.0,30), [0.1])#union(LinRange(0.0,0.2,10),LinRange(3.0,5.0,10),LinRange(10.0,20.0,10))
VList = union(LinRange(0.0,2.0,10), [0.5])#union(LinRange(0.0,2.0,25),[0.4,0.5,0.6])
X0 = read_inputs(fp1, fp2, LL)

@everywhere function solve_time_evolution(U::Float64, V::Float64, X0_in::Vector, L::Int, tspan, tsave, index, fpout, use_real)
    Q,P,LC,LK,LIm = MFDecoupling.gen_helpers(L)
    #X0, rhsf! = setup_calculation(X0, L; mode=:real)
    X0, rhsf! = setup_calculation(X0_in, L; mode= (use_real ? :real : :complex))
    p_0  = [L, U, V, 0.0, 0.0];

    linsolve = MFDecoupling.KrylovJL_GMRES()
    alg = MFDecoupling.AutoTsit5(MFDecoupling.KenCarp47(linsolve = linsolve, autodiff=use_real))
    prob = MFDecoupling.ODEProblem(rhsf!,X0,tspan,p_0)
    #idxs_list = union(collect(1:11),LIm .+ collect(1:11)) 
    idxs_list = collect(1:11)
    println("START $U/$V")
    flush(stdout)
    @time sol = MFDecoupling.solve(prob, alg; save_idxs=idxs_list, saveat=tsave, abstol=1e-9, reltol=1e-9);
    println("DONE with $U/$V")
    flush(stdout)
    jldopen(joinpath(fpout,"grid_res_$index.jld2"), "w") do f
        f["U"] = U
        f["V"] = V
        f["sol"] = sol[:,:]
        f["t"] = sol.t
    end
    return nothing #U, sol
end

wp = default_worker_pool()

futures = []
for (i,Ui) in enumerate(UList)
    for (j,Vi) in enumerate(VList)
        # println(j+(i-1)*length(VList))
        push!(futures, remotecall(solve_time_evolution, wp, Ui, Vi, X0, LL, tspan, tsave, j+(i-1)*length(VList), fpout, use_real))
    end
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
