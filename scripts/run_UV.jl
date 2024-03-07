using Distributed, ClusterManagers
using Hwloc

NNodes = 4
slurm_manager = SlurmManager(NNodes*96, ExponentialBackOff(n=20, first_delay=4,max_delay=10*4096, factor=2))
addprocs(slurm_manager, N = NNodes, topology=:all_to_all, lazy=true, exeflags=["--check-bounds=no"], job_file_loc = "./job_out")
println("addprocs done")
flush(stdout)

@everywhere using MFDecoupling
@everywhere using JLD2

fp = ARGS[1]
fpout = ARGS[2]
use_real = parse(Bool, ARGS[3])
ibatch = parse(Int,ARGS[4])
nbatch = parse(Int,ARGS[5])

fp1 = joinpath(fp,"CK_U0.1V1.0.dat") #ARGS[1]
fp2 = joinpath(fp,"rho_U0.1V1.0.dat")

#TODO: read from file
const LL::Int = 1000
const Uin::Float64 = 0.1
const Vin::Float64 = 1.0
const tmin::Float64 = 0.0
const tmax::Float64 = 500.0
tspan = (tmin,tmax)
tsave = LinRange(0.0,500.0,2000)

UList_1 = union(LinRange(0.0,5.0,60), [Uin])
VList_1 = LinRange(0.0,1.0,60)

# coarse
UList_2 = union(LinRange(0.0,12.0,50), [Uin])#union(LinRange(0.0,0.2,10),LinRange(3.0,5.0,10),LinRange(10.0,20.0,10))
VList_2 = union(LinRange(0.0,1.0,50), [Vin])#union(LinRange(0.0,2.0,25),[0.4,0.5,0.6])

# Upt
UList_3 = LinRange(3.0,5.0,40)
VList_3 = union(LinRange(0.0,1.0,30), [Vin])


UList = union(UList_1, UList_2, UList_3)
VList = union(VList_1, VList_2, VList_3)

X0 = read_inputs(fp1, fp2, LL)

@everywhere function solve_time_evolution(U::Float64, V::Float64, X0_in::Vector, L::Int, tspan, tsave, index, fpout, use_real)
    fname = joinpath(fpout,"grid_res_$index.jld2")
    if !isfile(fname)
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
        jldopen(fname, "w") do f
            f["U"] = U
            f["V"] = V
            f["sol"] = sol[:,:]
            f["t"] = sol.t
        end
    else
        println("File $index exists, skipping!")
    end
    return nothing #U, sol
end

wp = default_worker_pool()
UV_list = Base.product(UList,VList)
len_batches = trunc(Int,length(UV_list)/nbatch)
indices = ((ibatch*len_batches) + 1):((ibatch+1)*len_batches)
println("Running batch $ibatch/$nbatch with length $len_batches")
futures = []
for (i,el) in enumerate(UV_list)
    if i in indices
        U,V = el
    # println(j+(i-1)*length(VList))
    # println(i, " / ", U, " / ", V)
        push!(futures, remotecall(solve_time_evolution, wp, U, V, X0, LL, tspan, tsave, i, fpout, use_real))
    end
end

for fi in futures
    wait(fi)
end
