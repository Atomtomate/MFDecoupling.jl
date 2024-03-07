# ================================================================================ #
#                                  Run TE                                          #
# Description: Calculates equilibrium and time evolution for MFDecoupling.jl.      #
#              chain length `L` and initial and final values for `U` and `V` are   #
#              hardcoded in the script and need to be changed here.                #
# Arguments: - file path to data                                                   #
#            - rhs as real or complex (default should be 0,                        #
#                                       unless you need a specific solver)         #
#            - ibatch: index of batch. should be in [0, nbatch)                    #
#            - nbatch: number of batches, used to segment large lists of quenches  #
# Output: 2 Files: - eq_res.jld2 contains equilibirum results.                     #
#                  - te_res.jld2 contains time evolution results.                  #
# Author: Julian Stobbe                                                            #
# ================================================================================ #
using Distributed, ClusterManagers
using Hwloc

NNodes = 1
slurm_manager = SlurmManager(NNodes*96, ExponentialBackOff(n=20, first_delay=4,max_delay=10*4096, factor=2))
addprocs(slurm_manager, N = NNodes, lazy=true, exeflags=["--check-bounds=no"], job_file_loc = "./job_out")
println("addprocs done")
flush(stdout)

@everywhere using MFDecoupling
@everywhere using JLD2
@everywhere using Random
@everywhere using FileWatching
@everywhere lock_file = "lock.pid"

fp = ARGS[1]
use_real = parse(Bool, ARGS[2])
ibatch = parse(Int,ARGS[3])
nbatch = parse(Int,ARGS[4])


#TODO: read from file
const Vin::Float64 = 0.5
const tmin::Float64 = 0.0



LList = [1000] #[800, 1000, 1200] #[100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000]
const Uin_4::Float64 = 4.0
const Uin_5::Float64 = 3.5
const Uin_6::Float64 = 6.0
const Uin_7::Float64 = 1.0
const Vin::Float64 = 0.5

#60
UList_1 = union(LinRange(0.0,5.0,16), [Uin_4, Uin_5, Uin_6, Uin_7])
VList_1 = LinRange(0.0,1.0,16)

# coarse
#50
UList_2 = union(LinRange(0.0,12.0,15), [Uin_4, Uin_5, Uin_6, Uin_7])#union(LinRange(0.0,0.2,10),LinRange(3.0,5.0,10),LinRange(10.0,20.0,10))
VList_2 = union(LinRange(0.0,1.0,15), [Vin])#union(LinRange(0.0,2.0,25),[0.4,0.5,0.6])

# Upt
# #40/30
UList_3 = LinRange(3.0,5.0,14)
VList_3 = union(LinRange(0.0,1.0,13), [Vin])


UList = union(UList_1, UList_2, UList_3)
VList = union(VList_1, VList_2, VList_3)

@everywhere function solve_time_evolution(fp::String, Uin::Float64, Ufi::Float64, Vin::Float64, Vfi::Float64, L::Int, tspan, tsave, use_real)
    diff, ind = find_te_index(fp, Uin, Ufi, Vin, Vfi, L)
    if diff > 1e-8
        println("Starting time evolution for Uin = $Uin, Vin = $Vin, Ufi = $Ufi, Vfi = $Vfi, L = $L")
        X0_in, index = solve_eq(fp, Uin, Vin, L)
        pid_lock = joinpath(fp, "TE.pid")
        X0, rhsf! = setup_calculation(X0_in, L; mode= (use_real ? :real : :complex))
        p_0  = [L, Ufi, Vfi, 0.0, 0.0];
        linsolve = MFDecoupling.KrylovJL_GMRES()
        alg = MFDecoupling.AutoTsit5(MFDecoupling.KenCarp47(linsolve = linsolve, autodiff=use_real))
        prob = MFDecoupling.ODEProblem(rhsf!,X0,tspan,p_0)
        #idxs_list = union(collect(1:11),LIm .+ collect(1:11)) 
        idxs_list = collect(1:11)
        println("START $Ufi/$Vfi")
        flush(stdout)
        sol = MFDecoupling.solve(prob, alg; save_idxs=idxs_list, saveat=tsave, abstol=1e-9, reltol=1e-9);
        println("DONE with $Ufi/$Vfi")
        flush(stdout)
        fpout = joinpath(fp, "te_res.jld2")
        lock = FileWatching.Pidfile.mkpidlock(pid_lock)
        jldopen(fpout, "a+") do f
            content = newfile ? [] : (haskey(f,"content") ? f["content"] : [])
            push!(content, (Uin,Ufi,Vin,Vfi,L))
            new_index = length(content) + 1
            Base.delete!(f, "content")
            f["content"] = content
            f["$index/Uin"] = Uin
            f["$index/Vin"] = Vin
            f["$index/Ufi"] = Ufi
            f["$index/Vfi"] = Vfi
            f["$index/L"] = L
            f["$index/sol"] = sol[idxs_list,:]
            f["$index/t"] = sol.t
        end
        close(lock)
        println("Done with time evolution for Uin = $Uin, Vin = $Vin, Ufi = $Ufi, Vfi = $Vfi, L = $L")
    else
        println("Skipping time evolution for Uin = $Uin, Vin = $Vin, Ufi = $Ufi, Vfi = $Vfi, L = $L ! Found existing result at index $ind")
    end
    return nothing #U, sol
end

@everywhere function find_te_index(fp::String, Uin::Float64, Ufi::Float64, Vin::Float64, Vfi::Float64, L::Int)
    fp_eq = joinpath(fp, "te_res.jld2")
    pid_lock = joinpath(fp, "TE.pid")
    diff, ind = if isfile(fp_eq)
        lock = FileWatching.Pidfile.mkpidlock(pid_lock)
        jldopen(fp_eq, "a+") do f
            content = f["content"]
            findmin(x->abs(x[1]-Uin)+abs(x[2]-Ufi)+abs(x[3]-Vin)+abs(x[4]-Vfi)+abs(x[5]-L)/(100*L), content) 
        end
        close(lock)
    else
        Inf, NaN
    end
    return diff, ind
end
@everywhere function find_fg_init(fp::String, U::Float64, V::Float64, L::Int)
    fp_eq = joinpath(fp, "eq_res.jld2")
    pid_lock = joinpath(fp, "eq.pid")
    test = isfile(fp_eq)
    println("DBG3: ", fp_eq)
    println("DBG2: isfile $U/$V/$L = $test")
    diff, f_init, g_init, ind = if isfile(fp_eq)
        lock = FileWatching.Pidfile.mkpidlock(pid_lock)
        jldopen(fp_eq, "a+") do f
            content = f["content"]
            res, ind = findmin(x->abs(x[1]-U)+abs(x[2]-V)+abs(x[3]-L)/(100*L), content) 
            f_init = f["$ind/f_eq"]
            g_init = f["$ind/g_eq"]
            res, f_init, g_init, ind
        end
        close(lock)
    else
        Inf, 1.0, -1.0, NaN
    end
    return diff, f_init, g_init, ind
end

@everywhere function solve_eq(fp::String, U::Float64, V::Float64, L::Int)
    re = r"f(\h*)\=(\h*)(?<f>[+-]?(\d+([.]\d*)?([eE][+-]?\d+)?|[.]\d+([eE][+-]?\d+)?))(.*)g(\h*)=(\h*)(?<g>[+-]?(\d+([.]\d*)?([eE][+-]?\d+)?|[.]\d+([eE][+-]?\d+)?))"
    println("Starting equilibrium calculation for U = $U, V = $V, L = $L")
    fp_eq = joinpath(fp, "eq_res.jld2")
    pid_lock = joinpath(fp, "eq.pid")
    diff, f_init, g_init, ind = find_fg_init(fp, U, V, L) 
    println("DBG: U/V/L $U/$V/$L || diff = $diff")
    flush(stdout)
    rng = MersenneTwister();

    X, index = if diff < 1e-8
        lock = FileWatching.Pidfile.mkpidlock(pid_lock)
        jldopen(fp_eq, "r") do f
            f["$ind/X0"], ind
        end
        close(lock)
        println("Found equilibrium solution for U/V/L = $U/$V/$L")
        flush(stdout)
    else
        fileID = 100000*myid() + convert(Int, rand(rng, UInt16))
        while isfile(joinpath(fp,"rho0_$fileID.dat"))
            fileID = 100000*myid() + convert(Int, rand(rng, UInt16))
        end
        println("Starting equilibrium calculation for U/V/L = $U/$V/$L with id $fileID")
        flush(stdout)

        io = IOBuffer();
        cmd_p = `./main $U $V $g_init $f_init $fileID $L` 
        cmd = pipeline(cmd_p; stdout=io, stderr=devnull);
        run(cmd)
        str = String(take!(io))
        m = match(re, split(strip(str), "\n")[end-4])
        if isnothing(m)
            println("ERROR DURING EXECUTON OF U = $U, V = $V, L + $L")
        end
        fp1 = joinpath(fp,"CK_$fileID.dat") #ARGS[1]
        fp2 = joinpath(fp,"rho0_$fileID.dat")
        X0 = read_inputs(fp1, fp2, L)
        newfile = !isfile(fp_eq)

        lock = FileWatching.Pidfile.mkpidlock(pid_lock)
        new_index = jldopen(fp_eq, "a+") do f
            content = newfile ? [] : (haskey(f,"content") ? f["content"] : [])
            push!(content, (U,V,L))
            new_index = length(content) + 1
            println("adding EQ res at $new_index")
            Base.delete!(f, "content")
            f["content"] = content
            f["$new_index/X0"] = X0
            f["$new_index/f_eq"] = m[:f]
            f["$new_index/g_eq"] = m[:g]
            new_index
        end
        close(lock)
        rm(fp1) 
        rm(fp2)
        X0, new_index
    end
    println("Done with equilibrium calculation for U = $U, V = $V, L = $L")
    return X0, index
end

wp = default_worker_pool()
UV_List = Base.product(UList,VList)
len_batches = trunc(Int,length(UV_List)/nbatch)
indices = ((ibatch*len_batches) + 1):((ibatch+1)*len_batches)
ibatch == nbatch-1 && (indices = sort(union(indices, first(indices):length(UV_List))))
println("Running batch $ibatch/$nbatch with length $len_batches")
futures = []

for (iL, LL) in enumerate(LList)
    tmax::Float64 = LL/2.0
    tspan = (tmin,tmax)
    tsave = LinRange(0.0,tmax,trunc(Int,2*LL))
    if iL in indices
        for (i,el) in enumerate(UV_List)
            U,V = el
            push!(futures, remotecall(solve_time_evolution, wp, fp, U, U, V, V, LL, tspan, tsave, use_real))
            #solve_time_evolution(fp, U, V, LL, tspan, tsave, use_real)
        end
    end
end

for fi in futures
    wait(fi)
end
