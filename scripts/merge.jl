using Distributed
addprocs(90)

@everywhere using JLD2
@everywhere using MFDecoupling


read_dir = "data"
fn_out = ARGS[1]

# fn_out = "out_combined.jld2"
fn_in_list = readdir(read_dir, join=true, sort=true)
filter!(x->endswith(x, ".jld2"), fn_in_list)
# println(fn_in_list)
N = 2000 # 20000#10000 

sol_res = Array{ComplexF64,3}(undef, 11, N, length(fn_in_list))
U_res = Array{Float64,1}(undef, length(fn_in_list))
V_res = Array{Float64,1}(undef, length(fn_in_list))
t_list = Array{Float64,1}(undef, N)

@everywhere function read_file(fn)
    println("$fn")
    sol, U, V = jldopen(fn,"r") do f_in
        f_in["sol"], f_in["U"], f_in["V"]
    end
    return sol, U, V
end

wp = default_worker_pool()
futures = []



for fn_in in fn_in_list
        push!(futures, remotecall(read_file, wp, fn_in))
end

for (i,fi) in enumerate(futures)
    if i == 1
        t_list[:] = jldopen(fn_in_list[1],"r") do f_in
            f_in["t"]
        end
    end
    sol_i, U_i, V_i = try 
        fetch(fi)
    catch err
        println("file $(fn_in_list[i]) failed!")
    end
    sol_res[:,:,i] = sol_i
    U_res[i] = U_i
    V_res[i] = V_i
end

jldopen(fn_out, "w") do f_out
    f_out["solution"] = sol_res 
    f_out["U"] = U_res
    f_out["V"] = V_res
    f_out["t"] = t_list
end
