using MFDecoupling
using FourierAnalysis, Statistics
using Images
using JLD2

fn = ARGS[1]
fn_out = ARGS[2]

function find_lim(tx::AbstractVector, signal::AbstractVector; tx_start=tx[end-100], tx_end=490.0)
    start_ind = searchsortedfirst(tx, tx_start)    
    stop_ind = searchsortedfirst(tx, tx_end)    
    return find_lim(tx, signal, start_ind, stop_ind) 
end

function find_lim(tx::AbstractVector, signal::AbstractVector, start_ind::Int, stop_ind::Int)
    return mean(signal[start_ind:stop_ind])
end

function find_amplitude_lim(tx::AbstractVector, signal::AbstractVector, start_ind::Int, stop_ind::Int)
    ampl_lim = maximum(real(signal[start_ind:stop_ind])) - minimum(real(signal[start_ind:stop_ind]))
    return ampl_lim
end

function build_results(fn, fn_out)
    tx_start = 450.0
    tx_end   = 490.0
    UList = nothing
    VList = nothing
    tx = nothing
    f_up_lim = nothing
    rho11_lim = nothing
    f_ampl_lim = nothing
    rho11_ampl_lim = nothing

    jldopen(fn, "r") do f
        UList = f["U"]
        VList = f["V"]
        tx = f["t"]
        sol_raw = f["solution"]
        start_ind = searchsortedfirst(tx, tx_start)    
        stop_ind = searchsortedfirst(tx, tx_end)    
        println("interval: tx[$start_ind:$stop_ind]")

        f_up_lim = Vector{Float64}(undef, length(UList))
        rho11_lim = Vector{Float64}(undef, length(UList))
        f_ampl_lim = Vector{Float64}(undef, length(UList))
        rho11_ampl_lim = Vector{Float64}(undef, length(UList))

        for i in axes(sol_raw,3)
            f_up_lim[i] = real(find_lim(tx, sol_raw[3,:,i], start_ind, stop_ind) + find_lim(tx, sol_raw[7,:,i], start_ind, stop_ind))
            rho11_lim[i] = real(find_lim(tx, sol_raw[1,:,i], start_ind, stop_ind))
            f_ampl_lim[i] = find_amplitude_lim(tx, sol_raw[3,:,i] .+ sol_raw[7,:,i], start_ind, stop_ind)
            rho11_ampl_lim[i] = find_amplitude_lim(tx, sol_raw[1,:,i], start_ind, stop_ind)
        end
    end
    
    jldopen(fn_out, "w") do fout
        fout["UList"] = UList
        fout["VList"] = VList
        fout["tx"] = tx
        fout["f_up_lim"] = f_up_lim
        fout["rho11_lim"] = rho11_lim
        fout["f_ampl_lim"] = f_ampl_lim
        fout["rho11_ampl_lim"] = rho11_ampl_lim
    end
end

build_results(fn, fn_out)
