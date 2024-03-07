using Statistics, FourierAnalysis, Images

# Builds DataFrame from file obtained by merge.jl (after Uscan). 
# Expects .jld2 file (output of merge.jl) with data and path to file where dataframe will be stored.
# Example usage: julia build_df out_full.jld2 out_df.jld2

fn = ARGS[1]
dfFile = ARGS[2]

"""
     find_amplitude(tx, signal)

Estimation of amplitude, using local maxima/minima and linear interpolation.
Returns two arrays: `tx` (times of amplitudes) and `ampl` (amplitude).
"""
function find_amplitude(tx, signal)
    maxy_ix = findlocalmaxima(signal)
    miny_ix = findlocalminima(signal)
    if length(maxy_ix) != 0
        start = tx[miny_ix[1]] <  tx[maxy_ix[1]] ? 0 : 1
        ind = 2:length(maxy_ix)-2
        ampl = Vector{Float64}(undef, length(ind))
        t_ampl = tx[miny_ix[ind]]
        for (i,ii) in enumerate(ind)
            ampl[i] = mean(signal[maxy_ix[ii-1:ii]]) - signal[miny_ix[ii-start]]
        end
        return t_ampl, ampl
    else
        return [0.0], [0.0]
    end
end

function build_df(fn::String)
    df = DataFrame(U = Float64[], V = Float64[], t = Vector{Float64}[],        rho11 = Vector{ComplexF64}[], rho12 = Vector{ComplexF64}[],
                   rho13 = Vector{ComplexF64}[], rho14 = Vector{ComplexF64}[], rho22 = Vector{ComplexF64}[], rho23 = Vector{ComplexF64}[],
                   rho24 = Vector{ComplexF64}[], rho33= Vector{ComplexF64}[],  rho34 = Vector{ComplexF64}[], rho44 = Vector{ComplexF64}[], C01 = Vector{ComplexF64}[])
    if isfile(fn)
        rU, rV, rT, rS = jldopen(fn,"r") do ft
            ft["U"], ft["V"], ft["t"], ft["solution"]
        end
        for i in 1:length(rU)
            row = [rU[i], rV[i], rT, eachslice(rS[:,:,i],dims=1)...]
            push!(df,row)
        end
    elseif isdir(fn)
        fl = readdir(fn, join=true)
        filter!(x-> endswith(x,"jld2"), fl)
        for fi in fl
            rU, rV, rT, rS = jldopen(fi,"r") do fii
                tmp = keys(fii)[1]
                #V = haskey(fii["$tmp"], "V") ? fii["$tmp/V"] : NaN
                if haskey(fii, "U")
                    fii["U"], fii["V"], fii["t"], fii["sol"]
                else
                    fii["$tmp/U"], NaN, fii["$tmp/t"], fii["$tmp/sol"]
                end
            end
            row = [rU, rV, rT, eachslice(rS,dims=1)...]
            push!(df,row)
        end
    else
        error("$fn is not a result file or directory!")
    end
    df[!,:f_up] = df[:,:rho13] .+ df[:,:rho24]
    df[!,:f_do] = df[:,:rho12] .+ df[:,:rho34]
    
    # supplement mean value
    lim_ind = 300
    
    for cn in names(df)[4:end]
        df[!,Symbol(cn * "_lim")] = map(x-> mean(x[end-lim_ind:end]), df[:,Symbol(cn)])
    end

    ampl_cols = ["rho11", "f_up", "C01"]
    for cn in ampl_cols
        df[!,Symbol(cn * "_ampl")] = map(x->find_amplitude(x.t, real(x[cn])), eachrow(df[:,[:t, Symbol(cn)]]))
    end
    
    # calculate ffts of time series
    fft_cols = ["rho11", "rho12", "rho22", "rho24", "f_up", "C01"]
    fft_start_i = trunc(Int, length(df[1,:t]) / 10)
    tx_fft = df[1,:t][fft_start_i:end]
    freqs = collect(FourierAnalysis.fftfreq(length(tx_fft), 1.0/diff(tx_fft)[1]) |> FourierAnalysis.fftshift)
    df[!,:freqs] = [freqs for i in 1:size(df,1)]
    for cl in fft_cols
        df[!,Symbol(cl * "_fft")] = map(x-> FourierAnalysis.fft(x[1][fft_start_i:end] .- x[2]) ./ length(x[1]) |> FourierAnalysis.fftshift, eachrow(df[:,[Symbol(cl),Symbol(cl * "_lim")]]))
        df[!,Symbol(cl * "_fft_zm")] = map(x-> FourierAnalysis.fft(x[fft_start_i:end]) ./ length(x[fft_start_i:end]) |> FourierAnalysis.fftshift, df[:,Symbol(cl)])
    end

    return df
end

df = if false && isfile(dfFile)
    jldopen(dfFile) do ff
        ff["df"]
    end
else
    df_pt = build_df(fn_pt);
    df_coarse = build_df(fn_coarse);
    df_smallU = build_df(fn_smallU);
    
    # combine and remove duplicates
    df = combine(groupby(vcat(df_pt,df_coarse, df_smallU), [:U, :V]), last); 
    jldopen(dfFile, "w") do ff
        ff["df"] = df
    end
end
