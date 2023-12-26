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
    lim_ind = 1000
    var_window = trunc(Int, length(df[1,:t]) / 10)
    for cn in names(df)[4:end]
        df[!,Symbol(cn * "_lim")] = map(x-> mean(x[end-lim_ind:end]), df[:,Symbol(cn)])
        df[!,Symbol(cn * "_std")] = map(x-> [std(x[i-var_window:i]) for i in var_window+1:length(x)], df[:,Symbol(cn)])
    end

    # calculate ffts of time series
    fft_cols = ["rho11", "f_up", "C01"]
    for cl in fft_cols
        df[!,Symbol(cl * "_fft")] = map(x-> FourierAnalysis.fft(real(x[1] .- x[2])) ./ length(x[1]) |> FourierAnalysis.fftshift, eachrow(df[:,[Symbol(cl),Symbol(cl * "_lim")]]))
    end

    return df
end
df_pt= build_df(fn_pt);
df_coarse = build_df(fn_coarse);
df_smallU = build_df(fn_smallU);

# combine and remove duplicates
df = combine(groupby(vcat(df_pt,df_coarse, df_smallU), [:U, :V]), last); 

sort!(df_Vin, [:U])

VList = sort(unique(df.V))
