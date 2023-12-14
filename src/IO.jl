# ==================================================================================================== #
#                                              IO.jl                                                   #
# ---------------------------------------------------------------------------------------------------- #
#   Author          : Irakli Titvinidze, Julian Stobbe                                                 #
# ----------------------------------------- Description ---------------------------------------------- #
#   Input/Output operations for the module.                                                            #
# -------------------------------------------- TODO -------------------------------------------------- #
#   Documentation                                                                                      #
#   Cleanup                                                                                      #
# ==================================================================================================== #

"""
    read_inputs(fp1::String, fp2::String, L::Int)

Reads two files at locations `fp1`, `fp2`.... ?
"""
function read_inputs(fp1::String, fp2::String, L::Int)
    println("Reading File 1")
    CK= open(fp1, "r") do f
        readdlm(f)
    end;

    println("Reading File 2")
    ρ = open(fp2, "r") do f
        arr = readdlm(f)
    end;

    C=zeros(Float64, L+1, L+1)

    # println("L = $L")
    # println(size(CK))
    C[2:L+1,2:L+1]=deepcopy(CK[1:L,1:L]);
    C[1,2:L+1]=deepcopy(CK[2*L+1,1:L]+CK[2*L+2,1:L])
    C[2:L+1,1]=deepcopy(CK[1:L,2*L+1]+CK[1:L,2*L+2]);

    K=CK[L+1:2*L,1:L];

    lC = size(C,1)
    lK = size(K,1)
    lρ = size(ρ,1)

    X0 = ComplexF64[]

    for i in 1:lρ
        for j in i:lρ
            push!(X0, ρ[i,j]) 
        end
    end


    for j in 2:lC
        push!(X0, C[1,j])
    end

    for i in 2:lC
        for j in i:lC
            push!(X0, C[i,j])
        end
    end


    for i in 1:lK
        for j in (i+1):lK
            push!(X0, K[i,j])
        end
    end
    Q = UpperTriangular(Int[(L+1)*i - (i+1)*i/2 + j for i in 0:L, j in 0:L])
    P = UpperTriangular(Int[(i-1)*L - (i+1)*i/2 + j for i in 1:L, j in 1:L])
    # S = UpperTriangular(Int[(i-1)*4 - (i-1)*i/2 + j for i in 1:4, j in 1:4]);

    LC::Int     = floor(Int,(L+3)*L/2)    # Number of C
    LK::Int     = floor(Int,(L-1)*L/2)    # Number of K``
    return X0, Q, P, LC, LK
end

"""
    setup(fp1::String, fp2::String, L::Int; mode::Symbol=:real)

Sets up the calculation. 
Returns start vector `X0` and function `rhs` for the right hand side of differential equation. 
Mode can either be :real or :complex. The former uses a real vector representation (with twice as many vector entries) and is the (much fast) default.
"""
function setup_calculation(fp1::String, fp2::String, L::Int; mode=:real)

    X0, Q, P, LC, LK = read_inputs(fp1::String, fp2::String, L::Int)
    LIm::Int    =10+LC+LK
    Q = UpperTriangular(Int[(L+1)*i - (i+1)*i/2 + j for i in 1:L, j in 1:L])
    P = UpperTriangular(Int[(i-1)*L - (i+1)*i/2 + j for i in 1:L, j in 1:L])
    if mode == :real 
        X0_real = vcat(real(X0), imag(X0))
        rhs_f_r(dX::Vector, X::Vector, p::Vector, t::Float64)::Nothing = rhs_real!(dX, X, p, t, LC, LK, LIm, Q, P)
        return X0_real, rhs_f_r
    elseif mode == :complex
        rhs_f_c(dX::Vector, X::Vector, p::Vector, t::Float64)::Nothing = rhs!(dX, X, p, t, LC, LK, LIm, Q, P)
        return X0, rhs_f_c
    else
        error("Unkown mode $mode")
    end
end
