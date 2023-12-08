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
    read_inputs(fp1::String, fp2::String, LL::Int)

Reads two files at locations `fp1`, `fp2`.... ?
"""
function read_inputs(fp1::String, fp2::String, LL::Int)
    println("Reading File 1")
    CK= open(fp1, "r") do f
        readdlm(f)
    end;

    println("Reading File 2")
    ρ = open(fp2, "r") do f
        arr = readdlm(f)
    end;

    C=zeros(Float64, LL+1, LL+1)

    # println("LL = $LL")
    # println(size(CK))
    C[2:LL+1,2:LL+1]=deepcopy(CK[1:LL,1:LL]);
    C[1,2:LL+1]=deepcopy(CK[2*LL+1,1:LL]+CK[2*LL+2,1:LL])
    C[2:LL+1,1]=deepcopy(CK[1:LL,2*LL+1]+CK[1:LL,2*LL+2]);

    K=CK[LL+1:2*LL,1:LL];

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
    Q = UpperTriangular(Int[(LL+1)*i - (i+1)*i/2 + j for i in 0:LL, j in 0:LL])
    P = UpperTriangular(Int[(i-1)*LL - (i+1)*i/2 + j for i in 1:LL, j in 1:LL])
    S = UpperTriangular(Int[(i-1)*4 - (i-1)*i/2 + j for i in 1:4, j in 1:4]);

    LC::Int     = floor(Int,(LL+3)*LL/2)    # Number of C
    LK::Int     = floor(Int,(LL-1)*LL/2)    # Number of K``
    return X0, Q, P, S, LC, LK
end
