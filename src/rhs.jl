# ==================================================================================================== #
#                                             rhs.jl                                                   #
# ---------------------------------------------------------------------------------------------------- #
#   Author          : Irakli Titvinidze, Julian Stobbe                                                 #
# ----------------------------------------- Description ---------------------------------------------- #
#   Right hand side for time evolution differential equation                                           #
# -------------------------------------------- TODO -------------------------------------------------- #
#   Documentation                                                                                      #
#   Define index-lookup (maybe struct struct? Assigned to Julian!)                                     #
#   refactor with real arrays                                                                          #
#   run benchmarks (assigned to Julian)                                                                #
# ==================================================================================================== #


"""
    rhs!(dX::Vector, X::Vector, p, t)::Nothing

RHS for DGL.

Parameters 
-------
  - *p[1]*:  L,    size of the chain
  - *p[2]*:  U,    local Hubbard interaction on the impurity site
  - *p[3]*:  Vr,   coupling between the impurity and chain
  - *p[4]*:  μ_imp chemical potential for impurity
  - *p[5]*:  - μ_c   chemical potential for chain (onsite energy shift is includid in the chemical potential)
"""
function rhs_real!(dX::Vector, X::Vector, p::Vector, t::Float64, LC::Int, LK::Int, LIm::Int, Q::AbstractMatrix, P::AbstractMatrix)::Nothing
    
    J::Float64     = 1.0
    L::Int         = Int(p[1])
    U::Float64     = p[2]
    ϵ_imp::Float64 = -U/2     # In this version it has to be so!!!!
    Vr::Float64    = p[3]
    μ_imp::Float64 = -p[4]    # In this version it has to be so!!!!
    μ_c::Float64   = -p[5]    # In this version it has to be so!!!!
                              # Later on we should have only one μ


    """
    LC = (L+3)*L/2
    LK = (L-1)*L/2
    LIm= 10 + LC +LK

    Q[i,j] = (L+1)*i - (i+1)*i/2 + j    i <= j   (i=j= 0 is not considered)
            1, ..., LC
    P[i,j] = (i-1)*L - (i+1)*i/2 + j    0 < i < j
            1, ..., LK


    Matrix Elelements
    ρ:  [1,1], [1,2], [1,3], [1,4], [2,2], [2,3], [2,4], [3,3], [3,4], [4,4]
        Re  ->  [1, ..., 10]    and     Im  ->  [11 + LC + LK, ..., 20 + LC + LK]

        Re(ρ[1,1]) -> X[1],             Re(ρ[1,2]) -> X[2],             Re(ρ[1,3]) -> X[3],             Re(ρ[1,4]) -> X[4],
        Re(ρ[2,2]) -> X[5],             Re(ρ[2,3]) -> X[6],             Re(ρ[2,4]) -> X[7],
        Re(ρ[3,3]) -> X[8],             Re(ρ[3,4]) -> X[9],
        Re(ρ[3,3]) -> X[10]
        Im(ρ[1,1]) -> X[11+LC+LK],      Im(ρ[1,2]) -> X[12+LC+LK],      Im(ρ[1,3]) -> X[13+LC+LK],      Im(ρ[1,4]) -> X[14+LC+ LK],
        Im(ρ[2,2]) -> X[15+LC+LK],      Im(ρ[2,3]) -> X[16+LC+LK],      Im(ρ[2,4]) -> X[17+LC+LK],
        Im(ρ[3,3]) -> X[18+LC+LK],      Im(ρ[3,4]) -> X[19+LC+LK],
        Im(ρ[4,4]) -> X[20+LC+LK]

    C:  [0,1], [0,2], ..., [0,L], [1,1], ..., [1,L], [2,2], ..., [2,L], [3,3], ...., [L-1,L-1], [L-1,L], [L,L]
        Re  ->  [11, ..., 10 + LC]      and     Im  ->  [21 + LC + LK, ..., 20 + 2 LC + LK]

        Re(C[i,j]) -> X[10+Q[i,j]]
                    Re(C[0,1]) -> X[11],    Re(C[0,2]) -> X[12],    ...,
                    Re(C[1,1]) -> X[11+L],  ...,    Re(C[L,L]) -> X[10+LC]
        Im(C[i,j]) -> X[10+Q[i,j] + LIm]
                    Im(C[0,1]) -> X[11+LIm],      Im(C[0,2]) -> X[12+LIm],     ...,
                    Im(C[1,1]) -> X[11+L+LIm],    ...,    Im(C[L,L]) -> X[10+LC+LIm]

    K:  [1,2], ..., [1,L], [2,3], ..., [2,L], [3,4], ...., [3,L], [4,5], ..., [L-2,L-1], [L-2,L], [L-1,L]
        Re  ->  [11 + LC, ..., LIm]     and     Im  ->  [21 + 2 LC + LK, ..., 20 + 2 LC + 2 LK]

        Re(K[i,j]) -> X[10+LC+P[i,j]]
                    Re(K[1,2]) -> X[11+LC],     Re(K[1,3]) -> X[12+LC],     ...,
                    Re(K[2,3]) -> X[10+LC+L),   ...,    Re(K[L-1,L]) ->  X[LIm]
        Im(K[i,j]) -> X[10+LC+P[i,j]+LIm]
                    Im(K[1,2]) -> X[11+LC+LIm],     Im(K[1,3]) -> X[12+LC+LIm],     ...,
                    Im(K[2,3]) -> X[10+LC+L+LIm),   ...,    Im(K[L-1,L]) ->  X[2 LIm]
    """






    temp = -2 * imag(Vr*(X[11] - 1im*X[11+LIm]) * (X[2] + X[3] +1im*(X[2+LIm] + X[3+LIm])))
    dX[1] = temp
    dX[1+LIm] = 0
    temp = (U- μ_imp + ϵ_imp)*(X[2] + 1im*X[2+LIm]) + Vr * (X[11] + 1im*X[11+LIm]) * (X[5] + X[6] - X[1] + 1im*(X[5+LIm] - X[6+LIm] - X[1+LIm])) - Vr * (X[11] - 1im*X[11+LIm]) * (X[4] + 1im * X[4+LIm])
    dX[2] = imag(temp)
    dX[2+LIm] = -1 * real(temp)
    temp = (U- μ_imp + ϵ_imp)*(X[3] + 1im*X[3+LIm])  + Vr * (X[11] + 1im*X[11+LIm]) * (X[6] + X[8] - X[1] + 1im*(X[6+LIm] + X[8+LIm] - X[1+LIm])) - Vr * (X[11] - 1im*X[11+LIm]) * (X[4] + 1im * X[4+LIm])
    dX[3] = imag(temp)
    dX[3+LIm] = -1 * real(temp)
    temp = (U- 2*(μ_imp - ϵ_imp))*(X[4] + 1im * X[4+LIm]) + Vr * (X[11] + 1im*X[11+LIm]) * (X[7] + X[9] - X[2] -X[3] + 1im*(X[7+LIm] + X[9+LIm] - X[2+LIm] -X[3+LIm]) )
    dX[4] = imag(temp)
    dX[4+LIm] = -1 * real(temp)
    temp = 2 * imag(Vr * (X[11] - 1im*X[11+LIm]) *(X[2]  - X[7] + 1im*(X[2+LIm]  - X[7+LIm])))
    dX[5] = temp
    dX[5+LIm] = 0
    temp = Vr * (X[11] - 1im*X[11+LIm]) * (X[3] -X[7] + 1im*(X[3+LIm] -X[7+LIm])) + Vr * (X[11] + 1im*X[11+LIm]) * (X[9] - X[2] -1im*(X[9+LIm] - X[2+LIm]))
    dX[6] = imag(temp)
    dX[6+LIm] = -1 * real(temp)
    temp = Vr * (X[11] - 1im*X[11+LIm]) * (X[4] + 1im*X[4+LIm]) - (μ_imp - ϵ_imp) * (X[7] + 1im*X[7+LIm]) + Vr * (X[11] + 1im*X[11+LIm]) * (X[10] -X[5] - X[6] + 1im*(X[10+LIm] -X[5+LIm] - X[6+LIm]))
    dX[7] = imag(temp)
    dX[7+LIm] = -1 * real(temp)
    temp = 2 * imag(Vr * (X[11] - 1im*X[11+LIm]) *(X[3]  - X[9] +1im*(X[3+LIm]  - X[9+LIm])))
    dX[8] = temp
    dX[8+LIm] = 0
    temp = Vr * (X[11] - 1im*X[11+LIm]) * (X[4] + 1im*X[4+LIm]) - (μ_imp - ϵ_imp) * (X[9] + 1im*X[9+LIm]) + Vr * (X[11] + 1im*X[11+LIm]) * (X[10] - X[6] - X[8] + 1im*(X[10+LIm] + X[6+LIm] - X[8+LIm]))
    dX[9] = imag(temp)
    dX[9+LIm] = -1 * real(temp)
    temp =2 * imag(Vr * (X[11] - 1im*X[11+LIm]) *(X[7]  + X[9] +1im*(X[7+LIm]  + X[9+LIm])))
    dX[10] = temp
    dX[10+LIm] = 0

    temp = J * (X[12] + 1im*X[12+LIm]) - μ_c * (X[11] +1im*X[11+LIm]) - 2 * Vr * (X[3] + X[7]  + 1im*(X[3+LIm] + X[7+LIm])) * (X[10+Q[1,1]] + 1im*X[10+Q[1,1]+LIm]) + Vr * (X[3] + X[7]  + 1im*(X[3+LIm] + X[7+LIm]))
    dX[11] = imag(temp)
    dX[11+LIm] = -1 * real(temp)
    for j in 2:L-1
        temp = J * (X[9+j] + X[11+j] + 1im*(X[9+j+LIm] + X[11+j+LIm])) - μ_c * (X[10+j] + 1im*X[10+j+LIm]) - 2 * Vr * (X[3] + X[7]  + 1im*(X[3+LIm] + X[7+LIm])) * (X[10+Q[1,j]] + 1im*X[10+Q[1,j]+LIm]) + 2 * Vr * (X[3] + X[7]  - 1im*(X[3+LIm] + X[7+LIm])) * (X[10+LC+P[1,j]] + 1im*X[10+LC+P[1,j]+LIm])
        dX[10+j] =  imag(temp)
        dX[10+j+LIm] = -1 * real(temp)
    end
    temp = J*(X[9+L] + 1im*X[9+L+LIm]) - μ_c * (X[10+L] + 1im*X[10+L+LIm]) - 2 * Vr * (X[3] + X[7]  + 1im*(X[3+LIm] + X[7+LIm])) * (X[10+Q[1,L]] + 1im*X[10+Q[1,L]+LIm]) + 2 * Vr * (X[3] + X[7]  - 1im*(X[3+LIm] + X[7+LIm])) * (X[10+LC+P[1,L]] + 1im*X[10+LC+P[1,L]+LIm])
    dX[10+L] = imag(temp)
    dX[10+L+LIm] = -1 * real(temp)
    temp = -2 * imag(-J*(X[10+Q[1,2]] +1im*X[10+Q[1,2]+LIm]) + Vr*(X[3] + X[7]  - 1im*(X[3+LIm] + X[7+LIm]))*(X[11] +1im*X[11+LIm]))
    dX[10+Q[1,1]] = temp
    dX[10+Q[1,1]+LIm] = 0
    for j in 2 : L-1
        temp = J * (X[10+Q[2,j]] - X[10+Q[1,j-1]] - X[10+Q[1,j+1]] + 1im*(X[10+Q[2,j]+LIm] - X[10+Q[1,j-1]+LIm] - X[10+Q[1,j+1]+LIm])) + Vr*(X[3] + X[7]  - 1im*(X[3+LIm] + X[7+LIm])) * (X[10+j] + 1im*X[10+j+LIm])
        dX[10+Q[1,j]] = -1 * imag(temp)
        dX[10+Q[1,j]+LIm] = real(temp)
    end
    temp = J * (X[10+Q[2,L]] - X[10+Q[1,L-1]] + 1im*(X[10+Q[2,L]+LIm] - X[10+Q[1,L-1]+LIm])) + Vr*(X[3] + X[7]  - 1im*(X[3+LIm] + X[7+LIm])) * (X[10+L] + 1im*X[10+L+LIm])
    dX[10+Q[1,L]] = -1 * imag(temp)
    dX[10+Q[1,L]+LIm] = real(temp)
    for i in 2:(L-1)
        temp = -2 * J * imag(X[10+Q[i-1,i]] - X[10+Q[i,i+1]] + 1im*(X[10+Q[i-1,i]+LIm] - X[10+Q[i,i+1]+LIm]))
        dX[10+Q[i,i]] = temp
        dX[10+Q[i,i]+LIm] = 0
        for j in i+1 : L-1
            temp = J * (X[10+Q[i-1,j]] + X[10+Q[i+1,j]] - X[10+Q[i,j-1]] - X[10+Q[i,j+1]] + 1im*(X[10+Q[i-1,j]+LIm] + X[10+Q[i+1,j]+LIm] - X[10+Q[i,j-1]+LIm] - X[10+Q[i,j+1]+LIm]))
            dX[10+Q[i,j]] = -1 * imag(temp)
            dX[10+Q[i,j]+LIm] = real(temp)
        end
        temp = J * (X[10+Q[i-1,L]] + X[10+Q[i+1,L]] - X[10+Q[i,L-1]] + 1im*(X[10+Q[i-1,L]+LIm] + X[10+Q[i+1,L]+LIm] - X[10+Q[i,L-1]+LIm]))
        dX[10+Q[i,L]] = -1 * imag(temp)
        dX[10+Q[i,L]+LIm] = real(temp)
    end
    temp = -2 * J * imag(X[10+Q[L-1,L]] + 1im*X[10+Q[L-1,L]+LIm])
    dX[10+Q[L,L]] = temp
    dX[10+Q[L,L]+LIm] = 0

    temp = J * (X[10+LC+P[1,3]] + 1im*X[10+LC+P[1,3]+LIm]) - 2* μ_c * (X[10+LC+P[1,2]] + 1im*X[10+LC+P[1,2]+LIm]) + Vr*(X[3] + X[7]  + 1im*(X[3+LIm] + X[7+LIm])) * (X[12] + 1im*X[12+LIm])
    dX[10+LC+P[1,2]] = imag(temp)
    dX[10+LC+P[1,2]+LIm] = -1 * real(temp)
    for j in 3:L-1
        temp = J * (X[10+LC+P[2,j]] + X[10+LC+P[1,j-1]] + X[10+LC+P[1,j+1]] + 1im*(X[10+LC+P[2,j]+LIm] + X[10+LC+P[1,j-1]+LIm] + X[10+LC+P[1,j+1]+LIm])) - 2* μ_c*(X[10+LC+P[1,j]] + 1im*X[10+LC+P[1,j]+LIm]) + Vr*(X[3] + X[7]  + 1im*(X[3+LIm] + X[7+LIm]))*(X[10+j] + 1im*X[10+j+LIm])
        dX[10+LC+P[1,j]] = imag(temp)
        dX[10+LC+P[1,j]+LIm] = -1 * real(temp)
    end
    temp = J * (X[10+LC+P[2,L]] + X[10+LC+P[1,L-1]] + 1im*(X[10+LC+P[2,L]+LIm] + X[10+LC+P[1,L-1]+LIm])) - 2* μ_c*(X[10+LC+P[1,L]] + 1im*X[10+LC+P[1,L]+LIm]) + Vr*(X[3] + X[7]  + 1im*(X[3+LIm] + X[7+LIm]))*(X[10+L] + 1im*X[10+L+LIm])
    dX[10+LC+P[1,L]] = imag(temp)
    dX[10+LC+P[1,L]+LIm] = -1 * real(temp)
    for i in 2 : L-2
        temp = J * (X[10+LC+P[i-1,i+1]] + X[10+LC+P[i,i+2]] +1im*(X[10+LC+P[i-1,i+1]+LIm] + X[10+LC+P[i,i+2]+LIm])) - 2 * μ_c * (X[10+LC+P[i,i+1]] + 1im * X[10+LC+P[i,i+1]+LIm])
        dX[10+LC+P[i,i+1]] = imag(temp)
        dX[10+LC+P[i,i+1]+LIm] = -1 * real(temp)
        for j in i+2 : L-1
            temp = J * (X[10+LC+P[i-1,j]] + X[10+LC+P[i+1,j]] + X[10+LC+P[i,j-1]] + J*X[10+LC+P[i,j+1]] + 1im*(X[10+LC+P[i-1,j]+LIm] + X[10+LC+P[i+1,j]+LIm] + X[10+LC+P[i,j-1]+LIm] + J*X[10+LC+P[i,j+1]+LIm])) - 2* μ_c*(X[10+LC+P[i,j]] + 1im*X[10+LC+P[i,j]+LIm])
            dX[10+LC+P[i,j]] = imag(temp)
            dX[10+LC+P[i,j]+LIm] = -1 * real(temp)
        end
        temp = J * (X[10+LC+P[i-1,L]] + X[10+LC+P[i+1,L]] + X[10+LC+P[i,L-1]] + 1im*(X[10+LC+P[i-1,L]+LIm] + X[10+LC+P[i+1,L]+LIm] + X[10+LC+P[i,L-1]+LIm])) - 2* μ_c*(X[10+LC+P[i,L]] + 1im*X[10+LC+P[i,L]+LIm])
        dX[10+LC+P[i,L]] = imag(temp)
        dX[10+LC+P[i,L]+LIm] = -1 * real(temp)
    end
    temp = J * (X[10+LC+P[L-2,L]] + 1im*X[10+LC+P[L-2,L]+LIm]) - 2* μ_c * (X[10+LC+P[L-1,L]] +1im*X[10+LC+P[L-1,L]+LIm])
    dX[10+LC+P[L-1,L]] = imag(temp)
    dX[10+LC+P[L-1,L]+LIm] = -1 * real(temp)


    return nothing
end


"""
    rhs_real_test1!(dXin::Vector, Xin::Vector, p::Vector, t::Float64)::Nothing

Naive test for real version of [`rhs!`](@ref rhs!), [`rhs_real!`](@ref rhs_real!).
"""
function rhs_real_test1!(dXin::Vector, Xin::Vector, p::Vector, t::Float64)::Nothing

    J::Float64     = 1.0
    L::Int         = Int(p[1])
    U::Float64     = p[2]
    ϵ_imp::Float64 = -U/2     # In this version it has to be so!!!!
    Vr::Float64    = p[3]
    μ_imp::Float64 = -p[4]    # In this version it has to be so!!!!
    μ_c::Float64   = -p[5]    # In this version it has to be so!!!!
                              # Later on we should have only one μ
    LC::Int     = floor(Int,(L+3)*L/2)
    LK::Int     = floor(Int,(L-1)*L/2)
    LRe::Int    =10+LC+LK

    ##################################################################################################
    ####################### This Part has to be moved outside! I can't manage! #######################
    ##################################################################################################
        # TODO: Auslagern
    Q = UpperTriangular(Int[(L+1)*i - (i+1)*i/2 + j for i in 1:L, j in 1:L])
    P = UpperTriangular(Int[(i-1)*L - (i+1)*i/2 + j for i in 1:L, j in 1:L])
    #S = UpperTriangular(Int[(i-1)*4 - (i-1)*i/2 + j for i in 1:4, j in 1:4]);  We do not neet it !
    ##################################################################################################
    ##################################################################################################

    NH = trunc(Int, length(Xin)/2)
    X = Xin[1:NH] .+ Xin[NH+1:end] .* 1im
    dX = similar(X)

    dX[1] = -2 * imag(conj(Vr * X[11]) *( X[2]  + X[3]))
    dX[2] = -1im*( (U- μ_imp + ϵ_imp)*X[2] + Vr * X[11] * (X[5] + conj(X[6]) - X[1]) -conj(Vr * X[11]) * X[4])
    dX[3] = -1im*( (U- μ_imp + ϵ_imp)*X[3] + Vr * X[11] * (X[6] + X[8] - X[1]) -conj(Vr * X[11]) * X[4])
    dX[4] = -1im*( (U- 2*(μ_imp - ϵ_imp))*X[4] + Vr * X[11] * (X[7] + X[9] - X[2] -X[3]))
    dX[5] = 2 * imag(conj(Vr * X[11]) *( X[2]  - X[7]))
    dX[6] = -1im*(conj(Vr * X[11]) * (X[3] -X[7]) + Vr * X[11] * conj(X[9] - X[2]) )
    dX[7] = -1im*(conj(Vr * X[11]) * X[4] - (μ_imp - ϵ_imp) * X[7] + Vr * X[11] * (X[10] -X[5] - X[6]))
    dX[8] = 2 * imag(conj(Vr * X[11]) *( X[3]  - X[9]))
    dX[9] = -1im*(conj(Vr * X[11]) * X[4] - (μ_imp - ϵ_imp) * X[9] + Vr * X[11] * (X[10] -conj(X[6]) - X[8]))
    dX[10] =2 * imag(conj(Vr * X[11]) *( X[7]  + X[9]))

    dX[11]=-1im*(J*X[12] - μ_c*X[11] - 2*Vr*(X[3]+X[7])*X[10+Q[1,1]] + Vr*(X[3]+X[7]))
    for j in 2:L-1
        dX[10+j]=-1im*(J*X[9+j] + J*X[11+j] - μ_c*X[10+j]
                 - 2*Vr*(X[3]+X[7])*X[10+Q[1,j]] + 2*Vr*conj(X[3]+X[7])*X[10+LC+P[1,j]])
    end
    dX[10+L]=-1im*(J*X[9+L] - μ_c*X[10+L]
            - 2*Vr*(X[3]+X[7])*X[10+Q[1,L]] + 2*Vr*conj(X[3]+X[7])*X[10+LC+P[1,L]])
    dX[10+Q[1,1]]= -2 * imag(-J*X[10+Q[1,2]] + Vr*conj(X[3]+X[7])*X[11])
    for j in 2 : L-1
        dX[10+Q[1,j]]=1im*(J*X[10+Q[2,j]]-J*X[10+Q[1,j-1]]-J*X[10+Q[1,j+1]]+ Vr*conj(X[3]+X[7])*X[10+j])
    end
    dX[10+Q[1,L]]=1im*(J*X[10+Q[2,L]] - J*X[10+Q[1,L-1]] + Vr*conj(X[3]+X[7])*X[10+L])
    for i in 2 : L-1
        dX[10+Q[i,i]]=2*imag(-J*X[10+Q[i-1,i]] + X[10+Q[i,i+1]])
        for j in i+1 : L-1
            dX[10+Q[i,j]]=1im*(J*X[10+Q[i-1,j]]+J*X[10+Q[i+1,j]] - X[10+Q[i,j-1]] - X[10+Q[i,j+1]])
        end
        dX[10+Q[i,L]]=1im*(J*X[10+Q[i-1,L]] + J*X[10+Q[i+1,L]] - X[10+Q[i,L-1]])
    end
    dX[10+Q[L,L]]=-2*imag(J*X[10+Q[L-1,L]])



    dX[10+LC+P[1,2]]=-1im*(J*X[10+LC+P[1,3]] - 2* μ_c*X[10+LC+P[1,2]] + Vr*(X[3]+X[7])*X[12])
    for j in 3:L-1
        dX[10+LC+P[1,j]]=-1im*(J*X[10+LC+P[2,j]] + J*X[10+LC+P[1,j-1]] + J*X[10+LC+P[1,j+1]]
                          - 2* μ_c*X[10+LC+P[1,j]] + Vr*(X[3]+X[7])*X[10+j])
    end
    dX[10+LC+P[1,L]]=-1im*(J*X[10+LC+P[2,L]] + J*X[10+LC+P[1,L-1]]
                      - 2* μ_c*X[10+LC+P[1,L]] + Vr*(X[3]+X[7])*X[10+L])
    for i in 2 : L-2
        dX[10+LC+P[i,i+1]]=-1im*(J*X[10+LC+P[i-1,i+1]]+J*X[10+LC+P[i,i+2]] - 2* μ_c*X[10+LC+P[i,i+1]])
        for j in i+2 : L-1
            dX[10+LC+P[i,j]]=-1im*(J*X[10+LC+P[i-1,j]] + J*X[10+LC+P[i+1,j]] + J*X[10+LC+P[i,j-1]]
                              + J*X[10+LC+P[i,j+1]] - 2* μ_c*X[10+LC+P[i,j]])
        end
        dX[10+LC+P[i,L]]=-1im*(J*X[10+LC+P[i-1,L]] + J*X[10+LC+P[i+1,L]] + J*X[10+LC+P[i,L-1]] - 2* μ_c*X[10+LC+P[i,L]])
    end
    dX[10+LC+P[L-1,L]]=-1im*(J*X[10+LC+P[L-2,L]] - 2* μ_c*X[10+LC+P[L-1,L]])


    dXin[1:NH] = real(dX)
    dXin[NH+1:end] = imag(dX)


    return nothing
end
