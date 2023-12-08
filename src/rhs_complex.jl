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

    # wird zu:
    #XInd[11] = [11+1,11,3,7,10+Qrr,3,7]
    #dX[11]=-1im*(J*X[XInd[11][1]] + - μ_c*X[XInd[11][2]] - 2*Vr*(X[XInd[11][3]]+X[XInd[11][4]])*X[XInd[11][5]] + Vr*(X[XInd[11][6]]+X[XInd[11][7]]))


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


function rhs!(dX::Vector, X::Vector, p::Vector, t::Float64)::Nothing
    
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

    ##################################################################################################
    ####################### This Part has to be moved outside! I can't manage! #######################
    ##################################################################################################
        # TODO: Auslagern
    Q = UpperTriangular(Int[(L+1)*i - (i+1)*i/2 + j for i in 1:L, j in 1:L])
    P = UpperTriangular(Int[(i-1)*L - (i+1)*i/2 + j for i in 1:L, j in 1:L])
    #S = UpperTriangular(Int[(i-1)*4 - (i-1)*i/2 + j for i in 1:4, j in 1:4]);  We do not neet it !
    ##################################################################################################
    ##################################################################################################

    #  rho[1,1]
    dX[1] =-2 * imag(conj(Vr * X[11]) *( X[2]  + X[3]))
    #  rho[1,2]
    dX[2] = -1im*( (U- μ_imp + ϵ_imp)*X[2] + Vr * X[11] * (X[5] + conj(X[6]) - X[1]) -conj(Vr * X[11]) * X[4])
    #  rho[1,3]
    dX[3] = -1im*( (U- μ_imp + ϵ_imp)*X[3] + Vr * X[11] * (X[6] + X[8] - X[1]) -conj(Vr * X[11]) * X[4])
    #  rho[1,4]
    dX[4] = -1im*( (U- 2*(μ_imp - ϵ_imp))*X[4] + Vr * X[11] * (X[7] + X[9] - X[2] -X[3]))
    #  rho[2,2]
    dX[5] = 2 * imag(conj(Vr * X[11]) *( X[2]  - X[7]))
    #  rho[2,3]
    dX[6] = -1im*(conj(Vr * X[11]) * (X[3] -X[7]) + Vr * X[11] * conj(X[9] - X[2]) )
    #  rho[2,4]
    dX[7] = -1im*(conj(Vr * X[11]) * X[4] - (μ_imp - ϵ_imp) * X[7] + Vr * X[11] * (X[10] -X[5] - X[6]))
    #  rho[3,3]
    dX[8] =2 * imag(conj(Vr * X[11]) *( X[3]  - X[9]))
    #  rho[3,4]
    dX[9] = -1im*(conj(Vr * X[11]) * X[4] - (μ_imp - ϵ_imp) * X[9] + Vr * X[11] * (X[10] -conj(X[6]) - X[8]))
    #  rho[4,4]
    #dX[9] = -1im*(conj(Vr*X[11])*X[4]  + (-μ_imp + ϵ_imp )*X[9] +Vr*X[11]*real(X[10]) - conj(X[6])*Vr*X[11] - real(X[8])*Vr*X[11])
    dX[10] =2 * imag(conj(Vr * X[11]) *( X[7]  + X[9]))



    # C[0,1]
    dX[11]=-1im*(J*X[12] - μ_c*X[11] - 2*Vr*(X[3]+X[7])*X[10+Q[1,1]] + Vr*(X[3]+X[7]))
    # C[0,j] 1<j<L
    for j in 2:L-1
        dX[10+j]=-1im*(J*X[9+j] + J*X[11+j] - μ_c*X[10+j]
                 - 2*Vr*(X[3]+X[7])*X[10+Q[1,j]] + 2*Vr*conj(X[3]+X[7])*X[10+LC+P[1,j]])
    end
    # C[0,L]
    dX[10+L]=-1im*(J*X[9+L] - μ_c*X[10+L]
            - 2*Vr*(X[3]+X[7])*X[10+Q[1,L]] + 2*Vr*conj(X[3]+X[7])*X[10+LC+P[1,L]])




    # C[1,1]
    dX[10+Q[1,1]]= -2 * imag(-J*X[10+Q[1,2]] + Vr*conj(X[3]+X[7])*X[11])
    # C[1,j] 1<j<L
    for j in 2 : L-1
        dX[10+Q[1,j]]=1im*(J*X[10+Q[2,j]]-J*X[10+Q[1,j-1]]-J*X[10+Q[1,j+1]]+ Vr*conj(X[3]+X[7])*X[10+j])
    end

    # C[1,L]
    dX[10+Q[1,L]]=1im*(J*X[10+Q[2,L]] - J*X[10+Q[1,L-1]] + Vr*conj(X[3]+X[7])*X[10+L])



    for i in 2 : L-1
        # C[i,i]
        dX[10+Q[i,i]]=2*imag(-J*X[10+Q[i-1,i]] + X[10+Q[i,i+1]])
        # C[i,j] 1<i<j<L
        for j in i+1 : L-1
            dX[10+Q[i,j]]=1im*(J*X[10+Q[i-1,j]]+J*X[10+Q[i+1,j]] - X[10+Q[i,j-1]] - X[10+Q[i,j+1]])
        end
        # C[i,L] 1<i<L
        dX[10+Q[i,L]]=1im*(J*X[10+Q[i-1,L]] + J*X[10+Q[i+1,L]] - X[10+Q[i,L-1]])
    end
    # C[L,L]
    dX[10+Q[L,L]]=-2*imag(J*X[10+Q[L-1,L]])



    # K[1,2]
    dX[10+LC+P[1,2]]=-1im*(J*X[10+LC+P[1,3]] - 2* μ_c*X[10+LC+P[1,2]] + Vr*(X[3]+X[7])*X[12])
    # K[1,j]
    for j in 3:L-1
        dX[10+LC+P[1,j]]=-1im*(J*X[10+LC+P[2,j]] + J*X[10+LC+P[1,j-1]] + J*X[10+LC+P[1,j+1]]
                          - 2* μ_c*X[10+LC+P[1,j]] + Vr*(X[3]+X[7])*X[10+j])
    end
    # K[1,L]
    dX[10+LC+P[1,L]]=-1im*(J*X[10+LC+P[2,L]] + J*X[10+LC+P[1,L-1]]
                      - 2* μ_c*X[10+LC+P[1,L]] + Vr*(X[3]+X[7])*X[10+L])

    for i in 2 : L-2
        # K[i,i+1]  1<i<L-2
        dX[10+LC+P[i,i+1]]=-1im*(J*X[10+LC+P[i-1,i+1]]+J*X[10+LC+P[i,i+2]] - 2* μ_c*X[10+LC+P[i,i+1]])
        # K[i,j]    1<i<j<L
        for j in i+2 : L-1
            dX[10+LC+P[i,j]]=-1im*(J*X[10+LC+P[i-1,j]] + J*X[10+LC+P[i+1,j]] + J*X[10+LC+P[i,j-1]]
                              + J*X[10+LC+P[i,j+1]] - 2* μ_c*X[10+LC+P[i,j]])
        end
        # [i, L]   1<i
        dX[10+LC+P[i,L]]=-1im*(J*X[10+LC+P[i-1,L]] + J*X[10+LC+P[i+1,L]] + J*X[10+LC+P[i,L-1]] - 2* μ_c*X[10+LC+P[i,L]])
    end
    # [L-1, L]
    dX[10+LC+P[L-1,L]]=-1im*(J*X[10+LC+P[L-2,L]] - 2* μ_c*X[10+LC+P[L-1,L]])

    return nothing
end
