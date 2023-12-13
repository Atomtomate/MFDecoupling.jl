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

