using LinearAlgebra

# X0=Float64[]


# TODO: Auslagern
    

"""
    test!(dX::Vector, X::Vector, p, t)::Nothing

RHS for DGL.

Parameters 
-------
  - *p[1]*:  L,    size of the chain
  - *p[2]*:  U,    local Hubbard interaction on the impurity site
  - *p[3]*:  Vr,   coupling between the impurity and chain
  - *p[4]*:  μ_imp chemical potential for impurity
  - *p[5]*:  μ_c   chemical potential for chain (onsite energy shift is includid in the chemical potential)
"""
function test!(dX::Vector, X::Vector, p::Vector, t::Float64)::Nothing
    
    r::Int         = 1  # Must be one in this set up!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    J::Float64     = 1.0
    L::Int         = Int(p[1])
    U::Float64     = p[2]
    ϵ_imp::Float64 = U/2
    Vr::Float64    = p[3]
    μ_imp::Float64 = p[4]
    μ_c::Float64   = p[5]
    LC::Int     =floor(Int,(L+3)*L/2)
    LK::Int     = floor(Int,(L-1)*L/2)
    
    H_imp  = Hermitian([
             U  + 2(μ_imp-ϵ_imp)  Vr*X[10+r]         Vr*X[10+r]        0         ;
             conj(Vr*X[10+r])    μ_imp - ϵ_imp      0                 Vr*X[10+r];
             conj(Vr*X[10+r])    0                  μ_imp - ϵ_imp     Vr*X[10+r];
             0                   conj(Vr*X[10+r])   conj(Vr*X[10+r])  0     
             ])
    
    ρ  = Hermitian([X[1]         X[2]         X[3]         X[4]; 
                    conj(X[2])   X[5]         X[6]         X[7];
                    conj(X[3])   conj(X[6])   X[8]         X[9];
                    conj(X[4])   conj(X[7])   conj(X[9])   X[10]])
       
    # rho[i,j]
    dρ = -1im .* (H_imp * ρ - ρ * H_imp)
    dX[1:4] = dρ[1,:]
    dX[5:7] = dρ[2,2:end]
    dX[8:9] = dρ[3,3:end]
    dX[10]  = dρ[4,4]

    
    # C[0,r]
    Qrr=floor(Int, (L+1)*r - (r+1)*r/2 + r)
    dX[10+r]=-1im*(J*X[10+r+1] + μ_c*X[10+r] - 2*Vr*(X[3]+X[7])*X[10+Qrr] + Vr*(X[3]+X[7]))
        
    # C[0,j] 1<j<L
    for j in 2:L-1
        Qrj=floor(Int,(L+1)*r - (r+1)*r/2 + j)
        Prj=floor(Int,(r-1)*L - (r+1)*r/2 + j)
        dX[10+j]=-1im*(J*X[10+j-1] + J*X[10+j+1] + μ_c*X[10+j] 
                 - 2*Vr*(X[3]+X[7])*X[10+Qrj] + 2*Vr*conj(X[3]+X[7])*X[10+LC+Prj])
    end
    
    # C[0,L]

    QrL=floor(Int,(L+1)*r - (r+1)*r/2 + L)
    PrL=floor(Int,(r-1)*L - (r+1)*r/2 + L)
    dX[10+L]=-1im*(J*X[10+L-1] + μ_c*X[10+L] 
            - 2*Vr*(X[3]+X[7])*X[10+QrL] + 2*Vr*conj(X[3]+X[7])*X[10+LC+PrL])
     
    
    # if !all(dX[1:10+L] .≈ dX_test[1:10+L])
    #     error("test dX not equal to old!!!!!")
    # end

        
    # C[r,r]
    #Qrr=floor(Int,(L+1)*r - (r+1)*r/2 + r) already defined
    Qrrp1=floor(Int,(L+1)*r - (r+1)*r/2 + r+1)
    dX[10+Qrr]= -2 * imag(-J*X[10+Qrrp1] + Vr*conj(X[3]+X[7])*X[10+r])
    
    
    # C[r,j] 1<j<L
    for j in 2 : L-1
        Qrj=floor(Int,(L+1)*r - (r+1)*r/2 + j)
        Qrp1j=floor(Int,(L+1)*(r+1) - (r+1+1)*(r+1)/2 + j)
        Qrjm1=floor(Int,(L+1)*r - (r+1)*r/2 + j-1)
        Qrjp1=floor(Int,(L+1)*r - (r+1)*r/2 + j+1)
        dX[10+Qrj]=1im*(J*X[10+Qrp1j]-J*X[10+Qrjm1]-J*X[10+Qrjp1]+ Vr*conj(X[3]+X[7])*X[10+j]) 
    end
    
    # C[r,L]
    QrL=floor(Int,(L+1)*r - (r+1)*r/2 + L)
    Qrp1L=floor(Int,(L+1)*(r+1) - (r+1+1)*(r+1)/2 + L)
    QrLm1=floor(Int,(L+1)*r - (r+1)*r/2 + L-1)
    dX[10+QrL]=1im*(J*X[10+Qrp1L] - J*X[10+QrLm1] + Vr*conj(X[3]+X[7])*X[10+L]) 
    
    
    
    for i in 2 : L-1
        # C[i,i] 
        Qii=floor(Int,(L+1)*i - (i+1)*i/2 + i)
        Qim1i=floor(Int,(L+1)*(i-1) - (i-1+1)*(i-1)/2 + i)
        Qiip1=floor(Int,(L+1)*i - (i+1)*i/2 + i+1)
        dX[10+Qii]=2*imag(-J*X[10+Qim1i] + X[10+Qiip1]) 
               
        # C[i,j] 1<i<j<L
        for j in i+1 : L-1
            Qij=floor(Int,(L+1)*i - (i+1)*i/2 + j)    
            Qim1j=floor(Int,(L+1)*(i-1) - (i-1+1)*(i-1)/2 + j)
            Qip1j=floor(Int,(L+1)*(i+1) - (i+1+1)*(i+1)/2 + j)
            Qijm1=floor(Int,(L+1)*i - (i+1)*i/2 + j-1)    
            Qijp1=floor(Int,(L+1)*i - (i+1)*i/2 + j+1)    
            dX[10+Qij]=1im*(J*X[10+Qim1j]+J*X[10+Qip1j] - X[10+Qijm1] - X[10+Qijp1]) 
        end
        
        # C[i,L] 1<i<L
        QiL=floor(Int,(L+1)*i - (i+1)*i/2 + L)        
        Qim1L=floor(Int,(L+1)*(i-1) - (i-1+1)*(i-1)/2 + L)
        Qip1L=floor(Int,(L+1)*(i+1) - (i+1+1)*(i+1)/2 + L)
        QiLm1=floor(Int,(L+1)*i - (i+1)*i/2 + L-1)
        dX[10+QiL]=1im*(J*X[10+Qim1L] + J*X[10+Qip1L] - X[10+QiLm1]) 
    end
    
    # C[L,L]
    QLL=floor(Int,(L+1)*L - (L+1)*L/2 + L)
    QLm1L=floor(Int,(L+1)*(L-1) - (L-1+1)*(L-1)/2 + L)          
    dX[10+QLL]=-2*imag(J*X[10+QLm1L]) 
    

    
    # K[r,r+1]
    Prrp1=floor(Int,(r-1)*L - (r+1)*r/2 + r+1)
    Prrp1p1=floor(Int,(r-1)*L - (r+1)*r/2 + r+1+1) 
    dX[10+LC+Prrp1]=-1im*(J*X[10+LC+Prrp1p1] + 2*μ_c*X[10+LC+Prrp1] + Vr*(X[3]+X[7])*X[10+r+1]) 
    
    # K[r,j]
    for j in 3:L-1
        Prj=floor(Int,(r-1)*L - (r+1)*r/2 + j)
        Prp1j=floor(Int,(r+1-1)*L - (r+1+1)*(r+1)/2 + j)
        Prjm1=floor(Int,(r-1)*L - (r+1)*r/2 + j-1)
        Prjp1=floor(Int,(r-1)*L - (r+1)*r/2 + j+1)            
        dX[10+LC+Prj]=-1im*(J*X[10+LC+Prp1j] + J*X[10+LC+Prjm1] + J*X[10+LC+Prjp1] 
                          + 2*μ_c*X[10+LC+Prj] + Vr*(X[3]+X[7])*X[10+j])
    end
    
    # K[r,L]
    PrL=floor(Int,(r-1)*L - (r+1)*r/2 + L)
    Prp1L=floor(Int,(r+1-1)*L - (r+1+1)*(r+1)/2 + L)
    PrLm1=floor(Int,(r-1)*L - (r+1)*r/2 + L-1)                    
    dX[10+LC+PrL]=-1im*(J*X[10+LC+Prp1L] + J*X[10+LC+PrLm1] 
                      + 2*μ_c*X[10+LC+PrL] + Vr*(X[3]+X[7])*X[10+L]) 
    
    for i in 2 : L-2
        # K[i,i+1]  1<i<L-2    
        Piip1=floor(Int,(i-1)*L - (i+1)*i/2 + i+1)
        Pim1ip1=floor(Int,(i-1-1)*L - (i-1+1)*(i-1)/2 + i+1)    
        Piip1p1=floor(Int,(i-1)*L - (i+1)*i/2 + i+1+1)                
        dX[10+LC+Piip1]=-1im*(J*X[10+LC+Pim1ip1]+J*X[10+LC+Piip1p1]+2*μ_c*X[10+LC+Piip1])
        
        # K[i,j]    1<i<j<L
        for j in i+2 : L-1
            Pij=floor(Int,(i-1)*L - (i+1)*i/2 + j)
            Pim1j=floor(Int,(i-1-1)*L - (i-1+1)*(i-1)/2 + j)                
            Pip1j=floor(Int,(i+1-1)*L - (i+1+1)*(i+1)/2 + j)
            Pijm1=floor(Int,(i-1)*L - (i+1)*i/2 + j-1)
            Pijp1=floor(Int,(i-1)*L - (i+1)*i/2 + j+1)                
            dX[10+LC+Pij]=-1im*(J*X[10+LC+Pim1j] + J*X[10+LC+Pip1j] + J*X[10+LC+Pijm1]
                              + J*X[10+LC+Pijp1] + 2*μ_c*X[10+LC+Pij])     
        end
        
        # [i, L]   1<i
        PiL=floor(Int,(i-1)*L - (i+1)*i/2 + L)
        Pim1L=floor(Int,(i-1-1)*L - (i-1+1)*(i-1)/2 + L)                
        Pip1L=floor(Int,(i+1-1)*L - (i+1+1)*(i+1)/2 + L)
        PiLm1=floor(Int,(i-1)*L - (i+1)*i/2 + L-1)                     
        dX[10+LC+PiL]=-1im*(J*X[10+LC+Pim1L] + J*X[10+LC+Pip1L] + J*X[10+LC+PiLm1] + 2*μ_c*X[10+LC+PiL])         
    end
    
    # [L-1, L]  
    PLm1L=floor(Int,(L-1-1)*L - (L-1+1)*(L-1)/2 + L)
    PLm1m1L=floor(Int,(L-1-1-1)*L - (L-1-1+1)*(L-1-1)/2 + L)                        
    dX[10+LC+PLm1L]=-1im*(J*X[10+LC+PLm1m1L] + 2*μ_c*X[10+LC+PLm1L]) 

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
