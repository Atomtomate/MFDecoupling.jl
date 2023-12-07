using DifferentialEquations
using LinearAlgebra, DelimitedFiles
using Logging: global_logger
using TerminalLoggers: TerminalLogger
using ProgressMeter
using JLD2

global_logger(TerminalLogger());


fp1 = ARGS[1]
fp2 = ARGS[2]

LL = parse(Int, ARGS[3])
Uin = parse(Float64, ARGS[4]);
Vin = parse(Float64, ARGS[5]);
UU = parse(Float64, ARGS[6])
VV = parse(Float64, ARGS[7]);
tmin = parse(Float64, ARGS[8]);
tmax = parse(Float64, ARGS[9]);
outputf = ARGS[10]


tspan = (tmin,tmax)


println("Reading File 1")
CK= open(fp1, "r") do f
    readdlm(f)
end;

println("Reading File 2")
ρ = open(fp2, "r") do f
    arr = readdlm(f)
end;

C=zeros(Float64, LL+1, LL+1)

println("LL = $LL")
println(size(CK))
C[2:LL+1,2:LL+1]=deepcopy(CK[1:LL,1:LL]);
C[1,2:LL+1]=deepcopy(CK[2*LL+1,1:LL]+CK[2*LL+2,1:LL])
C[2:LL+1,1]=deepcopy(CK[1:LL,2*LL+1]+CK[1:LL,2*LL+2]);

K=CK[LL+1:2*LL,1:LL];

lC = size(C,1)
lK = size(K,1)
lρ = size(ρ,1)

# X0=Float64[]
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

dX = similar(X0);

# TODO: Auslagern
const Q = UpperTriangular(Int[(LL+1)*i - (i+1)*i/2 + j for i in 0:LL, j in 0:LL])
const P = UpperTriangular(Int[(i-1)*LL - (i+1)*i/2 + j for i in 1:LL, j in 1:LL])
const S = UpperTriangular(Int[(i-1)*4 - (i-1)*i/2 + j for i in 1:4, j in 1:4]);

const LC::Int     = floor(Int,(LL+3)*LL/2)    # Number of C
const LK::Int     = floor(Int,(LL-1)*LL/2)    # Number of K``
    

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
    
    
    H_imp  = Hermitian([
             UU + 2(μ_imp-ϵ_imp)  Vr*X[10+r]         Vr*X[10+r]        0         ; 
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




linsolve = KrylovJL_GMRES()
alg_expl = AutoTsit5(Rosenbrock23(autodiff=false))
    #SSPSDIRK2(autodiff=false)
alg_impl = AutoTsit5(ImplicitEuler(autodiff=false, linsolve = KrylovJL_GMRES()))
    #KenCarp47(linsolve = KrylovJL_GMRES(), autodiff=false)
    #

p_0  = [LL, UU, VV, 0.0, 0.0];

prob = ODEProblem{true}(test!,X0,tspan,p_0,
    progress = true,
    progress_steps = 1)

sol = solve(prob, alg_impl; save_everystep = true, abstol=1e-8, reltol=1e-8);

println("Errors:\n",sol.errors)
println(" =========================== ")
println("Algorithm Details:\n",sol.alg)
println(" =========================== ")
println("Solution stats:\n", sol.stats)

jldopen(outputf,"w") do f
  f["solution"] = sol.u
  f["t"] = sol.t
  f["Uin"] = Uin
  f["Vin"] = Vin
  f["UU"] = UU
  f["VV"]=VV
  f["L"]=LL
end

