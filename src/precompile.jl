# ==================================================================================================== #
#                                             rhs.jl                                                   #
# ---------------------------------------------------------------------------------------------------- #
#   Author          : Irakli Titvinidze, Julian Stobbe                                                 #
# ----------------------------------------- Description ---------------------------------------------- #
#   Precompilation workload for faster startup times.                                                  #
# -------------------------------------------- TODO -------------------------------------------------- #
#   Replace fp path by small test data                                                                 #
# ==================================================================================================== #

@setup_workload begin
    #TODO: local, very small, data here!!!!
    fp = "/scratch/projects/hhp00048/MFDecoupling/UScan/data_01"#ARGS[1]
    fpout = "tmp"
    fp1 = joinpath(fp,"CK_U0.1V0.5.dat") #ARGS[1]
    fp2 = joinpath(fp,"rho_U0.1V0.5.dat")

    const LL::Int = 1000
    Uin::Float64 = 0.1
    Vin::Float64 = 0.5
    Ufi::Float64 = 0.5
    Vfi::Float64 = 0.5
    tmin::Float64 = 0.0
    tmax::Float64 = 0.2
    tspan = (tmin,tmax)
    tsave = LinRange(0.1,0.2,2)



    X0 = read_inputs(fp1, fp2, LL)

    @compile_workload begin
        function solve_time_evolution(U::Float64, V::Float64, X0_in::Vector, L::Int, tspan, tsave, index, fpout, use_real)
            Q,P,LC,LK,LIm = MFDecoupling.gen_helpers(L)
            #X0, rhsf! = setup_calculation(X0, L; mode=:real)
            X0, rhsf! = setup_calculation(X0_in, L; mode= (use_real ? :real : :complex))
            p_0  = [L, U, V, 0.0, 0.0];

            linsolve = MFDecoupling.KrylovJL_GMRES()
            alg = MFDecoupling.AutoTsit5(MFDecoupling.KenCarp47(linsolve = linsolve, autodiff=use_real))
            prob = MFDecoupling.ODEProblem(rhsf!,X0,tspan,p_0)
            #idxs_list = union(collect(1:11),LIm .+ collect(1:11)) 
            idxs_list = collect(1:11)
            @time sol = MFDecoupling.solve(prob, alg; save_idxs=idxs_list, saveat=tsave, abstol=1e-9, reltol=1e-9);
            return nothing #U, sol
        end

        solve_time_evolution(Ufi, Vfi, X0, LL, tspan, tsave, 1, fpout, false)
    end
end
